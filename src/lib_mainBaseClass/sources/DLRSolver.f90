module DLRSolver_mod
   use mod_arrays
   use mod_nvtx
#ifndef NOACC
   use cudafor
#endif
   use mod_veclen

   use elem_qua
   use elem_hex
   use jacobian_oper
   use quadrature_rules
   use mesh_reader
   use inicond_reader
   use mass_matrix
   use mod_geom
   use mod_output
   use mod_period
   use time_integ
   use mod_analysis
   use mod_numerical_params
   use mod_time_ops
   use mod_fluid_viscosity
   use mod_postpro
   use mod_aver
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use CFDSolverPeriodicWithBoundaries_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: DLRSolver

      real(rp), public :: vo, M, delta, rho0, Re, to, po
      character(len=8), public :: tag, restart_step

   contains
      procedure, public :: fillBCTypes           => DLRSolver_fill_BC_Types
      procedure, public :: initializeParameters  => DLRSolver_initializeParameters
      procedure, public :: evalInitialConditions => DLRSolver_evalInitialConditions
      procedure, public :: initialBuffer => DLRSolver_initialBuffer
      procedure, public :: afterDt => DLRSolver_afterDt
   end type DLRSolver
contains

   subroutine DLRSolver_fill_BC_Types(this)
      class(DLRSolver), intent(inout) :: this

      bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! cylin
      bouCodes2BCType(2) = bc_type_non_slip_adiabatic ! wall
      bouCodes2BCType(3) = bc_type_far_field !in
      bouCodes2BCType(4) = bc_type_far_field !out
      bouCodes2BCType(5) = bc_type_unsteady_inlet !jet1
      bouCodes2BCType(6) = bc_type_unsteady_inlet !jet2
      bouCodes2BCType(7) = bc_type_slip_adiabatic !3d

   end subroutine DLRSolver_fill_BC_Types

   subroutine DLRSolver_initializeParameters(this)
      class(DLRSolver), intent(inout) :: this
      real(rp) :: mul, mur
      integer :: num_args, equal_pos, iarg
      character(len=64) :: arg, output_dir
      logical :: output_dir_exists

      ! get command line args, ie: mpirun -n 4 sod2d --tag=12 --restart_step=2500
      this%tag = "1"
      this%restart_step = ""
      num_args = command_argument_count()
      do iarg = 1, num_args
         call get_command_argument(iarg, arg)
         equal_pos = scan(adjustl(trim(arg)), "=")
         if (adjustl(trim(arg(:equal_pos-1))) .eq. "--tag") then
            this%tag = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--restart_step") then
            this%restart_step = trim(adjustl(arg(equal_pos+1:)))
         else
            stop "Unknown command line argument"
         end if
      end do

      ! create output dir if not existing
      output_dir = "./output_"//trim(adjustl(this%tag))//"/"
      inquire(file=trim(adjustl(output_dir)), exist=output_dir_exists)
      if (.not. output_dir_exists) call execute_command_line("mkdir -p "//trim(adjustl(output_dir)))

      ! cold start or restart
      if (this%restart_step == "") then
         this%loadResults = .false.
         this%continue_oldLogs = .false.
         this%load_step = 150000
      else
         this%loadResults = .true.
         this%continue_oldLogs = .true.
         read(this%restart_step, *) this%load_step
      end if

      ! sod paths
      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cylinder"

      write(this%results_h5_file_path,*) "./output_"//trim(adjustl(this%tag))//"/"
      write(this%results_h5_file_name,*) "results_"//trim(adjustl(this%tag))

      ! numerical params
      flag_les = 0
      flag_implicit = 0
      flag_les_ilsa=0
      implicit_solver = implicit_solver_bdf2_rk10
      flag_rk_order=4

      pseudo_cfl =0.95_rp
      pseudo_ftau= 4.0_rp
      maxIterNonLineal=1000
      tol=1e-3

      this%nstep = 90000001 !250001
      this%cfl_conv = 0.25_rp
      this%cfl_diff = 0.25_rp

      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1

      this%nleap = 20000 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 50  ! Saving interval, TODO: input
      this%nleapAVG = 20000000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 1.0_rp
      this%M  = 0.2_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp
      this%Re     =  100.0_rp

      mul = (this%rho0*this%delta*this%vo)/this%Re
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = mul/mur

      nscbc_u_inf = this%vo
      nscbc_p_inf = this%po
      nscbc_rho_inf = this%rho0
      nscbc_gamma_inf = this%gamma_gas
      nscbc_c_inf = sqrt(this%gamma_gas*this%po/this%rho0)
      nscbc_Rgas_inf = this%Rgas

      flag_buffer_on = .true.

      flag_buffer_on_east = .true.
      flag_buffer_e_min = 18.0_rp
      flag_buffer_e_size = 4.0_rp

      flag_buffer_on_west = .true.
      flag_buffer_w_min = 0.0_rp
      flag_buffer_w_size = 4.0_rp

   end subroutine DLRSolver_initializeParameters

   subroutine DLRSolver_evalInitialConditions(this)
      class(DLRSolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL,bcode
      logical :: readFiles

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo*(1.0-(((coordPar(iNodeL,2)-4.1_rp/2.0_rp)/4.1_rp)*2.0_rp)**2)
         u(iNodeL,2,2) = 0.0_rp
         u(iNodeL,3,2) = 0.0_rp
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%po/this%Rgas/this%to
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) = eta(iNodeL,2)
      end do
      !$acc end parallel loop

      !$acc kernels
      mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels
      call nvtxEndRange

   end subroutine DLRSolver_evalInitialConditions


   subroutine DLRSolver_initialBuffer(this)
      class(DLRSolver), intent(inout) :: this
      integer(4) :: iNodeL,bcode

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo*(1.0-(((coordPar(iNodeL,2)-4.1_rp/2.0_rp)/4.1_rp)*2.0_rp)**2)
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp
            if(bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
               bcode = bouCodesNodesPar(iNodeL) ! Boundary element code
               if (bcode == bc_type_unsteady_inlet .or. bcode == bc_type_non_slip_adiabatic) then
                  u_buffer(iNodeL,1) = 0.00_rp
                  u_buffer(iNodeL,2) = 0.00_rp
                  u_buffer(iNodeL,3) = 0.00_rp
                  u(iNodeL,1,2) = 0.00_rp
                  u(iNodeL,2,2) = 0.00_rp
                  u(iNodeL,3,2) = 0.00_rp
               end if
            end if
      end do
      !$acc end parallel loop

   end subroutine DLRSolver_initialBuffer

     subroutine DLRSolver_afterDt(this,istep)
      class(DLRSolver), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      integer(4) :: inode,bcode
      real(rp) :: x,y,z

      !$acc parallel loop
      do inode = 1,numNodesRankPar
         if(bouCodesNodesPar(inode) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(inode) ! Boundary element code
            if (bcode == bc_type_unsteady_inlet) then
               x = coordPar(inode,1)
               y = coordPar(inode,2)
               z = coordPar(inode,3)
               if(y>2.0_rp) then !top
                  u_buffer(inode,1) = real(0.06283185307179587,rp)*cos((2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2))) &
                                    - real(1.2217304763960306,rp))/real(0.05555555555555555,rp))*0.5_rp*cos(2*atan((y - 2.0_rp)/(x - 2.0_rp &
                                    + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2)))-real(1.2217304763960306,rp))
                  u_buffer(inode,2) = real(0.06283185307179587,rp)*cos((2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2))) &
                                    - real(1.2217304763960306,rp))/real(0.05555555555555555,rp))*0.5_rp*sin(2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp &
                                    + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2)))-real(1.2217304763960306,rp))
                  u_buffer(inode,3) = 0.00_rp
               else
                  u_buffer(inode,1) = real(0.06283185307179587,rp)*cos((2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2))) &
                                    + real(1.221730476396031,rp))/real(0.05555555555555555,rp))*0.5_rp*cos(2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp &
                                    + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2)))-real(1.221730476396031,rp))
                  u_buffer(inode,2) = real(0.06283185307179587,rp)*cos((2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2))) &
                                    + real(1.221730476396031,rp))/real(0.05555555555555555,rp))*0.5_rp*sin(2.0_rp*atan((y - 2.0_rp)/(x - 2.0_rp &
                                    + sqrt((x - 2.0_rp)**2+(y - 2.0_rp)**2)))-real(1.221730476396031,rp))
                  u_buffer(inode,3) = 0.00_rp
               end if
            end if
         end if
      end do
      !$acc end parallel loop
   end subroutine DLRSolver_afterDt

end module DLRSolver_mod
