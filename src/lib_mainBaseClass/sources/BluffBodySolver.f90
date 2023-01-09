module BluffBodySolver_mod
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
   use mod_constants
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BluffBodySolver

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po

   contains
      procedure, public :: fillBCTypes           =>BluffBodySolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolver_initializeParameters
      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditions
      procedure, public :: evalViscosityFactor=>BluffBodySolver_evalViscosityFactor
   end type BluffBodySolver
contains

   subroutine BluffBodySolver_fill_BC_Types(this)
      class(BluffBodySolver), intent(inout) :: this

      !bouCodes2BCType(1) = bc_type_inlet
      !bouCodes2BCType(2) = bc_type_inlet
      !bouCodes2BCType(3) = bc_type_inlet
      !bouCodes2BCType(4) = bc_type_inlet
      !bouCodes2BCType(5) = bc_type_non_slip_adiabatic

      bouCodes2BCType(1) = bc_type_non_slip_adiabatic
      bouCodes2BCType(2) = bc_type_inlet
   end subroutine BluffBodySolver_fill_BC_Types

   subroutine BluffBodySolver_initializeParameters(this)
      class(BluffBodySolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%gmsh_file_path,*) "./mesh/"
      write(this%gmsh_file_name,*) "cylin" 
      !write(this%gmsh_file_name,*) "naca" 

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cylin"
      !write(this%mesh_h5_file_name,*) "naca"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      this%isPeriodic = .true.
      this%loadMesh = .true.
      this%loadResults = .true.

      this%continue_oldLogs = .false.
      this%load_step = 50001

      this%nstep = 9000001 !250001
      this%cfl_conv = 2.0_rp !0.1_rp
      this%cfl_diff = 2.0_rp !0.1_rp

      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1

      this%nleap = 20000 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 50  ! Saving interval, TODO: input
      this%nleapAVG = 20000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 1.0_rp
      this%M  = 0.2_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp
      this%Re     =  10000.0_rp

      mul    = (this%rho0*this%delta*this%vo)/this%Re
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
      flag_buffer_x_min = 12.0_rp
      flag_buffer_x_max = 16.0_rp 

   end subroutine BluffBodySolver_initializeParameters

   subroutine BluffBodySolver_evalInitialConditions(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = 1.0_rp
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

   end subroutine BluffBodySolver_evalInitialConditions

   subroutine BluffBodySolver_evalViscosityFactor(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(4) :: iNodeL

      ! set out of the buffer zone
      ! remember that the mu_factor field has to be filled at least with the
      ! flag_mu_factor

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
   !     if(coordPar(iNodeL,1)<-12.0_rp) then
        !if(coordPar(iNodeL,1)<-8.0_rp) then
   !        mu_factor(iNodeL) = flag_mu_factor*100.0_rp
   !     end if
   !     if(coordPar(iNodeL,1)>12_rp) then
        !if(coordPar(iNodeL,1)>10_rp) then
   !        mu_factor(iNodeL) = flag_mu_factor*100.0_rp
   !     end if
   !     if(coordPar(iNodeL,2)<-12.0_rp) then
        !if(coordPar(iNodeL,2)<-10.0_rp) then
   !        mu_factor(iNodeL) = flag_mu_factor*100.0_rp
   !     end if
   !     if(coordPar(iNodeL,2)>12.0_rp) then
        !if(coordPar(iNodeL,2)>10.0_rp) then
   !        mu_factor(iNodeL) = flag_mu_factor*100.0_rp
   !     end if
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolver_evalViscosityFactor

end module BluffBodySolver_mod
