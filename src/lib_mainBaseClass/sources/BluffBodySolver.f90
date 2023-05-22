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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BluffBodySolver

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po

   contains
      procedure, public :: fillBCTypes           =>BluffBodySolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolver_initializeParameters
      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditions
      procedure, public :: initialBuffer => BluffBodySolver_initialBuffer
   end type BluffBodySolver
contains

   subroutine BluffBodySolver_fill_BC_Types(this)
      class(BluffBodySolver), intent(inout) :: this

      bouCodes2BCType(1) = bc_type_far_field
      bouCodes2BCType(2) = bc_type_far_field
      bouCodes2BCType(3) = bc_type_far_field
      bouCodes2BCType(4) = bc_type_far_field
      bouCodes2BCType(5) = bc_type_non_slip_adiabatic
      !$acc update device(bouCodes2BCType(:))

   end subroutine BluffBodySolver_fill_BC_Types

   subroutine BluffBodySolver_initializeParameters(this)
      class(BluffBodySolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cylin"
      !write(this%mesh_h5_file_name,*) "naca"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      ! numerical params
      flag_les = 1
      flag_implicit = 1
      flag_les_ilsa=0
      implicit_solver = implicit_solver_bdf2_rk10
      flag_rk_order=4

      pseudo_cfl =1.95_rp 
      pseudo_ftau= 4.5_rp  !pseudo_cfl*pseudo_ftau<10 !important
      maxIterNonLineal=300
      tol=1e-3

      this%cfl_conv = 100.0_rp 
      this%cfl_diff = 100.0_rp 

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 1000000 

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 50000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 50000
      this%loadRestartFile = .false.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%saveAvgFile = .false.
      this%loadAvgFile = .false.
      !----------------------------------------------

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 1.0_rp
      this%M  = 0.2_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp
      this%Re     =  3300.0_rp
      !this%Re     =  100000.0_rp

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

      !Witness points parameters
      this%have_witness          = .true.
      this%witness_inp_file_name = "witness.txt"
      this%witness_h5_file_name  = "resultwit.h5"
      this%leapwit               = 1
      this%nwit                  = 17986
      this%wit_save_u_i          = .true.
      this%wit_save_pr           = .true.
      this%wit_save_rho          = .true.
      this%continue_witness      = .false.     
 
      flag_buffer_on = .true.
     !cylinder
     flag_buffer_on_east = .true.
     flag_buffer_e_min = 100.0_rp
     flag_buffer_e_size = 10.0_rp

     flag_buffer_on_west = .true.
     flag_buffer_w_min = -15.0_rp
     flag_buffer_w_size = 10.0_rp

     flag_buffer_on_north = .true.
     flag_buffer_n_min = 15.0_rp
     flag_buffer_n_size = 10.0_rp

     flag_buffer_on_south = .true.
     flag_buffer_s_min = -15.0_rp
     flag_buffer_s_size = 10.0_rp

      !naca
     !flag_buffer_on_east = .true.
     !flag_buffer_e_min = 10.0_rp
     !flag_buffer_e_size = 5.0_rp 

     !flag_buffer_on_west = .true.
     !flag_buffer_w_min = -10.0_rp
     !flag_buffer_w_size = 2.5_rp 
     
     !flag_buffer_on_north = .true.
     !flag_buffer_n_min = 10.0_rp
     !flag_buffer_n_size = 2.5_rp 
     
     !flag_buffer_on_south = .true.
     !flag_buffer_s_min = -10.0_rp
     !flag_buffer_s_size = 2.5_rp 

   end subroutine BluffBodySolver_initializeParameters

   subroutine BluffBodySolver_evalInitialConditions(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo
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

   end subroutine BluffBodySolver_evalInitialConditions


   subroutine BluffBodySolver_initialBuffer(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolver_initialBuffer

end module BluffBodySolver_mod
