module BluffBodySolverIncomp_mod
   use mod_arrays
   use mod_nvtx
#ifndef NOACC
   use cudafor
#endif
   

   use elem_qua
   use elem_hex
   use jacobian_oper
   use quadrature_rules
   use mod_inicond_reader
   use mass_matrix
   use mod_geom
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
   use CFDSolverPeriodicWithBoundariesIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BluffBodySolverIncomp

      real(rp) , public  :: vo, M, delta, rho0, Re

   contains
      procedure, public :: fillBCTypes           =>BluffBodySolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => BluffBodySolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer => BluffBodySolverIncomp_initialBuffer
   end type BluffBodySolverIncomp
contains

   subroutine BluffBodySolverIncomp_fill_BC_Types(this)
      class(BluffBodySolverIncomp), intent(inout) :: this

      bouCodes2BCType(1) = bc_type_far_field
      bouCodes2BCType(2) = bc_type_outlet_incomp
      bouCodes2BCType(3) = bc_type_far_field
      bouCodes2BCType(4) = bc_type_far_field
      bouCodes2BCType(5) = bc_type_non_slip_adiabatic
      !$acc update device(bouCodes2BCType(:))

   end subroutine BluffBodySolverIncomp_fill_BC_Types

   subroutine BluffBodySolverIncomp_initializeParameters(this)
      class(BluffBodySolverIncomp), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cylin"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      ! numerical params
      flag_les = 1
      maxIter = 20
      tol=1e-3
   

      this%cfl_conv = 0.95_rp 
      this%cfl_diff = 0.95_rp 
      !flag_use_constant_dt = 1
      !this%dt = 2.5e-4
      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 1000000 

      this%save_logFile_first = 1
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 10000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%loadRestartFile = .false.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%saveAvgFile = .true.
      this%loadAvgFile = .false.
      !----------------------------------------------

      this%vo = 1.0_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%Re     =  10000.0_rp

      incomp_viscosity = (this%rho0*this%delta*this%vo)/this%Re
      flag_mu_factor = 1.0_rp

      nscbc_u_inf = this%vo
      nscbc_p_inf = 0.0_rp
      nscbc_rho_inf = this%rho0

      !Witness points parameters
      this%have_witness          = .false.
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
     flag_buffer_e_min = 40.0_rp
     flag_buffer_e_size = 10.0_rp

     flag_buffer_on_west = .true.
     flag_buffer_w_min = -20.0_rp
     flag_buffer_w_size = 10.0_rp

     flag_buffer_on_north = .true.
     flag_buffer_n_min = 20.0_rp
     flag_buffer_n_size = 10.0_rp

     flag_buffer_on_south = .true.
     flag_buffer_s_min = -20.0_rp
     flag_buffer_s_size = 10.0_rp

   end subroutine BluffBodySolverIncomp_initializeParameters

   subroutine BluffBodySolverIncomp_evalInitialConditions(this)
      class(BluffBodySolverIncomp), intent(inout) :: this
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
         pr(iNodeL,2) = 0.0_rp
         rho(iNodeL,2) = this%rho0
      end do
      !$acc end parallel loop

      !$acc kernels
      mu_e(:,:) = 0.0_rp 
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels
      call nvtxEndRange

   end subroutine BluffBodySolverIncomp_evalInitialConditions


   subroutine BluffBodySolverIncomp_initialBuffer(this)
      class(BluffBodySolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolverIncomp_initialBuffer

end module BluffBodySolverIncomp_mod
