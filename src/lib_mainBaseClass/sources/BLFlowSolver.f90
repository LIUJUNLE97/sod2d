module BLFlowSolver_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BLFlowSolver

      real(rp) , public  ::  M, d0, U0, rho0, Red0, Re, to, po, mu

   contains
      procedure, public :: fillBCTypes           => BLFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BLFlowSolver_initializeParameters
      procedure, public :: evalInitialConditions => BLFlowSolver_evalInitialConditions
      procedure, public :: initialBuffer => BLFlowSolver_initialBuffer
   end type BLFlowSolver
contains

   subroutine BLFlowSolver_fill_BC_Types(this)
      class(BLFlowSolver), intent(inout) :: this

      !bouCodes2BCType(1) = bc_type_non_slip_adiabatic
      bouCodes2BCType(1) = bc_type_slip_wall_model
      bouCodes2BCType(2) = bc_type_outlet

   end subroutine BLFlowSolver_fill_BC_Types

   subroutine BLFlowSolver_initializeParameters(this)
      class(BLFlowSolver), intent(inout) :: this
      real(rp) :: mur

      write(this%gmsh_file_path,*) "./mesh/"
      write(this%gmsh_file_name,*) "bl"

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "bl"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      this%isPeriodic = .true.
      this%loadMesh = .true.

      this%continue_oldLogs = .false.
      this%load_step = 20001

      this%nstep = 9000000 
      this%cfl_conv = 0.95_rp
      this%cfl_diff = 0.95_rp
      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1
      this%nleap = 20000 ! Saving interval, TODO: input
      this%nleap2 = 10  ! Saving interval, TODO: input
      this%nleapAVG = 20000

      this%Cp   = 1004.0_rp
      this%Prt  = 0.71_rp
      this%M    = 0.2_rp
      this%d0   = 1.0_rp
      this%U0   = 1.0_rp
      this%rho0 = 1.0_rp

      this%Red0  = 450.0_rp
      this%gamma_gas = 1.40_rp

      this%mu    = this%rho0*this%d0*this%U0/this%Red0

      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%U0*this%U0/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = this%mu/mur

      nscbc_u_inf = this%U0
      nscbc_p_inf = this%po
      nscbc_rho_inf = this%rho0
      nscbc_gamma_inf = this%gamma_gas
      nscbc_c_inf = sqrt(this%gamma_gas*this%po/this%rho0)
      nscbc_Rgas_inf = this%Rgas


      flag_buffer_on = .true.
      ! x outlet
      flag_buffer_on_east = .true.
      flag_buffer_e_min = 2800.0_rp
      flag_buffer_e_size = 200.0_rp 

   end subroutine BLFlowSolver_initializeParameters

   subroutine BLFlowSolver_initialBuffer(this)
      class(BLFlowSolver), intent(inout) :: this
      integer(4) :: iNodeL
      real(rp) :: yp

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2) 
         u_buffer(iNodeL,1) = (-1.0_rp/160000.0_rp)*yp*yp + yp*(1.0_rp/200.0_rp)
         u_buffer(iNodeL,2) = 0.0_rp
         u_buffer(iNodeL,3) = 0.0_rp
      end do
      !$acc end parallel loop

   end subroutine BLFlowSolver_initialBuffer

   subroutine BLFlowSolver_evalInitialConditions(this)
      class(BLFlowSolver), intent(inout) :: this
      integer(4) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, ti(3), yp,velo_aux1
      integer(4)   :: iLine,iNodeGSrl,auxCnt

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2) 
         u(iNodeL,1,2) = (-1.0_rp/160000.0_rp)*yp*yp + yp*(1.0_rp/200.0_rp)
         u(iNodeL,2,2) = 0.0_rp
         u(iNodeL,3,2) = 0.0_rp
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%rho0
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
   end subroutine BLFlowSolver_evalInitialConditions

end module BLFlowSolver_mod
