module TGVSolver_mod
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
   use CFDSolverPeriodic_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodic) :: TGVSolver

      real(rp) , public  ::  M, rho0,Re,to, po

   contains
      procedure, public :: initializeParameters  => TGVSolver_initializeParameters
      procedure, public :: evalInitialConditions => TGVSolver_evalInitialConditions
   end type TGVSolver
contains

   subroutine TGVSolver_initializeParameters(this)
      class(TGVSolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cube_per265"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      this%doGlobalAnalysis = .true.
      this%doTimerAnalysis = .true.

      this%loadResults = .false.
      this%continue_oldLogs = .false.
      this%load_step = 8001

      this%nstep = 5001

      this%cfl_conv = 0.5_rp
      this%cfl_diff = 0.5_rp
      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input

      this%nsaveAVG = 1!1000000
      this%nleap = 5000 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 250  ! Saving interval, TODO: input
      this%nleapAVG = 2000000000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%M  = 0.1_rp
      this%Re = 1600.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp

      mul    = (this%rho0*1.0_rp*1.0_rp)/this%Re
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = 1.0_rp*1.0_rp/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = mul/mur

      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas

   end subroutine TGVSolver_initializeParameters

   subroutine TGVSolver_evalInitialConditions(this)
      class(TGVSolver), intent(inout) :: this
      real(4) :: iniU(totalNumNodesSrl,ndime), iniRho(totalNumNodesSrl), iniP(totalNumNodesSrl)
      real(4) :: V0,L
      real(4) :: x,y,z
      integer(4) :: iNodeL,iNodeGSrl

      if(mpi_rank.eq.0) write(*,*) "--| TGV - Setting Initial Conditions..."

      V0 = 1.0_rp
      L  = 1.0_rp

      !acc parallel loop
      do iNodeL=1,numNodesRankPar
         x = coordPar(iNodeL,1)
         y = coordPar(iNodeL,2)
         z = coordPar(iNodeL,3)

         u(iNodeL,1,2) =  V0*sin(x/(L))*cos(y/(L))*cos(z/(L))
         u(iNodeL,2,2) = -V0*cos(x/(L))*sin(y/(L))*cos(z/(L))
         u(iNodeL,3,2) = 0.0

         pr(iNodeL,2)  = this%po+((1.0_rp*V0*V0)/(16.0_rp))*(cos(2.0_rp*x/L)+cos(2.0_rp*y/L))*(cos(2.0_rp*z/L)+2.0_rp)
         rho(iNodeL,2) = pr(iNodeL,2)/this%Rgas/this%to
      end do
      !acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
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

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop
   end subroutine TGVSolver_evalInitialConditions

end module TGVSolver_mod
