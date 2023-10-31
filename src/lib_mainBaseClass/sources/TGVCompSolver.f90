module TGVCompSolver_mod
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
   use CFDSolverPeriodic_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodic) :: TGVCompSolver

      real(rp) , public  ::  Ma, V0,L,T0,P0,Re,rho0

   contains
      procedure, public :: initializeParameters  => TGVCompSolver_initializeParameters
      procedure, public :: evalInitialConditions => TGVCompSolver_evalInitialConditions
   end type TGVCompSolver
contains

   subroutine TGVCompSolver_initializeParameters(this)
      class(TGVCompSolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cube"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      write(this%io_append_info,*) ""

      this%doGlobalAnalysis = .true.
      this%doTimerAnalysis = .true.
      this%saveInitialField = .false.

      flag_implicit = 1

      maxIter = 200
      tol = 1e-3

      this%maxPhysTime = 20.0_rp
      this%final_istep = 10000 
      this%cfl_conv = 0.8_rp
      this%cfl_diff = 1000.0_rp

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%loadRestartFile = .false.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.
      this%save_restartFile_first = 1
      this%save_restartFile_step = 500

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 500

      this%Ma = 1.25_rp
      this%L = 1.0_rp
      this%V0 = 1.0_rp
      this%Prt = 0.71_rp
      this%T0  = 1.0_rp
      this%Re = 1600.0_rp
      this%gamma_gas = 1.40_rp
      
      this%Rgas = 1.0_rp*1.0_rp/(1.4_rp*1.0_rp*1.25_rp*1.25_rp)
      this%Cp = this%gamma_gas*this%Rgas/(this%gamma_gas-1.0_rp)
      this%rho0 = (1.0_rp/(1.4_rp*1.25*1.25))/(this%Rgas*this%T0)

      mul    = (this%rho0*1.0_rp*1.0_rp)/this%Re
      this%P0 = 1.0_rp/(this%Ma*this%Ma*1.4_rp)
      mur = 0.000001458_rp*(this%T0**1.50_rp)/(this%T0+110.40_rp)
      flag_mu_factor = mul/mur

      nscbc_p_inf = this%P0
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas

   end subroutine TGVCompSolver_initializeParameters

   subroutine TGVCompSolver_evalInitialConditions(this)
      class(TGVCompSolver), intent(inout) :: this
      real(4) :: x,y,z
      integer(4) :: iNodeL

      if(mpi_rank.eq.0) write(*,*) "--| TGV - Setting Initial Conditions..."


      !$acc parallel loop
      do iNodeL=1,numNodesRankPar
         x = coordPar(iNodeL,1)
         y = coordPar(iNodeL,2)
         z = coordPar(iNodeL,3)

         u(iNodeL,1,2) =  this%V0*sin(x/(this%L))*cos(y/(this%L))*cos(z/(this%L))
         u(iNodeL,2,2) = -this%V0*cos(x/(this%L))*sin(y/(this%L))*cos(z/(this%L))
         u(iNodeL,3,2) = 0.0

         pr(iNodeL,2)  = this%P0+((1.0_rp*this%V0*this%V0)/(16.0_rp))*(cos(2.0_rp*x/this%L)+cos(2.0_rp*y/this%L))*(cos(2.0_rp*z/this%L)+2.0_rp)
         rho(iNodeL,2) = pr(iNodeL,2)/this%Rgas/this%T0
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) =  eta(iNodeL,2)
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
   end subroutine TGVCompSolver_evalInitialConditions

end module TGVCompSolver_mod
