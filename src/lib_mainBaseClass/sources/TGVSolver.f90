module TGVSolver_mod
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
      write(this%mesh_h5_file_name,*) "cube"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      write(this%io_append_info,*) ""

      this%doGlobalAnalysis = .true.
      this%doTimerAnalysis = .true.
      this%saveInitialField = .false.

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 500001
      this%maxPhysTime = 20.0_rp

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 1000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 1000
      this%loadRestartFile = .false.
      this%restartFile_to_load = 2 !1 or 2
      this%continue_oldLogs = .false.

      this%saveAvgFile = .false.
      this%loadAvgFile = .false.
      !----------------------------------------------

      ! numerical params
      flag_les = 0
      flag_implicit = 0

      maxIter = 200
      tol = 1e-3

      this%cfl_conv = 0.95_rp !0.5_rp
      this%cfl_diff = 100.0_rp !0.5_rp

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

      !Witness points parameters
      this%have_witness          = .false.
      this%witness_inp_file_name = "witness.txt"
      this%witness_h5_file_name  = "resultwit.h5"
      this%leapwit               = 1
      this%leapwitsave           = 20
      this%nwit                  = 12
      this%wit_save_u_i          = .true.
      this%wit_save_pr           = .true.
      this%wit_save_rho          = .true.
      this%continue_witness      = .false.

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

      !$acc parallel loop
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
   end subroutine TGVSolver_evalInitialConditions

end module TGVSolver_mod
