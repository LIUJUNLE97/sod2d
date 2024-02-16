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
      use json_module
      implicit none
      class(TGVCompSolver), intent(inout) :: this
      real(rp) :: mul, mur
      logical :: found
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      call json%get("save_logFile_first",this%save_logFile_first, found)
      call json%get("save_logFile_step",this%save_logFile_step, found)

      call json%get("save_resultsFile_first",this%save_resultsFile_first, found)
      call json%get("save_resultsFile_step" ,this%save_resultsFile_step, found)

      call json%get("save_restartFile_first",this%save_restartFile_first, found)
      call json%get("save_restartFile_step" ,this%save_restartFile_step, found)

      call json%get("loadRestartFile" ,this%loadRestartFile, found)
      call json%get("restartFile_to_load" ,this%restartFile_to_load, found)
      call json%get("continue_oldLogs" ,this%continue_oldLogs, found)
      call json%get("saveAvgFile" ,this%saveAvgFile, found)
      call json%get("loadAvgFile" ,this%loadAvgFile, found)

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found)

      call json%get("saveInitialField",this%saveInitialField, found)

      call json%get("mesh_h5_file_path",value, found)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found)
      write(this%mesh_h5_file_name,*) value
      
      call json%get("results_h5_file_path",value, found)
      write(this%results_h5_file_path,*) value
      call json%get("results_h5_file_name",value, found)
      write(this%results_h5_file_name,*) value

      call json%get("doGlobalAnalysis",this%doGlobalAnalysis, found)
      call json%get("doTimerAnalysis",this%doTimerAnalysis, found)

      call json%get("final_istep",this%final_istep, found)
      call json%get("maxPhysTime",this%maxPhysTime, found)

      call json%get("cfl_conv",this%cfl_conv, found)
      call json%get("cfl_diff",this%cfl_diff, found)

      call json%get("flag_implicit",flag_implicit, found)

      call json%get("maxIter",maxIter, found)
      call json%get("tol",tol, found)

      call json%get("M",this%Ma, found)
      call json%get("L",this%L, found)
      call json%get("v0",this%V0, found)
      call json%get("Prt",this%Prt, found)
      call json%get("T0",this%T0, found)
      call json%get("Re",this%Re, found)
      call json%get("gamma_gas",this%gamma_gas, found)

      ! fixed by the type of base class parameters
      this%Rgas = this%V0*this%V0/(this%gamma_gas*this%T0*this%Ma*this%Ma)
      this%Cp = this%gamma_gas*this%Rgas/(this%gamma_gas-1.0_rp)
      this%rho0 = (1.0_rp/(this%gamma_gas*this%Ma*this%Ma))/(this%Rgas*this%T0)

      mul    = (this%rho0*this%V0*this%L)/this%Re
      this%P0 = 1.0_rp/(this%Ma*this%Ma*this%gamma_gas)
      mur = 0.000001458_rp*(this%T0**1.50_rp)/(this%T0+110.40_rp)
      flag_mu_factor = mul/mur

      nscbc_p_inf = this%P0
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas

      call json%destroy()

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
