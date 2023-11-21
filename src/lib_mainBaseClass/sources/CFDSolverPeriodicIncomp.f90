module CFDSolverPeriodicIncomp_mod
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
   use time_integ_incomp
   use mod_analysis
   use mod_numerical_params
   use mod_time_ops
   use mod_fluid_viscosity
   use mod_postpro
   use mod_aver
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use CFDSolverBase_mod
   implicit none
   private

   type, public, extends(CFDSolverBase) :: CFDSolverPeriodicIncomp

   contains
      procedure, public :: callTimeIntegration              => CFDSolverPeriodicIncomp_callTimeIntegration
      procedure, public :: initializeDefaultParameters      => CFDSolverPeriodicIncomp_initializeDefaultParameters
      procedure, public :: evalInitialDt                    => CFDSolverPeriodicIncomp_evalInitialDt
      procedure, public :: evalDt                           => CFDSolverPeriodicIncomp_evalDt
      procedure, public :: initNSSolver                     => CFDSolverPeriodicIncomp_initNSSolver
      procedure, public :: endNSSolver                      => CFDSolverPeriodicIncomp_endNSSolver
   end type CFDSolverPeriodicIncomp
contains
   subroutine CFDSolverPeriodicIncomp_initNSSolver(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      
      call init_rk4_solver_incomp(numNodesRankPar) 

   end subroutine CFDSolverPeriodicIncomp_initNSSolver

   subroutine CFDSolverPeriodicIncomp_endNSSolver(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      
       call end_rk4_solver_incomp()

   end subroutine CFDSolverPeriodicIncomp_endNSSolver

   subroutine CFDSolverPeriodicIncomp_evalInitialDt(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this

      !$acc kernels
      csound(:) = 0.0_rp
      !$acc end kernels

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating initial dt..."
      call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
      if(mpi_rank.eq.0) write(111,*) "--| Initial time-step dt := ",this%dt,"s"

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverPeriodicIncomp_evalInitialDt

   subroutine CFDSolverPeriodicIncomp_evalDt(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      
      !$acc kernels
      csound(:) = 0.0_rp
      !$acc end kernels

      call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)

   end subroutine CFDSolverPeriodicIncomp_evalDt

   subroutine CFDSolverPeriodicIncomp_initializeDefaultParameters(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this

      write(this%mesh_h5_file_path,*) "./"
      write(this%mesh_h5_file_name,*) "meshFile"

      write(this%results_h5_file_path,*) "./"
      write(this%results_h5_file_name,*) "resultsFile"

      write(this%io_prepend_path,*) "./"
      write(this%io_append_info,*) ""

      this%time = 0.0_rp
      this%initial_istep = 1
      this%maxPhysTime = 1.0e6_rp

      !--------------------------------------------------------------------------
      this%save_logFile_first = 1
      this%save_logFile_step = 100
      !--------------------------------------------------------------------------
      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%restartFile_to_load = 1
      this%restartFileCnt = 1
      !--------------------------------------------------------------------------
      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 100000
      !--------------------------------------------------------------------------
      this%initial_avgTime = 0.0_rp
      this%elapsed_avgTime = 0.0_rp
      !--------------------------------------------------------------------------

      this%loadRestartFile    = .false.
      this%saveAvgFile        = .false.
      this%loadAvgFile        = .false.
      this%continue_oldLogs   = .false.
      this%doGlobalAnalysis   = .false.
      this%doTimerAnalysis    = .false.
      this%isFreshStart       = .true.
      this%saveInitialField   = .false.
      this%saveSurfaceResults = .false.
      this%isWallModelOn      = .false.
      this%isSymmetryOn=.false.
      !@JORDI: discuss which other parameters can be set as default....

      this%numNodeScalarFields2save    = 0
      this%numNodeVectorFields2save    = 0
      this%numElemGpScalarFields2save  = 0
      this%save_scalarField_rho        = .false.
      this%save_scalarField_muFluid    = .true.
      this%save_scalarField_pressure   = .true.
      this%save_scalarField_energy     = .false.
      this%save_scalarField_entropy    = .true.
      this%save_scalarField_csound     = .false.
      this%save_scalarField_machno     = .false.
      this%save_scalarField_divU       = .false.
      this%save_scalarField_qcrit      = .true.
      this%save_scalarField_muSgs      = .true.
      this%save_scalarField_muEnvit    = .true.
      this%save_vectorField_vel        = .true.
      this%save_vectorField_gradRho    = .false.
      this%save_vectorField_curlU      = .true.

      this%numAvgNodeScalarFields2save    = 0
      this%numAvgNodeVectorFields2save    = 0
      this%numAvgElemGpScalarFields2save  = 0
      this%save_avgScalarField_rho     = .false.
      this%save_avgScalarField_pr      = .true.
      this%save_avgScalarField_mueff   = .true.
      this%save_avgVectorField_vel     = .true.
      this%save_avgVectorField_ve2     = .true.
      this%save_avgVectorField_vex     = .true.
      this%save_avgVectorField_vtw     = .true.

      flag_real_diff=1
      flag_diff_suth=0
      flag_walave = .true.

   end subroutine CFDSolverPeriodicIncomp_initializeDefaultParameters

   subroutine CFDSolverPeriodicIncomp_callTimeIntegration(this,istep)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      integer(4)                    , intent(in)    :: istep

      this%noBoundaries = .true.
      call ab_main_incomp(istep,this%save_logFile_next,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbnNodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor,mue_l)

   end subroutine CFDSolverPeriodicIncomp_callTimeIntegration

end module CFDSolverPeriodicIncomp_mod
