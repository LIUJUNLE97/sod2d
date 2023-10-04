module ABlFlowSolverIncomp_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: ABlFlowSolverIncomp

      real(rp) , public  :: rough,vinf,Lhub,rho0,Lz,ustar

   contains
      procedure, public :: fillBCTypes           => ABlFlowSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => ABlFlowSolverIncomp_initializeParameters
      procedure, public :: initializeSourceTerms => ABlFlowSolverIncomp_initializeSourceTerms
      procedure, public :: evalInitialConditions => ABlFlowSolverIncomp_evalInitialConditions
   end type ABlFlowSolverIncomp
contains

   subroutine ABlFlowSolverIncomp_fill_BC_Types(this)
      class(ABlFlowSolverIncomp), intent(inout) :: this

      bouCodes2BCType(1) = bc_type_slip_wall_model
      bouCodes2BCType(2) = bc_type_slip_adiabatic
      !$acc update device(bouCodes2BCType(:))

   end subroutine ABlFlowSolverIncomp_fill_BC_Types

   subroutine ABlFlowSolverIncomp_initializeSourceTerms(this)
      class(ABlFlowSolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL

      allocate(source_term(numNodesRankPar,ndime))
      
      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) =  this%rho0*this%ustar**2/this%Lz
         source_term(iNodeL,2) = 0.00_rp
         source_term(iNodeL,3) = 0.00_rp
      end do
      !$acc end parallel loop

   end subroutine ABlFlowSolverIncomp_initializeSourceTerms

   subroutine ABLFlowSolverIncomp_initializeParameters(this)
      class(ABlFlowSolverIncomp), intent(inout) :: this
      real(rp) :: mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "abl"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 10000001

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 5000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 5000
      this%loadRestartFile = .false.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%saveAvgFile = .true.
      this%loadAvgFile = .false.

      this%saveSurfaceResults = .false.
      !----------------------------------------------

      ! numerical params
      flag_les = 1
      flag_type_wmles = wmles_type_abl

      this%cfl_conv = 0.9_rp 
      this%cfl_diff = 0.9_rp
      
      this%rho0   = 1.0_rp !1.2
      this%Lz = 1500.0_rp
      this%Lhub = 90.0_rp
      this%vinf = 8.0_rp
      this%rough = 0.1682_rp
      this%ustar = this%vinf*0.41_rp/log(1.0_rp+this%Lhub/this%rough)
      incomp_viscosity = 1.81e-5

      flag_mu_factor = 1.0_rp

      nscbc_p_inf = 0.0_rp

      maxIter = 200
      tol = 1e-3

      flag_fs_fix_pressure = .false.
      
      period_walave   = 200.0_rp
      this%initial_avgTime = 3600.0_rp

   end subroutine ABlFlowSolverIncomp_initializeParameters

   subroutine ABlFlowSolverIncomp_evalInitialConditions(this)
      class(ABlFlowSolverIncomp), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, rti(3), zp,velo_aux1
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      logical :: readFiles
      character(512) :: initialField_filePath

      readFiles = .false.

      if(readFiles) then
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         write(initialField_filePath,*) ""
         call read_veloc_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = 0.0_rp 
            rho(iNodeL,2) = this%rho0            

            rho(iNodeL,3) = rho(iNodeL,2)
            pr(iNodeL,3) =  pr(iNodeL,2)
         end do
         !$acc end parallel loop
      else
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         auxCnt = 1
         !!!!$acc parallel loop 
         serialLoop : do iLine = 1,totalNumNodesSrl
            call random_number(rti)
            if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
               iNodeL = matGidSrlOrdered(auxCnt,1)
               auxCnt=auxCnt+1

               zp = coordPar(iNodeL,3)
               velo =  this%ustar*log(1.0_rp+zp/this%rough)/0.41_rp

               u(iNodeL,1,2) = velo*(1.0_rp + 0.1_rp*(rti(1) -0.5_rp))
               u(iNodeL,2,2) = velo*(0.1_rp*(rti(2) -0.5_rp))
               u(iNodeL,3,2) = velo*(0.1_rp*(rti(3) -0.5_rp))
            end if
            if(auxCnt.gt.numNodesRankPar) then
               exit serialLoop
            end if
         end do serialLoop
         !!!!$acc end parallel loop

         !$acc update device(u(:,:,:))

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = 0.0_rp 
            rho(iNodeL,2) = this%rho0            
         end do
         !$acc end parallel loop
      end if

      !$acc kernels
      mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      zo(:) = this%rough
      !$acc end kernels
      call nvtxEndRange

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

      !
      ! Initialize exponential averaging for wall law 
      !
      call nvtxStartRange("Wall Average init")
      if(flag_walave .eqv. .true.) then
         !$acc kernels
         walave_u(:,:) = u(:,:,2)
         !$acc end kernels
      end if
      call nvtxEndRange

   end subroutine ABlFlowSolverIncomp_evalInitialConditions

end module ABlFlowSolverIncomp_mod
