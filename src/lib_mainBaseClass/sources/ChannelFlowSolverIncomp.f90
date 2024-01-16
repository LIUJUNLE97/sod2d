module ChannelFlowSolverIncomp_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: ChannelFlowSolverIncomp

      real(rp) , public  :: vo, delta, rho0, Retau, Re, utau, mu

   contains
      procedure, public :: fillBCTypes           => ChannelFlowSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => ChannelFlowSolverIncomp_initializeParameters
      procedure, public :: initializeSourceTerms => ChannelFlowSolverIncomp_initializeSourceTerms
      procedure, public :: evalInitialConditions => ChannelFlowSolverIncomp_evalInitialConditions
   end type ChannelFlowSolverIncomp
contains

   subroutine ChannelFlowSolverIncomp_fill_BC_Types(this)
      class(ChannelFlowSolverIncomp), intent(inout) :: this

      !bouCodes2BCType(1) = bc_type_slip_wall_model
      bouCodes2BCType(1) = bc_type_non_slip_adiabatic
      !$acc update device(bouCodes2BCType(:))

   end subroutine ChannelFlowSolverIncomp_fill_BC_Types

   subroutine ChannelFlowSolverIncomp_initializeSourceTerms(this)
      class(ChannelFlowSolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL

      allocate(source_term(numNodesRankPar,ndime))
      !$acc enter data create(source_term(:,:))

      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) = (this%utau*this%utau*this%rho0/this%delta)
         source_term(iNodeL,2) = 0.00_rp
         source_term(iNodeL,3) = 0.00_rp

         !just a momentary trick
         pr(iNodeL,2) = 0.0_rp 
         rho(iNodeL,2) = this%rho0            

         rho(iNodeL,3) = rho(iNodeL,2)
         pr(iNodeL,3) =  pr(iNodeL,2)
      end do
      !$acc end parallel loop

      !$acc kernels
      mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      !$acc end kernels
   end subroutine ChannelFlowSolverIncomp_initializeSourceTerms

   subroutine ChannelFlowSolverIncomp_initializeParameters(this)
      class(ChannelFlowSolverIncomp), intent(inout) :: this
      real(rp) :: mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "channel"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 10000001

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 10000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%loadRestartFile = .true.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%saveAvgFile = .true.
      this%loadAvgFile = .false.

      this%saveSurfaceResults = .false.
      !----------------------------------------------

      ! numerical params
      flag_les = 1

      !this%cfl_conv = 0.9_rp 
      !this%cfl_diff = 0.9_rp
      flag_use_constant_dt = 1
      this%dt = 5.0e-4
      flag_cg_prec_bdc = .false.

      
      this%vo = 1.0_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%Retau  = 950.0_rp

      this%Re     = exp((1.0_rp/0.88_rp)*log(this%Retau/0.09_rp))
      this%mu    = (this%rho0*2.0_rp*this%delta*this%vo)/this%Re
      this%utau   = (this%Retau*this%mu)/(this%delta*this%rho0)

      incomp_viscosity = this%mu
      flag_mu_factor = 1.0_rp

      nscbc_p_inf = 0.0_rp

      maxIter = 20
      tol = 1e-3

      flag_fs_fix_pressure = .false.
      
      period_walave   = 200.0_rp

   end subroutine ChannelFlowSolverIncomp_initializeParameters

   subroutine ChannelFlowSolverIncomp_evalInitialConditions(this)
      class(ChannelFlowSolverIncomp), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, rti(3), yp,velo_aux1
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
               if(coordPar(iNodeL,2)<this%delta) then
                  yp = coordPar(iNodeL,2)*this%utau*this%rho0/this%mu
               else
                  yp = abs(coordPar(iNodeL,2)-2.0_rp*this%delta)*this%utau*this%rho0/this%mu
               end if

               velo = this%utau*((1.0_rp/0.41_rp)*log(1.0_rp+0.41_rp*yp)+7.8_rp*(1.0_rp-exp(-yp/11.0_rp)-(yp/11.0_rp)*exp(-yp/3.0_rp))) 

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

   end subroutine ChannelFlowSolverIncomp_evalInitialConditions

end module ChannelFlowSolverIncomp_mod
