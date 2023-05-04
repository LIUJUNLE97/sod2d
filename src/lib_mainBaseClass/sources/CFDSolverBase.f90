module mod_arrays
      use mod_numerical_params

      implicit none

      ! main allocatable arrays
      ! integer ---------------------------------------------------
      integer(4), allocatable :: lelpn(:),point2elem(:),bouCodes2BCType(:)
      integer(4), allocatable :: atoIJ(:),atoIJK(:),listHEX08(:,:),lnbn(:,:),invAtoIJK(:,:,:),gmshAtoI(:),gmshAtoJ(:),gmshAtoK(:),lnbnNodes(:)

      ! real ------------------------------------------------------
      real(rp), allocatable :: normalsAtNodes(:,:)
      real(rp), allocatable :: helem(:),helem_l(:,:)
      real(rp), allocatable :: xgp(:,:), wgp(:), xgp_b(:,:), wgp_b(:)
      real(rp), allocatable :: Ngp(:,:), dNgp(:,:,:), Ngp_b(:,:), dNgp_b(:,:,:)
      real(rp), allocatable :: Ngp_l(:,:), dNgp_l(:,:,:),dlxigp_ip(:,:,:)
      real(rp), allocatable :: Je(:,:), He(:,:,:,:), bou_norm(:,:)
      real(rp), allocatable :: gpvol(:,:,:), gradRho(:,:), curlU(:,:), divU(:), Qcrit(:)
      real(rp), allocatable :: u(:,:,:),q(:,:,:),rho(:,:),pr(:,:),E(:,:),Tem(:,:),e_int(:,:),csound(:),eta(:,:),machno(:),tauw(:,:)
      real(rp), allocatable :: Ml(:)
      real(rp), allocatable :: mu_e(:,:),mu_fluid(:),mu_sgs(:,:),mu_factor(:)
      real(rp), allocatable :: source_term(:,:)
      real(rp), allocatable :: acurho(:), acupre(:), acuvel(:,:), acuve2(:,:), acumueff(:),acuvex(:,:),acutw(:,:)
      real(rp), allocatable :: avrho(:), avpre(:), avvel(:,:), avve2(:,:), avmueff(:),avvex(:,:),avtw(:,:)
      real(rp), allocatable :: kres(:),etot(:),au(:,:),ax1(:),ax2(:),ax3(:)
      real(rp), allocatable :: Fpr(:,:), Ftau(:,:)
      real(rp), allocatable :: u_buffer(:,:)

end module mod_arrays

module CFDSolverBase_mod
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
      use mod_comms
      use mod_comms_boundaries
      use mod_smartredis
   implicit none
   private

   type, public :: CFDSolverBase

      ! main integer parameters
      !NOTES @JORDI: -> it would be nice if nelem, npoin, nper... etc dissapear from here!
      integer(4), public :: nstep, nper, numCodes,currentNonLinealIter
      integer(4), public :: nsave, nleap
      integer(4), public :: nsave2, nleap2
      integer(4), public :: nsaveAVG, nleapAVG
      integer(4), public :: counter, npoin_w
      integer(4), public :: load_step, initial_istep

      ! main logical parameters
      logical, public    :: doGlobalAnalysis=.false.,isFreshStart=.true.,doTimerAnalysis=.false.
      logical, public    :: loadResults=.false.,continue_oldLogs=.false.,saveInitialField=.true.,isWallModelOn=.false.
      logical, public    :: useIntInComms=.false.,useRealInComms=.false.

      ! main char variables
      character(512) :: log_file_name
      character(512) :: mesh_h5_file_path,mesh_h5_file_name
      character(512) :: results_h5_file_path,results_h5_file_name

      ! main real parameters
      real(rp) , public                   :: cfl_conv, cfl_diff, acutim
      real(rp) , public                   :: leviCivi(3,3,3), surfArea, EK, VolTot, eps_D, eps_S, eps_T, maxmachno
      real(rp) , public                   :: dt, Cp, Rgas, gamma_gas,Prt,tleap,time, maxPhysTime
      logical  , public                   :: noBoundaries

   contains
      procedure, public :: printDt => CFDSolverBase_printDt
      procedure, public :: printAll => CFDSolverBase_printAll
      procedure, public :: run => CFDSolverBase_run
      procedure, public :: initializeDefaultParameters => CFDSolverBase_initializeDefaultParameters
      procedure, public :: initializeParameters => CFDSolverBase_initializeParameters
      procedure, public :: initializeSourceTerms => CFDSolverBase_initializeSourceTerms
      procedure, public :: openMesh => CFDSolverBase_openMesh
      procedure, public :: evalCharLength => CFDSolverBase_evalCharLength
      procedure, public :: boundaryFacesToNodes => CFDSolverBase_boundaryFacesToNodes
      procedure, public :: normalFacesToNodes => CFDSolverBase_normalFacesToNodes
      procedure, public :: fillBCTypes => CFDSolverBase_fill_BC_Types
      procedure, public :: allocateVariables => CFDSolverBase_allocateVariables
      procedure, public :: evalOrLoadInitialConditions => CFDSolverBase_evalOrLoadInitialConditions
      procedure, public :: evalInitialConditions => CFDSolverBase_evalInitialConditions
      procedure, public :: evalInitialViscosity =>CFDSolverBase_evalInitialViscosity
      procedure, public :: evalViscosityFactor=>CFDSolverBase_evalViscosityFactor
      procedure, public :: evalInitialDt =>CFDSolverBase_evalInitialDt
      procedure, public :: evalShapeFunctions =>CFDSolverBase_evalShapeFunctions
      procedure, public :: evalBoundaryNormals =>CFDSolverBase_evalBoundaryNormals
      procedure, public :: evalJacobians =>CFDSolverBase_evalJacobians
      procedure, public :: evalAtoIJKInverse =>CFDSolverBase_evalAtoIJKInverse
      procedure, public :: eval_elemPerNode_and_nearBoundaryNode =>CFDSolverBase_eval_elemPerNode_and_nearBoundaryNode
      procedure, public :: evalMass=>CFDSolverBase_evalMass
      procedure, public :: evalFirstOutput =>CFDSolverBase_evalFirstOutput
      procedure, public :: evalTimeIteration =>CFDSolverBase_evalTimeIteration
      procedure, public :: callTimeIntegration =>CFDSolverBase_callTimeIntegration
      procedure, public :: saveAverages =>CFDSolverBase_saveAverages
      procedure, public :: savePosprocessingFields =>CFDSolverBase_savePosprocessingFields
      procedure, public :: afterDt =>CFDSolverBase_afterDt

      procedure, public :: initialBuffer =>CFDSolverBase_initialBuffer

      procedure :: open_log_file
      procedure :: close_log_file
      procedure :: open_analysis_files
      procedure :: close_analysis_files
      procedure :: flush_log_file
      procedure :: eval_vars_after_load_hdf5_resultsFile
      procedure :: eval_initial_mu_sgs
      procedure :: checkIfWallModelOn
   end type CFDSolverBase
contains

   subroutine CFDSolverBase_printDt(this)
      class(CFDSolverBase), intent(inout) :: this

      write(*,*) " Dt ",this%dt
   end subroutine CFDSolverBase_printDt

   subroutine CFDSolverBase_printAll(this)
      class(CFDSolverBase), intent(inout) :: this

      call this%printDt()
   end subroutine CFDSolverBase_printAll

   subroutine CFDSolverBase_initializeSourceTerms(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iNodeL

      allocate(source_term(numNodesRankPar,ndime))
      !$acc kernels
      source_term(:,:) = 0.00_rp
      !$acc end kernels

   end subroutine CFDSolverBase_initializeSourceTerms

   subroutine CFDSolverBase_initializeDefaultParameters(this)
      class(CFDSolverBase), intent(inout) :: this

      write(this%mesh_h5_file_path,*) "./"
      write(this%mesh_h5_file_name,*) "meshFile"

      write(this%results_h5_file_path,*) "./"
      write(this%results_h5_file_name,*) "resultsFile"

      this%time = 0.0_rp
      this%load_step = 0
      this%initial_istep = 1
      this%maxPhysTime = 1.0e6_rp

      this%loadResults = .false.
      this%continue_oldLogs= .false.
      this%doGlobalAnalysis = .false.
      this%doTimerAnalysis = .false.
      this%isFreshStart=.true.
      this%saveInitialField=.true.
      this%isWallModelOn=.false.
      !@JORDI: discuss which other parameters can be set as default....

   end subroutine CFDSolverBase_initializeDefaultParameters

   subroutine CFDSolverBase_initializeParameters(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_initializeParameters

   subroutine CFDSolverBase_openMesh(this)
      class(CFDSolverBase), intent(inout) :: this

      call nvtxStartRange("Open mesh")
      ! init hdf5
      !call init_hdf5_interface(this%mesh_h5_file_path,this%mesh_h5_file_name,this%results_h5_file_path,this%results_h5_file_name)
      call set_hdf5_meshFile_name(this%mesh_h5_file_path,this%mesh_h5_file_name,mpi_size)
      call set_hdf5_baseResultsFile_name(this%results_h5_file_path,this%results_h5_file_name,this%mesh_h5_file_name,mpi_size)

      this%useIntInComms=.true.
      this%useRealInComms=.true.

      call load_hdf5_meshfile()
      ! init comms
      call init_comms(this%useIntInComms,this%useRealInComms)
      !----- init comms boundaries
      call init_comms_bnd(this%useIntInComms,this%useRealInComms)

      call nvtxEndRange

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_openMesh

   subroutine CFDSolverBase_fill_BC_Types(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_fill_BC_Types

   subroutine checkIfWallModelOn(this)
      class(CFDSolverBase), intent(inout) :: this
      integer :: iBound,bcCode,auxBoundCnt

      this%isWallModelOn = .false.
      do iBound = 1,numBoundCodes
         if(bouCodes2BCType(iBound) .eq. bc_type_slip_wall_model) then
            this%isWallModelOn = .true.
            if(mpi_rank.eq.0) write(111,*) "--| Wall-Model activated in Boundary id",iBound
         end if
      end do

      numBoundsWMRankPar = 0
      do iBound = 1,numBoundsRankPar
         bcCode = bouCodes2BCType(bouCodesPar(iBound))
         if(bcCode .eq. bc_type_slip_wall_model) then
            numBoundsWMRankPar = numBoundsWMRankPar + 1
         end if
      end do
      !write(*,*) '[',mpi_rank,'] numBoundsWMRankPar',numBoundsWMRankPar

      allocate(listBoundsWallModel(numBoundsWMRankPar))
      auxBoundCnt = 0
      do iBound = 1,numBoundsRankPar
         bcCode = bouCodes2BCType(bouCodesPar(iBound))
         if(bcCode .eq. bc_type_slip_wall_model) then
            auxBoundCnt = auxBoundCnt + 1
            listBoundsWallModel(auxBoundCnt) = iBound
         end if
      end do

   end subroutine checkIfWallModelOn

   subroutine CFDSolverBase_normalFacesToNodes(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), allocatable    :: aux1(:)
      integer(4) :: iNodeL,iBound,ipbou,iElem,jgaus,kgaus,idime,iAux
      real(rp) :: aux(3), normaux,sig

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating Normals at Nodes for Wall-Model"

      allocate(normalsAtNodes(numNodesRankPar,ndime))

      !$acc kernels
      normalsAtNodes(:,:) = 0.0_rp
      !$acc end kernels

      !$acc parallel loop gang
      do iAux = 1,numBoundsWMRankPar
         iBound = listBoundsWallModel(iAux)
         iElem = point2elem(boundPar(iBound,npbou)) ! I use an internal face node to be sure is the correct element
         jgaus = connecParWork(iElem,nnode)         ! internal node
         !$acc loop vector private(aux)
         do ipbou = 1,npbou
            kgaus = boundPar(iBound,ipbou) ! node at the boundary
            sig=1.0_rp
            aux(1) = boundNormalPar(iBound,(ipbou-1)*ndime+1)
            aux(2) = boundNormalPar(iBound,(ipbou-1)*ndime+2)
            aux(3) = boundNormalPar(iBound,(ipbou-1)*ndime+3)
            normaux = sqrt(dot_product(aux,aux))
            if(dot_product(coordPar(jgaus,:)-coordPar(kgaus,:), aux(:)) .lt. 0.0_rp ) then
               sig=-1.0_rp
            end if
            !$acc loop seq
            do idime = 1,ndime
               aux(idime) = aux(idime)*sig/normaux
            end do
            normalsAtNodes(kgaus,1) = normalsAtNodes(kgaus,1) + aux(1)
            normalsAtNodes(kgaus,2) = normalsAtNodes(kgaus,2) + aux(2)
            normalsAtNodes(kgaus,3) = normalsAtNodes(kgaus,3) + aux(3)
         end do
      end do
      !$acc end parallel loop

      if(mpi_size.ge.2) then
         call mpi_halo_boundary_atomic_update_real(normalsAtNodes(:,1))
         call mpi_halo_boundary_atomic_update_real(normalsAtNodes(:,2))
         call mpi_halo_boundary_atomic_update_real(normalsAtNodes(:,3))
      end if

      !$acc parallel loop  private(aux)
      do iNodeL = 1,numNodesRankPar
         aux(1) = normalsAtNodes(iNodeL,1)
         aux(2) = normalsAtNodes(iNodeL,2)
         aux(3) = normalsAtNodes(iNodeL,3)
         normaux = sqrt(dot_product(aux,aux))

         if(normaux .gt. 1e-10) then
            normalsAtNodes(iNodeL,1) = aux(1)/normaux
            normalsAtNodes(iNodeL,2) = aux(2)/normaux
            normalsAtNodes(iNodeL,3) = aux(3)/normaux
         end if
      end do
      !$acc end parallel loop

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_normalFacesToNodes

   subroutine CFDSolverBase_boundaryFacesToNodes(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), allocatable    :: aux1(:)
      integer(4) :: iNodeL,iBound,ipbou,ielem,jgaus,kgaus,idime
      real(rp) :: aux(3), normaux,sig

      allocate(bouCodesNodesPar(numNodesRankPar))
      allocate(aux1(numNodesRankPar))
      allocate(bouCodes2BCType(numBoundCodes))

      bouCodes2BCType(:) = 0

      call this%fillBCTypes()

      !$acc kernels
      aux1(:) = max_num_bou_codes
      bouCodesNodesPar(:) =  max_num_bou_codes
      !$acc end kernels

      !$acc parallel loop gang
      do iBound = 1,numBoundsRankPar
         !$acc loop vector
         do ipbou = 1,npbou
            aux1(boundPar(iBound,ipbou)) = min(aux1(boundPar(iBound,ipbou)),bouCodes2BCType(bouCodesPar(iBound)))
         end do
      end do
      !$acc end parallel loop

      if((isMeshBoundaries).and.(mpi_size.ge.2)) then
         call mpi_halo_min_boundary_update_int_iSendiRcv(aux1)
      end if

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(aux1(iNodeL) .lt. max_num_bou_codes) then
            bouCodesNodesPar(iNodeL) = aux1(iNodeL)
         end if
      end do
      !$acc end parallel loop

      call this%checkIfWallModelOn()

      deallocate(aux1)

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_boundaryFacesToNodes

   subroutine CFDSolverBase_evalCharLength(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp) :: he_aux
      integer(4) :: iElem

      call nvtxStartRange("Elem size compute")
      allocate(helem(numElemsRankPar))
      do iElem = 1,numElemsRankPar
         call char_length(iElem,numElemsRankPar,numNodesRankPar,connecParOrig,coordPar,he_aux)
         helem(iElem) = he_aux
      end do
      call nvtxEndRange

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalCharLength

   subroutine CFDSolverBase_allocateVariables(this)
      class(CFDSolverBase), intent(inout) :: this

      if(mpi_rank.eq.0) write(111,*) "--| ALLOCATING MAIN VARIABLES"
      call nvtxStartRange("Allocate main vars")
      !
      ! Last rank is for prediction-advance related to entropy viscosity,
      ! where 1 is prediction, 2 is final value
      !
      allocate(u(numNodesRankPar,ndime,2))  ! Velocity
      allocate(q(numNodesRankPar,ndime,3))  ! momentum
      allocate(rho(numNodesRankPar,3))      ! Density
      allocate(pr(numNodesRankPar,2))       ! Pressure
      allocate(E(numNodesRankPar,3))        ! Total Energy
      allocate(Tem(numNodesRankPar,2))      ! Temperature
      allocate(e_int(numNodesRankPar,2))    ! Internal Energy
      allocate(eta(numNodesRankPar,3))      ! entropy
      allocate(csound(numNodesRankPar))     ! Speed of sound
      allocate(machno(numNodesRankPar))     ! Speed of sound
      allocate(mu_fluid(numNodesRankPar))   ! Fluid viscosity
      allocate(mu_factor(numNodesRankPar))  ! Fluid viscosity
      allocate(mu_e(numElemsRankPar,ngaus))  ! Elemental viscosity
      allocate(mu_sgs(numElemsRankPar,ngaus))! SGS viscosity
      allocate(u_buffer(numNodesRankPar,ndime))  ! momentum at the buffer

      allocate(tauw(numNodesRankPar,ndime))  ! momentum at the buffer

      !$acc kernels
      u(:,:,:) = 0.0_rp
      q(:,:,:) = 0.0_rp
      rho(:,:) = 0.0_rp
      pr(:,:) = 0.0_rp
      E(:,:) = 0.0_rp
      Tem(:,:) = 0.0_rp
      e_int(:,:) = 0.0_rp
      eta(:,:) = 0.0_rp
      csound(:) = 0.0_rp
      machno(:) = 0.0_rp
      mu_fluid(:) = 0.0_rp
      mu_factor(:) = 1.0_rp
      mu_e(:,:) = 0.0_rp
      mu_sgs(:,:) = 0.0_rp

      u_buffer(:,:) = 0.0_rp
      tauw(:,:) = 0.0_rp
      !$acc end kernels

      !ilsa
      allocate(kres(numNodesRankPar))
      allocate(etot(numNodesRankPar))
      allocate(au(numNodesRankPar,ndime))
      allocate(ax1(numNodesRankPar))
      allocate(ax2(numNodesRankPar))
      allocate(ax3(numNodesRankPar))
      !$acc kernels
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      au(:,:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      !$acc end kernels

      !boundary
      if(numBoundCodes .ge. 1) then
         allocate(Fpr(ndime,numBoundCodes))
         allocate(Ftau(ndime,numBoundCodes))
         !$acc kernels
         Fpr(:,:) = 0.0_rp
         Ftau(:,:) = 0.0_rp
         !$acc end kernels
      end if

      !*********************************************************************!
      ! Derivative-related fields                                           !
      !*********************************************************************!
      allocate(gradRho(numNodesRankPar,ndime))
      allocate(curlU(numNodesRankPar,ndime))
      allocate(divU(numNodesRankPar))
      allocate(Qcrit(numNodesRankPar))

      !$acc kernels
      gradRho(:,:) = 0.0_rp
      curlU(:,:) = 0.0_rp
      divU(:) = 0.0_rp
      Qcrit(:) = 0.0_rp
      !$acc end kernels

      !*********************************************************************!
      ! Allocate accumulators for averaging process                         !
      !*********************************************************************!
      allocate(acurho(numNodesRankPar))
      allocate(acupre(numNodesRankPar))
      allocate(acumueff(numNodesRankPar))
      allocate(acuvel(numNodesRankPar,ndime))
      allocate(acuve2(numNodesRankPar,ndime))
      allocate(acuvex(numNodesRankPar,ndime))
      allocate(acutw(numNodesRankPar,ndime))

      allocate(avrho(numNodesRankPar))
      allocate(avpre(numNodesRankPar))
      allocate(avmueff(numNodesRankPar))
      allocate(avvel(numNodesRankPar,ndime))
      allocate(avve2(numNodesRankPar,ndime))
      allocate(avvex(numNodesRankPar,ndime))
      allocate(avtw(numNodesRankPar,ndime))

      !$acc kernels
      acurho(:) = 0.0_rp
      acupre(:) = 0.0_rp
      acumueff(:) = 0.0_rp
      acuvel(:,:) = 0.0_rp
      acuve2(:,:) = 0.0_rp
      acuvex(:,:) = 0.0_rp
      acutw(:,:) = 0.0_rp

      avrho(:) = 0.0_rp
      avpre(:) = 0.0_rp
      avmueff(:) = 0.0_rp
      avvel(:,:) = 0.0_rp
      avve2(:,:) = 0.0_rp
      avvex(:,:) = 0.0_rp
      avtw(:,:) = 0.0_rp
      !$acc end kernels
      this%acutim = 0.0_rp

      call nvtxEndRange

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_allocateVariables

   subroutine CFDSolverBase_evalOrLoadInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

      if(this%loadResults) then
         if(mpi_rank.eq.0) write(111,*) "--| Loading results load_step ",this%load_step
         call load_hdf5_resultsFile(this%load_step,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,mu_sgs)
         if(mpi_rank.eq.0) write(111,*) "   --| Loaded results for load_step",this%load_step,"time",this%time
         call this%eval_vars_after_load_hdf5_resultsFile()

         if(this%continue_oldLogs) then
            this%initial_istep = this%load_step
            this%nsave = this%load_step+this%nleap
            this%nsaveAVG = this%load_step+this%nleapAVG
            this%nsave2 = this%load_step+this%nleap2
            this%isFreshStart = .false.
         else
            this%time = 0.0_rp
         end if
      else
         if(mpi_rank.eq.0) write(111,*) "--| Evaluating Initial Conditions..."
         call this%evalInitialConditions()
      end if

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalOrLoadInitialConditions

   subroutine CFDSolverBase_evalInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_evalInitialConditions

   subroutine CFDSolverBase_evalInitialViscosity(this)
      class(CFDSolverBase), intent(inout) :: this

        if (flag_real_diff == 1) then
           if (flag_diff_suth == 0) then
              call constant_viscosity(numNodesRankPar,0.000055_rp,mu_fluid)
           else
              call sutherland_viscosity(numNodesRankPar,Tem(:,2),mu_factor,mu_fluid)
           end if
        else if (flag_real_diff == 0) then
           !$acc kernels
           mu_fluid(:) = 0.0_rp
           !$acc end kernels
        else
           if(mpi_rank.eq.0) write(111,*) "--| DIFFUSION FLAG MUST BE EITHER 0 OR 1, NOT: ",flag_real_diff
           stop 1
        end if

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalInitialViscosity

   subroutine CFDSolverBase_evalViscosityFactor(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iNodeL

      ! set out of the buffer zone
      ! remember that the mu_factor field has to we filled at least with the
      ! flag_mu_factor

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

   end subroutine CFDSolverBase_evalViscosityFactor

   subroutine eval_initial_mu_sgs(this)
      class(CFDSolverBase), intent(inout) :: this

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating initial mu_sgs..."

      call nvtxStartRange("MU_SGS")
      if(flag_les_ilsa == 1) then
         this%dt = 1.0_rp !To avoid 0.0 division inside sgs_ilsa_visc calc
         call sgs_ilsa_visc(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,Ngp,dNgp,He,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,this%dt,rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3)
      else
         call sgs_visc(numElemsRankPar,numNodesRankPar,connecParWork,Ngp,dNgp,He,gpvol,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),Ml,mu_sgs)
      end if
      call nvtxEndRange

   end subroutine eval_initial_mu_sgs

   subroutine CFDSolverBase_evalInitialDt(this)
      class(CFDSolverBase), intent(inout) :: this

      !*********************************************************************!
      ! Compute initial time-step size                                      !
      !*********************************************************************!
#if 0
      if((this%loadResults .eqv. .false.).and.(flag_les == 1)) then
         !if not loading results and the case is turbulent, an initial mu_sgs is calculated to properly eval first time-step dt
         call this%eval_initial_mu_sgs()
      end if
#endif
      if(mpi_rank.eq.0) write(111,*) "--| Evaluating initial dt..."
      if (flag_real_diff == 1) then
         call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,rho(:,2))
      else
         call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
      end if
      if(mpi_rank.eq.0) write(111,*) "--| Initial time-step dt := ",this%dt,"s"

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalInitialDt

   subroutine CFDSolverBase_evalShapeFunctions(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp)                    :: s, t, z
      integer(4) :: igaus

      !*********************************************************************!
      ! Generate GLL table                                                  !
      !*********************************************************************!
      if(mpi_rank.eq.0) write(111,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."
      call nvtxStartRange("Gaussian Quadrature")

      !*********************************************************
      !           Allocating required arrays!
      allocate(atoIJ(16))
      allocate(atoIJK(64))
      !allocate(vtk_atoIJK(64))
      allocate(listHEX08((porder**ndime),2**ndime))

      allocate(xgp(ngaus,ndime))
      allocate(wgp(ngaus))
      allocate(xgp_b(npbou,ndime-1))
      allocate(wgp_b(npbou))

      allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))
      allocate(Ngp_l(ngaus,nnode),dNgp_l(ndime,nnode,ngaus))
      allocate(Ngp_b(npbou,npbou),dNgp_b(ndime-1,npbou,npbou))
      allocate(dlxigp_ip(ngaus,ndime,porder+1))
      !*********************************************************

      call set_hex64_lists(atoIJK,listHEX08)
      call set_qua16_lists(atoIJ)

      if(mpi_rank.eq.0) write(111,*) "  --| Generating Gauss-Lobatto-Legendre table..."
      call GaussLobattoLegendre_hex(atoIJK,xgp,wgp)
      call GaussLobattoLegendre_qua(atoIJ,xgp_b,wgp_b)

      call nvtxEndRange

      !*********************************************************************!
      ! Generate N and dN for all GP                                        !
      !*********************************************************************!

      ! TODO: Allow for more element types

      if(mpi_rank.eq.0) write(111,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."
      call nvtxStartRange("N and dN")

      do igaus = 1,ngaus
         s = xgp(igaus,1)
         t = xgp(igaus,2)
         z = xgp(igaus,3)
         call hex64(s,t,z,atoIJK,Ngp(igaus,:),dNgp(:,:,igaus),Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
      end do
      !
      ! Compute N andd dN for boundary elements
      !
      do igaus = 1,npbou
         s = xgp_b(igaus,1)
         t = xgp_b(igaus,2)
         call qua16(s,t,atoIJ,Ngp_b(igaus,:),dNgp_b(:,:,igaus))
      end do

      call nvtxEndRange

      !
      ! Compute Levi-Civita tensor
      !
      this%leviCivi = 0.0_rp
      this%leviCivi(2,3,1) =  1.0_rp
      this%leviCivi(3,2,1) = -1.0_rp
      this%leviCivi(1,3,2) = -1.0_rp
      this%leviCivi(3,1,2) =  1.0_rp
      this%leviCivi(1,2,3) =  1.0_rp
      this%leviCivi(2,1,3) = -1.0_rp

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalShapeFunctions

   subroutine CFDSolverBase_evalBoundaryNormals(this)
      class(CFDSolverBase), intent(inout) :: this

      if (isMeshBoundaries) then
         if(mpi_rank.eq.0) write(111,*) "--| COMPUTING BOUNDARY ELEMENT NORMALS"
         allocate(boundNormalPar(numBoundsRankPar,ndime*npbou))
         call nvtxStartRange("Bou normals")
         call boundary_normals(numNodesRankPar,numBoundsRankPar,boundParOrig,this%leviCivi,coordPar,dNgp_b,boundNormalPar)
         call nvtxEndRange
      end if

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalBoundaryNormals

   subroutine CFDSolverBase_evalJacobians(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ielem, igaus,iElemG
      real(8) :: vol_rank, vol_tot_d

      !*********************************************************************!
      ! Generate Jacobian related information                               !
      !*********************************************************************!

      if(mpi_rank.eq.0) write(111,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

      call nvtxStartRange("Jacobian info")
      allocate(He(ndime,ndime,ngaus,numElemsRankPar))
      allocate(gpvol(1,ngaus,numElemsRankPar))
      call elem_jacobian(numElemsRankPar,numNodesRankPar,connecParOrig,coordPar,dNgp,wgp,gpvol,He)
      call  nvtxEndRange

      vol_rank  = 0.0
      vol_tot_d = 0.0
      do ielem = 1,numElemsRankPar
         do igaus = 1,ngaus
            vol_rank = vol_rank+gpvol(1,igaus,ielem)
         end do
      end do

      call MPI_Allreduce(vol_rank,vol_tot_d,1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

      this%VolTot = real(vol_tot_d,rp)

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
      if(mpi_rank.eq.0) write(111,*) '--| DOMAIN VOLUME := ',this%VolTot

   end subroutine CFDSolverBase_evalJacobians

   subroutine CFDSolverBase_evalAtoIJKInverse(this)
      class(CFDSolverBase), intent(inout) :: this

      allocate(invAtoIJK(porder+1,porder+1,porder+1))
      allocate(gmshAtoI(nnode))
      allocate(gmshAtoJ(nnode))
      allocate(gmshAtoK(nnode))
      call atioIJKInverse(atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)

   end subroutine CFDSolverBase_evalAtoIJKInverse

   subroutine CFDSolverBase_eval_elemPerNode_and_nearBoundaryNode(this)
      class(CFDSolverBase), intent(inout) :: this

      !*********************************************************************!
      ! Compute list of elements per node (connectivity index)              !
      !*********************************************************************!
      ! evaluate near boundaries for the inlets and outlets
      !not the best place Oriol!
      if(mpi_rank.eq.0) write(111,*) "--| Doing near boundary calculations..."
      allocate(lelpn(numNodesRankPar))
      allocate(point2elem(numNodesRankPar))
      if(mpi_rank.eq.0) write(111,*) '  --| Evaluating point2elem array...'
      call elemPerNode(numElemsRankPar,numNodesRankPar,connecParWork,lelpn,point2elem)

      if(mpi_rank.eq.0) write(111,*) '  --| Evaluating lnbn & lnbnNodes arrays...'
      allocate(lnbn(numBoundsRankPar,npbou))
      allocate(lnbnNodes(numNodesRankPar))
      call nearBoundaryNode(numElemsRankPar,numNodesRankPar,numBoundsRankPar,connecParWork,coordPar,boundPar,bouCodesNodesPar,point2elem,atoIJK,lnbn,lnbnNodes)

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_eval_elemPerNode_and_nearBoundaryNode

   subroutine CFDSolverBase_evalMass(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iElem

      !*********************************************************************!
      ! Compute mass matrix (Lumped and Consistent) and set solver type     !
      !*********************************************************************!

      if(mpi_rank.eq.0) write(111,*) '--| COMPUTING LUMPED MASS MATRIX...'
      call nvtxStartRange("Lumped mass compute")
      allocate(Ml(numNodesRankPar))
      call lumped_mass_spectral(numElemsRankPar,numNodesRankPar,connecParWork,gpvol,Ml)
      call nvtxEndRange

      !charecteristic length for spectral elements for the entropy
      !stablisation
      allocate(helem_l(numElemsRankPar,nnode))
      do iElem = 1,numElemsRankPar
         call char_length_spectral(iElem,numElemsRankPar,numNodesRankPar,connecParOrig,coordPar,Ml,helem_l)
      end do

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

   end subroutine CFDSolverBase_evalMass

   subroutine CFDSolverBase_evalFirstOutput(this)
      class(CFDSolverBase), intent(inout) :: this
      character(500) :: tmpname
      integer(4) :: iCode

      if(this%saveInitialField) then
         call this%savePosprocessingFields(0)
      end if

      !*********************************************************************!
      ! Compute surface forces and area                                                                !
      !*********************************************************************!
      if (isMeshBoundaries) then
         do iCode = 1,numBoundCodes
            call nvtxStartRange("Surface info")
            call surfInfo(0,0.0_rp,numElemsRankPar,numNodesRankPar,numBoundsRankPar,iCode,connecParWork,boundPar,point2elem,&
               bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
               mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(:,iCode),Ftau(:,iCode))
            call nvtxEndRange
         end do
      end if

      call compute_fieldDerivs(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)

      call volAvg_EK(numElemsRankPar,numNodesRankPar,connecParWork,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
      call visc_dissipationRate(numElemsRankPar,numNodesRankPar,connecParWork,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
      call maxMach(numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,machno,this%maxmachno)
      call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
      if(mpi_rank.eq.0) then
         write(111,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
         write(111,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
         20 format(6(F16.8,2X))
      end if

   end subroutine CFDSolverBase_evalFirstOutput

   subroutine CFDSolverBase_callTimeIntegration(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      if(mpi_rank.eq.0) write(111,*) " Time integration should be overwritted"
      stop 1

   end subroutine CFDSolverBase_callTimeIntegration

   subroutine CFDSolverBase_saveAverages(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in)   :: istep

      call eval_average_window(isMeshPeriodic,numNodesRankPar,numElemsRankPar,acuvel,acuve2,acuvex,acurho,acupre,acumueff,acutw,this%acutim,&
											avvel,avve2,avvex,avrho,avpre,avmueff,avtw,nPerRankPar,masSlaRankPar)

      call save_hdf5_avgResultsFile(istep,avvel,avve2,avvex,avrho,avpre,avmueff,avtw)

   end subroutine CFDSolverBase_saveAverages

   subroutine CFDSolverBase_afterDt(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in) :: istep

   end subroutine CFDSolverBase_afterDt

   subroutine CFDSolverBase_savePosprocessingFields(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in) :: istep

      call save_hdf5_resultsFile(istep,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),eta(:,2),csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs)

   end subroutine CFDSolverBase_savePosprocessingFields

   subroutine CFDSolverBase_initialBuffer(this)
      class(CFDSolverBase), intent(inout) :: this

      !$acc kernels
      u_buffer(:,1) = nscbc_u_inf
      !$acc end kernels

   end subroutine CFDSolverBase_initialBuffer

   subroutine CFDSolverBase_evalTimeIteration(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: icode,counter,istep,flag_emac,flag_predic,inonLineal
      character(4) :: timeStep
      real(8) :: iStepTimeRank,iStepTimeMax,iStepEndTime,iStepStartTime,iStepAvgTime
      real(rp) :: inv_iStep, aux_rho(numNodesRankPar),aux_E(numNodesRankPar),aux_eta(numNodesRankPar),aux_q(numNodesRankPar,ndime),aux_pseudo_cfl
      real(rp) :: aux_envit(numElemsRankPar,nnode)
      logical :: do__iteration

      counter = 1
      flag_emac = 0
      flag_predic=0

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

      call nvtxStartRange("Start RK4")
      if(mpi_rank.eq.0) then
         write(*,*) 'Strarting evalTimeItarion! All info will be written in the log file: ',this%log_file_name
         write(111,*) 'Doing evalTimeItarion. Ini step:',this%initial_istep,'| End step:',this%nstep
      end if
      ! periodic with boundaries
      do istep = this%initial_istep,this%nstep
         !if (istep==this%nsave.and.mpi_rank.eq.0) write(111,*) '   --| STEP: ', istep
         !
         ! Prediction
         !
         flag_predic = 1
         call nvtxStartRange("Init pred "//timeStep,istep)
         !$acc kernels
         rho(:,1) = rho(:,2)
         u(:,:,1) = u(:,:,2)
         q(:,:,1) = q(:,:,2)
         pr(:,1) = pr(:,2)
         E(:,1) = E(:,2)
         Tem(:,1) = Tem(:,2)
         e_int(:,1) = e_int(:,2)
         eta(:,1) = eta(:,2)
         !$acc end kernels
         call nvtxEndRange

         ! nvtx range for full RK
         ! write(timeStep,'(i4)') istep
         call nvtxStartRange("RK4 step "//timeStep,istep)
         !
         ! Advance with entropy viscosity
         !
         flag_predic = 0

         if(flag_implicit == 1) then
            !$acc kernels
            aux_rho(:) = rho(:,2)
            aux_E(:) = E(:,2)
            aux_q(:,:) = q(:,:,2)
            aux_eta(:) = eta(:,2)
            aux_envit(:,:) = mu_e(:,:)
            !$acc end kernels
            aux_pseudo_cfl = pseudo_cfl
         end if

         do__iteration = .true.
         inonLineal = 1
         do while(do__iteration .eqv. .true.)
            if(flag_implicit == 1) then
               !$acc kernels
               rho(:,2) = aux_rho(:)
               E(:,2) = aux_E(:)
               q(:,:,2) = aux_q(:,:)
               eta(:,2) = aux_eta(:)
               mu_e(:,:) = aux_envit(:,:)
               !$acc end kernels
            else
               do__iteration = .false.
            end if
            if(this%doTimerAnalysis) iStepStartTime = MPI_Wtime()
            call this%callTimeIntegration(istep)
            if(flag_implicit == 1) then
               if((this%currentNonLinealIter .gt. maxIterNonLineal) .and. (inonLineal .lt. 4)) then
                  inonLineal = inonLineal + 1
                  pseudo_cfl = pseudo_cfl*0.5_rp
                  if(mpi_rank.eq.0) write(111,*)"(WARRNING)  non lineal iteration failed in time ",istep," new pseudo cfl ",pseudo_cfl," non lineal tries ",inonLineal
                  call flush(111)
               else
                  do__iteration = .false.
               end if
            end if
         end do

         if(flag_implicit == 1) then
            !$acc kernels
            rho(:,3) = rho(:,1)
            E(:,3) = E(:,1)
            q(:,:,3) = q(:,:,1)
            eta(:,3) = eta(:,1)
            !$acc end kernels
            pseudo_cfl = aux_pseudo_cfl
         end if

         if(this%doTimerAnalysis) then
            iStepEndTime = MPI_Wtime()
            iStepTimeRank = iStepEndTime - iStepStartTime
            call MPI_Allreduce(iStepTimeRank,iStepTimeMax,1,mpi_datatype_real8,MPI_MAX,MPI_COMM_WORLD,mpi_err)
            inv_iStep = 1.0_rp/real(istep)
            iStepAvgTime = (iStepAvgTime*(istep-1)+iStepTimeMax)*inv_iStep

            if((mpi_rank.eq.0).and.(istep==this%nsave2)) then
               write(123,*) istep,iStepTimeMax,iStepAvgTime
               call flush(123)
            end if
         end if

         this%time = this%time+this%dt

         if (istep == this%nsave2 .and. (this%doGlobalAnalysis)) then
            call volAvg_EK(numElemsRankPar,numNodesRankPar,connecParWork,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
            call visc_dissipationRate(numElemsRankPar,numNodesRankPar,connecParWork,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
            call maxMach(numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,machno,this%maxmachno)
            call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
            if(mpi_rank.eq.0) then
               write(111,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
               write(111,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
               20 format(6(F16.8,2X))
               call flush(666)
            end if
         end if

         if (flag_real_diff == 1) then
            call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,rho(:,2))
         else
            call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
         end if

         call nvtxEndRange

         !
         ! Update the accumulators for averaging
         !
         call nvtxStartRange("Accumulate"//timeStep,istep)
         call favre_average(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,this%dt,rho,u,pr, &
                            mu_fluid,mu_e,mu_sgs,tauw,this%acutim,acurho,acupre,acuvel,acuve2,acuvex,acumueff,acutw)
         call nvtxEndRange

         if (istep == this%nsave2) then
            if(istep==this%nsave2.and.mpi_rank.eq.0) write(111,*) "step",istep,"time:",this%time,"s (dt",this%dt,"s)"

            if (isMeshBoundaries) then
               do icode = 1,numBoundCodes!this%numCodes
                  call nvtxStartRange("Surface info")
                  call surfInfo(istep,this%time,numElemsRankPar,numNodesRankPar,numBoundsRankPar,icode,connecParWork,boundPar,point2elem, &
                     bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
                     mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(:,iCode),Ftau(:,iCode))

                  call nvtxEndRange
                  if(mpi_rank.eq.0) call flush(888+icode)
               end do
            end if
         end if

         !
         ! Output the averages after some steps (user defined)
         !
         if (istep == this%nsaveAVG) then
            call nvtxStartRange("Output AVG"//timeStep,istep)
            call this%saveAverages(istep)
            this%nsaveAVG = this%nsaveAVG+this%nleapAVG
            call nvtxEndRange
         end if

         !
         ! Call VTK output
         !
         if (istep == this%nsave) then
            if (mpi_rank.eq.0) write(111,*) ' -Saving file step: ',istep
            call nvtxStartRange("Output "//timeStep,istep)
            call compute_fieldDerivs(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)
            call this%savePosprocessingFields(istep)
            this%nsave = this%nsave+this%nleap
            call nvtxEndRange
         end if

         call this%afterDt(istep)

         if(istep==this%nsave2) then
            this%nsave2 = this%nsave2+this%nleap2
            if(mpi_rank.eq.0) call flush(111)
         end if

         counter = counter+1

         ! End simulation when physical time is reached (user defined)
         if (this%time .ge. this%maxPhysTime) then
            write(111,*) "--| Time integration finished at step: ",istep,"| time: ",this%time
            ! TODO: check if we want to save the last step
            exit
         end if
      end do
      call nvtxEndRange
   end subroutine CFDSolverBase_evalTimeIteration

   subroutine open_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      character(len=1024) :: aux_string_mpisize
      integer :: iCode

      if(mpi_rank.eq.0) then
         write(aux_string_mpisize,'(I0)') mpi_size

         this%log_file_name = 'sod2d_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.log'
         if(this%continue_oldLogs) then
            open(unit=111,file=this%log_file_name,status='old',position='append')
         else
            open(unit=111,file=this%log_file_name,status='replace')
         end if
      end if

      if(mpi_rank.eq.0) then
         write(111,*) "--| Flags defined in the current case:"
         write(111,*) "    flag_real_diff: ",           flag_real_diff
         write(111,*) "    flag_diff_suth: ",           flag_diff_suth
         write(111,*) "    flag_rk_order: ",            flag_rk_order
         write(111,*) "    flag_les: ",                 flag_les
         write(111,*) "    flag_les_ilsa: ",            flag_les_ilsa
         write(111,*) "    flag_solver_type: ",         flag_solver_type
         write(111,*) "    flag_spectralElem: ",        flag_spectralElem
         write(111,*) "    flag_normalise_entropy: ",   flag_normalise_entropy
         write(111,*) "--------------------------------------"
         write(111,*) "    ce: ",      ce
         write(111,*) "    cmax: ",    cmax
         write(111,*) "    c_sgs: ",   c_sgs
         write(111,*) "    cglob: ",   cglob
         write(111,*) "    c_rho: ",   c_rho
         write(111,*) "    c_ener: ",  c_ener
         write(111,*) "    stau: ",    stau
         write(111,*) "    T_ilsa: ",  T_ilsa
         write(111,*) "--------------------------------------"
      end if

   end subroutine open_log_file

   subroutine open_analysis_files(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      character(len=1024) :: filenameAnalysis,fileNameTimer,filenameBound,aux_string_mpisize,aux_string_code
      integer :: iCode

      if(mpi_rank.eq.0) then
         write(aux_string_mpisize,'(I0)') mpi_size

         if(this%doGlobalAnalysis) then
            filenameAnalysis = 'analysis_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.dat'
            if(this%continue_oldLogs) then
               open(unit=666,file=filenameAnalysis,status='old',position='append')
            else
               open(unit=666,file=filenameAnalysis,status='replace')
            end if
         end if

         if(this%doTimerAnalysis) then
            fileNameTimer = 'timer_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.log'
            open(unit=123,file=fileNameTimer,status='replace')
            write(123,*) "iter iteTime iteTimeAvg"
         end if

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               write(aux_string_code,'(I0)') iCode
               filenameBound = 'surf_code_'//trim(aux_string_code)//'-'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.dat'
               if(this%continue_oldLogs) then
                  open(unit=888+iCode,form='formatted',file=filenameBound,status='old',position='append')
               else
                  open(unit=888+iCode,form='formatted',file=filenameBound,status='replace')
                  write(888+iCode,60) "ITER", "TIME", "AREA", "FPRE_X", "FPRE_Y", "FPRE_Z", "FTAU_X", "FTAU_Y", "FTAU_Z"
                  60 format(9(3X,A,5X))
               end if
            end do
         end if
      end if

   end subroutine open_analysis_files

   subroutine close_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         close(unit=111)
      end if
   end subroutine close_log_file

   subroutine close_analysis_files(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         if(this%doGlobalAnalysis) close(unit=666)
         if(this%doTimerAnalysis)  close(unit=123)

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               close(unit=888+iCode)
            end do
         end if
      end if
   end subroutine close_analysis_files

   subroutine flush_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         flush(unit=111)

         if(this%doGlobalAnalysis) then
            flush(unit=666)
         end if

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               flush(unit=888+iCode)
            end do
         end if
      end if

   end subroutine flush_log_file

   subroutine eval_vars_after_load_hdf5_resultsFile(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iNodeL,idime
      real(rp) :: umag

      !values loaded -> rho,u,pr,E,mu_e,mu_sgs

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         e_int(iNodeL,2) = E(iNodeL,2)/rho(iNodeL,2) - 0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))
         !e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
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
         umag = 0.0_rp
         !$acc loop seq
         do idime = 1,ndime
            umag = umag + u(iNodeL,idime,2)**2
         end do
         umag = sqrt(umag)
         machno(iNodeL) = umag/csound(iNodeL)
         !machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

      !$acc kernels
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels

   end subroutine eval_vars_after_load_hdf5_resultsFile

   subroutine CFDSolverBase_run(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this

      ! Init MPI
      call init_mpi()

      ! Init HDF5 interface
      call init_hdf5_interface()

      ! Main simulation parameters
      call this%initializeDefaultParameters()
      call this%initializeParameters()

      ! Define vector length to be used
      call define_veclen()

      ! Open log file
      call this%open_log_file()

      ! read the mesh
      call this%openMesh()

      ! Open analysis files
      call this%open_analysis_files

      ! Eval shape Functions
      call this%evalShapeFunctions()

      ! Allocate variables
      call this%allocateVariables()

      ! Eval or load initial conditions
      call this%evalOrLoadInitialConditions()

      ! Init of the source terms
      call this%initializeSourceTerms()

      ! Eval  viscosty factor
      call this%evalViscosityFactor()

      ! Eval initial viscosty
      call this%evalInitialViscosity()

      ! Compute characteristic size of the elements
      call this%evalCharLength()

      ! Eval boundary element normal
      call this%evalBoundaryNormals()

      ! Eval Jacobian information
      call this%evalJacobians()

      ! Eval AtoIJK inverse
      call this%evalAtoIJKInverse()

      ! Eval BoundaryFacesToNodes
      call  this%boundaryFacesToNodes()

      ! Eval list Elems per Node and Near Boundary Node
      call this%eval_elemPerNode_and_nearBoundaryNode()

      ! Eval mass
      call this%evalMass()

      ! Eval first output
      if(this%isFreshStart) call this%evalFirstOutput()
      call this%flush_log_file()

      ! Eval Normal Faces To Nodes
      if(this%isWallModelOn) call  this%normalFacesToNodes()

      ! Eval initial time step
      call this%evalInitialDt()

      call this%initialBuffer()

      call this%flush_log_file()

      ! Do the time iteration
      call this%evalTimeIteration()

      call this%close_log_file()
      call this%close_analysis_files()

      ! End hdf5 interface
      call end_hdf5_interface()

      ! End comms
      call end_comms()
      call end_comms_bnd()

      ! End MPI
      call end_mpi()

   end subroutine CFDSolverBase_run

end module CFDSolverBase_mod
