module mod_arrays
        use mod_constants

      implicit none

      ! main allocatable arrays
      !integer(4), allocatable  :: connecVTK(:,:)!, connec(:,:), bound(:,:), ldof(:), lbnodes(:), bou_codes(:,:)
      !integer(4), allocatable  :: masSla(:,:), connec_orig(:,:), bound_orig(:,:), lpoin_w(:)
      integer(4), allocatable  :: lelpn(:),point2elem(:)
      integer(4), allocatable  :: atoIJ(:),atoIJK(:),listHEX08(:,:),lnbn(:,:),invAtoIJK(:,:,:),gmshAtoI(:),gmshAtoJ(:),gmshAtoK(:)
!      integer(4), allocatable  :: atoIJ(:), atoIJK(:), vtk_atoIJK(:), listHEX08(:,:), connecLINEAR(:,:),lnbn(:,:),invAtoIJK(:,:,:),gmshAtoI(:),gmshAtoJ(:),gmshAtoK(:)
!      real(rp), allocatable    :: coord(:,:), helem(:),helem_l(:,:)
      real(rp), allocatable    :: helem(:),helem_l(:,:)
      real(rp), allocatable    :: xgp(:,:), wgp(:), xgp_b(:,:), wgp_b(:)
      real(rp), allocatable    :: Ngp(:,:), dNgp(:,:,:), Ngp_b(:,:), dNgp_b(:,:,:)
      real(rp), allocatable    :: Ngp_l(:,:), dNgp_l(:,:,:),dlxigp_ip(:,:,:)
      real(rp), allocatable    :: Je(:,:), He(:,:,:,:), bou_norm(:,:)
      real(rp) , allocatable   :: gpvol(:,:,:), gradRho(:,:), curlU(:,:), divU(:), Qcrit(:)
      real(rp), allocatable    :: u(:,:,:), q(:,:,:), rho(:,:), pr(:,:), E(:,:), Tem(:,:), e_int(:,:), csound(:), eta(:,:), machno(:)
      real(rp), allocatable    :: Ml(:)
      real(rp), allocatable    :: mu_e(:,:),mu_fluid(:),mu_sgs(:,:),mu_factor(:)
      real(rp), allocatable    :: source_term(:)
      real(rp), allocatable    :: acurho(:), acupre(:), acuvel(:,:), acuve2(:,:), acumueff(:)
      real(rp), allocatable    :: avrho(:), avpre(:), avvel(:,:), avve2(:,:), avmueff(:)
      real(rp), allocatable    :: kres(:),etot(:),au(:,:),ax1(:),ax2(:),ax3(:)
      real(rp), allocatable    :: Fpr(:,:), Ftau(:,:)

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
      use mod_constants
      use mod_time_ops
      use mod_fluid_viscosity
      use mod_postpro
      use mod_aver
      use mod_mpi
      use mod_mpi_mesh
      use mod_hdf5
      use mod_comms
   implicit none
   private

   type, public :: CFDSolverBase

      ! main integer parameters
      !NOTES @JORDI: -> it would be nice if nelem, npoin, nper... etc dissapear from here!
      integer(4), public :: nstep, nper, numCodes 
      !integer(4), public  :: nelem, npoin, nboun, nbcodes
      !integer(4), public  :: ndof, nbnodes
      integer(4), public :: ppow
      integer(4), public :: nsave, nleap
      integer(4), public :: nsave2, nleap2
      integer(4), public :: nsaveAVG, nleapAVG
      integer(4), public :: counter, npoin_w
      integer(4), public :: load_step, initial_istep
      logical, public    :: isPeriodic=.false.,loadMesh=.false.,doGlobalAnalysis=.false.,isFreshStart=.true.
      logical, public    :: loadResults=.false.,continue_oldLogs=.false.,interpInitialResults=.false. 
      logical, public    :: useIntInComms=.false.,useFloatInComms=.false.,useDoubleInComms=.false.

      ! main char variables
      character(512) :: log_file_name
      character(512) :: gmsh_file_path,gmsh_file_name
      character(512) :: mesh_h5_file_path,mesh_h5_file_name
      character(512) :: results_h5_file_path,results_h5_file_name

      ! main real parameters
      real(rp) , public                   :: cfl_conv, cfl_diff, acutim
      real(rp) , public                   :: leviCivi(3,3,3), surfArea, EK, VolTot, eps_D, eps_S, eps_T, maxmachno
      real(rp) , public                   :: dt, Cp, Rgas, gamma_gas,Prt,tleap,time

   contains
      procedure, public :: printDt => CFDSolverBase_printDt
      procedure, public :: printAll => CFDSolverBase_printAll
      procedure, public :: run => CFDSolverBase_run
      procedure, public :: initializeDefaultParameters => CFDSolverBase_initializeDefaultParameters
      procedure, public :: initializeParameters => CFDSolverBase_initializeParameters
      procedure, public :: initializeSourceTerms => CFDSolverBase_initializeSourceTerms
      procedure, public :: openMesh => CFDSolverBase_openMesh
      procedure, public :: evalCharLength => CFDSolverBase_evalCharLength
      !procedure, public :: splitBoundary => CFDSolverBase_splitBoundary
      procedure, public :: allocateVariables => CFDSolverBase_allocateVariables
      procedure, public :: evalOrLoadInitialConditions => CFDSolverBase_evalOrLoadInitialConditions
      procedure, public :: evalInitialConditions => CFDSolverBase_evalInitialConditions
      procedure, public :: evalInitialViscosity =>CFDSolverBase_evalInitialViscosity
      procedure, public :: evalInitialDt =>CFDSolverBase_evalInitialDt
      procedure, public :: evalShapeFunctions =>CFDSolverBase_evalShapeFunctions
      procedure, public :: interpolateInitialConditions =>CFDSolverBase_interpolateInitialConditions
      procedure, public :: evalBoundaryNormals =>CFDSolverBase_evalBoundaryNormals
      procedure, public :: evalJacobians =>CFDSolverBase_evalJacobians
      !procedure, public :: evalVTKconnectivity =>CFDSolverBase_evalVTKconnectivity
      procedure, public :: evalAtoIJKInverse =>CFDSolverBase_evalAtoIJKInverse
      procedure, public :: evalPeriodic =>CFDSolverBase_evalPeriodic
      procedure, public :: evalMass=>CFDSolverBase_evalMass
      procedure, public :: evalFirstOutput =>CFDSolverBase_evalFirstOutput
      procedure, public :: evalTimeIteration =>CFDSolverBase_evalTimeIteration
      procedure, public :: callTimeIntegration =>CFDSolverBase_callTimeIntegration
      procedure, public :: saveAverages =>CFDSolverBase_saveAverages
      procedure, public :: savePosprocessingFields =>CFDSolverBase_savePosprocessingFields
      procedure, public :: afterDt =>CFDSolverBase_afterDt

      procedure :: open_log_file
      procedure :: close_log_file
      procedure :: eval_vars_after_load_hdf5_resultsFile
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

        allocate(source_term(ndime))
        source_term(1) = 0.00_rp
        source_term(2) = 0.00_rp
        source_term(3) = 0.00_rp

   end subroutine CFDSolverBase_initializeSourceTerms

   subroutine CFDSolverBase_initializeDefaultParameters(this)
      class(CFDSolverBase), intent(inout) :: this

      write(this%gmsh_file_path,*) "./mesh/"
      write(this%gmsh_file_name,*) "meshName" 

      write(this%mesh_h5_file_path,*) "./"
      write(this%mesh_h5_file_name,*) "meshFile"

      write(this%results_h5_file_path,*) "./"
      write(this%results_h5_file_name,*) "resultsFile"

      this%time = 0.0_rp
      this%load_step = 0
      this%initial_istep = 1

      this%isPeriodic = .false.
      this%loadMesh = .false.
      this%loadResults = .false.
      this%continue_oldLogs= .false.
      this%doGlobalAnalysis = .false.
      this%interpInitialResults = .false.
      this%isFreshStart=.true.
      !@JORDI: discuss which other parameters can be set as default....

   end subroutine CFDSolverBase_initializeDefaultParameters

   subroutine CFDSolverBase_initializeParameters(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_initializeParameters

   subroutine CFDSolverBase_openMesh(this)
      class(CFDSolverBase), intent(inout) :: this

      call nvtxStartRange("Open mesh")
#if 1
!---------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) 
      !write(111,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(this%gmsh_file_path))," !"
      !call read_alya_mesh_files(this%gmsh_file_path,this%gmsh_file_name,this%isPeriodic)

      if(this%loadMesh) then
         call load_hdf5_meshfile(this%mesh_h5_file_path,this%mesh_h5_file_name)
      else
         call read_alyaMesh_part_and_create_hdf5Mesh(this%gmsh_file_path,this%gmsh_file_name,this%isPeriodic,&
                                                    this%mesh_h5_file_path,this%mesh_h5_file_name)
      end if
!---------------------------------------------------------------------------------------------------------------
#else
      call read_dims(this%file_path,this%file_name,this%npoin,this%nelem,this%nboun)
      allocate(connec(this%nelem,nnode))
      if (this%nboun .ne. 0) then
         allocate(bound(this%nboun,npbou))
         allocate(bou_codes(this%nboun,2))
         allocate(bou_norm(this%nboun,ndime*npbou))
         call read_fixbou(this%file_path,this%file_name,this%nboun,this%nbcodes,bou_codes)
         this%numCodes = maxval(bou_codes(:,2))
         allocate(Fpr(this%numCodes,ndime))
         allocate(Ftau(this%numCodes,ndime))
         write(111,*) "--| TOTAL BOUNDARY CODES :", this%numCodes
      end if
      allocate(coord(this%npoin,ndime))
      call read_geo_dat(this%file_path,this%file_name,this%npoin,this%nelem,this%nboun,connec,bound,coord)
      if (this%isPeriodic == 1) then
         allocate(masSla(this%nper,2))
         allocate(connec_orig(this%nelem,nnode))
         if (this%nboun .ne. 0) then
            allocate(bound_orig(this%nboun,npbou))
         end if
         call read_periodic(this%file_path,this%file_name,this%nper,masSla)
      end if
#endif
      call nvtxEndRange

   end subroutine CFDSolverBase_openMesh
#if 0
   subroutine CFDSolverBase_splitBoundary(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), allocatable    :: aux1(:)
      integer(4) :: iNodeL,iboun,ipbou,idof,ibnodes

      !@TODO: REVIEW FUNC
      !NOT REQUIRED FOR FULL PERIODIC
      if (this%nboun .ne. 0) then
         write(111,*) "--| SPLITTING BOUNDARY NODES FROM DOFs..."
         call nvtxStartRange("Bnodes split")
         allocate(aux1(numNodesRankPar))

         !
         ! Fill aux1 with all nodes in order
         !
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            aux1(iNodeL) = iNodeL
         end do
         !$acc end parallel loop

         !
         ! If node is on boundary, zero corresponding aux1 entry
         !
         !$acc parallel loop gang 
         do iboun = 1,this%nboun
            !$acc loop vector
            do ipbou = 1,npbou
               aux1(bound(iboun,ipbou)) = 0
            end do
         end do
         !$acc end parallel loop

         !
         ! Determine how many nodes are boundary nodes
         !
         this%ndof = 0
         do iNodeL = 1,numNodesRankPar
            if (aux1(iNodeL) == 0) then
               this%ndof = this%ndof+1
            end if
         end do

         this%nbnodes = this%ndof    ! Nodes on boundaries
         this%ndof = this%npoin-this%ndof ! Free nodes
         write(111,*) '--| TOTAL FREE NODES := ',this%ndof
         write(111,*) '--| TOTAL BOUNDARY NODES := ',this%nbnodes

         allocate(ldof(this%ndof))
         allocate(lbnodes(this%nbnodes))

         !
         ! Split aux1 into the 2 lists
         !
         idof = 0    ! Counter for free nodes
         ibnodes = 0 ! Counter for boundary nodes
         !$acc parallel loop reduction(+:idof,ibnodes)
         do iNodeL = 1,numNodesRankPar
            if (aux1(iNodeL) == 0) then
               ibnodes = ibnodes+1
               lbnodes(ibnodes) = iNodeL
            else
               idof = idof+1
               ldof(idof) = aux1(iNodeL)
            end if
         end do
         !$acc end parallel loop

         deallocate(aux1)
         call nvtxEndRange
      end if
   end subroutine CFDSolverBase_splitBoundary
#endif
   subroutine CFDSolverBase_evalCharLength(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp) :: he_aux
      integer(4) :: iElem

      call nvtxStartRange("Elem size compute")
      allocate(helem(numElemsInRank))
      do iElem = 1,numElemsInRank
         call char_length(iElem,numElemsInRank,numNodesRankPar,connecParOrig,coordPar,he_aux)
         helem(iElem) = he_aux
      end do
      call nvtxEndRange

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
      allocate(q(numNodesRankPar,ndime,2))  ! momentum
      allocate(rho(numNodesRankPar,2))      ! Density
      allocate(pr(numNodesRankPar,2))       ! Pressure
      allocate(E(numNodesRankPar,2))        ! Total Energy
      allocate(Tem(numNodesRankPar,2))      ! Temperature
      allocate(e_int(numNodesRankPar,2))    ! Internal Energy
      allocate(eta(numNodesRankPar,2))      ! entropy
      allocate(csound(numNodesRankPar))     ! Speed of sound
      allocate(machno(numNodesRankPar))     ! Speed of sound
      allocate(mu_fluid(numNodesRankPar))   ! Fluid viscosity
      allocate(mu_factor(numNodesRankPar))  ! Fluid viscosity
      allocate(mu_e(numElemsInRank,ngaus))  ! Elemental viscosity
      allocate(mu_sgs(numElemsInRank,ngaus))! SGS viscosity
      !ilsa
      allocate(kres(numNodesRankPar))
      allocate(etot(numNodesRankPar))
      allocate(au(numNodesRankPar,ndime))
      allocate(ax1(numNodesRankPar))
      allocate(ax2(numNodesRankPar))
      allocate(ax3(numNodesRankPar))

      !boundary
      if(numBoundCodes .ge. 1) then
         allocate(Fpr(numBoundCodes,ndime))
         allocate(Ftau(numBoundCodes,ndime))
      end if

      !*********************************************************************!
      ! Allocate accumulators for averaging process                         !
      !*********************************************************************!
      allocate(acurho(numNodesRankPar))
      allocate(acupre(numNodesRankPar))
      allocate(acumueff(numNodesRankPar))
      allocate(acuvel(numNodesRankPar,ndime))
      allocate(acuve2(numNodesRankPar,ndime))
      allocate(avrho(numNodesRankPar))
      allocate(avpre(numNodesRankPar))
      allocate(avmueff(numNodesRankPar))
      allocate(avvel(numNodesRankPar,ndime))
      allocate(avve2(numNodesRankPar,ndime))

      !*********************************************************************!
      ! Derivative-related fields                                           !
      !*********************************************************************!
      allocate(gradRho(numNodesRankPar,ndime))
      allocate(curlU(numNodesRankPar,ndime))
      allocate(divU(numNodesRankPar))
      allocate(Qcrit(numNodesRankPar))

      !$acc kernels
      acurho(:) = 0.0_rp
      acupre(:) = 0.0_rp
      acumueff(:) = 0.0_rp
      acuvel(:,:) = 0.00_rp
      acuve2(:,:) = 0.00_rp
      acumueff(:) = 0.00_rp
      avrho(:) = 0.0_rp
      avpre(:) = 0.0_rp
      avmueff(:) = 0.0_rp
      avvel(:,:) = 0.00_rp
      avve2(:,:) = 0.00_rp
      avmueff(:) = 0.00_rp
      !$acc end kernels
      this%acutim = 0.0_rp

      call nvtxEndRange

   end subroutine CFDSolverBase_allocateVariables

   subroutine CFDSolverBase_evalOrLoadInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

      if(this%loadResults) then
         call load_hdf5_resultsFile(this%load_step,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,mu_sgs)
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
         call this%evalInitialConditions()
      end if

   end subroutine CFDSolverBase_evalOrLoadInitialConditions

   subroutine CFDSolverBase_evalInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_evalInitialConditions

   subroutine CFDSolverBase_evalInitialViscosity(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ipoin

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
           STOP(1)
        end if

   end subroutine CFDSolverBase_evalInitialViscosity

   subroutine CFDSolverBase_evalInitialDt(this)
      class(CFDSolverBase), intent(inout) :: this

      !*********************************************************************!
      ! Compute initial time-step size                                      !
      !*********************************************************************!

      if (flag_real_diff == 1) then
         call adapt_dt_cfl(numElemsInRank,numNodesRankPar,connecParOrig,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,rho(:,2))
         if(mpi_rank.eq.0) write(111,*) "--| TIME STEP SIZE dt := ",this%dt,"s"
      else
         call adapt_dt_cfl(numElemsInRank,numNodesRankPar,connecParOrig,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
         if(mpi_rank.eq.0) write(111,*) "--| TIME STEP SIZE dt := ",this%dt,"s"
      end if

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

      if(mpi_rank.eq.0) write(111,*) "  --| GENERATING CHEBYSHEV TABLE..."
      call chebyshev_hex(atoIJK,xgp,wgp)
      call chebyshev_qua(atoIJ,xgp_b,wgp_b)

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

   end subroutine CFDSolverBase_evalShapeFunctions

   subroutine CFDSolverBase_interpolateInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp),    allocatable    :: aux_1(:,:), aux_2(:)
      integer(4) :: ielem,inode,igaus,idime

      !*********************************************************************!
      ! Adjust element nodes and variables if spectral 
      ! element type being used                                             !
      !*********************************************************************!

      !----------------------------------------------------------------------------
      !     COORDS
      if(not(this%loadMesh)) then
         if(mpi_rank.eq.0) write(*,*) "--| Interpolating nodes coordinates..."
         allocate(aux_1(numNodesRankPar,ndime))
         aux_1(:,:) = coordPar(:,:)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               do idime = 1,ndime
                  call var_interpolate(aux_1(connecParOrig(ielem,:),idime),Ngp_l(inode,:),coordPar(connecParOrig(ielem,inode),idime))
               end do
            end do
         end do
         deallocate(aux_1)
         call overwrite_coordinates_hdf5()
      end if
      !----------------------------------------------------------------------------
      if(this%interpInitialResults) then
         if(mpi_rank.eq.0) write(*,*) "--| Interpolating initial results..."
         allocate(aux_1(numNodesRankPar,ndime))
         aux_1(:,:) = u(:,:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               do idime = 1,ndime
                  call var_interpolate(aux_1(connecParOrig(ielem,:),idime),Ngp_l(inode,:),u(connecParOrig(ielem,inode),idime,2))
               end do
            end do
         end do
         aux_1(:,:) = q(:,:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               do idime = 1,ndime
                  call var_interpolate(aux_1(connecParOrig(ielem,:),idime),Ngp_l(inode,:),q(connecParOrig(ielem,inode),idime,2))
               end do
            end do
         end do
         deallocate(aux_1)
         allocate(aux_2(numNodesRankPar))
         aux_2(:) = rho(:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),rho(connecParOrig(ielem,inode),2))
            end do
         end do
         aux_2(:) = pr(:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),pr(connecParOrig(ielem,inode),2))
            end do
         end do
         aux_2(:) = E(:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),E(connecParOrig(ielem,inode),2))
            end do
         end do
         aux_2(:) = Tem(:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),Tem(connecParOrig(ielem,inode),2))
            end do
         end do
         aux_2(:) = e_int(:,2)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),e_int(connecParOrig(ielem,inode),2))
            end do
         end do
         aux_2(:) = csound(:)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),csound(connecParOrig(ielem,inode)))
            end do
         end do
         aux_2(:) = machno(:)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),machno(connecParOrig(ielem,inode)))
            end do
         end do
         aux_2(:) = mu_fluid(:)
         do ielem = 1,numElemsInRank
            do inode = (2**ndime)+1,nnode
               call var_interpolate(aux_2(connecParOrig(ielem,:)),Ngp_l(inode,:),mu_fluid(connecParOrig(ielem,inode)))
            end do
         end do
         deallocate(aux_2)
      end if

   end subroutine CFDSolverBase_interpolateInitialConditions

   subroutine CFDSolverBase_evalBoundaryNormals(this)
      class(CFDSolverBase), intent(inout) :: this

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
      if (isMeshBoundaries) then
         if(mpi_rank.eq.0) write(111,*) "--| COMPUTING BOUNDARY ELEMENT NORMALS"
         allocate(boundNormalPar(numBoundsRankPar,ndime*npbou))
         call nvtxStartRange("Bou normals")
         call boundary_normals(numNodesRankPar,numBoundsRankPar,boundPar,this%leviCivi,coordPar,dNgp_b,boundNormalPar)
         call nvtxEndRange
      end if
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
      allocate(He(ndime,ndime,ngaus,numElemsInRank))
      allocate(gpvol(1,ngaus,numElemsInRank))
      call elem_jacobian(numElemsInRank,numNodesRankPar,connecParOrig,coordPar,dNgp,wgp,gpvol,He) 
      call  nvtxEndRange

      vol_rank  = 0.0
      vol_tot_d = 0.0
      do ielem = 1,numElemsInRank
         do igaus = 1,ngaus
            vol_rank = vol_rank+gpvol(1,igaus,ielem)
         end do
      end do

      call MPI_Allreduce(vol_rank,vol_tot_d,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,mpi_err)

      this%VolTot = real(vol_tot_d,rp) 

      if(mpi_rank.eq.0) write(111,*) '--| DOMAIN VOLUME := ',this%VolTot

   end subroutine CFDSolverBase_evalJacobians

   subroutine CFDSolverBase_evalVTKconnectivity(this)
      class(CFDSolverBase), intent(inout) :: this

      !allocate(connecVTK(numElemsInRank,nnode))
      !allocate(connecLINEAR(numElemsInRank*(porder**ndime),2**ndime))
      !call create_connecVTK(numElemsInRank,connecParOrig,atoIJK,vtk_atoIJK,connecVTK)
      !call linearMeshOutput(numElemsInRank,connecParOrig,listHEX08,connecLINEAR)

   end subroutine CFDSolverBase_evalVTKconnectivity
   
   subroutine CFDSolverBase_evalAtoIJKInverse(this)
      class(CFDSolverBase), intent(inout) :: this

      allocate(invAtoIJK(porder+1,porder+1,porder+1))
      allocate(gmshAtoI(nnode))
      allocate(gmshAtoJ(nnode))
      allocate(gmshAtoK(nnode))
      call atioIJKInverse(atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)

   end subroutine CFDSolverBase_evalAtoIJKInverse

   subroutine CFDSolverBase_evalPeriodic(this)
      class(CFDSolverBase), intent(inout) :: this
      !integer(4) ipoin

      !REVISAR FUNCIONS PREVIES SI NO FAIG EL WORKING LIST AQUI! SOBRETOT PERQUE EM CANVIA EL CONNEC!
      !call create_working_lists() !located in mod_mpi_mesh, rethink name of the func...
#if 0
      if (this%isPeriodic) then
         if (this%nboun .eq. 0) then
            call periodic_ops(numElemsInRank,numNodesRankPar,this%nboun,this%npoin_w,this%nper, &
               lpoin_w,connecPar,connecParOrig,masSla)
         else
            call periodic_ops(numElemsInRank,numNodesRankPar,this%nboun,this%npoin_w,this%nper, &
               lpoin_w,connecPar,connecParOrig,masSla,bound,bound_orig)
         end if
      else
         this%npoin_w = this%npoin
         allocate(lpoin_w(this%npoin_w)) ! All nodes are working nodes
         !$acc parallel loop
         do ipoin = 1,this%npoin_w
            lpoin_w(ipoin) = ipoin
         end do
         !$acc end parallel loop
      end if
#endif
      !*********************************************************************!
      ! Compute list of elements per node (connectivity index)              !
      !*********************************************************************!
      ! evaluate near boundaries for the inlets and outlets
      !not the best place Oriol!
      allocate(lelpn(numNodesRankPar))
      allocate(point2elem(numNodesRankPar))
      if(mpi_rank.eq.0) write(111,*) '--| POINT 2 ELEM begin'
      call elemPerNode(numElemsInRank,numNodesRankPar,connecParWork,lelpn,point2elem)

      if(mpi_rank.eq.0) write(111,*) '--| POINT 2 ELEM done'
      allocate(lnbn(numBoundsRankPar,npbou))
      call nearBoundaryNode(numElemsInRank,numNodesRankPar,numBoundsRankPar,connecParWork,coordPar,boundPar,point2elem,atoIJK,lnbn)

   end subroutine CFDSolverBase_evalPeriodic

   subroutine CFDSolverBase_evalMass(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iElem

      !*********************************************************************!
      ! Compute mass matrix (Lumped and Consistent) and set solver type     !
      !*********************************************************************!

      if(mpi_rank.eq.0) write(111,*) '--| COMPUTING LUMPED MASS MATRIX...'
      call nvtxStartRange("Lumped mass compute")
      allocate(Ml(numNodesRankPar))
      call lumped_mass_spectral(numElemsInRank,numNodesRankPar,connecParWork,gpvol,Ml)
      call nvtxEndRange

      !charecteristic length for spectral elements for the entropy
      !stablisation
      allocate(helem_l(numElemsInRank,nnode))
      do iElem = 1,numElemsInRank
         call char_length_spectral(iElem,numElemsInRank,numNodesRankPar,connecParWork,coordPar,Ml,helem_l)
      end do
   end subroutine CFDSolverBase_evalMass

   subroutine CFDSolverBase_evalFirstOutput(this)
      class(CFDSolverBase), intent(inout) :: this
      character(500) :: tmpname
      integer(4) :: iCode

      !*********************************************************************!
      ! Compute surface forces and area                                                                !
      !*********************************************************************!
      if (isMeshBoundaries) then
         do iCode = 1,numBoundCodes
            call nvtxStartRange("Surface info")
            call surfInfo(0,0.0_rp,numElemsInRank,numNodesRankPar,numBoundsRankPar,iCode,connecParWork,boundPar,point2elem,&
               bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
               mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(iCode,:),Ftau(iCode,:))
            call nvtxEndRange
         end do
      end if

      call compute_fieldDerivs(numElemsInRank,numNodesRankPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)

      call volAvg_EK(numElemsInRank,numNodesRankPar,connecParWork,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
      call visc_dissipationRate(numElemsInRank,numNodesRankPar,connecParWork,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
      call maxMach(numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,machno,this%maxmachno)
      call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
      if(mpi_rank.eq.0) then
         write(111,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
         write(111,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
         20 format(6(F16.8,2X))
      end if

   end subroutine CFDSolverBase_evalFirstOutput

   subroutine CFDSolverBase_callTimeIntegration(this)
      class(CFDSolverBase), intent(inout) :: this
      if(mpi_rank.eq.0) write(111,*) " Time integration should be overwritted"
      STOP(1)

   end subroutine CFDSolverBase_callTimeIntegration

   subroutine CFDSolverBase_saveAverages(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      if(mpi_rank.eq.0) write(111,*) " Save Averages should be overwritted"
      STOP(1)

   end subroutine CFDSolverBase_saveAverages

   subroutine CFDSolverBase_afterDt(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

   end subroutine CFDSolverBase_afterDt

   subroutine CFDSolverBase_savePosprocessingFields(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      if(mpi_rank.eq.0) write(111,*) " Save Posprocessing Fields should be overwritted"
      STOP(1)

   end subroutine CFDSolverBase_savePosprocessingFields

   subroutine CFDSolverBase_evalTimeIteration(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: icode,counter,istep,ppow=1,flag_emac,flag_predic
      character(4) :: timeStep

      counter = 1
      flag_emac = 0
      flag_predic=0

      call nvtxStartRange("Start RK4")
      if(mpi_rank.eq.0) write(*,*) ' Doing evalTimeItarion. Ini step:',this%initial_istep,'| end step:',this%nstep
      ! periodic with boundaries
      do istep = this%initial_istep,this%nstep
         if (istep==this%nsave.and.mpi_rank.eq.0) write(111,*) '   --| STEP: ', istep
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

         call this%callTimeIntegration()

         this%time = this%time+this%dt

         if (istep == this%nsave2 .and. (this%doGlobalAnalysis)) then
            call volAvg_EK(numElemsInRank,numNodesRankPar,connecParWork,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
            call visc_dissipationRate(numElemsInRank,numNodesRankPar,connecParWork,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
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
            call adapt_dt_cfl(numElemsInRank,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,rho(:,2))
            if(istep==this%nsave2.and.mpi_rank.eq.0) write(111,*) "DT := ",this%dt,"s time := ",this%time,"s"
         else
            call adapt_dt_cfl(numElemsInRank,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
            if(istep==this%nsave2.and.mpi_rank.eq.0) write(111,*) "DT := ",this%dt,"s time := ",this%time,"s"
         end if

         call nvtxEndRange

         !
         ! Update the accumulators for averaging
         !
         call nvtxStartRange("Accumulate"//timeStep,istep)
         call favre_average(numElemsInRank,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,this%dt,rho,u,pr, &
            mu_fluid,mu_e,mu_sgs,this%acutim,acurho,acupre,acuvel,acuve2,acumueff)
         call nvtxEndRange

         if (istep == this%nsave2) then
            if (isMeshBoundaries) then
               do icode = 1,numBoundCodes!this%numCodes
                  call nvtxStartRange("Surface info")
                  call surfInfo(istep,this%time,numElemsInRank,numNodesRankPar,numBoundsRankPar,icode,connecParWork,boundPar,point2elem, &
                     bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
                     mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(icode,:),Ftau(icode,:))

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
            call nvtxStartRange("Output "//timeStep,istep)
            call compute_fieldDerivs(numElemsInRank,numNodesRankPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)
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

      end do
      call nvtxEndRange
   end subroutine CFDSolverBase_evalTimeIteration

   subroutine open_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      character(len=1024) :: filename,filenameAnalysis,filenameBound,aux_string_mpisize,aux_string_code
      integer :: iCode

      if(mpi_rank.eq.0) then
         write(aux_string_mpisize,'(I0)') mpi_size

         filename = 'sod2d_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.log'
         if(this%continue_oldLogs) then
            open(unit=111,file=filename,status='old',position='append')
         else
            open(unit=111,file=filename,status='replace')
         end if

         if(this%doGlobalAnalysis) then
            filenameAnalysis = 'analysis_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.dat'
            if(this%continue_oldLogs) then
               open(unit=666,file=filenameAnalysis,status='old',position='append')
            else
               open(unit=666,file=filenameAnalysis,status='replace')
            end if
         end if

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               write(aux_string_code,'(I0)') iCode
               filenameBound = 'surf_code_'//trim(aux_string_code)//'-'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//'.dat'
               if(this%continue_oldLogs) then
                  open(unit=888+iCode,form='formatted',file=filenameBound,status='old',position='append')
               else
                  open(unit=888+iCode,form='formatted',file=filenameBound,status='replace')
               end if
               write(888+iCode,60) "ITER", "TIME", "AREA", "FPRE_X", "FPRE_Y", "FPRE_Z", "FTAU_X", "FTAU_Y", "FTAU_Z"
               60 format(9(3X,A,5X))
            end do
         end if
      end if

   end subroutine open_log_file

   subroutine close_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         close(unit=111)

         if(this%doGlobalAnalysis) then
            close(unit=666)
         end if

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               close(unit=888+iCode)
            end do
         end if
      end if

   end subroutine close_log_file

   subroutine eval_vars_after_load_hdf5_resultsFile(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iNodeL,idime
      real(rp) :: umag

      !values loaded -> rho,u,pr,E

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         e_int(iNodeL,2) = E(iNodeL,2)/rho(iNodeL,2) - 0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))
         !e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
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
      !mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      !mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels

      this%interpInitialResults = .false.

   end subroutine eval_vars_after_load_hdf5_resultsFile

   subroutine CFDSolverBase_run(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this

      ! Init MPI
      call init_mpi()

        ! Main simulation parameters

        call this%initializeDefaultParameters()         
        call this%initializeParameters()

        ! Init of the source terms

        call this%initializeSourceTerms()

        ! Define vector length to be used 
        
        call define_veclen()

      ! init hdf5

      call init_hdf5_interface(this%mesh_h5_file_path,this%mesh_h5_file_name,this%results_h5_file_path,this%results_h5_file_name)

      ! read the mesh

        call this%openMesh()

        ! init comms
        this%useIntInComms=.true.
        this%useFloatInComms=.true.
        this%useDoubleInComms=.false.
        call init_comms(this%useIntInComms,this%useFloatInComms,this%useDoubleInComms)

         ! Open log file
         call this%open_log_file()

        ! Compute characteristic size of the elements

        call this%evalCharLength()

        ! Splitting boundary nodes
        ! now this is done in the parallel mesh!
        !call this%splitBoundary()

        ! Allocate variables

        call this%allocateVariables()

        ! Eval or load initial conditions

        call this%evalOrLoadInitialConditions()
        !call this%evalInitialConditions()

        ! Eval initial viscosty

        call this%evalInitialViscosity()

        ! Eval initial time step

        call this%evalInitialDt()

        ! Eval shape Functions

        call this%evalShapeFunctions()

        ! Interpolate initial conditions

        call this%interpolateInitialConditions()

        ! Eval boundary element normal

        call this%evalBoundaryNormals()

        ! Eval Jacobian information

        call this%evalJacobians()

        ! Eval vtk connectivity

        !call this%evalVTKconnectivity()

        ! Eval AtoIJK inverse

        call this%evalAtoIJKInverse()

        ! Eval periodic

        call this%evalPeriodic()

        ! Eval mass 

        call this%evalMass()

        ! Eval first output
        if(this%isFreshStart) call this%evalFirstOutput()

        ! Do the time iteration

        call this%evalTimeIteration()

      call this%close_log_file()

      ! End hdf5 interface
      call end_hdf5_interface()

      ! End comms
      call end_comms()

      ! End MPI      
      call end_mpi()

   end subroutine CFDSolverBase_run

end module CFDSolverBase_mod
