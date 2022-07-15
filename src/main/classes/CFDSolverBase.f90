
module CFDSolverBase_mod
        use mod_nvtx
        use cudafor
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
   implicit none
   private

   type, public :: CFDSolverBase

      ! main integer parameters
      integer(4), public         :: nstep, nper, numCodes
      integer(4), public         :: nelem, npoin, nboun, nbcodes
      integer(4), public         :: ndof, nbnodes
      integer(4), public         :: ppow
      integer(4), public         :: nsave, nleap
      integer(4), public         :: nsave2, nleap2
      integer(4), public         :: nsaveAVG, nleapAVG
      integer(4), public         :: counter, npoin_w
      integer(4), public         :: isPeriodic, getForces,globalAnalysis

      ! main char variables
      character(500), public    :: file_path
      character(500), public    :: file_name
       

      ! main allocatable arrays
      integer(4), allocatable, public     :: connec(:,:), connecVTK(:,:), bound(:,:), ldof(:), lbnodes(:), bou_codes(:,:)
      integer(4), allocatable, public     :: masSla(:,:), connec_orig(:,:), bound_orig(:,:), lelpn(:),point2elem(:)
      integer(4), allocatable, public     :: lpoin_w(:), atoIJ(:), atoIJK(:), vtk_atoIJK(:), listHEX08(:,:), connecLINEAR(:,:),lnbn(:,:),invAtoIJK(:,:,:),gmshAtoI(:),gmshAtoJ(:),gmshAtoK(:)
      real(rp),    allocatable, public    :: coord(:,:), helem(:),helem_l(:,:)
      real(rp),    allocatable, public    :: xgp(:,:), wgp(:), xgp_b(:,:), wgp_b(:)
      real(rp),    allocatable, public    :: Ngp(:,:), dNgp(:,:,:), Ngp_b(:,:), dNgp_b(:,:,:)
      real(rp),    allocatable, public    :: Ngp_l(:,:), dNgp_l(:,:,:),dlxigp_ip(:,:,:)
      real(rp),    allocatable, public    :: Je(:,:), He(:,:,:,:), bou_norm(:,:)
      real(rp),    allocatable, public    :: gpvol(:,:,:), gradRho(:,:), curlU(:,:), divU(:), Qcrit(:)
      real(rp),    allocatable, public    :: u(:,:,:), q(:,:,:), rho(:,:), pr(:,:), E(:,:), Tem(:,:), e_int(:,:), csound(:), eta(:,:), machno(:)
      real(rp),    allocatable, public    :: Ml(:)
      real(rp),    allocatable, public    :: mu_e(:,:),mu_fluid(:),mu_sgs(:,:),mu_factor(:)
      real(rp),    allocatable, public    :: source_term(:)
      real(rp),    allocatable, public    :: acurho(:), acupre(:), acuvel(:,:), acuve2(:,:), acumueff(:)
      real(rp),    allocatable, public    :: kres(:),etot(:),au(:,:),ax1(:),ax2(:),ax3(:)
      real(rp),    allocatable, public    :: Fpr(:), Ftau(:)

      ! main real parameters
      real(rp) , public                   :: cfl_conv, cfl_diff, acutim
      real(rp) , public                   :: leviCivi(3,3,3)
      real(rp) , public                   :: dt

   contains
      procedure, public :: printDt => CFDSolverBase_printDt
      procedure, public :: printAll => CFDSolverBase_printAll
      procedure, public :: run => CFDSolverBase_run
      procedure, public :: initializeParameters => CFDSolverBase_initializeParameters
      procedure, public :: initializeSourceTerms => CFDSolverBase_initializeSourceTerms
      procedure, public :: readMesh => CFDSolverBase_readMesh
      procedure, public :: evalCharLength => CFDSolverBase_evalCharLength
      procedure, public :: splitBoundary => CFDSolverBase_splitBoundary
      procedure, public :: allocateVariables => CFDSolverBase_allocateVariables
      procedure, public :: evalInitialConditions => CFDSolverBase_evalInitialConditions
      procedure, public :: evalInitialViscosity =>CFDSolverBase_evalInitialViscosity
      procedure, public :: evalInitialDt =>CFDSolverBase_evalInitialDt
      procedure, public :: evalShapeFunctions =>CFDSolverBase_evalShapeFunctions
      procedure, public :: interpolateInitialConditions =>CFDSolverBase_interpolateInitialConditions
      procedure, public :: evalBoundaryNormals =>CFDSolverBase_evalBoundaryNormals
      procedure, public :: evalJacobians =>CFDSolverBase_evalJacobians
      procedure, public :: evalVTKconnectivity =>CFDSolverBase_evalVTKconnectivity
      procedure, public :: evalAtoIJKInverse =>CFDSolverBase_evalAtoIJKInverse
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

        allocate(this%source_term(ndime))
        this%source_term(1) = 0.00_rp
        this%source_term(2) = 0.00_rp
        this%source_term(3) = 0.00_rp

   end subroutine CFDSolverBase_initializeSourceTerms

   subroutine CFDSolverBase_initializeParameters(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_initializeParameters

   subroutine CFDSolverBase_readMesh(this)
      class(CFDSolverBase), intent(inout) :: this

      write(this%file_path,*) "./mesh/"
      write(1,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(this%file_path))," !"
      call nvtxStartRange("Read mesh")

      call read_dims(this%file_path,this%file_name,this%npoin,this%nelem,this%nboun)
      allocate(this%connec(this%nelem,nnode))
      if (this%nboun .ne. 0) then
         allocate(this%bound(this%nboun,npbou))
         allocate(this%bou_codes(this%nboun,2))
         allocate(this%bou_norm(this%nboun,ndime*npbou))
         call read_fixbou(this%file_path,this%file_name,this%nboun,this%nbcodes,this%bou_codes)
         this%numCodes = maxval(this%bou_codes(:,2))
         write(1,*) "--| TOTAL BOUNDARY CODES :", this%numCodes
      end if
      allocate(this%coord(this%npoin,ndime))
      call read_geo_dat(this%file_path,this%file_name,this%npoin,this%nelem,this%nboun,this%connec,this%bound,this%coord)
      if (this%isPeriodic == 1) then
         allocate(this%masSla(this%nper,2))
         allocate(this%connec_orig(this%nelem,nnode))
         if (this%nboun .ne. 0) then
            allocate(this%bound_orig(this%nboun,npbou))
         end if
         call read_periodic(this%file_path,this%file_name,this%nper,this%masSla)
      end if
      call nvtxEndRange

   end subroutine CFDSolverBase_readMesh

   subroutine CFDSolverBase_splitBoundary(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), allocatable    :: aux1(:)
      integer(4) :: ipoin,iboun,ipbou,idof,ibnodes

      if (this%nboun .ne. 0) then
         write(1,*) "--| SPLITTING BOUNDARY NODES FROM DOFs..."
         call nvtxStartRange("Bnodes split")
         allocate(aux1(this%npoin))

         !
         ! Fill aux1 with all nodes in order
         !
         !$acc parallel loop
         do ipoin = 1,this%npoin
            aux1(ipoin) = ipoin
         end do
         !$acc end parallel loop

         !
         ! If node is on boundary, zero corresponding aux1 entry
         !
         !$acc parallel loop gang 
         do iboun = 1,this%nboun
            !$acc loop vector
            do ipbou = 1,npbou
               aux1(this%bound(iboun,ipbou)) = 0
            end do
         end do
         !$acc end parallel loop

         !
         ! Determine how many nodes are boundary nodes
         !
         this%ndof = 0
         do ipoin = 1,this%npoin
            if (aux1(ipoin) == 0) then
               this%ndof = this%ndof+1
            end if
         end do

         this%nbnodes = this%ndof    ! Nodes on boundaries
         this%ndof = this%npoin-this%ndof ! Free nodes
         write(1,*) '--| TOTAL FREE NODES := ',this%ndof
         write(1,*) '--| TOTAL BOUNDARY NODES := ',this%nbnodes

         allocate(this%ldof(this%ndof))
         allocate(this%lbnodes(this%nbnodes))

         !
         ! Split aux1 into the 2 lists
         !
         idof = 0    ! Counter for free nodes
         ibnodes = 0 ! Counter for boundary nodes
         !$acc parallel loop reduction(+:idof,ibnodes)
         do ipoin = 1,this%npoin
            if (aux1(ipoin) == 0) then
               ibnodes = ibnodes+1
               this%lbnodes(ibnodes) = ipoin
            else
               idof = idof+1
               this%ldof(idof) = aux1(ipoin)
            end if
         end do
         !$acc end parallel loop

         deallocate(aux1)
         call nvtxEndRange
      end if
   end subroutine CFDSolverBase_splitBoundary

   subroutine CFDSolverBase_evalCharLength(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp) :: he_aux
      integer(4) :: ielem

      call nvtxStartRange("Elem size compute")
      allocate(this%helem(this%nelem))
      do ielem = 1,this%nelem
         call char_length(ielem,this%nelem,this%npoin,this%connec,this%coord,he_aux)
         this%helem(ielem) = he_aux
      end do
      call nvtxEndRange

   end subroutine CFDSolverBase_evalCharLength

   subroutine CFDSolverBase_allocateVariables(this)
      class(CFDSolverBase), intent(inout) :: this

      WRITE(1,*) "--| ALLOCATING MAIN VARIABLES"
      call nvtxStartRange("Allocate main vars")
      !
      ! Last rank is for prediction-advance related to entropy viscosity,
      ! where 1 is prediction, 2 is final value
      !
      allocate(this%u(this%npoin,ndime,2))  ! Velocity
      allocate(this%q(this%npoin,ndime,2))  ! momentum
      allocate(this%rho(this%npoin,2))      ! Density
      allocate(this%pr(this%npoin,2))       ! Pressure
      allocate(this%E(this%npoin,2))        ! Total Energy
      allocate(this%Tem(this%npoin,2))      ! Temperature
      allocate(this%e_int(this%npoin,2))    ! Internal Energy
      allocate(this%eta(this%npoin,2))      ! entropy
      allocate(this%csound(this%npoin))     ! Speed of sound
      allocate(this%machno(this%npoin))     ! Speed of sound
      allocate(this%mu_fluid(this%npoin))   ! Fluid viscosity
      allocate(this%mu_factor(this%npoin))   ! Fluid viscosity
      allocate(this%mu_e(this%nelem,ngaus))  ! Elemental viscosity
      allocate(this%mu_sgs(this%nelem,ngaus))! SGS viscosity
      !ilsa
      allocate(this%kres(this%npoin))
      allocate(this%etot(this%npoin))
      allocate(this%au(this%npoin,ndime))
      allocate(this%ax1(this%npoin))
      allocate(this%ax2(this%npoin))
      allocate(this%ax3(this%npoin))
      call nvtxEndRange

   end subroutine CFDSolverBase_allocateVariables

   subroutine CFDSolverBase_evalInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_evalInitialConditions

   subroutine CFDSolverBase_evalInitialViscosity(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ipoin

        !$acc parallel loop
        do ipoin = 1,this%npoin
           this%mu_factor(ipoin) = flag_mu_factor
        end do
        !$acc end parallel loop

        if (flag_real_diff == 1) then
           if (flag_diff_suth == 0) then
              call constant_viscosity(this%npoin,0.000055_rp,this%mu_fluid)
           else
              call sutherland_viscosity(this%npoin,this%Tem(:,2),this%mu_factor,this%mu_fluid)
           end if
        else if (flag_real_diff == 0) then
           !$acc kernels
           this%mu_fluid(:) = 0.0_rp
           !$acc end kernels
        else
           write(1,*) "--| DIFFUSION FLAG MUST BE EITHER 0 OR 1, NOT: ",flag_real_diff
           STOP(1)
        end if

   end subroutine CFDSolverBase_evalInitialViscosity

   subroutine CFDSolverBase_evalInitialDt(this)
      class(CFDSolverBase), intent(inout) :: this

      !*********************************************************************!
      ! Compute initial time-step size                                      !
      !*********************************************************************!

      if (flag_real_diff == 1) then
         call adapt_dt_cfl(this%nelem,this%npoin,this%connec,this%helem,this%u(:,:,2),this%csound,this%cfl_conv,this%dt,this%cfl_diff,this%mu_fluid,this%mu_sgs,this%rho(:,2))
         write(1,*) "--| TIME STEP SIZE dt := ",this%dt,"s"
      else
         call adapt_dt_cfl(this%nelem,this%npoin,this%connec,this%helem,this%u(:,:,2),this%csound,this%cfl_conv,this%dt)
         write(1,*) "--| TIME STEP SIZE dt := ",this%dt,"s"
      end if

   end subroutine CFDSolverBase_evalInitialDt

   subroutine CFDSolverBase_evalShapeFunctions(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp)                    :: s, t, z
      integer(4) :: igaus

      !*********************************************************************!
      ! Generate GLL table                                                  !
      !*********************************************************************!

      allocate(this%atoIJ(16))
      allocate(this%atoIJK(64))
      allocate(this%vtk_atoIJK(64))
      call hex64(1.0_rp,1.0_rp,1.0_rp,this%atoIJK,this%vtk_atoIJK)
      call qua16(1.0_rp,1.0_rp,this%atoIJ)

      write(1,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."

      call nvtxStartRange("Gaussian Quadrature")

      allocate(this%xgp(ngaus,ndime))
      allocate(this%wgp(ngaus))
      allocate(this%xgp_b(npbou,ndime-1))
      allocate(this%wgp_b(npbou))

      write(1,*) "  --| GENERATING CHEBYSHEV TABLE..."
      call chebyshev_hex(this%atoIJK,this%xgp,this%wgp)
      call chebyshev_qua(this%atoIJ,this%xgp_b,this%wgp_b)

      call nvtxEndRange

      !*********************************************************************!
      ! Generate N and dN for all GP                                        !
      !*********************************************************************!

      ! TODO: Allow for more element types

      write(1,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."
      call nvtxStartRange("N and dN")

      allocate(this%Ngp(ngaus,nnode),this%dNgp(ndime,nnode,ngaus))
      allocate(this%Ngp_l(ngaus,nnode),this%dNgp_l(ndime,nnode,ngaus))
      allocate(this%Ngp_b(npbou,npbou),this%dNgp_b(ndime-1,npbou,npbou))
      allocate(this%dlxigp_ip(ngaus,ndime,porder+1))

      do igaus = 1,ngaus
         s = this%xgp(igaus,1)
         t = this%xgp(igaus,2)
         z = this%xgp(igaus,3)
         allocate(this%listHEX08((porder**ndime),2**ndime))
         call hex64(s,t,z,this%atoIJK,this%vtk_atoIJK,this%listHEX08,this%Ngp(igaus,:),this%dNgp(:,:,igaus),this%Ngp_l(igaus,:),this%dNgp_l(:,:,igaus),this%dlxigp_ip(igaus,:,:))
      end do
      !
      ! Compute N andd dN for boundary elements
      !
      do igaus = 1,npbou
         s = this%xgp_b(igaus,1)
         t = this%xgp_b(igaus,2)
         call qua16(s,t,this%atoIJ,this%Ngp_b(igaus,:),this%dNgp_b(:,:,igaus))
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

      allocate(aux_1(this%npoin,ndime))
      aux_1(:,:) = this%coord(:,:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(this%connec(ielem,:),idime),this%Ngp_l(inode,:),this%coord(this%connec(ielem,inode),idime))
            end do
         end do
      end do
      aux_1(:,:) = this%u(:,:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(this%connec(ielem,:),idime),this%Ngp_l(inode,:),this%u(this%connec(ielem,inode),idime,2))
            end do
         end do
      end do
      aux_1(:,:) = this%q(:,:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(this%connec(ielem,:),idime),this%Ngp_l(inode,:),this%q(this%connec(ielem,inode),idime,2))
            end do
         end do
      end do
      deallocate(aux_1)
      allocate(aux_2(this%npoin))
      aux_2(:) = this%rho(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%rho(this%connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = this%pr(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%pr(this%connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = this%E(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%E(this%connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = this%Tem(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%Tem(this%connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = this%e_int(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%e_int(this%connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = this%csound(:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%csound(this%connec(ielem,inode)))
         end do
      end do
      aux_2(:) = this%machno(:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%machno(this%connec(ielem,inode)))
         end do
      end do
      aux_2(:) = this%mu_fluid(:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(this%connec(ielem,:)),this%Ngp_l(inode,:),this%mu_fluid(this%connec(ielem,inode)))
         end do
      end do
      deallocate(aux_2)

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
      if (this%nboun .ne. 0) then
         allocate(this%Fpr(ndime))
         allocate(this%Ftau(ndime))
         write(1,*) "--| COMPUTING BOUNDARY ELEMENT NORMALS"
         call nvtxStartRange("Bou normals")
         call boundary_normals(this%npoin,this%nboun,this%bound,this%leviCivi,this%coord,this%dNgp_b,this%bou_norm)
         call nvtxEndRange
      end if
   end subroutine CFDSolverBase_evalBoundaryNormals

   subroutine CFDSolverBase_evalJacobians(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ielem, igaus
      real(rp) :: VolTot

      !*********************************************************************!
      ! Generate Jacobian related information                               !
      !*********************************************************************!

      write(1,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

      call nvtxStartRange("Jacobian info")
      allocate(this%He(ndime,ndime,ngaus,this%nelem))
      allocate(this%gpvol(1,ngaus,this%nelem))
      call elem_jacobian(this%nelem,this%npoin,this%connec,this%coord,this%dNgp,this%wgp,this%gpvol,this%He)
      call  nvtxEndRange
      VolTot = 0.0_rp
      do ielem = 1,this%nelem
         do igaus = 1,ngaus
            VolTot = VolTot+this%gpvol(1,igaus,ielem)
         end do
      end do
      write(1,*) '--| DOMAIN VOLUME := ',VolTot

   end subroutine CFDSolverBase_evalJacobians

   subroutine CFDSolverBase_evalVTKconnectivity(this)
      class(CFDSolverBase), intent(inout) :: this

      allocate(this%connecVTK(this%nelem,nnode))
      allocate(this%connecLINEAR(this%nelem*(porder**ndime),2**ndime))
      call create_connecVTK(this%nelem,this%connec,this%atoIJK,this%vtk_atoIJK,this%connecVTK)
      call linearMeshOutput(this%nelem,this%connec,this%listHEX08,this%connecLINEAR)

   end subroutine CFDSolverBase_evalVTKconnectivity
   
   subroutine CFDSolverBase_evalAtoIJKInverse(this)
      class(CFDSolverBase), intent(inout) :: this

      allocate(this%invAtoIJK(porder+1,porder+1,porder+1))
      allocate(this%gmshAtoI(nnode))
      allocate(this%gmshAtoJ(nnode))
      allocate(this%gmshAtoK(nnode))
      call atioIJKInverse(this%atoIJK,this%invAtoIJK,this%gmshAtoI,this%gmshAtoJ,this%gmshAtoK)

   end subroutine CFDSolverBase_evalAtoIJKInverse

   subroutine CFDSolverBase_run(this)
      class(CFDSolverBase), intent(inout) :: this

        open(unit=1,file="sod2d.log",status="replace")

        ! Main simulation parameters
        
        call this%initializeParameters()

        ! Init of the source terms

        call this%initializeSourceTerms()

        ! Define vector length to be used 
        
        call define_veclen()

        ! Read the mesh

        call this%readMesh()

        ! Compute characteristic size of the elements

        call this%evalCharLength()

        ! Splitting boundary nodes

        call this%splitBoundary()

        ! Allocate variables                                                  !

        call this%allocateVariables()

        ! Eval initial conditions

        call this%evalInitialConditions()

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

        call this%evalVTKconnectivity()

        ! Eval AtoIJK inverse

        call this%evalAtoIJKInverse()

   end subroutine CFDSolverBase_run

end module CFDSolverBase_mod
