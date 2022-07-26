module mod_arrays
        use mod_constants

      implicit none

      ! main allocatable arrays
      integer(4), allocatable  :: connec(:,:), connecVTK(:,:), bound(:,:), ldof(:), lbnodes(:), bou_codes(:,:)
      integer(4), allocatable  :: masSla(:,:), connec_orig(:,:), bound_orig(:,:), lelpn(:),point2elem(:)
      integer(4), allocatable  :: lpoin_w(:), atoIJ(:), atoIJK(:), vtk_atoIJK(:), listHEX08(:,:), connecLINEAR(:,:),lnbn(:,:),invAtoIJK(:,:,:),gmshAtoI(:),gmshAtoJ(:),gmshAtoK(:)
      real(rp), allocatable    :: coord(:,:), helem(:),helem_l(:,:)
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
      real(rp), allocatable    :: kres(:),etot(:),au(:,:),ax1(:),ax2(:),ax3(:)
      real(rp), allocatable    :: Fpr(:), Ftau(:)

end module mod_arrays

module CFDSolverBase_mod
        use mod_arrays
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
      integer(4), public         :: isPeriodic=0,globalAnalysis=0

      ! main char variables
      character(500), public    :: file_path
      character(500), public    :: file_name
       


      ! main real parameters
      real(rp) , public                   :: cfl_conv, cfl_diff, acutim
      real(rp) , public                   :: leviCivi(3,3,3), surfArea, EK, VolTot, eps_D, eps_S, eps_T, maxmachno
      real(rp) , public                   :: dt, Cp, Rgas, gamma_gas,Prt,tleap,time

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
      procedure, public :: evalPeriodic =>CFDSolverBase_evalPeriodic
      procedure, public :: evalMass=>CFDSolverBase_evalMass
      procedure, public :: evalFirstOutput =>CFDSolverBase_evalFirstOutput
      procedure, public :: evalTimeIteration =>CFDSolverBase_evalTimeIteration
      procedure, public :: callTimeIntegration =>CFDSolverBase_callTimeIntegration
      procedure, public :: saveAverages =>CFDSolverBase_saveAverages
      procedure, public :: savePosprocessingFields =>CFDSolverBase_savePosprocessingFields
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

   subroutine CFDSolverBase_initializeParameters(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_initializeParameters

   subroutine CFDSolverBase_readMesh(this)
      class(CFDSolverBase), intent(inout) :: this

      write(this%file_path,*) "./mesh/"
      write(1,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(this%file_path))," !"
      call nvtxStartRange("Read mesh")

      call read_dims(this%file_path,this%file_name,this%npoin,this%nelem,this%nboun)
      allocate(connec(this%nelem,nnode))
      if (this%nboun .ne. 0) then
         allocate(bound(this%nboun,npbou))
         allocate(bou_codes(this%nboun,2))
         allocate(bou_norm(this%nboun,ndime*npbou))
         call read_fixbou(this%file_path,this%file_name,this%nboun,this%nbcodes,bou_codes)
         this%numCodes = maxval(bou_codes(:,2))
         write(1,*) "--| TOTAL BOUNDARY CODES :", this%numCodes
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
               aux1(bound(iboun,ipbou)) = 0
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

         allocate(ldof(this%ndof))
         allocate(lbnodes(this%nbnodes))

         !
         ! Split aux1 into the 2 lists
         !
         idof = 0    ! Counter for free nodes
         ibnodes = 0 ! Counter for boundary nodes
         !$acc parallel loop reduction(+:idof,ibnodes)
         do ipoin = 1,this%npoin
            if (aux1(ipoin) == 0) then
               ibnodes = ibnodes+1
               lbnodes(ibnodes) = ipoin
            else
               idof = idof+1
               ldof(idof) = aux1(ipoin)
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
      allocate(helem(this%nelem))
      do ielem = 1,this%nelem
         call char_length(ielem,this%nelem,this%npoin,connec,coord,he_aux)
         helem(ielem) = he_aux
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
      allocate(u(this%npoin,ndime,2))  ! Velocity
      allocate(q(this%npoin,ndime,2))  ! momentum
      allocate(rho(this%npoin,2))      ! Density
      allocate(pr(this%npoin,2))       ! Pressure
      allocate(E(this%npoin,2))        ! Total Energy
      allocate(Tem(this%npoin,2))      ! Temperature
      allocate(e_int(this%npoin,2))    ! Internal Energy
      allocate(eta(this%npoin,2))      ! entropy
      allocate(csound(this%npoin))     ! Speed of sound
      allocate(machno(this%npoin))     ! Speed of sound
      allocate(mu_fluid(this%npoin))   ! Fluid viscosity
      allocate(mu_factor(this%npoin))   ! Fluid viscosity
      allocate(mu_e(this%nelem,ngaus))  ! Elemental viscosity
      allocate(mu_sgs(this%nelem,ngaus))! SGS viscosity
      !ilsa
      allocate(kres(this%npoin))
      allocate(etot(this%npoin))
      allocate(au(this%npoin,ndime))
      allocate(ax1(this%npoin))
      allocate(ax2(this%npoin))
      allocate(ax3(this%npoin))

      !*********************************************************************!
      ! Allocate accumulators for averaging process                         !
      !*********************************************************************!
      allocate(acurho(this%npoin))
      allocate(acupre(this%npoin))
      allocate(acumueff(this%npoin))
      allocate(acuvel(this%npoin,ndime))
      allocate(acuve2(this%npoin,ndime))

      !$acc kernels
      acurho(:) = 0.0_rp
      acupre(:) = 0.0_rp
      acumueff(:) = 0.0_rp
      acuvel(:,:) = 0.00_rp
      acuve2(:,:) = 0.00_rp
      acumueff(:) = 0.00_rp
      !$acc end kernels
      this%acutim = 0.0_rp

      call nvtxEndRange

   end subroutine CFDSolverBase_allocateVariables

   subroutine CFDSolverBase_evalInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_evalInitialConditions

   subroutine CFDSolverBase_evalInitialViscosity(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ipoin

        if (flag_real_diff == 1) then
           if (flag_diff_suth == 0) then
              call constant_viscosity(this%npoin,0.000055_rp,mu_fluid)
           else
              call sutherland_viscosity(this%npoin,Tem(:,2),mu_factor,mu_fluid)
           end if
        else if (flag_real_diff == 0) then
           !$acc kernels
           mu_fluid(:) = 0.0_rp
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
         call adapt_dt_cfl(this%nelem,this%npoin,connec,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,rho(:,2))
         write(1,*) "--| TIME STEP SIZE dt := ",this%dt,"s"
      else
         call adapt_dt_cfl(this%nelem,this%npoin,connec,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
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

      allocate(atoIJ(16))
      allocate(atoIJK(64))
      allocate(vtk_atoIJK(64))
      call hex64(1.0_rp,1.0_rp,1.0_rp,atoIJK,vtk_atoIJK)
      call qua16(1.0_rp,1.0_rp,atoIJ)

      write(1,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."

      call nvtxStartRange("Gaussian Quadrature")

      allocate(xgp(ngaus,ndime))
      allocate(wgp(ngaus))
      allocate(xgp_b(npbou,ndime-1))
      allocate(wgp_b(npbou))

      write(1,*) "  --| GENERATING CHEBYSHEV TABLE..."
      call chebyshev_hex(atoIJK,xgp,wgp)
      call chebyshev_qua(atoIJ,xgp_b,wgp_b)

      call nvtxEndRange

      !*********************************************************************!
      ! Generate N and dN for all GP                                        !
      !*********************************************************************!

      ! TODO: Allow for more element types

      write(1,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."
      call nvtxStartRange("N and dN")

      allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))
      allocate(Ngp_l(ngaus,nnode),dNgp_l(ndime,nnode,ngaus))
      allocate(Ngp_b(npbou,npbou),dNgp_b(ndime-1,npbou,npbou))
      allocate(dlxigp_ip(ngaus,ndime,porder+1))

      do igaus = 1,ngaus
         s = xgp(igaus,1)
         t = xgp(igaus,2)
         z = xgp(igaus,3)
         allocate(listHEX08((porder**ndime),2**ndime))
         call hex64(s,t,z,atoIJK,vtk_atoIJK,listHEX08,Ngp(igaus,:),dNgp(:,:,igaus),Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
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

      allocate(aux_1(this%npoin,ndime))
      aux_1(:,:) = coord(:,:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(connec(ielem,:),idime),Ngp_l(inode,:),coord(connec(ielem,inode),idime))
            end do
         end do
      end do
      aux_1(:,:) = u(:,:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(connec(ielem,:),idime),Ngp_l(inode,:),u(connec(ielem,inode),idime,2))
            end do
         end do
      end do
      aux_1(:,:) = q(:,:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(connec(ielem,:),idime),Ngp_l(inode,:),q(connec(ielem,inode),idime,2))
            end do
         end do
      end do
      deallocate(aux_1)
      allocate(aux_2(this%npoin))
      aux_2(:) = rho(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),rho(connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = pr(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),pr(connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = E(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),E(connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = Tem(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),Tem(connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = e_int(:,2)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),e_int(connec(ielem,inode),2))
         end do
      end do
      aux_2(:) = csound(:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),csound(connec(ielem,inode)))
         end do
      end do
      aux_2(:) = machno(:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),machno(connec(ielem,inode)))
         end do
      end do
      aux_2(:) = mu_fluid(:)
      do ielem = 1,this%nelem
         do inode = (2**ndime)+1,nnode
            call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),mu_fluid(connec(ielem,inode)))
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
         allocate(Fpr(ndime))
         allocate(Ftau(ndime))
         write(1,*) "--| COMPUTING BOUNDARY ELEMENT NORMALS"
         call nvtxStartRange("Bou normals")
         call boundary_normals(this%npoin,this%nboun,bound,this%leviCivi,coord,dNgp_b,bou_norm)
         call nvtxEndRange
      end if
   end subroutine CFDSolverBase_evalBoundaryNormals

   subroutine CFDSolverBase_evalJacobians(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ielem, igaus

      !*********************************************************************!
      ! Generate Jacobian related information                               !
      !*********************************************************************!

      write(1,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

      call nvtxStartRange("Jacobian info")
      allocate(He(ndime,ndime,ngaus,this%nelem))
      allocate(gpvol(1,ngaus,this%nelem))
      call elem_jacobian(this%nelem,this%npoin,connec,coord,dNgp,wgp,gpvol,He)
      call  nvtxEndRange
      this%VolTot = 0.0_rp
      do ielem = 1,this%nelem
         do igaus = 1,ngaus
            this%VolTot = this%VolTot+gpvol(1,igaus,ielem)
         end do
      end do
      write(1,*) '--| DOMAIN VOLUME := ',this%VolTot

   end subroutine CFDSolverBase_evalJacobians

   subroutine CFDSolverBase_evalVTKconnectivity(this)
      class(CFDSolverBase), intent(inout) :: this

      allocate(connecVTK(this%nelem,nnode))
      allocate(connecLINEAR(this%nelem*(porder**ndime),2**ndime))
      call create_connecVTK(this%nelem,connec,atoIJK,vtk_atoIJK,connecVTK)
      call linearMeshOutput(this%nelem,connec,listHEX08,connecLINEAR)

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
      integer(4) ipoin

      if (this%isPeriodic == 1) then
         if (this%nboun .eq. 0) then
            call periodic_ops(this%nelem,this%npoin,this%nboun,this%npoin_w,this%nper, &
               lpoin_w,connec,connec_orig,masSla)
         else
            call periodic_ops(this%nelem,this%npoin,this%nboun,this%npoin_w,this%nper, &
               lpoin_w,connec,connec_orig,masSla,bound,bound_orig)
         end if
      else if (this%isPeriodic == 0) then
         this%npoin_w = this%npoin
         allocate(lpoin_w(this%npoin_w)) ! All nodes are working nodes
         !$acc parallel loop
         do ipoin = 1,this%npoin_w
            lpoin_w(ipoin) = ipoin
         end do
         !$acc end parallel loop
      end if
      !*********************************************************************!
      ! Compute list of elements per node (connectivity index)              !
      !*********************************************************************!
      ! evaluate near boundaries for the inlets and outlets
      !not the best place Oriol!
      allocate(lelpn(this%npoin))
      allocate(point2elem(this%npoin))
      write(1,*) '--| POINT 2 ELEM begin'
      call elemPerNode(this%nelem,this%npoin,connec,lelpn,point2elem)

      write(1,*) '--| POINT 2 ELEM done'

      allocate(lnbn(this%nboun,npbou))
      call nearBoundaryNode(this%nelem,this%npoin,this%nboun,connec,coord,bound,point2elem,atoIJK,lnbn)

   end subroutine CFDSolverBase_evalPeriodic

   subroutine CFDSolverBase_evalMass(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ielem

      !*********************************************************************!
      ! Compute mass matrix (Lumped and Consistent) and set solver type     !
      !*********************************************************************!

      write(1,*) '--| COMPUTING LUMPED MASS MATRIX...'
      call nvtxStartRange("Lumped mass compute")
      allocate(Ml(this%npoin))
      call lumped_mass_spectral(this%nelem,this%npoin,connec,gpvol,Ml)
      call nvtxEndRange

      !charecteristic length for spectral elements for the entropy
      !stablisation
      allocate(helem_l(this%nelem,nnode))
      do ielem = 1,this%nelem
         call char_length_spectral(ielem,this%nelem,this%npoin,connec,coord,Ml,helem_l)
      end do
   end subroutine CFDSolverBase_evalMass

   subroutine CFDSolverBase_evalFirstOutput(this)
      class(CFDSolverBase), intent(inout) :: this
      character(500)             :: tmpname
      integer(4) :: icode

      !*********************************************************************!
      ! Compute surface forces and area                                                                !
      !*********************************************************************!
      if (this%nboun .ne. 0) then
         do icode = 1,this%numCodes
            write(tmpname,'("surfcode_",i0,".dat")') icode
            open(unit=888+icode,form='formatted',file=tmpname,status='replace')
            write(888+icode,60) "ITER", "TIME", "AREA", "FPRE_X", "FPRE_Y", "FPRE_Z", "FTAU_X", "FTAU_Y", "FTAU_Z"
            60            format(9(3X,A,5X))
            call nvtxStartRange("Surface info")
            call surfInfo(0,0.0_rp,this%nelem,this%npoin,this%nboun,icode,connec,bound,point2elem, &
               bou_codes,bou_norm,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coord, &
               mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr,Ftau)
            call nvtxEndRange
         end do
      end if

      !*********************************************************************!
      ! Compute derivative-related fields and produce the 1st output        !
      !*********************************************************************!
      allocate(gradRho(this%npoin,ndime))
      allocate(curlU(this%npoin,ndime))
      allocate(divU(this%npoin))
      allocate(Qcrit(this%npoin))

      call compute_fieldDerivs(this%nelem,this%npoin,connec,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)

      open(unit=666,file="analysis.dat",status="replace")
      call volAvg_EK(this%nelem,this%npoin,connec,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
      call visc_dissipationRate(this%nelem,this%npoin,connec,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
      call maxMach(this%npoin,this%npoin_w,lpoin_w,machno,this%maxmachno)
      call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
      write(1,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
      write(1,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
20      format(6(F16.8,2X))

   end subroutine CFDSolverBase_evalFirstOutput

   subroutine CFDSolverBase_callTimeIntegration(this)
      class(CFDSolverBase), intent(inout) :: this

      write(1,*) " Time integration should be overwritted"
      STOP(1)

   end subroutine CFDSolverBase_callTimeIntegration

   subroutine CFDSolverBase_saveAverages(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      write(1,*) " Save Averages should be overwritted"
      STOP(1)

   end subroutine CFDSolverBase_saveAverages

   subroutine CFDSolverBase_savePosprocessingFields(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      write(1,*) " Save Posprocessing Fields should be overwritted"
      STOP(1)

   end subroutine CFDSolverBase_savePosprocessingFields

   subroutine CFDSolverBase_evalTimeIteration(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: icode, counter,istep, ppow=1, flag_emac, flag_predic
      character(4) :: timeStep


      counter = 1
      flag_emac = 0
      flag_predic=0

      call nvtxStartRange("Start RK4")

      ! periodic with boundaries
      do istep = 1,this%nstep
         if (istep == this%nsave)write(1,*) '   --| STEP: ', istep
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
         
         if (istep == this%nsave2 .and. (this%globalAnalysis .eq. 1)) then
            call volAvg_EK(this%nelem,this%npoin,connec,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
            call visc_dissipationRate(this%nelem,this%npoin,connec,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
            call maxMach(this%npoin,this%npoin_w,lpoin_w,machno,this%maxmachno)
            call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
            write(1,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
            write(1,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
20      format(6(F16.8,2X))
         end if

         if (flag_real_diff == 1) then
            call adapt_dt_cfl(this%nelem,this%npoin,connec,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,rho(:,2))
            if (istep == this%nsave2)write(1,*) "DT := ",this%dt,"s time := ",this%time,"s"
         else
            call adapt_dt_cfl(this%nelem,this%npoin,connec,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
            if (istep == this%nsave2)write(1,*) "DT := ",this%dt,"s time := ",this%time,"s"
         end if

         call nvtxEndRange

         !
         ! Update the accumulators for averaging
         !
         call nvtxStartRange("Accumulate"//timeStep,istep)
         call favre_average(this%nelem,this%npoin,this%npoin_w,lpoin_w,connec,this%dt,rho,u,pr, &
            mu_fluid,mu_e,mu_sgs,this%acutim,acurho,acupre,acuvel,acuve2,acumueff)
         call nvtxEndRange

         if (istep == this%nsave2) then
            if (this%nboun .ne. 0) then
               do icode = 1,this%numCodes
                  call nvtxStartRange("Surface info")
                  call surfInfo(istep,this%time,this%nelem,this%npoin,this%nboun,icode,connec,bound,point2elem, &
                     bou_codes,bou_norm,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coord, &
                     mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr,Ftau)

                  call nvtxEndRange
                  call flush(888+icode)
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
            call compute_fieldDerivs(this%nelem,this%npoin,connec,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)
            call this%savePosprocessingFields(istep)
            this%nsave = this%nsave+this%nleap
            call nvtxEndRange
         end if

         if(istep==this%nsave2) then
            this%nsave2 = this%nsave2+this%nleap2
            call flush(1)
         end if

         counter = counter+1

      end do
      call nvtxEndRange
   end subroutine CFDSolverBase_evalTimeIteration

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

        ! Eval periodic

        call this%evalPeriodic()

        ! Eval mass 

        call this%evalMass()

        this%time = 0.0_rp

        ! Eval first output

        call this%evalFirstOutput()


        ! Do the time iteration

        call this%evalTimeIteration()

   end subroutine CFDSolverBase_run

end module CFDSolverBase_mod
