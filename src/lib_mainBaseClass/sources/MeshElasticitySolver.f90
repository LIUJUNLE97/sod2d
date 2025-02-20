module MeshElasticitySolver_mod
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
   use mod_InSitu
   use mod_saveFields
   use CFDSolverBase_mod
   use mod_solver_meshElasticity
   implicit none
   private

   type, public, extends(CFDSolverBase) :: MeshElasticitySolver

     real(rp) , public  :: E_young
     real(rp) , public  :: nu_poisson
     
     real(rp), allocatable :: quality_e(:)!,:)
     
     real(rp), allocatable :: coord_input_safe(:,:)
     logical :: is_imposed_displacement = .false.
     real(rp), allocatable :: imposed_displacement(:,:)
     
     real(rp) ,public:: factor_deformation=1.0_rp ! to analytically increase deformation in tests
     
   contains
      procedure, public :: fillBCTypes           => MeshElasticitySolver_fill_BC_Types
      procedure, public :: initializeParameters  => MeshElasticitySolver_initializeParameters
      procedure, public :: run                   => MeshElasticitySolver_run
      procedure, public :: initialBuffer         => imposedDisplacement_elasticitySolverBuffer
      procedure, public :: initializeDefaultSaveFields => CFDSolverBase_initializeDefaultSaveFields_elasticity

      procedure, public :: computeQuality
      procedure, public :: assessElasticityParameters
   end type MeshElasticitySolver
contains
  

  subroutine CFDSolverBase_initializeDefaultSaveFields_elasticity(this)
     class(MeshElasticitySolver), intent(inout) :: this

     call initializeDefaultSaveFieldsElasticity()

  end subroutine CFDSolverBase_initializeDefaultSaveFields_elasticity

  subroutine MeshElasticitySolver_fill_BC_Types(this)
     class(MeshElasticitySolver), intent(inout) :: this

    call this%readJSONBCTypes() 

  end subroutine MeshElasticitySolver_fill_BC_Types

  subroutine MeshElasticitySolver_initializeParameters(this)
     use json_module
     implicit none
     class(MeshElasticitySolver), intent(inout) :: this
     logical :: found, found_aux = .false.
     type(json_file) :: json
     character(len=:) , allocatable :: value

     call json%initialize()
     call json%load_file(json_filename)

     ! get(label,target,is found?, default value)

     call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
     write(this%mesh_h5_file_path,*) value
     call json%get("mesh_h5_file_name",value, found,"channel"); call this%checkFound(found,found_aux)
     write(this%mesh_h5_file_name,*) value
    
     call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
     call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)   
     
     call json%get("E",this%E_young, found,10.0_rp); call this%checkFound(found,found_aux)
     call json%get("nu",this%nu_poisson, found,0.4_rp); call this%checkFound(found,found_aux)  

     call json%get("saveInitialField",this%saveInitialField, found,.true.); call this%checkFound(found,found_aux)
     !call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)
     
     call json%destroy()

     if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'

  end subroutine MeshElasticitySolver_initializeParameters

  subroutine MeshElasticitySolver_run(this)
    !
    implicit none
    class(MeshElasticitySolver), intent(inout) :: this
    !
    real(rp)   :: minQ,maxQ
    !
    ! Init MPI
    call init_mpi()

    ! Init HDF5 interface
    call init_hdf5_interface()

    ! Init Save Fields vars and arrays
    call init_saveFields()

    ! Main simulation parameters
    call this%initializeDefaultParameters()
    call this%initializeParameters()

    call this%optimizeParameters()

    call read_json_saveFields(json_filename)

    ! Open log file
    call this%open_log_file()

    ! read the mesh
    call this%openMesh()

    ! Init hdf5 auxiliar saving arrays
    call init_hdf5_auxiliar_saving_arrays()

    ! Open analysis files
    call this%open_analysis_files

    ! Eval shape Functions
    call this%evalShapeFunctions()

    ! Allocate variables
    call this%allocateVariables()

    ! Setting fields to be saved
    call setFields2Save(rho(:,2),mu_fluid,pr(:,2),E(:,2),eta(:,2),csound,machno,divU,qcrit,Tem(:,2),&
      u(:,:,2),gradRho,curlU,mu_sgs,mu_e,&
      avrho,avpre,avpre2,avmueff,avvel,avve2,avvex,avtw)!,quality_e=quality_e)
    
    ! Eval or load initial conditions
    call this%evalOrLoadInitialConditions()

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

    call this%set_mappedFaces_linkingNodes()   

    ! Eval mass
    print*,'MOVE MASS TO BEFORE CON GRAD'
    call this%evalMass()

    ! Eval first output
    if(this%isFreshStart) call this%evalFirstOutput()
    call this%flush_log_file()

    !call  this%normalFacesToNodes()
    !

    !
    print*,'bouCodes as rho...'
   rho(:,2) = real(bouCodesNodesPar,rp)
    !
    print*,'commented analysis of elasticity parameters'
    !call this%assessElasticityParameters()
    !
    !
    if(.false.) then ! do_curveFromPrescribedDisplacement
      print*,'curving mesh now with analytical boun, to straighten it later and curve it again ... :-)'
      print*,'once it works, test in real testcases'
      !
      call this%initialBuffer()

      if (this%noBoundaries .eqv. .false.) then
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
      end if
      !if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u(:,:,2),u_buffer)
      !
      print*,'Quality before elasticity'
      call this%computeQuality(minQ,maxQ)
      print*,'    minQ: ',minQ
      !
      call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
         numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
         gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
         this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
         bouCodesNodesPar,normalsAtNodes,u_buffer)
      !
      coordPar = coordPar+u(:,:,2)
      print*,'Quality after elasticity'
      call this%computeQuality(minQ,maxQ)
      print*,'    minQ: ',minQ
      !
      call this%saveInstResultsFiles(1)    
    end if
    !
    if(.true.) then ! do_curveInteriorMesh
      u(:,:,2) = 0.0_rp
      call this%saveInstResultsFiles(0)
      print*,'1- Input curved mesh quality'
      call this%computeQuality(minQ,maxQ)
      print*,'       minQ: ',minQ,'    maxQ: ',maxQ
      
      print*,'2- Save input (boundary curved) coordinates'
      call save_input_coordinates(numNodesRankPar,ndime,coordPar,this%coord_input_safe)
      
      print*,'3- Straighten mesh (coordpar)'
      call compute_straight_mesh(numNodesRankPar,ndime,coordPar,numElemsRankPar,nnode,connecParWork,this%coord_input_safe)
      u(:,:,2) = coordPar-this%coord_input_safe ! -> displacement to see in paraview the straight mesh
      call this%saveInstResultsFiles(1)
      u(:,:,2) = 0.0_rp
      print*,'   Quality straight-sided mesh'
      call this%computeQuality(minQ,maxQ)
      print*,'       minQ: ',minQ,'    maxQ: ',maxQ!,' ---> should be all one?!?!?'
      
      print*,'4- Compute displacement of the boundary'
      call compute_displacement_straight_mesh(numNodesRankPar,ndime,coordPar,this%coord_input_safe,bouCodesNodesPar,&
        this%is_imposed_displacement,this%imposed_displacement)
        
      print*,'5- Impose elasticity boundary conditions'
      call this%initialBuffer()
      if (this%noBoundaries .eqv. .false.) then
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
      end if
      
      print*,'6- Call conjGrad_meshElasticity to compute displacements with linear elasticity'
      call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
         numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
         gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
         this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
         bouCodesNodesPar,normalsAtNodes,u_buffer)
      !
      coordPar = coordPar+u(:,:,2)
      print*,'   Quality after elasticity-based curved mesh'
      call this%computeQuality(minQ,maxQ)
      print*,'       minQ: ',minQ,'    maxQ: ',maxQ
      u(:,:,2) = coordPar-this%coord_input_safe ! -> displacement to see in paraview the new curved mesh
      call this%saveInstResultsFiles(2)
    end if
    !

    call this%close_log_file()
    call this%close_analysis_files()

    ! Deallocate the variables
    call this%deallocateVariables()
    
    if(allocated(this%coord_input_safe)) deallocate(this%coord_input_safe)
    if(allocated(this%imposed_displacement)) deallocate(this%imposed_displacement)

    ! End hdf5 auxiliar saving arrays
    call end_hdf5_auxiliar_saving_arrays()

    ! End hdf5 interface
    call end_hdf5_interface()

    ! Finalize InSitu   - bettre done before end_comms
    call end_InSitu()

    ! End comms
    call end_comms()
    call end_comms_bnd()

    ! End MPI
    call end_mpi()
    !
    print*,'IMPLEMENTED THINGS: '
    print*,' - Mesh curving from a "bad" input curved mesh'
    print*,' - Mesh curving after imposing displacement (ALE)'
    print*,' FUTURE: '
    print*,' - Add an option in the json to read the displacement of the ALE vs do mesh curving?'
    print*,' - Optimization?'
    print*,' - Think about doing my own stuff in qualities to allow optimization'
    print*,'TODO solve nasty thing:'
    print*,' - I am doing something nasty... created a link in the meshElasticitySOlver folder to the mod_meshquality'
    print*,' - When Lucas rearranges folders in sod, remove this and do properly'
    
  end subroutine MeshElasticitySolver_run
  !
  !
  !
  subroutine save_input_coordinates(npoin,ndime,coordinates,coord_input_safe)

     implicit none
     !class(MeshElasticitySolver), intent(inout) :: this

     integer(4),             intent(in)    :: npoin, ndime
     real(rp),               intent(in)    :: coordinates(npoin,ndime)
     real(rp), allocatable,  intent(inout) :: coord_input_safe(:,:) !npoin,ndime)
     integer(4) :: inode
     
     !integer(4),             intent(in)    :: npoin, nboun,  bou_codes_nodes(npoin)
     !real(rp), allocatable,  intent(inout) :: coord_input_safe(npoin,ndime)
     !integer(4)                 :: iboun,bcode,ipbou,inode,idime,iBoundNode
    
!      print*,'Getting boundaries bc_type_non_slip_adiabatic, but could use any boundary tag.. '
    
     allocate(coord_input_safe(npoin,ndime))
    
     !$acc parallel loop  
     do inode = 1,npoin
       coord_input_safe(inode,:) = coordinates(inode,:)
!         if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
!            bcode = bou_codes_nodes(inode) ! Boundary element code
!            if (bcode == bc_type_non_slip_adiabatic ) then ! non_slip wall adiabatic -> could be anything
!               x_bou(inode,1) = 0.0_rp
!               x_bou(inode,2) = 0.0_rp
!               x_bou(inode,3) = 0.0_rp
!            end if
!         end if ! This guy
     end do
     !$acc end parallel loop

  end subroutine save_input_coordinates
  !
  !
  !
  subroutine compute_straight_mesh(npoin,ndime,coords,nelem,nnode,connec,coord_input_safe)
!     use mod_maths!, only:
    use mod_arrays, only:  xgp ! high-order nodes
    implicit none

    integer(4),  intent(in)    :: npoin, ndime, nelem, nnode
    real(rp),    intent(inout) :: coords(npoin,ndime)
    integer(4),  intent(in   ) :: connec(nelem,nnode)
    real(rp),    intent(in   ) :: coord_input_safe(npoin,ndime)
    integer(4) :: ielem, inode, idime,i,j,k,id_vertex,theNode,theVertex
    
    real(rp) :: Nx,Ny,Nz,N(nnode,8),xnode_lin(3),minus_plus_one(2)
    
    ! isoparametric mapping: construct linear shape functions
    minus_plus_one(1) = -1.0_rp
    minus_plus_one(2) =  1.0_rp
    id_vertex = 0
    do inode = 1,nnode !->eval linear shape functions on master coordinate "inode" of the high-order elem
      ! compute shape function of linear element
      do k=1,2
        Nz = ( 1.0_rp + (minus_plus_one(k) *xgp(inode,3)) )/2.0_rp
        do j=1,2
          Ny = ( 1.0_rp + (minus_plus_one(j) *xgp(inode,2)) )/2.0_rp
          do i=1,2
            Nx = ( 1.0_rp + (minus_plus_one(i) *xgp(inode,1)) )/2.0_rp
            
            ! eval shape fun on master coordinates xgp
            id_vertex = invAtoIJK(i,j,k)
            N(inode,id_vertex) = Nx*Ny*Nz
          
          end do !i
        end do !j
      end do !k
    end do !inode
    
    !$acc parallel loop  
    do ielem = 1,nelem
      do inode = 1,nnode
        !for each element and node, we now compute the straight position
      
        theNode = connec(ielem,inode)
        coords(theNode,:) = 0.0_rp
        
        do id_vertex=1,8
          theVertex = connec(ielem,id_vertex)
          xnode_lin = coord_input_safe(theVertex,:)
        
          do idime =1,3
            coords(theNode,idime) = coords(theNode,idime) + xnode_lin(idime)*N(inode,id_vertex)
          end do !idime
        end do

!         print*,'theNode:',theNode,'   ielem: ',ielem,' inode: ',inode
!         print*,coords(theNode,:)
!         print*,coord_input_safe(theNode,:)
        
      end do!inode
    end do!ielem
    !$acc end parallel loop
    !
  end subroutine compute_straight_mesh
  !
  !
  !
  subroutine compute_displacement_straight_mesh(npoin,ndime,coords,coord_input_safe,bou_codes_nodes,&
    is_imposed_displacement,imposed_displacement)
    implicit none
    !
    integer(4),  intent(in)    :: npoin, ndime
    real(rp),    intent(in   ) :: coords(npoin,ndime)
    real(rp),    intent(in   ) :: coord_input_safe(npoin,ndime)
    integer(4),  intent(in   ) :: bou_codes_nodes(npoin)
    logical, intent(inout) :: is_imposed_displacement
    real(rp),allocatable, intent(inout) :: imposed_displacement(:,:)
    integer(4) :: inode, bcode
    integer(4) :: iBound,iElem, innerNodeL,bndNodeL,ipbou
    !
    is_imposed_displacement = .true.
    allocate(imposed_displacement(npoin,ndime))
    !$acc parallel loop
    do inode = 1,npoin
      imposed_displacement(inode,:) = 0.0_rp
      if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
         bcode = bou_codes_nodes(inode) ! Boundary element code
         if (bcode == bc_type_non_slip_adiabatic ) then ! non_slip wall adiabatic -> could be anything
           imposed_displacement(inode,:) = coord_input_safe(inode,:)-coords(inode,:)
         end if
      end if
    end do
    !$acc end parallel loop
    
!     rho(:,2)=0.0_rp
!     do iBound = 1,numBoundsRankPar
!       iElem = point2elem(boundPar(iBound,npbou)) ! I use an internal face node to be sure is the correct element
!       innerNodeL = connecParWork(iElem,nnode)         ! internal node
!       !$acc loop vector private(aux)
!       do ipbou = 1,npbou
!          bndNodeL = boundPar(iBound,ipbou) ! node at the boundary
!          imposed_displacement(bndNodeL,:) = coord_input_safe(bndNodeL,:)-coords(bndNodeL,:)
!          rho(bndNodeL,2) = 1.0_rp
!        end do
!      end do
    !
  end subroutine compute_displacement_straight_mesh
  !
  !
  !
  subroutine assessElasticityParameters(this)
    implicit none
    class(MeshElasticitySolver), intent(inout) :: this
   
    integer(4) :: iyoung,ipoisson,num_young,num_poisson,ideformation
    real(rp)   :: fact_young, fact_poisson, ini_young, ini_poisson, end_poisson
    real(rp)   :: minQ,maxQ
    character(len=2) :: str_name
   
    real(rp) :: E_safe, nu_safe
   
    nu_safe = this%nu_poisson
    E_safe  = this%E_young
   
    !"E":10, "nu":0.4,
    num_young   = 3!8
    num_poisson = 3!10
    ini_young   = 0.001_rp !0.0001_rp
    fact_young  = 10.0_rp
    ini_poisson = 0.1_rp ! 0.05_rp
    end_poisson = 0.49_rp
    fact_poisson= (end_poisson-ini_poisson)/num_poisson
   
    do ideformation=1,5
      this%factor_deformation = real(ideformation,rp)!1.0_rp+real(ideformation,rp)
     
      this%factor_deformation = this%factor_deformation / porder
     
      print*,'this%factor_deformation: ',this%factor_deformation
     
      ! HERE WE CAN DO THE MESH MAGIC
      call this%initialBuffer()

      if (this%noBoundaries .eqv. .false.) then
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
      end if
     
      write(str_name, '(I2)') ideformation 
      open(unit=666, file="paramQual_"//TRIM(str_name)//'.txt', status="replace", action="write")
      write(666,*) "      YOUNG            POISSON              MIN_Q             MAX_Q"
      print*, "                      YOUNG            POISSON              MIN_Q             MAX_Q"
      do ipoisson = 0,num_poisson
        do iyoung = 0,num_young
          this%nu_poisson = ini_poisson + fact_poisson*ipoisson
          this%E_young    = ini_young * fact_young**iyoung
          call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
             numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
             gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
             this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
             bouCodesNodesPar,normalsAtNodes,u_buffer)

          call this%computeQuality(minQ,maxQ)

          write(666, *) this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ

          print*,iyoung + (num_young+1)*ipoisson,' -> ',  this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ

          u(:,:,2) = 0.0_rp
        end do
      end do
      close(666)
    end do
   
    this%nu_poisson = nu_safe
    this%E_young = E_safe 
   
  end subroutine assessElasticityParameters
  !
  !
  !
  subroutine computeQuality(this,minQ,maxQ)
   
    use mod_mesh_quality, only: eval_ElemQuality
   
    class(MeshElasticitySolver), intent(inout) :: this
    real(rp),intent(out)   :: minQ, maxQ
    integer(4) :: ielem
    real(8)    :: quality_vec(2)

    integer(4) :: id_quality = 2 ! 1 for aniso, 2 for isso cube ideal 

    minQ = 1.0_rp
    maxQ = 0.0_rp
    do ielem = 1,numElemsRankPar
      !
      !numWorkingNodesRankPar or numNodesRankPar?
      call eval_ElemQuality(size(connecParWork,2),ngaus,numNodesRankPar,numElemsRankPar,&
        ielem,real(coordPar+u(:,:,2),8),connecParWork,real(dNgp,8),real(wgp,8),quality_vec)
      !
      minQ = min(minQ,real(quality_vec(id_quality),rp))
      maxQ = max(maxQ,real(quality_vec(id_quality),rp))
      !
    end do
    !print*,'Min q: ',minQ
    !print*,'Max q: ',maxQ
    !print*,' YoungPoissonMinMax ',this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ
    if(minQ<0) minQ = 0.0_rp
    !     open(unit=666, file="paramQual.txt", status="old", action="write")
    !     !open(unit=iunit, file="output.txt", status="old", action="write", iostat=i)
    !     write(666, *) this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ
    !     close(666)
   
  end subroutine computeQuality
  !
  !
  !
  subroutine imposedDisplacement_elasticitySolverBuffer(this)
    !
    class(MeshElasticitySolver), intent(inout) :: this
    !
    integer(4) :: iNodeL, bcode
    real(rp) :: x0,x1,xpoin,ypoin,zpoin,pertx,perty,pertz,blend_bou
    real(rp) :: factor_sincos
    !
    if(this%is_imposed_displacement) then
      !
      !$acc parallel loop
      do iNodeL=1,numNodesRankPar
        u_buffer(iNodeL,:) = this%imposed_displacement(iNodeL,:)
      end do
      !$acc end parallel loop
      ! 
    else
      !
      factor_sincos = this%factor_deformation
      print*,'factor_sincos: ',factor_sincos,' ---------------------------------------------------------'

      x0=minval(coordPar(:,1)) !it is a cube: assumed!
      x1=maxval(coordPar(:,1)) !it is a cube: assumed!
      
      ! TOY ANALYTICA DISPLACEMENT 
      !$acc parallel loop
      do iNodeL=1,numNodesRankPar
        u_buffer(iNodeL,:) = 0.0_rp
     
        if(coordPar(iNodeL,3)<1e-14) then ! lower cube boundary, z=0
          xpoin = coordPar(iNodeL,1)
          ypoin = coordPar(iNodeL,2)
          zpoin = coordPar(iNodeL,3)
       
          pertx = (xpoin-x0)*(x1-xpoin)/((x1-x0)/2.0_rp)**2.0_rp
          perty = (ypoin-x0)*(x1-ypoin)/((x1-x0)/2.0_rp)**2.0_rp
          pertz = (pertx*perty)
          u_buffer(iNodeL,3) = -pertz/25.0_rp     !quadratic displacement
       
          blend_bou = pertz**2.0_rp
       
          pertx =  sin(xpoin)
          perty =  sin(ypoin)
          pertz = (pertx*perty)
          u_buffer(iNodeL,3) = -(pertz*blend_bou) *factor_sincos     !sinusoidal displacement
        end if
      end do
      !$acc end parallel loop
      ! 
    end if
    !
  end subroutine imposedDisplacement_elasticitySolverBuffer
  !
  !
  !
end module MeshElasticitySolver_mod
