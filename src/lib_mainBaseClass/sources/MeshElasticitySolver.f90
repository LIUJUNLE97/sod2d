module mod_arrays_me
  use mod_constants

  implicit none

  real(rp), allocatable :: quality_e(:)
  real(rp), allocatable :: coord_input_safe(:,:)
  real(rp), allocatable :: imposed_displacement(:,:)
  real(rp), allocatable :: metric(:,:,:)
  real(rp), allocatable :: N_lin(:,:)

end module mod_arrays_me

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
   use mod_arrays_me
   implicit none
   private
   
   type, public, extends(CFDSolverBase) :: MeshElasticitySolver
          
     real(rp) , public  :: E_young    = 10.0_rp
     real(rp) , public  :: nu_poisson = 0.4_rp                    
     
     logical :: is_imposed_displacement = .false.               
     logical :: saveNewCoords = .false.       
     
     real(rp), public:: factor_deformation=1.0_rp ! to analytically increase deformation in stress tests
     
   contains
      procedure, public :: fillBCTypes           => MeshElasticitySolver_fill_BC_Types
      procedure, public :: initializeParameters  => MeshElasticitySolver_initializeParameters
      procedure, public :: run                   => MeshElasticitySolver_run
      procedure, public :: initialBuffer         => imposedDisplacement_elasticitySolverBuffer
      procedure, public :: initializeDefaultSaveFields => CFDSolverBase_initializeDefaultSaveFields_elasticity
      
      procedure, public :: computeQuality
      procedure, public :: assessElasticityParameters_forDifferentDefs
      procedure, public :: assessBestElasticityParameters
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
    
     call json%get("save_logFile_step",this%save_logFile_step, found, 10); call this%checkFound(found,found_aux)

     call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
     call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)   
     
     call json%get("E",this%E_young, found,10.0_rp); call this%checkFound(found,found_aux)
     call json%get("nu",this%nu_poisson, found,0.4_rp); call this%checkFound(found,found_aux)  

     call json%get("saveInitialField",this%saveInitialField, found,.true.); call this%checkFound(found,found_aux)
     !call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)

     call json%get("saveNewCoords",this%saveNewCoords, found,.false.); call this%checkFound(found,found_aux)
     
     call this%readJSONMeshElasticityTypes()

     call json%destroy()

     if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'

  end subroutine MeshElasticitySolver_initializeParameters

  subroutine MeshElasticitySolver_run(this)
    !
    implicit none
    class(MeshElasticitySolver), intent(inout) :: this
    !
    real(rp)   :: minQ,maxQ,factor_backtrack
    real(rp)   :: minQTot
    integer(4):: numInv,numLow,numBacktracks
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
    call this%evalMass()

    ! Eval first output
    if(this%isFreshStart) call this%evalFirstOutput()
    call this%flush_log_file()

    !call  this%normalFacesToNodes()
    !
    !
    !print*,'Exporting bouCodes as rho...'
    !$acc kernels
    rho(:,2) = real(bouCodesNodesPar,rp)
    !$acc end kernels
    !
    
    !this%elasticity_problemType = elasticity_fromALE
    
    !
    !print*,'commented analysis of elasticity parameters'
    !call this%assessElasticityParameters_forDifferentDefs()
    !
    !
    if( elasticity_problemType == elasticity_fromALE ) then ! do_curveFromPrescribedDisplacement
      if(mpi_rank.eq.0) write(*,*) '  --| Imposing a prescribed displacement in dirichlet boundaries...'
      !
      call this%initialBuffer()

      if (this%noBoundaries .eqv. .false.) then
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
      end if
      !if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u(:,:,2),u_buffer)
      !
      if(mpi_rank.eq.0) write(*,*) '  --| Quality before elasticity'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      !call MPI_Reduce(minQ,minQTot,1,mpi_datatype_real,MPI_MIN,0,app_comm,mpi_err)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ
      !
      if(mpi_rank.eq.0) write(*,*) '  --| conjGrad_meshElasticity...'
      call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
         numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
         gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
         this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
         bouCodesNodesPar,normalsAtNodes,u_buffer)
      !
      !$acc kernels
      coordPar = coordPar+u(:,:,2)
      !$acc end kernels
      
      if(mpi_rank.eq.0) write(*,*) '  --| Quality after elasticity'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      !call MPI_Reduce(minQ,minQTot,1,mpi_datatype_real,MPI_MIN,0,app_comm,mpi_err)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ
      !
      call this%saveInstResultsFiles(1)    
    end if
    !
    if( elasticity_problemType == elasticity_fromBouCurving) then ! do_curveInteriorMesh
      if(mpi_rank.eq.0) write(*,*) '  --| Elasticity for boundary mesh curving...'
      
      !$acc kernels
      u(:,:,2) = 0.0_rp
      !$acc end kernels
      if(mpi_rank.eq.0) write(*,*) '  --| Input curved mesh quality'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ
      call this%saveInstResultsFiles(0)
      
      call save_input_coordinates(numNodesRankPar,ndime,coordPar,coord_input_safe)
      
      if(mpi_rank.eq.0) write(*,*) '  --| Straighten mesh (coordpar)'
      call compute_straight_mesh(numNodesRankPar,ndime,coordPar,numElemsRankPar,nnode,connecParWork,coord_input_safe)
      if(mpi_rank.eq.0) write(*,*) '  --| Quality straight-sided mesh'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ

      !$acc kernels
      u(:,:,2) = coordPar-coord_input_safe ! -> displacement to see in paraview the straight mesh
      !$acc end kernels
      call this%saveInstResultsFiles(1)
      
      !print*,'- Compute displacement of the boundary'
      call compute_displacement_straight_mesh(numNodesRankPar,ndime,coordPar,coord_input_safe,bouCodesNodesPar,&
        this%is_imposed_displacement,imposed_displacement)
      !
      call this%initialBuffer()
      if (this%noBoundaries .eqv. .false.) then
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
      end if
      !
      !if(mpi_rank.eq.0) write(*,*) '  --| Assess Best Elasticity Parameters'
      !call this%assessBestElasticityParameters()
      !
      if(mpi_rank.eq.0) write(*,*) '  --| conjGrad_meshElasticity...'
      call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
         numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
         gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
         this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
         bouCodesNodesPar,normalsAtNodes,u_buffer)
      !
      !$acc kernels
      coordPar = coordPar+u(:,:,2)
      !$acc end kernels
      if(mpi_rank.eq.0) write(*,*) '  --| Quality after elasticity-based curved mesh'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ

      !$acc kernels
      u(:,:,2) = coordPar-coord_input_safe ! -> displacement to see in paraview the new curved mesh
      !$acc end kernels
      call this%saveInstResultsFiles(2)
      
      if(this%saveNewCoords) call save_coordinates_hdf5(this%meshFile_h5_full_name)
    end if
    !
    if( elasticity_problemType == elasticity_fromMetric ) then ! do_curveFromMetric
      print*,'STOP!!! elasticity_fromMetric STILL NOT IN PRODUCTION'
      stop
      
      if(mpi_rank.eq.0) print*,'use metric in the domain to drive elasticity'
      this%is_imposed_displacement = .true.
      allocate(imposed_displacement(numNodesRankPar,ndime))
      !$acc enter data create(imposed_displacement(:,:))
      
      !$acc kernels
      imposed_displacement = 0.0_rp
      !$acc end kernels
      !
      call this%initialBuffer()

      if (this%noBoundaries .eqv. .false.) then
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
      end if
      !if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u(:,:,2),u_buffer)
      !
      if(mpi_rank.eq.0) print*,'Quality before elasticity'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      !call MPI_Reduce(real(minQ,8),minQTot,1,mpi_datatype_real8,MPI_MIN,0,app_comm,mpi_err)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ
      !
      call computeAnalyticalMetric(numNodesRankPar,ndime,coordPar,metric)
      
      call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
         numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
         gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
         this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
         bouCodesNodesPar,normalsAtNodes,u_buffer,metric)
      !
      call save_input_coordinates(numNodesRankPar,ndime,coordPar,coord_input_safe)
      
      !u(:,:,2) = -u(:,:,2) ! did we solve the problem the other way around?
      if(mpi_rank.eq.0) print*,'Enlarging displacement artificially to emulate optimization'
      u(:,:,2) = u(:,:,2) * 5
      coordPar = coord_input_safe + u(:,:,2)
      if(mpi_rank.eq.0) print*,'Quality after elasticity'
      call this%computeQuality(minQ,maxQ,numInv,numLow)
      !call MPI_Reduce(real(minQ,8),minQTot,1,mpi_datatype_real8,MPI_MIN,0,app_comm,mpi_err)
      if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ
      
      coordPar = coord_input_safe 
      call this%saveInstResultsFiles(1)
      
      factor_backtrack = 0.9_rp
      numBacktracks = 0
      do while(minQ.le.0) 
        numBacktracks = numBacktracks+1
        u(:,:,2) = u(:,:,2)*factor_backtrack
        coordPar = coord_input_safe + u(:,:,2)
        if(mpi_rank.eq.0) print*,'  Quality after backtrack: ',numBacktracks
        call this%computeQuality(minQ,maxQ,numInv,numLow)
        !call MPI_Reduce(real(minQ,8),minQTot,1,mpi_datatype_real8,MPI_MIN,0,app_comm,mpi_err)
        if(mpi_rank.eq.0) write(*,*) '  --| minQ: ',minQ
        if(mpi_rank.eq.0) print*,'     minMaxUx: ',minval(u(:,1,2)),' / ',maxval(u(:,1,2))
        if(mpi_rank.eq.0) print*,'     minMaxUy: ',minval(u(:,2,2)),' / ',maxval(u(:,2,2))
        if(mpi_rank.eq.0) print*,'     minMaxUz: ',minval(u(:,3,2)),' / ',maxval(u(:,3,2))
      end do
      !
      coordPar = coord_input_safe
      call this%saveInstResultsFiles(2)
      
    end if
    !
    call this%close_log_file()
    call this%close_analysis_files()

    ! Deallocate the variables
    call this%deallocateVariables()
    
    if(allocated(coord_input_safe)) deallocate(coord_input_safe)
    if(allocated(imposed_displacement)) deallocate(imposed_displacement)

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
    if(mpi_rank.eq.0) then
      print*,'IMPLEMENTED THINGS: '
      print*,' - Mesh curving from a "bad" input curved mesh'
      print*,' - Mesh curving after imposing displacement (ALE)'
      print*,' - r-adaptation though an analytical metric (isotropic)'
      print*,' FUTURE: '
      print*,' - Optimization?'
      print*,'TODO solve nasty thing:'
      print*,' - I am doing something nasty... created a link in the meshElasticitySOlver folder to the mod_meshquality'
      print*,' - When Lucas rearranges folders in sod, remove this and do properly'
    end if
    
  end subroutine MeshElasticitySolver_run
  !
  !
  !
  subroutine compute_straight_mesh(npoin,ndime,coords,nelem,nnode,connec,coord_input_safe)
      use mod_arrays, only:  xgp ! high-order nodes
      implicit none
  
      integer(4),  intent(in)    :: npoin, ndime, nelem, nnode
      real(rp),    intent(inout) :: coords(npoin,ndime)
      integer(4),  intent(in   ) :: connec(nelem,nnode)
      real(rp),    intent(in   ) :: coord_input_safe(npoin,ndime)
      integer(4) :: ielem, inode, idime,i,j,k,id_vertex,theNode,theVertex
      
      real(rp) :: Nx,Ny,Nz,xnode_lin(3),minus_plus_one(2),coordLocal(ndime)
      
      allocate(N_lin(nnode,nnode))
      !$acc enter data create(N_lin(:,:))
  
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
              N_lin(inode,id_vertex) = Nx*Ny*Nz
            end do !i
          end do !j
        end do !k
      end do !inode
      
      !$acc update device(N_lin(:,:))
  
      !$acc parallel loop private(coordLocal)
      do ielem = 1,nelem
        !$acc loop seq
        do inode = 1,nnode ! avoid the vertices of the hex
          !For each element and node, we now compute the straight position
          
          !Contributions of the straight (linear) hexahedral element (loop on linear vertices):
          coordLocal(:) = 0.0_rp
          !$acc loop seq
          do k=1,2
            !$acc loop seq
            do j=1,2
              !$acc loop seq
              do i=1,2
                id_vertex = invAtoIJK(i,j,k)
                theVertex = connec(ielem,id_vertex)
                xnode_lin = coord_input_safe(theVertex,:)

                !$acc loop seq
                do idime =1,3
                  coordLocal(idime) = coordLocal(idime) + xnode_lin(idime)*N_lin(inode,id_vertex)
                end do !idime
              end do
            end do
          end do  
  
          theNode = connec(ielem,inode)
          !$acc loop seq        
          do idime =1,ndime
            !$acc atomic write
            coords(theNode, idime) = coordLocal(idime)
            !$acc end atomic 
          end do
         
        end do!inode
      end do!ielem
      !$acc end parallel loop
      !
    end subroutine compute_straight_mesh
  !
  !
  !
  subroutine computeQuality(this,minQ,maxQ,countInvalid,countLowQ)
    !
    class(MeshElasticitySolver), intent(inout) :: this
    real(rp),   intent(out) :: minQ, maxQ
    integer(4), intent(out) :: countInvalid,countLowQ 
    integer(4) :: ielem
    real(rp)   :: coordElem(size(connecParWork,2),ndime)
    real(rp)   :: quality(numElemsRankPar),distortion(numElemsRankPar)
    ! variables from eval_ElemQuality_simple subroutine
    integer(4) :: igaus, idime, jdime
    real(rp) :: eta_g, eta_elem, volume, modulus, gpvolIdeal, detIdeal, detIdealCube
    real(rp) :: elemJ(ndime, ndime), idealCubeJ(ndime, ndime)
    ! variables from shape quality subroutine
    real(rp) :: S(ndime,ndime), StS(ndime,ndime), sigma, Sf, detS
    real(rp),parameter :: d=3.0_rp
    !
    !$acc update host(coordPar(:,:))

    idealCubeJ = 0.0_rp
    do idime = 1, ndime
        idealCubeJ(idime, idime) = 1.0_rp
    end do
    detIdealCube = 1.0_rp
    !$acc enter data copyin(idealCubeJ,detIdealCube)

    !$acc parallel loop gang  private(coordElem,eta_elem,volume)
    do ielem = 1,numElemsRankPar
      !
      coordElem = coordPar(connecParWork(ielem,:),:)
      !
      eta_elem = 0.0_rp
      volume   = 0.0_rp
      !$acc loop vector private(elemJ,S,detS,sigma,StS,Sf,eta_g,gpvolIdeal) reduction(+:eta_elem, volume)
      do igaus = 1, ngaus
          elemJ(:, :) = 0.0_rp
          !$acc loop seq
          do idime = 1, ndime
              !$acc loop seq
              do jdime = 1, ndime
                  elemJ(jdime, idime) = dot_product(dNgp(idime, :, igaus), coordElem(:, jdime))
              end do
          end do
  
          ! SHAPE MEASURE
          S     = matmul(elemJ, idealCubeJ)
          detS  = S(1,1)*(S(2,2)*S(3,3) - S(2,3)*S(3,2)) &
                - S(1,2)*(S(2,1)*S(3,3) - S(2,3)*S(3,1)) &
                + S(1,3)*(S(2,1)*S(3,2) - S(2,2)*S(3,1))
          sigma = (detS + abs(detS))/2
          StS   = matmul(transpose(S), S)
          Sf    = StS(1,1) + StS(2,2) + StS(3,3)
          eta_g = Sf/(d*sigma**(2.0d0/d)) 
          gpvolIdeal = detIdealCube * wgp(igaus)

          eta_elem   = eta_elem + eta_g * eta_g * gpvolIdeal
          volume     = volume + 1.0_rp * 1.0_rp * gpvolIdeal
      end do
  
      distortion(ielem) = sqrt(eta_elem) / sqrt(volume)
      quality(ielem)    = 1.0_rp / distortion(ielem)
      modulus = modulo(quality(ielem), 1.0_rp)
      if (int(modulus) .ne. 0) then
          quality(ielem) = -1.0_rp
          distortion(ielem) = 1.0e10_rp
      end if
      !
      mu_e(ielem,:) = quality(ielem) ! for postprocessing
    end do
    !$acc end parallel loop

    !$acc update device(mu_e,quality,distortion)
    
    countInvalid = 0
    countLowQ    = 0
    minQ         =  1.0d30
    maxQ         = -1.0d30
    !$acc parallel loop reduction(+:countInvalid,countLowQ) reduction(min:minQ) reduction(max:maxQ)
    do ielem = 1, numElemsRankPar
        if (quality(ielem)    < 0.0d0 ) countInvalid = countInvalid + 1
        if (distortion(ielem) > 1.0d6 ) countLowQ    = countLowQ + 1
        minQ = min(minQ, quality(ielem))
        maxQ = max(maxQ, quality(ielem))
    end do
    !
    call MPI_Reduce(countInvalid,countInvalid,1,mpi_datatype_int4,MPI_SUM,0,app_comm,mpi_err)
    call MPI_Reduce(countLowQ   ,countLowQ   ,1,mpi_datatype_int4,MPI_SUM,0,app_comm,mpi_err)
    call MPI_Reduce(minQ        ,minQ        ,1,mpi_datatype_real,MPI_MIN,0,app_comm,mpi_err)
    call MPI_Reduce(maxQ        ,maxQ        ,1,mpi_datatype_real,MPI_MAX,0,app_comm,mpi_err)
    !
    if(minQ<0) minQ = 0.0_rp        
    !
  end subroutine computeQuality
  !
  !
  !
  subroutine computeQualitykk(this,minQ,maxQ,countInvalid,countLowQ)
   
    !use mod_quality, only: eval_ElemQuality_simple
   
    class(MeshElasticitySolver), intent(inout) :: this
    real(rp),   intent(out) :: minQ, maxQ
    integer(4), intent(out) :: countInvalid,countLowQ 
    integer(4) :: ielem,mnode
    real(rp)   :: coordElem(size(connecParWork,2),ndime)
    real(rp)   :: quality(numElemsRankPar),distortion(numElemsRankPar)
    !
    mnode = size(connecParWork,2)

    !$acc update host(coordPar(:,:))

    print*,"numElemsRankPar: ",numElemsRankPar,' switched to test how it works'

    ! !$acc data copyin(dNgp, wgp, coordPar, connecParWork)  ! Copy dNgp, wgp, coordPar, and connecParWork to the device
    !$acc parallel loop private(coordElem)
    do ielem = 1,1e4!numElemsRankPar
      !
      coordElem = coordPar(connecParWork(ielem,:),:)
      !call eval_ElemQuality_simple(mnode,ngaus,coordElem,dNgp,wgp,quality(ielem),distortion(ielem))
      !
      mu_e(ielem,:) = quality(ielem)
    end do
    !$acc end parallel loop

    !$acc update device(mu_e(:,:))
    
    countInvalid = 0
    !$acc parallel loop reduction(+:countInvalid)
    do ielem = 1,numElemsRankPar
        if (quality(ielem) < 0.0d0) countInvalid = countInvalid + 1
    end do
    countLowQ = 0
    !$acc parallel loop reduction(+:countLowQ)
    do ielem = 1,numElemsRankPar
        if (distortion(ielem) > 1.0d10) countLowQ = countLowQ + 1
    end do
    minQ = 1.0d30
    !$acc parallel loop reduction(min:minQ)
    do ielem = 1,numElemsRankPar
        minQ = min(minQ, quality(ielem))
    end do
    maxQ = -1.0d30
    !$acc parallel loop reduction(max:maxQ)
    do ielem = 1,numElemsRankPar
        maxQ = max(maxQ, quality(ielem))
    end do
    !
    call MPI_Reduce(countInvalid,countInvalid,1,mpi_datatype_int4,MPI_SUM,0,app_comm,mpi_err)
    call MPI_Reduce(countLowQ   ,countLowQ   ,1,mpi_datatype_int4,MPI_SUM,0,app_comm,mpi_err)
    call MPI_Reduce(minQ        ,minQ        ,1,mpi_datatype_real,MPI_MIN,0,app_comm,mpi_err)
    call MPI_Reduce(maxQ        ,maxQ        ,1,mpi_datatype_real,MPI_MAX,0,app_comm,mpi_err)
    !
    if(minQ<0) minQ = 0.0_rp        
            !
    print*,'qelems: ',quality(1:10)
    print*,minQ,maxQ,countLowQ,countInvalid
    !
  end subroutine computeQualitykk
  !
  !
  !
  subroutine computeQuality_old(this,minQ,maxQ,countInvalid,countLowQ)
   
    use mod_quality, only: eval_ElemQuality
   
    class(MeshElasticitySolver), intent(inout) :: this
    real(rp),   intent(out) :: minQ, maxQ
    integer(4), intent(out) :: countInvalid,countLowQ 
    integer(4) :: ielem
    real(8)    :: quality_vec(2)

    integer(4) :: id_quality = 1 ! 1 for aniso, 2 for isso cube ideal 
    
    real(8)   :: q_gauss(ngaus)
    !
    countInvalid = 0
    countLowQ = 0
    minQ = 1.0_rp
    maxQ = 0.0_rp

    !$acc update host(coordPar(:,:))

    do ielem = 1,numElemsRankPar
      !
      call eval_ElemQuality(size(connecParWork,2),ngaus,numNodesRankPar,numElemsRankPar,&
        ielem,real(coordPar,8),connecParWork,real(dNgp,8),real(wgp,8),&
        quality_vec,q_gauss)
      !
      minQ = min(minQ,real(quality_vec(id_quality),rp))
      maxQ = max(maxQ,real(quality_vec(id_quality),rp))
      !
      mu_e(ielem,:) = real(q_gauss,4)
      !
      if(real(quality_vec(id_quality),rp)<0) then
        countInvalid = countInvalid+1
      end if
      if(minval(q_gauss)<0) then
        countLowQ = countLowQ+1
      end if
    end do

    !$acc update device(mu_e(:,:))

    !
    if(minQ<0) minQ = 0.0_rp
    !print*,'   Num invalid: ',countInvalid,'     Num lowQ: ',countLowQ
    !
  end subroutine computeQuality_old
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
    !$acc enter data create(imposed_displacement(:,:))
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
    
    ! I COULD DO THE FOLLOWING TO NOT CHECK ANY BC IN PARTICULAR, AND JUST GET ALL BOUNDARIES
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
     !$acc enter data create(coord_input_safe(:,:))
    
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
     !$acc update host(coord_input_safe(:,:))


  end subroutine save_input_coordinates
  !
  !
  !
  subroutine computeAnalyticalMetric(npoin,ndime,coordinates,metric)

     implicit none
     !class(MeshElasticitySolver), intent(inout) :: this

     integer(4),             intent(in)    :: npoin, ndime
     real(rp),               intent(in)    :: coordinates(npoin,ndime)
     real(rp), allocatable,  intent(inout) :: metric(:,:,:)
     
     integer(4) :: inode
     real(rp) :: t_aux,hdesired,hz_maz,hz_min
     real(rp) :: eigenx,eigeny,eigenz
     real(rp) :: h_ref, coords_aux(3),V(3,3),Mnode(3,3)
     
     
     h_ref = 6.28319_rp/8.0_rp ! for the cube mesh...
     
     allocate(metric(ndime,ndime,npoin))
      !$acc enter data create(metric(:,:,:))
    
     !$acc parallel loop  
     do inode = 1,npoin
       !!! ---- boundary layer ------
!        t_aux = coordinates(inode,3)/6.28319_rp
!        hz_maz = 2.0_rp !riemanian
!        hz_min = 0.1_rp !riemanian
!        !t_aux = t_aux**2 ! quadratic transition from min (aniso) to max (iso)
!        ! elem size in between [0.1 1] according to taux (height, normalized)
!        hdesired = hz_maz*t_aux  + hz_min*(1-t_aux)
!
!        hdesired = hdesired/h_ref
!
!        eigenx = 1.0_rp
!        eigeny = 1.0_rp
!        !eigenz = 1.0_rp/hdesired**2
!        eigenz = hdesired!**2
!
!        metric(:,:,inode) = 0.0_rp
!        metric(1,1,inode) = eigenx
!        metric(2,2,inode) = eigeny
!        metric(3,3,inode) = eigenz
       
       !!! ---- radial sizing ------
!        hz_maz = 2.0_rp !riemanian
!        hz_min = 0.1_rp !riemanian
!        t_aux = sqrt(sum((coordinates(inode,:)-3.14)**2)) ! radius on center of cube
!        t_aux = t_aux/3.14 ! t\in[0,1]
!        t_aux = t_aux*t_aux
!        hdesired = hz_maz*t_aux  + hz_min*(1-t_aux)
!
!        !eigenx = 1.0_rp/hdesired**2
!        eigenx = hdesired!**2
!        eigeny = eigenx
!        eigenz = eigenx
!
!        metric(:,:,inode) = 0.0_rp
!        metric(1,1,inode) = eigenx
!        metric(2,2,inode) = eigeny
!        metric(3,3,inode) = eigenz
       
       !!! ---- boundary layer Diagonal------
       coords_aux = 2*coords_aux/6.28319_rp
       coords_aux = coords_aux-1
       t_aux = abs(coords_aux(1) + coords_aux(2))/2!sqrt(2.0_rp)
       hz_maz = 2.0_rp !riemanian
       hz_min = 0.1_rp !riemanian
       t_aux = t_aux**2 ! quadratic transition from min (aniso) to max (iso)
       hdesired = hz_maz*t_aux  + hz_min*(1-t_aux)
       
       hdesired = hdesired/h_ref
       
       eigenx = 1.0_rp
       eigeny = 1.0_rp
       eigenz = hdesired!**2
       
       Mnode(:,:) = 0.0_rp
       Mnode(1,1) = eigenx
       Mnode(2,2) = eigeny
       Mnode(3,3) = eigenz
       
       V(1,:)=(/1.0_rp,-1.0_rp,0.0_rp/)/sqrt(2.0_rp)
       V(2,:)=(/0.0_rp,0.0_rp,1.0_rp/)
       V(3,:)=(/1.0_rp,1.0_rp,0.0_rp/)/sqrt(2.0_rp)

       Mnode = matmul(Mnode,transpose(V))
       metric(:,:,inode) = matmul(V,Mnode)
       
       Mnode(:,:) = 0.0_rp
       Mnode(1,1) = eigenz
       Mnode(2,2) = eigenz
       Mnode(3,3) = eigenz
       metric(:,:,inode) = Mnode ! I am using V to change gradU of basis
       
!        print*;'check that I am imposing the sizing correctly... imagina que lestic cagant aqui...'
!        stop 1
     end do
     !$acc end parallel loop

  end subroutine computeAnalyticalMetric
  !
  !
  !
  subroutine assessElasticityParameters_forDifferentDefs(this)
    implicit none
    class(MeshElasticitySolver), intent(inout) :: this
   
    integer(4) :: iyoung,ipoisson,num_young,num_poisson,ideformation,numInv,numLow
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

          call this%computeQuality(minQ,maxQ,numInv,numLow)

          write(666, *) this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ

          print*,iyoung + (num_young+1)*ipoisson,' -> ',  this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ

          u(:,:,2) = 0.0_rp
        end do
      end do
      close(666)
    end do
   
    this%nu_poisson = nu_safe
    this%E_young = E_safe 
   
  end subroutine assessElasticityParameters_forDifferentDefs
  !
  !
  !
  subroutine assessBestElasticityParameters(this)
    implicit none
    class(MeshElasticitySolver), intent(inout) :: this
   
    integer(4) :: iyoung,ipoisson,num_young,num_poisson,numInv,numLow
    real(rp)   :: fact_young, fact_poisson, ini_young, ini_poisson, end_poisson
    real(rp)   :: minQ,maxQ
    character(len=2) :: str_name
   
    real(rp) :: E_safe, nu_safe, E_best, nu_best, q_best = 0.0_rp
   
    nu_safe = this%nu_poisson
    E_safe  = this%E_young
    nu_best = nu_safe 
    E_best  = E_safe  
   
    !"E":10, "nu":0.4,
    num_young   = 3!8
    num_poisson = 3!10
    ini_young   = 0.01_rp !0.0001_rp
    fact_young  = 10.0_rp ! 10.0_rp
    ini_poisson = 0.139_rp!0.1_rp ! 0.05_rp
    end_poisson = 0.49_rp
    fact_poisson= (end_poisson-ini_poisson)/num_poisson
   
    ! HERE WE CAN DO THE MESH MAGIC
!     call this%initialBuffer()
!
!     if (this%noBoundaries .eqv. .false.) then
!        call temporary_bc_routine_dirichlet_prim_meshElasticity(&
!          numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)
!     end if
    !
    if(mpi_rank.eq.0) then
      print*,'Range Poisson: [',ini_poisson,',',ini_poisson + fact_poisson*num_poisson,'] -> num: ', num_poisson+1
      print*,'Range Young: [',ini_young,',',ini_young * fact_young**num_young,']  -> num: ',num_young+1
      str_name = "";
      open(unit=666, file="paramQual"//TRIM(str_name)//'.txt', status="replace", action="write")
      write(666,*) "      YOUNG            POISSON              MIN_Q             MAX_Q"
      print*, "                      YOUNG            POISSON              MIN_Q             MAX_Q"
    end if
    do ipoisson = 0,num_poisson
      do iyoung = 0,num_young
        
        this%nu_poisson = ini_poisson + fact_poisson*ipoisson
        this%E_young    = ini_young * fact_young**iyoung
        
        call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
           numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
           gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
           this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final
           bouCodesNodesPar,normalsAtNodes,u_buffer)
        !
        coordPar = coordPar+u(:,:,2)
        call this%computeQuality(minQ,maxQ,numInv,numLow)
        
        if(minQ>q_best) then
          q_best = minQ
          nu_best = this%nu_poisson
          E_best  = this%E_young
        end if

        if(mpi_rank.eq.0) then
          write(666, *) this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ
          print*,iyoung + (num_young+1)*ipoisson,' -> ',  this%E_young,' ',this%nu_poisson ,' ',minQ,' ',maxQ,&
          ' ',numInv ,' ',numLow
        end if

        coordPar = coordPar-u(:,:,2)
        u(:,:,2) = 0.0_rp
      end do
    end do
    if(mpi_rank.eq.0) then
      close(666)
    end if
     
    print*,'     -> Best parameters: ',nu_best,'  ',E_best
    this%nu_poisson = nu_best
    this%E_young = E_best
   
  end subroutine assessBestElasticityParameters
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
        u_buffer(iNodeL,:) = imposed_displacement(iNodeL,:)
      end do
      !$acc end parallel loop
      ! 
    else
      !
      if(mpi_rank.eq.0)  print*,'IMPOSING AN ANALYTICAL DISPLACEMENT (SHOULD BE READ FROM SOMEWHERE INSTEAD)'
      
      factor_sincos = this%factor_deformation
      !print*,'factor_sincos: ',factor_sincos,' ---------------------------------------------------------'

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
