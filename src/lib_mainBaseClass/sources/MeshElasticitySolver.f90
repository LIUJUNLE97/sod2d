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

   contains
      procedure, public :: fillBCTypes           => MeshElasticitySolver_fill_BC_Types
      procedure, public :: initializeParameters  => MeshElasticitySolver_initializeParameters
      procedure, public :: run                   => MeshElasticitySolver_run
      procedure, public :: initialBuffer         => imposedDisplacement_elasticitySolverBuffer
   end type MeshElasticitySolver
contains

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

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'

   end subroutine MeshElasticitySolver_initializeParameters

   subroutine MeshElasticitySolver_run(this)
      implicit none
      class(MeshElasticitySolver), intent(inout) :: this

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
      if(flag_use_species .eqv. .true.) then
         call setFields2Save(rho(:,2),mu_fluid,pr(:,2),E(:,2),eta(:,2),csound,machno,divU,qcrit,Tem(:,2),&
                             u(:,:,2),gradRho,curlU,mu_sgs,mu_e,&
                             avrho,avpre,avpre2,avmueff,avvel,avve2,avvex,avtw,Yk(:,:,2),mu_e_Yk)
      else
         call setFields2Save(rho(:,2),mu_fluid,pr(:,2),E(:,2),eta(:,2),csound,machno,divU,qcrit,Tem(:,2),&
                             u(:,:,2),gradRho,curlU,mu_sgs,mu_e,&
                             avrho,avpre,avpre2,avmueff,avvel,avve2,avvex,avtw)
      end if

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
      print*,'this%isFreshStart: ',this%isFreshStart
      if(this%isFreshStart) call this%evalFirstOutput()
      call this%flush_log_file()

      call  this%normalFacesToNodes()

      call this%initialBuffer()

      print*,'this%initialBuffer():   min         /        max'
      print*,'    ',minval(u_buffer(:,1)),' / ',maxval(u_buffer(:,1))
      print*,'    ',minval(u_buffer(:,2)),' / ',maxval(u_buffer(:,2))
      print*,'    ',minval(u_buffer(:,3)),' / ',maxval(u_buffer(:,3))
      
      print*,'noBoundaries: ', this%noBoundaries
      if (this%noBoundaries .eqv. .false.) then
        !if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet_incomp(u(:,:,2))
!          call temporary_bc_routine_dirichlet_prim_meshElasticity(&
!            numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,2),u_buffer)
         print*,' i changed dirichlet bc from u(:,:,2) to u(:,:,1)'
         call temporary_bc_routine_dirichlet_prim_meshElasticity(&
           numNodesRankPar,numBoundsRankPar,bouCodesNodesPar,lbnodesPar,normalsAtNodes,u(:,:,1),u_buffer)

         print*,'u(:,:,1):   min         /        max'
         print*,'    ',minval(u(:,1,1)),' / ',maxval(u(:,1,1))
         print*,'    ',minval(u(:,2,1)),' / ',maxval(u(:,2,1))
         print*,'    ',minval(u(:,3,1)),' / ',maxval(u(:,3,1))
         print*,'u(:,:,2):   min         /        max'
         print*,'    ',minval(u(:,1,2)),' / ',maxval(u(:,1,2))
         print*,'    ',minval(u(:,2,2)),' / ',maxval(u(:,2,2))
         print*,'    ',minval(u(:,3,2)),' / ',maxval(u(:,3,2))
      end if
      !if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u(:,:,2),u_buffer)
      
      call conjGrad_meshElasticity(1,this%save_logFile_next,this%noBoundaries,numElemsRankPar,numNodesRankPar,&
         numWorkingNodesRankPar,numBoundsRankPar,connecParWork,workingNodesPar,invAtoIJK,&
         gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,&
         this%nu_poisson,this%E_young,u(:,:,1),u(:,:,2), &!u1 condicion inicial u2 terme font y solucio final 
         bouCodesNodesPar,normalsAtNodes,u_buffer) 

      call this%saveInstResultsFiles(1)    

      call this%close_log_file()
      call this%close_analysis_files()

      ! Deallocate the variables
      call this%deallocateVariables()

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

   end subroutine MeshElasticitySolver_run


   subroutine imposedDisplacement_elasticitySolverBuffer(this)
      class(MeshElasticitySolver), intent(inout) :: this
      integer(4) :: iNodeL, bcode
      
      real(rp) :: x0,x1,xpoin,ypoin,zpoin,pertx,perty,pertz

      !$acc parallel loop
!       do iNodeL = 1,numNodesRankPar
!             u_buffer(iNodeL,1) = this%vo*cos(this%aoa*v_pi/180.0_rp)
!             u_buffer(iNodeL,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)
!             u_buffer(iNodeL,3) = 0.0_rp
!       end do
      
      x0=minval(coordPar(:,1))
      x1=maxval(coordPar(:,1))
      do iNodeL=1,numNodesRankPar
        u_buffer(iNodeL,:) = 0.0_rp
        
        if(coordPar(iNodeL,3)<1e-14) then
          xpoin = coordPar(iNodeL,1)
          ypoin = coordPar(iNodeL,2)
          zpoin = coordPar(iNodeL,3)
          
          pertx = (xpoin-x0)*(x1-xpoin)
          perty = (ypoin-x0)*(x1-ypoin)
          pertz = (pertx*perty)/50.0_rp
          u_buffer(iNodeL,3) = -pertz*2
          
          pertx =  sin(xpoin)
          perty =  sin(ypoin)
          pertz = (pertx*perty)
          u_buffer(iNodeL,3) = -pertz
          
        end if
        
      end do
      !$acc end parallel loop

   end subroutine imposedDisplacement_elasticitySolverBuffer


end module MeshElasticitySolver_mod
