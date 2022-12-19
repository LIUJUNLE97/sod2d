module mod_meshConversorTool
   use mod_constants
   use mod_mpi
   use mod_utils
   use mod_gmsh_indices
   use iso_c_binding
   implicit none

!-----------------------------------   
#define _CHECK_ 1
#define int_size 4
!-----------------------------------   
#ifndef GEMPAINTERFACE
   interface
      subroutine gempa_do_partition(numElemsInRank,numRanksToPart,x,y,z,weights,part) bind(c)
        import c_int, c_double
        integer(kind=c_int), intent(in), value :: numElemsInRank
        integer(kind=c_int), intent(in), value :: numRanksToPart
        real(kind=c_double), intent(in) :: x(numElemsInRank)
        real(kind=c_double), intent(in) :: y(numElemsInRank)
        real(kind=c_double), intent(in) :: z(numElemsInRank)
        integer(kind=c_int), intent(in) :: weights(numElemsInRank)
        integer(kind=c_int), intent(out) :: part(numElemsInRank)
      end subroutine gempa_do_partition
   end interface
#endif
!-----------------------------------   

type :: vector_int
   integer(int_size),dimension(:), allocatable :: elems
end type vector_int

type :: matrix_int
   integer(int_size),dimension(:,:), allocatable :: elems
end type matrix_int

type :: vector_rp
   real(rp),dimension(:), allocatable :: elems
end type vector_rp

type :: matrix_rp
   real(rp),dimension(:,:), allocatable :: elems
end type matrix_rp

type :: jagged_vector_int
   type(vector_int),dimension(:), allocatable :: vector
end type jagged_vector_int

type :: jagged_matrix_int
   type(matrix_int),dimension(:), allocatable :: matrix
end type jagged_matrix_int

type :: jagged_vector_rp
   type(vector_rp),dimension(:), allocatable :: vector
end type jagged_vector_rp

type :: jagged_matrix_rp
   type(matrix_rp),dimension(:), allocatable :: matrix
end type jagged_matrix_rp

! ################################################################################################
! ----------------- VARS for new Par mesh FORMAT -------------------------------------------------
! ################################################################################################

integer(int_size) :: totalNumElements,totalNumNodesPar,totalNumNodesSrl,totalNumBoundsSrl
integer(int_size), allocatable :: numElemsRankPar(:), rankElemStart(:), rankElemEnd(:)
integer(int_size), allocatable :: rankNodeStart(:), rankNodeEnd(:),numNodesRankPar(:)
integer(int_size), allocatable :: numBoundsRankPar(:)

type(jagged_vector_int) :: elemGid_jv, globalIdSrl_jv, globalIdPar_jv

type(jagged_vector_int) :: connecVTK_jv
type(jagged_matrix_int) :: connecParOrig_jm,connecParWork_jm

type(jagged_matrix_rp) :: coordPar_jm

type(jagged_vector_int) :: workingNodesPar_jv
integer(int_size),allocatable :: numWorkingNodesRankPar(:)

type(jagged_vector_int) :: masSlaRankPar_jv
integer(int_size),allocatable :: nPerRankPar(:)

integer(int_size) :: numBoundCodes
type(jagged_vector_int) :: bouCodesPar
type(jagged_matrix_int) :: bouPar

integer(int_size) :: mshRankInMpiRankStart,mshRankInMpiRankEnd,numMshRanksInMpiRank
integer(int_size), allocatable :: mshRanksInMpiRank(:),mapMshRankToMpiRank(:)

! ################################################################################################
! ------------------------ VARS for MPI COMMS ----------------------------------------------------
! ################################################################################################

type(jagged_matrix_int) :: matrixCommScheme_jm
type(jagged_vector_int) :: commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv
integer(int_size),allocatable :: numNodesToComm(:),numRanksWithComms(:)




! ################################################################################################
integer, parameter :: fileId_geo = 99, fileId_bou = 98, fileId_dims = 97
! ################################################################################################
!     FOR DEBUG
character(128) :: file_name_deb,aux_string_mpirank_deb,aux_string_mshrank_deb
character(*), parameter :: fmt_csv = '(1x,*(g0,","))'
! ################################################################################################

contains

   subroutine read_gmsh_files_and_do_partitioning_in_parallel(file_path,file_name,isPeriodic,numMshRanks2Part)
      implicit none
      character(500), intent(in)  :: file_path, file_name
      logical, intent(in) :: isPeriodic
      integer(4),intent(in) :: numMshRanks2Part
      integer(4) :: numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numPerNodes
      integer(4) :: numElemsMpiRank,numNodesMpiRank,numBoundsMpiRank
      integer, allocatable :: listElemsMpiRank(:),listNodesMpiRank(:),listElemsBoundsMpiRank(:,:),connecMpiRank(:,:),boundFacesMpiRank(:,:)
      real(rp), allocatable :: coordNodesMpiRank(:,:)
      integer,allocatable :: numBoundaryNodes(:),numMpiBoundaryNodes(:),numMshBoundaryNodes(:),numInnerNodes(:)

      type(jagged_vector_int) :: listNodesMshRank_jv,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,boundaryNodes_jv,mshBoundaryNodes_jv,mpiBoundaryNodes_jv
      type(jagged_matrix_int) :: connecMshRank_jm,listElemsBoundsMshRank_jm,boundFacesMshRank_jm
      type(jagged_matrix_rp) :: coordMshRank_jm
      integer,dimension(0:numMshRanks2Part-1) :: iNodeStartPar

      integer(4) :: iMshRank,mshRank

      call read_dims_file(file_path,file_name,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh)

      !totalNumElements  = numElems
      !totalNumNodesSrl  = numNodes
      !totalNumBoundsSrl = numBoundFaces

      call read_geo_dat_file_in_parallel(file_path,file_name,isPeriodic,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,&
                        numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank,&
                        numBoundsMpiRank,listElemsBoundsMpiRank,boundFacesMpiRank,connecMpiRank,coordNodesMpiRank)

      if(isPeriodic) then
         write(*,*) 'TODO line 135 for isPeriodic'
         !isMeshPeriodic = .true.
         !tinc que cridar/fer funcio que llegeixi les boundaries paraleles
      end if

      call distribute_ranks2Part_in_mpiRank(numMshRanks2Part,mshRankInMpiRankStart,mshRankInMpiRankEnd,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank)

      call do_element_partitioning_gempa_in_parallel(file_path,file_name,isPeriodic,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numMshRanks2Part,&
                        numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank,connecMpiRank,coordNodesMpiRank,&
                        listNodesMshRank_jv,connecMshRank_jm,listElemsBoundsMshRank_jm,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,boundFacesMshRank_jm,coordMshRank_jm)

      !-----------------------------------------------------------------------
      allocate(boundaryNodes_jv%vector(numMshRanksInMpiRank))
      allocate(mpiBoundaryNodes_jv%vector(numMshRanksInMpiRank))
      allocate(mshBoundaryNodes_jv%vector(numMshRanksInMpiRank))

      allocate(numBoundaryNodes(numMshRanksInMpiRank))
      allocate(numMpiBoundaryNodes(numMshRanksInMpiRank))
      allocate(numMshBoundaryNodes(numMshRanksInMpiRank))
      allocate(numInnerNodes(numMshRanksInMpiRank))
      !-----------------------------------------------------------------------

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank) 

         call get_rankPartitionBoundaryNodes_in_parallel(mshRank,numElemsRankPar(iMshRank),numNodesRankPar(iMshRank),numBoundsRankPar(iMshRank),&
                  elemGid_jv%vector(iMshRank)%elems,listNodesMshRank_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,&
                  coordMshRank_jm%matrix(iMshRank)%elems,listElemsBoundsMshRank_jm%matrix(iMshRank)%elems,listBoundFacesMshRank_jv%vector(iMshRank)%elems,&
                  boundFacesMshRank_jm%matrix(iMshRank)%elems,numBoundaryNodes(iMshRank),numMshBoundaryNodes(iMshRank),numMpiBoundaryNodes(iMshRank),&
                  numInnerNodes(iMshRank),boundaryNodes_jv%vector(iMshRank)%elems,mpiBoundaryNodes_jv%vector(iMshRank)%elems,mshBoundaryNodes_jv%vector(iMshRank)%elems)
      end do

      call define_parallelNodePartitioning(numMshRanks2Part,iNodeStartPar)

      !--------------------------------------------------------------------------------------

      allocate(globalIdSrl_jv%vector(numMshRanksInMpiRank))
      allocate(globalIdPar_jv%vector(numMshRanksInMpiRank))
      allocate(connecVTK_jv%vector(numMshRanksInMpiRank))
      allocate(connecParOrig_jm%matrix(numMshRanksInMpiRank))
      allocate(connecParWork_jm%matrix(numMshRanksInMpiRank))
      allocate(coordPar_jm%matrix(numMshRanksInMpiRank))
      allocate(numWorkingNodesRankPar(numMshRanksInMpiRank))
      allocate(workingNodesPar_jv%vector(numMshRanksInMpiRank))

      allocate(nPerRankPar(numMshRanksInMpiRank))
      allocate(masSlaRankPar_jv%vector(numMshRanksInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)

         allocate(globalIdSrl_jv%vector(iMshRank)%elems(numNodesRankPar(iMshRank)))
         allocate(globalIdPar_jv%vector(iMshRank)%elems(numNodesRankPar(iMshRank)))

         allocate(connecVTK_jv%vector(iMshRank)%elems(numElemsRankPar(iMshRank)*nnode))
         allocate(connecParOrig_jm%matrix(iMshRank)%elems(numElemsRankPar(iMshRank),nnode))
         allocate(connecParWork_jm%matrix(iMshRank)%elems(numElemsRankPar(iMshRank),nnode))

         allocate(coordPar_jm%matrix(iMshRank)%elems(numNodesRankPar(iMshRank),3))

         !$acc kernels
         globalIdSrl_jv%vector(iMshRank)%elems(:) = -1
         globalIdPar_jv%vector(iMshRank)%elems(:) = -1
         connecVTK_jv%vector(iMshRank)%elems(:) = 0
         connecParOrig_jm%matrix(iMshRank)%elems(:,:) = 0
         connecParWork_jm%matrix(iMshRank)%elems(:,:) = 0
         coordPar_jm%matrix(iMshRank)%elems(:,:) = 0.0_rp
         !$acc end kernels

         !reordering nodes
         call reorder_nodes_in_mshRank(mshRank,numMshRanks2Part,numElemsRankPar(iMshRank),numNodesRankPar(iMshRank),numBoundsRankPar(iMshRank),&
                           elemGid_jv%vector(iMshRank)%elems,listNodesMshRank_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,iNodeStartPar,&
                           globalIdSrl_jv%vector(iMshRank)%elems,globalIdPar_jv%vector(iMshRank)%elems,&
                           connecVTK_jv%vector(iMshRank)%elems,connecParOrig_jm%matrix(iMshRank)%elems,&
                           boundFacesMshRank_jm%matrix(iMshRank)%elems)

         call set_nodesCoordinates(mshRank,numElemsRankPar(iMshRank),numNodesRankPar(iMshRank),globalIdSrl_jv%vector(iMshRank)%elems,&
            listNodesMshRank_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems,&
            coordPar_jm%matrix(iMshRank)%elems)

         if(isPeriodic) then
            write(*,*) 'TODO line 219 for isPeriodic -> create_masSla? maybe can be done before now, but must be done before calling create_working_lists_parallel'
            !call create_masSla_parallel()
         else
            nPerRankPar(iMshRank) = 0
            allocate(masSlaRankPar_jv%vector(numMshRanksInMpiRank)%elems(nPerRankPar(iMshRank)))
         end if
#if 1 
         call create_working_lists_parallel(isPeriodic,mshRank,numElemsRankPar(iMshRank),numNodesRankPar(iMshRank),globalIdSrl_jv%vector(iMshRank)%elems,&
               connecParOrig_jm%matrix(iMshRank)%elems,nPerRankPar(iMshRank),masSlaRankPar_jv%vector(iMshRank)%elems,&
               connecParWork_jm%matrix(iMshRank)%elems,numWorkingNodesRankPar(iMshRank),workingNodesPar_jv%vector(iMshRank)%elems)
#endif

      end do


      call generate_mpi_comm_scheme_parallel(numMshRanks2Part,numMpiBoundaryNodes,mpiBoundaryNodes_jv)

      !remaining things to do
      !-------------------------------------------------------------------------------------------

      !call define_boundPar_and_boundCodesPar_parallel()






      !call define_numDoFRankPar_and_numBoundaryNodesRankPar()
      !integer(int_size) :: ndofRankPar, numBoundaryNodesRankPar


      !-------------------------------------------------------------------------------------------

      deallocate(listElemsMpiRank)
      deallocate(listNodesMpiRank)
      !deallocate(listElemsBoundsMpiRank)
      deallocate(connecMpiRank)
      deallocate(coordNodesMpiRank)


      !call do_node_partitioning_and_connectivity()
      !!--------------
      !INSIDE do_node_partitioning-----------------------------------------
      !call get_rankPartitionBoundaryNodes(boundaryNodes)
      !call define_parallelNodePartitioning(iNodeStartPar,iNodeEndPar)
      !call define_mpi_boundaries_inPar(boundaryNodes,vecSharedBN_full)
      !deallocate(boundaryNodes)
      !call reorder_nodes_in_proc(iNodeStartPar)
      !call generate_mpi_comm_scheme(vecSharedBN_full)
      !deallocate(vecSharedBN_full)
      !INSIDE do_node_partitioning-----------------------------------------

      !!------------------------------------------------

      !call create_nodesCoordinates()

      !if(isMeshPeriodic) then
      !   call create_masSla_parallel()
      !end if

      !call create_working_lists()

      !call create_hdf5_meshFile_from_tool()

   end subroutine read_gmsh_files_and_do_partitioning_in_parallel
!---------------------------------------------------------------------------------------------------------------------------
   subroutine read_dims_file(file_path,file_name,numElems,numNodes,numBounds)
      implicit none     
      character(500), intent(in) :: file_path, file_name
      integer(4), intent(out) :: numElems,numNodes,numBounds
      character(500) :: file_type, line  

      write(file_type,*) ".dims.dat"

      if(mpi_rank.eq.0)write(*,*) "--| Reading dims.dat file..."
      open(unit=fileId_dims,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      read(fileId_dims,*) line, numNodes
      if(mpi_rank.eq.0)write(*,*) "--| Nodes on mesh: ",numNodes
      read(fileId_dims,*) line, numElems
      if(mpi_rank.eq.0)write(*,*) "--| Elemensts on mesh: ",numElems
      read(fileId_dims,*) line, numBounds
      if(mpi_rank.eq.0)write(*,*) "--| Boundary faces on mesh : ",numBounds
      if(mpi_rank.eq.0)write(*,*) "--| End of dims.dat file!"
      close(unit=fileId_dims)
    
   end subroutine read_dims_file
!---------------------------------------------------------------------------------------------------------------------------
   subroutine open_geo_dat_file(file_path,file_name)
      implicit none
      character(500), intent(in) :: file_path, file_name
      character(128) :: file_type_geo

      write(file_type_geo,*) ".geo.dat"

      if(mpi_rank.eq.0) write(*,*) "--| Opening geo.dat file..."
      open(unit=fileId_geo,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type_geo)),status="old")
   end subroutine open_geo_dat_file

   subroutine close_geo_dat_file()
      implicit none

      if(mpi_rank.eq.0) write(*,*) "--| Closing geo.dat file"
      close(unit=fileId_geo)
   end subroutine close_geo_dat_file
!---------------------------------------------------------------------------------------------------------------------------
   subroutine open_fix_bou_file(file_path,file_name)
      implicit none
      character(500), intent(in) :: file_path, file_name
      character(128) :: file_type_bou

      write(file_type_bou,*) ".fix.bou"

      if(mpi_rank.eq.0) write(*,*) "--| Opening fix.bou file..."
      open(unit=fileId_bou,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type_bou)),status="old")
   end subroutine open_fix_bou_file

   subroutine close_fix_bou_file()
      implicit none

      if(mpi_rank.eq.0) write(*,*) "--| Closing fix.bou file"
      close(unit=fileId_bou)
   end subroutine close_fix_bou_file
!---------------------------------------------------------------------------------------------------------------------------

   subroutine read_element_section_geo_file(numElemsSrl)
      implicit none
      integer(4),intent(in) :: numElemsSrl
      integer(4) :: iline

      ! Nodes/element section, blank for now
      !------------------------------------------------------------------------------------------
      read(fileId_geo,*) ! Section header
      do iline = 1,numElemsSrl
         read(fileId_geo,*)
      end do
      read(fileId_geo,*) ! Section ender
      !------------------------------------------------------------------------------------------
   end subroutine read_element_section_geo_file

   subroutine read_elem_connec_and_nodes_coords_from_geo_file(numElemsSrl,numNodesSrl,numElemsInRank,listElemsRank,numNodesInRank,connecRank,listNodesRank,coordNodesRank)
      implicit none
      integer,intent(in) :: numElemsSrl,numNodesSrl,numElemsInRank,listElemsRank(numElemsInRank)
      integer,intent(out) :: numNodesInRank,connecRank(numElemsInRank,nnode)
      integer,allocatable,intent(inout) :: listNodesRank(:)
      real(rp),allocatable,intent(inout) :: coordNodesRank(:,:)
      integer(4) :: auxConnec(nnode+1)
      integer(4) :: iLine,iElem,iNode,iNodeG,elemCnt,nodeCnt,iAux,jAux
      real(8) :: x,y,z
      character(2000) :: line
      integer(4) :: nextElemToRead,nextNodeToRead,prevNodeG

      integer,allocatable :: rawNodeListRank(:)

      if(mpi_rank.eq.0)write(*,*) "--| Reading element table..."

      !------------------------------------------------------------------------------------------
      read(fileId_geo,*) ! Section header
      elemCnt=0
      nextElemToRead=listElemsRank(elemCnt+1)
      do iElem=1,numElemsSrl
         !write(*,*) '(',mpi_rank,')iElem',iElem
         read(fileId_geo,'(a)') line
         if(iElem.eq.nextElemToRead) then
            read(line,*) (auxConnec(iNode), iNode=1,nnode+1)
            elemCnt=elemCnt+1

            !write(*,*) 'connecRank[',mpi_rank,'](',elemCnt,')'
            connecRank(elemCnt,1:nnode) = auxConnec(2:nnode+1)
            
            if(elemCnt.lt.numElemsInRank) nextElemToRead=listElemsRank(elemCnt+1)
            !if(mpi_rank.eq.0) write(*,*) iElem,connecRank(elemCnt,:)
         end if
      end do
      read(fileId_geo,*) ! Section ender

      !------------------------------------------------------------------------------------------

      allocate(rawNodeListRank(numElemsInRank*nnode))
      !$acc kernels
      rawNodeListRank(:) = 0
      !$acc end kernels

      ! add all the iNodeGs of this mpirank
      nodeCnt=0
      do iAux = 1,numElemsInRank
         do jAux = 1,nnode
            iNodeG = connecRank(iAux,jAux)
            nodeCnt=nodeCnt+1
            rawNodeListRank(nodeCnt) = iNodeG
         end do
      end do

      !sorting the nodes id
      call quicksort_array_int(rawNodeListRank)

      prevNodeG=0
      nodeCnt=0
      do iAux = 1,(numElemsInRank*nnode)
         if((rawNodeListRank(iAux)).ne.(prevNodeG)) then
            nodeCnt=nodeCnt+1
            prevNodeG = rawNodeListRank(iAux)
         end if
      end do
      numNodesInRank = nodeCnt
      !write(*,*) 'newMethod.[',mpi_rank,'] numNodesInRank ',numNodesInRank

      allocate(listNodesRank(numNodesInRank))
      allocate(coordNodesRank(numNodesInRank,ndime))

      !$acc kernels
      coordNodesRank(:,:) = 0
      listNodesRank(:) = 0
      !$acc end kernels

      prevNodeG=0
      nodeCnt=0
      do iAux = 1,(numElemsInRank*nnode)
         if((rawNodeListRank(iAux)).ne.(prevNodeG)) then
            nodeCnt=nodeCnt+1
            listNodesRank(nodeCnt) = rawNodeListRank(iAux)
            prevNodeG = rawNodeListRank(iAux)
         end if
      end do

      ! Nodal coordinates section
      !------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0)write(*,*) "--| Reading coordinates..."
      if (ndime .ne. 3) then
         write(*,*) 'Crashing in read_elem_connec_and_nodes_coords_from_geo_file! ndime is not 3! EXIT!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      nodeCnt=0
      nextNodeToRead=listNodesRank(nodeCnt+1)
      read(fileId_geo,*) ! Section header
      do iLine =1,numNodesSrl
         read(fileId_geo,*) iNode,x,y,z !iLine is equal to iNode (or should be!)
         !write(*,*) 'iNode',iNode,x,y,z
         if(iNode.eq.nextNodeToRead) then
            nodeCnt=nodeCnt+1
            coordNodesRank(nodeCnt,1)=real(x,rp)
            coordNodesRank(nodeCnt,2)=real(y,rp) 
            coordNodesRank(nodeCnt,3)=real(z,rp)

            if(nodeCnt.lt.numNodesInRank) nextNodeToRead=listNodesRank(nodeCnt+1)
         end if
      end do
      read(fileId_geo,*) ! Section ender

   end subroutine read_elem_connec_and_nodes_coords_from_geo_file

   subroutine read_boundaries(numBoundFacesSrl,numElemsInRank,listElemsInRank,connecInRank,numBoundsInRank,listElemsBoundsInRank,listBoundFacesInRank,boundFacesCodesInRank,boundFacesInRank)
      implicit none
      integer,intent(in) :: numBoundFacesSrl,numElemsInRank,listElemsInRank(numElemsInRank),connecInRank(numElemsInRank,nnode)
      integer(4), intent(out) :: numBoundsInRank
      integer, allocatable, intent(inout) :: listElemsBoundsInRank(:,:),listBoundFacesInRank(:),boundFacesCodesInRank(:),boundFacesInRank(:,:)
      integer(4) :: iBoundNodes(npbou),bou_aux(npbou+1),bou_code_aux
      character(2000) :: line
      integer(4),allocatable :: auxBoundFacesInRank(:,:)
      integer(4) :: iLine,iNode,iBound,iBoundBou,iElem,ind_gmsh,iElemG,iAux,jAux
      integer(4) :: iNodeG_inFace,iNodeG_inElem,nodeCnt
      integer :: f_iNodeG(4),e_iNodeG(8) !for the vertex of squares of the face and the element
      logical :: vertexFound

      read(fileId_geo,*) line! Section header
      read(fileId_bou,*) line! Section header

      !for the moment I will use this auxiliar array, but if we see that is allocating too much memory
      !maybe we can try another strategy, without allocating an array of numBoundFaceSrl
      !but doing a double loop
      !1st loop: determine the bounds that are in the msh rank
      !2nd loop: store only the nodes in the rank
      allocate(auxBoundFacesInRank(numBoundFacesSrl,2+npbou))
      auxBoundFacesInRank(:,:) = 0
      allocate(listElemsBoundsInRank(numElemsInRank,maxBoundsPerElem))
      listElemsBoundsInRank(:,:) = 0

      numBoundsInRank=0
      do iLine = 1,numBoundFacesSrl
         !reading boundary nodes-------------------------------
         read(fileId_geo,'(a)') line
         read(line,*) (bou_aux(iNode), iNode=1,npbou+1)
         iBound = bou_aux(1) !iLine & iBound must be the same
         iBoundNodes = bou_aux(2:npbou+1)
         !-----------------------------------------------------

         !reading boundary code--------------------------------
         read(fileId_bou,'(a)') line
         read(line,*) iBoundBou,bou_code_aux
         !-----------------------------------------------------

         if(iBound.ne.iBoundBou) then
            write(*,*) 'Error in read_boundaries: iBound',iBound,' not equal to iBoundBou',iBoundBou
            call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
         end if
         !write(*,*) 'iBound',iBound,'iBoundBou',iBoundBou,'bouCode',bou_code_aux,'iBoundNodes',iBoundNodes

         !fill the corners of the face to check
         do iAux=1,4
            ind_gmsh = gmsh2ij(posFaceVertices(iAux))
            f_iNodeG(iAux) = iBoundNodes(ind_gmsh)
         end do
         !write(*,*) '[',mpi_rank,']iBound',iBound,' -> f_iNodeG ',f_iNodeG(:)

         !Search to which element this face belongs
         do iElem=1,numElemsInRank
            nodeCnt=0
            iElemG=listElemsInRank(iElem)
            !fill the corners of the element 
            do iAux=1,8
               ind_gmsh = gmsh2ijk(posElemVertices(iAux))
               e_iNodeG(iAux) = connecInRank(iElem,ind_gmsh)
            end do

            !isBoundInElem=.false.
            fLoop: do iAux=1,4
               vertexFound=.false.
               iNodeG_inFace = f_iNodeG(iAux)
               eLoop: do jAux=1,8
                  iNodeG_inElem = e_iNodeG(jAux)
                  if(iNodeG_inFace .eq. iNodeG_inElem) then
                     nodeCnt=nodeCnt+1
                     vertexFound=.true.
                     exit eLoop
                  end if
               end do eLoop              
               if(not(vertexFound)) exit fLoop
            end do fLoop     
            if(nodeCnt.ge.4) then
               !isBoundInElem=.true.
               numBoundsInRank=numBoundsInRank+1
               !write(*,*) '[',mpi_rank,']iBound',iBound,' -> elem ',iElemG
               auxBoundFacesInRank(numBoundsInRank,1)=iBound
               auxBoundFacesInRank(numBoundsInRank,2)=bou_code_aux
               auxBoundFacesInRank(numBoundsInRank,3:npbou+2)=iBoundNodes(1:npbou)

               addLoop : do iAux = 1,maxBoundsPerElem
                  if(listElemsBoundsInRank(iElem,iAux).eq.0) then
                     !write(*,*) '[',mpi_rank,'] adding iBound',iBound,'to elem',iElemG
                     listElemsBoundsInRank(iElem,iAux) = iBound

                     exit addLoop
                  end if
               end do addLoop
            end if

         end do
      end do

      read(fileId_geo,*) ! Section ender
      read(fileId_bou,*) ! Section ender

      write(*,*) '[',mpi_rank,'] numBoundsInRank',numBoundsInRank
      !do iElem=1,numElemsInRank
      !   if(listElemsBoundsInRank(iElem,1).ne.0) then
      !      write(*,*) 'i',iElem,'iG',listElemsInRank(iElem),'bounds',listElemsBoundsInRank(iElem,:)
      !   end if
      !end do

      allocate(listBoundFacesInRank(numBoundsInRank))
      allocate(boundFacesCodesInRank(numBoundsInRank))
      allocate(boundFacesInRank(numBoundsInRank,npbou))

      do iBound=1,numBoundsInRank
         listBoundFacesInRank(iBound)  = auxBoundFacesInRank(iBound,1)
         boundFacesCodesInRank(iBound) = auxBoundFacesInRank(iBound,2)
         boundFacesInRank(iBound,:)    = auxBoundFacesInRank(iBound,3:npbou+2)
         !write(*,*) '[',mpi_rank,']bfmpirank(',iBound,')',boundFacesInRank(iBound,:)
      end do

      deallocate(auxBoundFacesInRank)

   end subroutine read_boundaries


   subroutine read_geo_dat_file_in_parallel(file_path,file_name,isPeriodic,numElemsSrl,numNodesSrl,numBoundFacesSrl,&
                                            numElemsMpiRank,listElemsMpiRank,&
                                            numNodesMpiRank,listNodesMpiRank,numBoundsMpiRank,listElemsBoundsMpiRank,&
                                            boundFacesMpiRank,connecMpiRank,coordNodesMpiRank)
      implicit none
      character(500), intent(in) :: file_path, file_name
      logical, intent(in)        :: isPeriodic
      integer(4), intent(in)     :: numElemsSrl,numNodesSrl,numBoundFacesSrl
      integer(4), intent(inout)  :: numElemsMpiRank,numNodesMpiRank,numBoundsMpiRank
      integer, allocatable, intent(inout) :: listElemsMpiRank(:),listNodesMpiRank(:),listElemsBoundsMpiRank(:,:),connecMpiRank(:,:),boundFacesMpiRank(:,:)
      real(rp), allocatable, intent(inout) :: coordNodesMpiRank(:,:)
      integer(4)                 :: iElemMpiRankStart,iElemMpiRankEnd
      integer(4)                 :: auxCnt,iAux,jAux
      integer(4)                 :: nodeCnt
      integer(4)                 :: iline,iNode,iBound,iBoundBou,iElem,ind_gmsh,iNodeG,iElemG
      integer(4)                 :: iNodeG_inFace,iNodeG_inElem
      integer(4)                 :: iBoundNodes(npbou),aux(nnode+1),bou_aux(npbou+1),bou_code_aux
      character(2000)            :: line
      real(8)                    :: x,y,z
      logical                    :: vertexFound
      integer(4), allocatable    :: auxBoundFacesMpiRank(:,:)
      integer                    :: f_iNodeG(4),e_iNodeG(8) !for the vertex of squares of the face and the element

      !----------------------------------------------------------------------------------------------

      call do_element_partitioning_serial(numElemsSrl,iElemMpiRankStart,iElemMpiRankEnd,numElemsMpiRank)

      !----------------------------------------------------------------------------------------------

      ! Connectivity table section
      !------------------------------------------------------------------------------------------
      allocate(connecMpiRank(numElemsMpiRank,nnode))
      allocate(listElemsMpiRank(numElemsMpiRank))

      !$acc kernels
      connecMpiRank(:,:) = 0
      listElemsMpiRank(:) = 0
      !$acc end kernels

      auxCnt=0
      do iElem = iElemMpiRankStart,iElemMpiRankEnd
         auxCnt=auxCnt+1
         listElemsMpiRank(auxCnt) = iElem
      end do

      call open_geo_dat_file(file_path,file_name)
      call read_element_section_geo_file(numElemsSrl)
      call read_elem_connec_and_nodes_coords_from_geo_file(numElemsSrl,numNodesSrl,numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,connecMpiRank,listNodesMpiRank,coordNodesMpiRank)

#if 0
      if(mpi_rank.eq.0)write(*,*) "--| Reading element table..."
      read(fileId_geo,*) ! Section header
      do iline = 1,(iElemProcStart-1)
         read(fileId_geo,*)
      end do
      auxCnt=0
      do iline = iElemProcStart,iElemProcEnd
         read(fileId_geo,'(a)') line
         read(line,*) (aux(inode), inode=1,nnode+1)
         auxCnt=auxCnt+1
         listElemsMpiRank(auxCnt) = iline
         connecMpiRank(auxCnt,1:nnode) = aux(2:nnode+1)
         !if(mpi_rank.eq.2) write(*,*) 'connecMpiRank(',auxCnt,'): ', connecMpiRank(auxCnt,:)
      end do
      do iline = (iElemProcEnd+1),numElemsSrl
         read(fileId_geo,*)
      end do
      read(fileId_geo,*) ! Section ender
      !------------------------------------------------------------------------------------------

      allocate(isNodeMpiRank(numNodesSrl))
      !$acc kernels
      isNodeMpiRank(:) = .false.
      !$acc end kernels

      nodeCnt=0
      do iAux = 1,numElemsMpiRank
         do jAux = 1,nnode
            iNode = connecMpiRank(iAux,jAux)
            if(not(isNodeMpiRank(iNode))) then
               isNodeMpiRank(iNode) = .true.
               nodeCnt=nodeCnt+1
            end if
         end do
      end do
      numNodesMpiRank = nodeCnt
      !write(*,*) '1.[',mpi_rank,'] nodeCnt ',nodeCnt

      ! Nodal coordinates section
      !------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0)write(*,*) "--| Reading coordinates..."
      if (ndime .ne. 3) then
         write(*,*) 'Crashing in read_geo_dat_file_in_parallel! ndime is not 3! EXIT!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      allocate(listNodesMpiRank(numNodesMpiRank))
      allocate(coordNodesMpiRank(numNodesMpiRank,ndime))

      !$acc kernels
      coordNodesMpiRank(:,:) = 0
      listNodesMpiRank(:) = 0
      !$acc end kernels

      auxCnt=0
      read(fileId_geo,*) ! Section header
      do iline =1,numNodesSrl
         read(fileId_geo,*) iNode,x,y,z !iline is equal to iNode (or should be!)
         if(isNodeMpiRank(iNode)) then
               auxCnt=auxCnt+1
               listNodesMpiRank(auxCnt) = iNode
               coordNodesMpiRank(auxCnt,1)=real(x,rp)
               coordNodesMpiRank(auxCnt,2)=real(y,rp) 
               coordNodesMpiRank(auxCnt,3)=real(z,rp)
               !if(mpi_rank.eq.2) write(*,*) '2.[',mpi_rank,'] lNiR(',auxCnt,')',iNode
         end if
      end do
      read(fileId_geo,*) ! Section ender
#endif
      ! Boundary nodes section
      !------------------------------------------------------------------------------------------

      if(isPeriodic) then
         if(mpi_rank.eq.0)write(*,*) "--| Shit for periodic...."
         write(*,*) 'ABORTING, SHIT FOR PERIODIC YET TO BE IMPLEMENTED!'
         call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
         !call open_fix_bou_file(file_path,file_name)
         !call read_boundaries(numBoundFacesSrl,numElemsMpiRank,listElemsMpiRank,connecMpiRank,numBoundsMpiRank,listElemsBoundsMpiRank,boundFacesMpiRank)

#if 0
      read(fileId_geo,*) line! Section header
      read(fileId_bou,*) line! Section header

      allocate(auxBoundFacesMpiRank(numBoundFacesSrl,2+npbou))
      auxBoundFacesMpiRank(:,:) = 0
      allocate(listElemsBoundsMpiRank(numElemsMpiRank,maxBoundsPerElem))
      listElemsBoundsMpiRank(:,:) = 0

      numBoundsMpiRank=0
      do iline = 1,numBoundFacesSrl
         !reading boundary nodes-------------------------------
         read(fileId_geo,'(a)') line
         read(line,*) (bou_aux(inode), inode=1,npbou+1)
         iBound = bou_aux(1) 
         iBoundNodes = bou_aux(2:npbou+1)
         !-----------------------------------------------------

         !reading boundary code--------------------------------
         read(fileId_bou,'(a)') line
         read(line,*) iBoundBou,bou_code_aux
         !-----------------------------------------------------

         if(iBound.ne.iBoundBou) then
            write(*,*) 'Error in read_geo_dat_file_in_parallel: iBound',iBound,' not equal to iBoundBou',iBoundBou
            call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
         end if
         !write(*,*) 'iBound',iBound,'iBoundBou',iBoundBou,'bouCode',bou_code_aux,'iBoundNodes',iBoundNodes

         !fill the corners of the face to check
         do iAux=1,4
            ind_gmsh = gmsh2ij(posFaceVertices(iAux))
            f_iNodeG(iAux) = iBoundNodes(ind_gmsh)
         end do
         !write(*,*) '[',mpi_rank,']iBound',iBound,' -> f_iNodeG ',f_iNodeG(:)

         !Search to which element this face belongs
         do iElem=1,numElemsMpiRank
            nodeCnt=0
            iElemG=listElemsMpiRank(iElem)
            !fill the corners of the element 
            do iAux=1,8
               ind_gmsh = gmsh2ijk(posElemVertices(iAux))
               e_iNodeG(iAux) = connecMpiRank(iElem,ind_gmsh)
            end do

            !isBoundInElem=.false.
            fLoop: do iAux=1,4
               vertexFound=.false.
               iNodeG_inFace = f_iNodeG(iAux)
               eLoop: do jAux=1,8
                  iNodeG_inElem = e_iNodeG(jAux)
                  if(iNodeG_inFace .eq. iNodeG_inElem) then
                     nodeCnt=nodeCnt+1
                     vertexFound=.true.
                     exit eLoop
                  end if
               end do eLoop              
               if(not(vertexFound)) exit fLoop
            end do fLoop     
            if(nodeCnt.ge.4) then
               !isBoundInElem=.true.
               numBoundsMpiRank=numBoundsMpiRank+1
               !write(*,*) '[',mpi_rank,']iBound',iBound,' -> elem ',iElemG
               auxBoundFacesMpiRank(numBoundsMpiRank,1)=iBound
               auxBoundFacesMpiRank(numBoundsMpiRank,2)=iElemG
               auxBoundFacesMpiRank(numBoundsMpiRank,3:npbou+2)=iBoundNodes(1:npbou)

               addLoop : do iAux = 1,maxBoundsPerElem
                  if(listElemsBoundsMpiRank(iElem,iAux).eq.0) then
                     !write(*,*) '[',mpi_rank,'] adding iBound',iBound,'to elem',iElemG
                     listElemsBoundsMpiRank(iElem,iAux) = iBound

                     exit addLoop
                  end if
               end do addLoop
            end if

         end do
         !boundGMSH(iline,1:npbou) = bou_aux(2:npbou+1)
      end do

      read(fileId_geo,*) ! Section ender
      read(fileId_bou,*) ! Section ender

      write(*,*) '[',mpi_rank,'] numBoundsMpiRank',numBoundsMpiRank
      !do iElem=1,numElemsMpiRank
      !   if(listElemsBoundsMpiRank(iElem,1).ne.0) then
      !      write(*,*) 'i',iElem,'iG',listElemsMpiRank(iElem),'bounds',listElemsBoundsMpiRank(iElem,:)
      !   end if
      !end do

      call MPI_Allreduce(numBoundsMpiRank,aux_sum_NBP,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,mpi_err)

      if(aux_sum_NBP.ne.numBoundFacesSrl) then
         write(*,*) 'Error in read_geo_dat_file_in_parallel: aux_sum_NBP',aux_sum_NBP,' not equal to numBoundFacesSrl',numBoundFacesSrl
         call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
      end if

      allocate(boundFacesMpiRank(numBoundsMpiRank,2+npbou))

      do iBound=1,numBoundsMpiRank
         boundFacesMpiRank(iBound,:)=auxBoundFacesMpiRank(iBound,:)
         !write(*,*) '[',mpi_rank,']bfmpirank(',iBound,')',boundFacesMpiRank(iBound,:)
      end do

      deallocate(auxBoundFacesMpiRank)
#endif

         !call close_fix_bou_file()
      end if

      call close_geo_dat_file()

      !------------------------------------------------------------------------------------------

      if(mpi_rank.eq.0)write(*,*) "--| End of reading GEO.DAT file in parallel!"
      if(mpi_rank.eq.0)write(*,*) "--------------------------------------------"
    
   end subroutine read_geo_dat_file_in_parallel

   subroutine do_element_partitioning_gempa_in_parallel(file_path,file_name,isPeriodic,numElemsSrl,numNodesSrl,numBoundFacesSrl,&
                  numMshRanks2Part,numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank,connecMpiRank,coordNodesMpiRank,&
                  listNodesMshRank_jv,connecMshRank_jm,listElemsBoundsMshRank_jm,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,boundFacesMshRank_jm,coordMshRank_jm)
      implicit none
      character(500), intent(in) :: file_path, file_name
      logical,intent(in) :: isPeriodic
      integer(4),intent(in) :: numElemsSrl,numNodesSrl,numBoundFacesSrl
      integer(4),intent(in) :: numElemsMpiRank,numNodesMpiRank,numMshRanks2Part
      integer,intent(in) :: listElemsMpiRank(numElemsMpiRank),listNodesMpiRank(numNodesMpiRank),connecMpiRank(numElemsMpiRank,nnode)
      real(rp),intent(in) :: coordNodesMpiRank(numNodesMpiRank,ndime)
      type(jagged_vector_int),intent(out) :: listNodesMshRank_jv,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv
      type(jagged_matrix_int),intent(out) :: connecMshRank_jm,listElemsBoundsMshRank_jm,boundFacesMshRank_jm
      type(jagged_matrix_rp),intent(out) :: coordMshRank_jm

      real(8), allocatable, dimension(:) :: x,y,z
      integer, allocatable :: listElems2Part(:),weightElems2Part(:),unfoldedElems(:)
      integer, allocatable :: linkedElems(:,:),elemPart(:,:),elemPartAux(:,:)

      integer :: ii,jj,kk,m,iElem,iNode,iElemG,iNodeG,iBound,iPos
      integer :: mshRank,mpiRank,iMshRank,iMpiRank
      integer :: numElems2Part,numValues2get,aux_numElemsMshRank,aux_sum_NBP,aux_sum_NB_mpiRank

      real(8) :: x_a,y_a,z_a
      integer,dimension(0:numMshRanks2Part-1) :: numElems2mshRankInMpiRank,numElemsInMshRank
      integer,dimension(0:numMshRanks2Part-1,0:mpi_size-1) :: matNumElems2mshRankInMpiRank

      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size,target_displacement
      !logical, allocatable :: auxMatrixNumNodesInMpiRank(:,:)

      if(mpi_rank.eq.0) write(*,*) '--| Doing element partitioning GEMPA in parallel...'

      !1. obtain the list element 2 part using gempa
      if(isPeriodic) then
         !TODO
         !call get_listElems2Par_Periodic(listElems2Part,weightElems2Par,linkedElems)
      else
         call get_listElems2Par_in_parallel(numElemsMpiRank,listElemsMpiRank,listElems2Part,weightElems2Part,numElems2Part)
      end if

      allocate(elemPart(numElems2Part,4))
      allocate(x(numElems2Part))
      allocate(y(numElems2Part))
      allocate(z(numElems2Part))

      do iElem=1,numElems2Part
         elemPart(iElem,1) = listElems2Part(iElem)
         elemPart(iElem,3) = weightElems2Part(iElem)
         elemPart(iElem,4) = iElem

         x_a=0.0_rp
         y_a=0.0_rp
         z_a=0.0_rp
         do ii=1,8
            m = gmsh2ijk(posElemVertices(ii))
            iNodeG = connecMpiRank(iElem,m)

            iPos = binarySearch_int_i(listNodesMpiRank,iNodeG)

            x_a = x_a + real(coordNodesMpiRank(iPos,1),8)
            y_a = y_a + real(coordNodesMpiRank(iPos,2),8)
            z_a = z_a + real(coordNodesMpiRank(iPos,3),8)
         end do

         x(iElem) = x_a/8.0d0
         y(iElem) = y_a/8.0d0
         z(iElem) = z_a/8.0d0
         !write(*,*) '[',mpi_rank,']iElem',iElem,'iElemG',elemPart(iElem,1),'w', elemPart(iElem,3),'x',x(iElem),'y',y(iElem),'z',z(iElem)
      end do

#ifndef GEMPAINTERFACE
      !---- CALLING GEMPA in PARALLEL-----------------------
      call gempa_do_partition(numElems2Part,numMshRanks2Part,x,y,z,elemPart(:,3),elemPart(:,2))
      !---------------------------------------------------------------------
#else
      write(*,*) 'FATAL ERROR! Trying to call do_element_partitioning_gempa() for program compiled with flag GEMPAINTERFACE!'
      call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
#endif

      if(isPeriodic) then
      ! if parallel, now we have to put the rank to the 'slave' elements who were not included in the partitioning process
      ! since they were linked to a 'master' element
         !TODO
      end if

      ! @now we have to info all the ranks about which are the elements that they will 'store'
      ! since we have done the partitioning in paralel...

      numElems2mshRankInMpiRank(:)=0
      !!!$acc parallel loop
      do iElem=1,numElems2Part
         !iElemG   = elemPart(iElem,1)
         mshRank = elemPart(iElem,2)-1 !elemPart came with rank=1:mpi_size
         !!!$acc atomic update
         numElems2mshRankInMpiRank(mshRank) = numElems2mshRankInMpiRank(mshRank) + 1
         !!!$acc end atomic
      end do
      !!!$acc end parallel loop

      !write(*,*) 'numElems2mshRankInMpiRank[',mpi_rank,'] ',numElems2mshRankInMpiRank(:)

#if 0
      allocate(auxMatrixNumNodesInMpiRank(numNodesMpiRank,0:numMshRanks2Part-1))
      auxMatrixNumNodesInMpiRank = .false.

      do iElem=1,numElems2Part
         iElemG   = elemPart(iElem,1)
         mshRank = elemPart(iElem,2)-1 !elemPart came with rank=1:mpi_size

         do ii=1,nnode
            iNodeG = connecMpiRank(iElem,ii)
            iPos = binarySearch_int_i(listNodesMpiRank,iNodeG)
            auxMatrixNumNodesInMpiRank(iPos,iMshRank) = .true.
         end do

      end do

      numNodes2mshRankInMpiRank(:)=0
      do iMshRank=0,numMshRanks2Part-1
         do iNode=1,numNodesMpiRank
            if(auxMatrixNumNodesInMpiRank(iNode,iMshRank)) then
               numNodes2mshRankInMpiRank(iMshRank) = numNodes2mshRankInMpiRank(iMshRank) + 1
            end if
         end do
      end do

      write(*,*) '[',mpi_rank,']numNodes2mshRankInMpiRank',numNodes2mshRankInMpiRank(:)

      deallocate(auxMatrixNumNodesInMpiRank)
#endif
#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'elemPartition_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,iElem'
      do ii=1,numElems2Part
         iElemG   = elemPart(ii,1)
         mshRank = elemPart(ii,2)-1 
         iElem    = elemPart(ii,4)

         x_a=0.
         y_a=0.
         z_a=0.
         do jj=1,8
            m = gmsh2ijk(posElemVertices(jj))
            iNodeG = connecMpiRank(iElem,m)

            iPos = binarySearch_int_i(listNodesMpiRank,iNodeG)

            x_a = x_a + real(coordNodesMpiRank(iPos,1),8)
            y_a = y_a + real(coordNodesMpiRank(iPos,2),8)
            z_a = z_a + real(coordNodesMpiRank(iPos,3),8)
         end do
         x_a = x_a/8.d0
         y_a = y_a/8.d0
         z_a = z_a/8.d0

         write(1,fmt_csv) x_a,y_a,z_a,iElemG,mshRank,iElem
      end do
      close(1)
#endif


      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(numElems2mshRankInMpiRank,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMpiRank=0,mpi_size-1
         call MPI_Get(matNumElems2mshRankInMpiRank(:,iMpiRank),numMshRanks2Part,MPI_INTEGER,iMpiRank,target_displacement,&
                                               numMshRanks2Part,MPI_INTEGER,window_id,mpi_err)
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      call quicksort_matrix_int(elemPart,2)

      ii=1
      do mshRank=0,numMshRanks2Part-1
         jj=numElems2mshRankInMpiRank(mshRank)
         if(jj.ne.0) then
            kk=ii+jj-1
            call quicksort_matrix_int(elemPart,1,ii,kk)
            ii=kk+1
         endif
      end do

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'elemPartitionOrdered_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,iElem'
      do ii=1,numElems2Part
         iElemG  = elemPart(ii,1)
         mshRank = elemPart(ii,2)-1
         iElem   = elemPart(ii,4)

         x_a=0.
         y_a=0.
         z_a=0.
         do jj=1,8
            m = gmsh2ijk(posElemVertices(jj))
            iNodeG = connecMpiRank(iElem,m)

            iPos = binarySearch_int_i(listNodesMpiRank,iNodeG)

            x_a = x_a + real(coordNodesMpiRank(iPos,1),8)
            y_a = y_a + real(coordNodesMpiRank(iPos,2),8)
            z_a = z_a + real(coordNodesMpiRank(iPos,3),8)
         end do
         x_a = x_a/8.d0
         y_a = y_a/8.d0
         z_a = z_a/8.d0

         write(1,fmt_csv) x_a,y_a,z_a,iElemG,mshRank,iElem
      end do
      close(1)
#endif

      deallocate(x)
      deallocate(y)
      deallocate(z)

      !--------------------------------------------------------------------------------------

      numElemsInMshRank(:)=0

      do mpiRank=0,mpi_size-1
         do mshRank=0,numMshRanks2Part-1
            numElemsInMshRank(mshRank) = numElemsInMshRank(mshRank) + matNumElems2mshRankInMpiRank(mshRank,mpiRank)
         end do
      end do

      !write(*,*) 'numElemsInMshRank[',mpi_rank,'] ',numElemsInMshRank(:)
      !if(mpi_rank.eq.0) then
      !   do mshRank=0,numMshRanks2Part-1
      !      write(*,*) 'mshRank',mshRank,'numElemsInRankP',numElemsInMshRank(mshRank)
      !   end do
      !end if

      !call distribute_ranks2Part_in_mpiRank(numMshRanks2Part,iMshRankInMpiRankStart,iMshRankInMpiRankEnd,numMshRanksInMpiRank,mshRanksInMpiRank)

      !maxNumElemsInMpiRank = 0
      !do iMshRank=1,numMshRanks2Part
      !   maxNumElemsInMpiRank=max(numElemsInMshRank(iMshRank),maxNumElemsInMpiRank)
      !end do
      !write(*,*) '[',mpi_rank,']numMshRanksInMpiRank',numMshRanksInMpiRank,'maxNumElemsInRank',maxNumElemsInRank

      !allocate(elemGid(maxNumElemsInMpiRank,numMshRanksInMpiRank))
      !$acc kernels
      !elemGid(:,:) = -1
      !$acc end kernels

      allocate(elemGid_jv%vector(numMshRanksInMpiRank))
      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         allocate(elemGid_jv%vector(iMshRank)%elems(numElemsInMshRank(mshRank)))
         !$acc kernels
         elemGid_jv%vector(iMshRank)%elems(:) = -1
         !$acc end kernels
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*numElems2Part

      call MPI_Win_create(elemPart(:,1),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         iPos=1
         mshRank=mshRanksInMpiRank(iMshRank)

         do mpiRank=0,mpi_size-1
            numValues2get=matNumElems2mshRankInMpiRank(mshRank,mpiRank)

            if(numValues2get.ne.0) then
               target_displacement=0
               do jj=0,mshRank-1
                  target_displacement=target_displacement+matNumElems2mshRankInMpiRank(jj,mpiRank)
               end do

               !write(*,*) 'mpiRank',mpiRank,'mshRank',mshRank,'iPos',iPos,'nv2g',numValues2get,'td',target_displacement
               !call MPI_Get(elemGid(iPos,iMshRank),numValues2get,MPI_INTEGER,mpiRank,target_displacement,&
               !           numValues2get,MPI_INTEGER,window_id,mpi_err)
               call MPI_Get(elemGid_jv%vector(iMshRank)%elems(iPos),numValues2get,MPI_INTEGER,mpiRank,target_displacement,&
                          numValues2get,MPI_INTEGER,window_id,mpi_err)

               iPos = iPos + numValues2get
            end if

         end do
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         write(aux_string_mshrank_deb,'(I0)') mshRank
         file_name_deb = 'elemGid_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
         open(1,file=file_name_deb)
         do iElem=1,numElemsInMshRank(iMshRank)
            !write(1,fmt_csv) iElem,elemGid(iElem,iMshRank)
            write(1,fmt_csv) iElem,elemGid_jv%vector(iMshRank)%elems(iElem)
         end do
         close(1)
      end do
#endif
!---------------------------------------------------------------------------------

      allocate(numElemsRankPar(numMshRanksInMpiRank))
      allocate(rankElemStart(numMshRanksInMpiRank))
      allocate(rankElemEnd(numMshRanksInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank= mshRanksInMpiRank(iMshRank)
         numElemsRankPar(iMshRank) = numElemsInMshRank(mshRank)
      end do
      !write(*,*) 'numElemsRankPar[',mpi_rank,']',numElemsRankPar(:)

      rankElemStart(:)=1
      do iMshRank=1,numMshRanksInMpiRank
         mshRank= mshRanksInMpiRank(iMshRank)
         do jj=0,mshRank-1
            rankElemStart(iMshRank)=rankElemStart(iMshRank)+numElemsInMshRank(jj)
         end do
         rankElemEnd(iMshRank) = rankElemStart(iMshRank) + numElemsInMshRank(jj)-1
      end do
      !write(*,*) 'numElemsRankPar[',mpi_rank,']',numElemsRankPar(:),'eStart',rankElemStart(:),'eEnd ',rankElemEnd(:)

!------------------------------------------------------------------------------

      allocate(numNodesRankPar(numMshRanksInMpiRank))
      allocate(numBoundsRankPar(numMshRanksInMpiRank))
      allocate(listNodesMshRank_jv%vector(numMshRanksInMpiRank))
      allocate(connecMshRank_jm%matrix(numMshRanksInMpiRank))!,boundMshRank_jm,bou_codesMshRank_jm
      allocate(coordMshRank_jm%matrix(numMshRanksInMpiRank))
      
      if(numBoundFacesSrl.ne.0) then
         !ull, aixo potser canviarho si ara TOTES son periodiques...
         !see below..
         allocate(listBoundFacesMshRank_jv%vector(numMshRanksInMpiRank))
         allocate(boundFacesCodesMshRank_jv%vector(numMshRanksInMpiRank))
         allocate(listElemsBoundsMshRank_jm%matrix(numMshRanksInMpiRank))
         allocate(boundFacesMshRank_jm%matrix(numMshRanksInMpiRank))
      end if

      numNodesRankPar(:)  = 0
      numBoundsRankPar(:) = 0

      !ara puc aplicar la nova idea!
      do iMshRank=1,numMshRanksInMpiRank
         aux_numElemsMshRank=numElemsRankPar(iMshRank)
         mshRank= mshRanksInMpiRank(iMshRank)
         write(*,*) '#Reading stuff for mshRank',mshRank,'in mpi_rank',mpi_rank,'(iMshRank',iMshRank,')','aux_numElemsMshRank',aux_numElemsMshRank

         allocate(connecMshRank_jm%matrix(iMshRank)%elems(aux_numElemsMshRank,nnode))

         call open_geo_dat_file(file_path,file_name)
         call read_element_section_geo_file(numElemsSrl)
         call read_elem_connec_and_nodes_coords_from_geo_file(numElemsSrl,numNodesSrl,aux_numElemsMshRank,elemGid_jv%vector(iMshRank)%elems,numNodesRankPar(iMshRank),&
         connecMshRank_jm%matrix(iMshRank)%elems,listNodesMshRank_jv%vector(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems)
         
         if(numBoundFacesSrl.ne.0) then
            !ull, aixo potser canviarho si ara TOTES son periodiques...
            !quan implementi el periodic tenirho en compte ja que potser no tinc cap boundary com a tal
            !si no que totes son periodiques... de moment ho deixo

            call open_fix_bou_file(file_path,file_name)
            call read_boundaries(numBoundFacesSrl,aux_numElemsMshRank,elemGid_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,&
               numBoundsRankPar(iMshRank),listElemsBoundsMshRank_jm%matrix(iMshRank)%elems,&
               listBoundFacesMshRank_jv%vector(iMshRank)%elems,boundFacesCodesMshRank_jv%vector(iMshRank)%elems,boundFacesMshRank_jm%matrix(iMshRank)%elems)

            call close_fix_bou_file()
         end if

         call close_geo_dat_file()

      end do

      !fer el check aqui de que s'han llegit correctament totes les boundaries
      if(numBoundFacesSrl.ne.0) then
         aux_sum_NB_mpiRank=0
         do iMshRank=1,numMshRanksInMpiRank
            aux_sum_NB_mpiRank = aux_sum_NB_mpiRank + numBoundsRankPar(iMshRank)
         end do

         call MPI_Allreduce(aux_sum_NB_mpiRank,aux_sum_NBP,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,mpi_err)

         if(aux_sum_NBP.ne.numBoundFacesSrl) then
            write(*,*) '[',mpi_rank,']Error in do_element_partitioning_gempa_in_parallel(..) after read_boundaries(..): aux_sum_NBP',aux_sum_NBP,' not equal to  numBoundFacesSrl',numBoundFacesSrl
            call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
         end if
      end if

#if 0
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      !now we need to comm the corresponding connectivty,nodes,coords and boundaries for the elems we send to each proc!

      !listElemsMpiRank(numElemsMpiRank)
      !listNodesMpiRank(numNodesMpiRank)
      !connecMpiRank(numElemsMpiRank,nnode)
      !listElemsBoundsMpiRank(numElemsMpiRank,maxBoundsPerElem)
      !coordNodesMpiRank(numNodesMpiRank,ndime)

      !1. so first: define the connec that each mpiRank need to send to each corresponding mshRank



      !allocate(aux_connecMpiRank(numElems2Part*nnode))
      aux_array_size = numElems2Part*nnode
      allocate(aux_array2snd_int(aux_array_size))
      allocate(aux_array2rcv_int(maxNumElemsInMpiRank*nnode,numMshRanksInMpiRank))
      allocate(connecMshRank(maxNumElemsInMpiRank,nnode,numMshRanksInMpiRank))

      !fill the aux_array
      iPos=0
      do ii=1,numElems2Part
         iElem = elemPart(ii,4)
         do jj=1,nnode
            iPos=iPos+1
            iNodeG = connecMpiRank(iElem,jj)
            aux_array2snd_int(iPos)=iNodeG
         end do
      end do

!-----------------------------------------------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*aux_array_size

      call MPI_Win_create(aux_array2snd_int,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      do ii=1,numMshRanksInMpiRank
         iPos=1
         iMshRank=mshRanksInMpiRank(ii)
         do iMpiRank=0,mpi_size-1
            numValues2get=(nnode*matNumElems2mshRankInMpiRank(iMshRank,iMpiRank))

            if(numValues2get.ne.0) then
               target_displacement=0
               do jj=0,iMshRank-1
                  target_displacement=target_displacement+(nnode*matNumElems2mshRankInMpiRank(jj,iMpiRank))
               end do

               !write(*,*) 'iMpiRank',iMpiRank,'iMshRank',iMshRank,'iPos',iPos,'nv2g',numValues2get,'td',target_displacement
               call MPI_Get(aux_array2rcv_int(iPos,ii),numValues2get,MPI_INTEGER,iMpiRank,target_displacement,&
                         numValues2get,MPI_INTEGER,window_id,mpi_err)

               iPos = iPos + numValues2get
            end if
         end do

      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)

      do ii=1,numMshRanksInMpiRank
         iMshRank=mshRanksInMpiRank(ii)
         jj=0
         do iElem=1,numElemsInMshRank(iMshRank) 
            do iNode=1,nnode
               jj=jj+1
               connecMshRank(iElem,iNode,ii) = aux_array2rcv_int(jj,ii)
            end do
         end do
      end do
!---------------------------------------------------------------------------------
#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      do ii=1,numMshRanksInMpiRank
         iMshRank=mshRanksInMpiRank(ii)
         write(aux_string_mshrank_deb,'(I0)') iMshRank
         file_name_deb = 'connec_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
         open(1,file=file_name_deb)
         do iElem=1,numElemsInMshRank(iMshRank)
            write(1,fmt_csv) elemGid(iElem,ii),connecMshRank(iElem,:,ii)
         end do
         close(1)
      end do
#endif
!---------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------

      deallocate(aux_array2snd_int)
      deallocate(aux_array2rcv_int)

      aux_array_size = numElems2Part*maxBoundsPerElem
      allocate(aux_array2snd_int(aux_array_size))
      allocate(aux_array2rcv_int(maxNumElemsInMpiRank*maxBoundsPerElem,numMshRanksInMpiRank))
      allocate(listElemsBoundsMshRank(maxNumElemsInMpiRank,maxBoundsPerElem,numMshRanksInMpiRank))
      
      !fill the aux_array
      iPos=0
      do ii=1,numElems2Part
         iElem = elemPart(ii,4)
         do jj=1,maxBoundsPerElem
            iPos=iPos+1
            iBound = listElemsBoundsMpiRank(iElem,jj)
            aux_array2snd_int(iPos)=iBound
         end do
      end do

!-----------------------------------------------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*aux_array_size

      call MPI_Win_create(aux_array2snd_int,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      do ii=1,numMshRanksInMpiRank
         iPos=1
         iMshRank=mshRanksInMpiRank(ii)
         do iMpiRank=0,mpi_size-1
            numValues2get=(maxBoundsPerElem*matNumElems2mshRankInMpiRank(iMshRank,iMpiRank))

            if(numValues2get.ne.0) then
               target_displacement=0
               do jj=0,iMshRank-1
                  target_displacement=target_displacement+(maxBoundsPerElem*matNumElems2mshRankInMpiRank(jj,iMpiRank))
               end do

               !write(*,*) 'iMpiRank',iMpiRank,'iMshRank',iMshRank,'iPos',iPos,'nv2g',numValues2get,'td',target_displacement
               call MPI_Get(aux_array2rcv_int(iPos,ii),numValues2get,MPI_INTEGER,iMpiRank,target_displacement,&
                         numValues2get,MPI_INTEGER,window_id,mpi_err)

               iPos = iPos + numValues2get
            end if
         end do

      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)

      do ii=1,numMshRanksInMpiRank
         iMshRank=mshRanksInMpiRank(ii)
         jj=0
         do iElem=1,numElemsInMshRank(iMshRank) 
            do iBound=1,maxBoundsPerElem
               jj=jj+1
               listElemsBoundsMshRank(iElem,iBound,ii) = aux_array2rcv_int(jj,ii)
            end do
         end do
      end do

!---------------------------------------------------------------------------------
#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      do ii=1,numMshRanksInMpiRank
         iMshRank=mshRanksInMpiRank(ii)
         write(aux_string_mshrank_deb,'(I0)') iMshRank
         file_name_deb = 'listElemsBounds_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
         open(1,file=file_name_deb)
         do iElem=1,numElemsInMshRank(iMshRank)
            write(1,fmt_csv) elemGid(iElem,ii),listElemsBoundsMshRank(iElem,:,ii)
         end do
         close(1)
      end do
#endif
!---------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------

      !2. then identify the required nodes and its coords

      !3. the same bor the bounds linked to the elems

      !4. em puc basar en el esquema comuniacio dels elements elemPart(:,1)... i si aprofito per crear varies finestres i ferho tot? :o
 #endif  

      deallocate(listElems2Part)
      deallocate(elemPart)


!------------------------------------------------------------------------------

   end subroutine do_element_partitioning_gempa_in_parallel
!------------------------------------------------------------------------------------------------

   subroutine do_element_partitioning_serial(numElems2Par,iElemStart,iElemEnd,iElemsInRank)
      implicit none
      integer, intent(in) :: numElems2Par
      integer, intent(out) :: iElemStart,iElemEnd,iElemsInRank

      !TODO REVISAR EL PARTICIONAMENT DELS ELEMENTS PERQUE LA DISTRIBUCIO QUE FAIG
      !Al iElemsInRank pot donar problemes quan es un valor molt petit!
      !---- number of elements and range this process will write
      iElemsInRank = (numElems2Par + mpi_size - 1) / mpi_size
      iElemStart = iElemsInRank * mpi_rank + 1
      iElemEnd   = iElemsInRank * (mpi_rank + 1)
      if (iElemEnd .gt. numElems2Par) iElemEnd = numElems2Par
      iElemsInRank = iElemEnd - iElemStart + 1

      !write(*,*) '#rank ',mpi_rank,' iElemsInRank ',iElemsInRank,' iElemS ',iElemStart,' iElemE ',iElemEnd

   end subroutine do_element_partitioning_serial

   subroutine get_listElems2Par_in_parallel(numElemsMpiRank,listElemsMpiRank,listElems2Part,weightElems2Part,numElems2Part)
      implicit none
      integer(4), intent(in)  :: numElemsMpiRank,listElemsMpiRank(numElemsMpiRank)
      integer, allocatable, intent(out) :: listElems2Part(:),weightElems2Part(:)
      integer, intent(out) :: numElems2Part
      integer :: iElem, iElemG

      allocate(listElems2Part(numElemsMpiRank))
      allocate(weightElems2Part(numElemsMpiRank))
      !$acc parallel loop
      do iElem=1,numElemsMpiRank
         listElems2Part(iElem) = listElemsMpiRank(iElem)
         weightElems2Part(iElem) = 1
      end do
      !$acc end parallel loop

      numElems2Part = numElemsMpiRank

   end subroutine get_listElems2Par_in_parallel

   subroutine distribute_ranks2Part_in_mpiRank(numMshRanks2Part,iMshRankStart,iMshRankEnd,numRanksMpiRank,ranksMpiRank,mapRanksToMpiRank)
      implicit none
      integer, intent(in) :: numMshRanks2Part
      integer, intent(out) :: iMshRankStart,iMshRankEnd,numRanksMpiRank
      integer, allocatable, intent(out) :: ranksMpiRank(:),mapRanksToMpiRank(:)
      integer, dimension(0:mpi_size-1) :: vecRanksMpiRank
      integer :: iMpiRank,iMshRankCnt,ii

      if(numMshRanks2Part.ge.mpi_size) then
         !---- number of elements and range this process will write
         call distribution_algorithm(numMshRanks2Part,mpi_size,vecRanksMpiRank)

         numRanksMpiRank = vecRanksMpiRank(mpi_rank)
         iMshRankStart=0
         do iMpiRank=0,(mpi_rank-1) !find the iMshRankStart
            iMshRankStart = iMshRankStart + vecRanksMpiRank(iMpiRank)
         end do
         iMshRankEnd = iMshRankStart + (numRanksMpiRank - 1)

         allocate(ranksMpiRank(numRanksMpiRank))
         do iMpiRank=1,numRanksMpiRank
            ranksMpiRank(iMpiRank)= iMshRankStart + (iMpiRank-1)
         end do

         allocate(mapRanksToMpiRank(numMshRanks2Part))
         iMshRankCnt=0
         do iMpiRank=0,(mpi_size-1)
            do ii=1,vecRanksMpiRank(iMpiRank)
               iMshRankCnt=iMshRankCnt+1
               mapRanksToMpiRank(iMshRankCnt)=iMpiRank
            end do
         end do

      else
         write(*,*) 'numMshRanks2Part<mpi_size ... haig de pensar en aquest cas...'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      !write(*,*) '#rank',mpi_rank,'numRanksMpiRank',numRanksMpiRank,'ranksInP',ranksMpiRank(:),'mapRanksToMpiRank',mapRanksToMpiRank(:)

   end subroutine distribute_ranks2Part_in_mpiRank

   subroutine get_rankPartitionBoundaryNodes_in_parallel(mshRank,numElemsInRank,numNodesInRank,numBoundsInRank,elemGid,listNodesInRank,connecInRank,coordInRank,&
               listElemsBoundsInRank,listBoundFacesInRank,boundFacesInRank,numBoundNodesRankPar,numMshBoundNodesRankPar,numMpiBoundNodesRankPar,numInnerNodesRankPar,boundaryNodes,mpiBoundaryNodes,mshBoundaryNodes)
      implicit none
      integer, intent(in) :: mshRank,numElemsInRank,numNodesInRank,numBoundsInRank,elemGid(numElemsInRank),listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode)
      real(rp),intent(in) :: coordInRank(numNodesInRank,3)
      integer, intent(in) :: listElemsBoundsInRank(numElemsInRank,maxBoundsPerElem),listBoundFacesInRank(numBoundsInRank),boundFacesInRank(numBoundsInRank,npbou)
      
      integer, intent(out) :: numBoundNodesRankPar,numMshBoundNodesRankPar,numMpiBoundNodesRankPar,numInnerNodesRankPar
      integer, allocatable, intent(out) :: boundaryNodes(:),mpiBoundaryNodes(:),mshBoundaryNodes(:)
      
      integer :: nodeOwnedCnt(numNodesInRank)
      logical :: nodeInBoundary(numNodesInRank),nodeInMshBoundary(numNodesInRank),nodeInMpiBoundary(numNodesInRank)
      integer :: i,ind,ind_f,iElemL,iElemG,iNode,iNodeG,iNodeG_inFace,indexNode,indexNode_inFace
      integer :: auxCnt,mpiAuxCnt,mshAuxCnt
      integer :: iBound,numBoundFacesToCheck,boundFacesToCheck(maxBoundsPerElem)
      integer, parameter :: checkFacePos = 6
      !-------------------------------------------

      !$acc kernels
      nodeOwnedCnt(:)=0
      nodeInBoundary(:)=.false.
      nodeInMshBoundary(:)=.false.
      nodeInMpiBoundary(:)=.false.
      !listBoundsInRank(:)=boundFacesInRank(:,1)
      !$acc end kernels

      do iElemL=1,numElemsInRank
         do ind = 1,nnode
            !iNodeG = connecGMSH(iElemG,gmsh2ijk(ind))
            !nodeOwned(iNodeG) = nodeOwned(iNodeG) + 1
            iNodeG = connecInRank(iElemL,ind)
            indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
            nodeOwnedCnt(indexNode) = nodeOwnedCnt(indexNode) + 1
         end do
      end do

      !check if face is in boundary

      do iElemL=1,numElemsInRank

         numBoundFacesToCheck=0
         do iBound=1,maxBoundsPerElem
            if(listElemsBoundsInRank(iElemL,iBound).ne.0) then 
               numBoundFacesToCheck=numBoundFacesToCheck+1
               boundFacesToCheck(numBoundFacesToCheck) = listElemsBoundsInRank(iElemL,iBound)
            end if
         end do

         call new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceFront2ijk, listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceLeft2ijk,  listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceTop2ijk,   listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceBack2ijk,  listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceRight2ijk, listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceBottom2ijk,listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)

      end do

      numBoundNodesRankPar=0
      numMshBoundNodesRankPar=0
      numMpiBoundNodesRankPar=0
      numInnerNodesRankPar=0

      do iNode=1,numNodesInRank
         if(nodeInBoundary(iNode)) then !node is boundary
            numBoundNodesRankPar=numBoundNodesRankPar+1
            if(nodeInMshBoundary(iNode)) numMshBoundNodesRankPar=numMshBoundNodesRankPar+1
            if(nodeInMpiBoundary(iNode)) numMpiBoundNodesRankPar=numMpiBoundNodesRankPar+1
         else !node is inner
            numInnerNodesRankPar=numInnerNodesRankPar+1
         end if
      end do

      write(*,*) '1.#rank[',mpi_rank,'] numNodesInRank',numNodesInRank,'bN',numBoundNodesRankPar,'bMshN',numMshBoundNodesRankPar,'bMpiN',numMpiBoundNodesRankPar,'iN',numInnerNodesRankPar
      allocate(boundaryNodes(numBoundNodesRankPar))
      allocate(mpiBoundaryNodes(numMpiBoundNodesRankPar))
      allocate(mshBoundaryNodes(numMshBoundNodesRankPar))

      auxCnt=0
      mpiAuxCnt=0
      mshAuxCnt=0
      do iNode=1,numNodesInRank
         iNodeG=listNodesInRank(iNode)
         if(nodeInBoundary(iNode)) then !node is boundary
            auxCnt=auxCnt+1
            boundaryNodes(auxCnt) = iNodeG
            if(nodeInMshBoundary(iNode)) then
               mshAuxCnt=mshAuxCnt+1
               mshBoundaryNodes(mshAuxCnt) = iNodeG
            end if
            if(nodeInMpiBoundary(iNode)) then
               mpiAuxCnt=mpiAuxCnt+1
               mpiBoundaryNodes(mpiAuxCnt) = iNodeG
            end if
            !write(*,*) '1.#rank[',mpi_rank,']mshRank(',mshRank,')iNodeG',iNodeG,'auxCnt',auxCnt
         end if
      end do

     !i=0
     !do iNodeG=1,totalNumNodesSrl
     !   if(nodeInBoundary(iNodeG).eq.1) then !node is boundary
     !      i=i+1
     !      boundaryNodes(i) = iNodeG
     !   end if
     !end do
      !write(*,*) '#rank[',mpi_rank,']  nodesInRank ',numNodesRankPar, " bN ", numBNodesRankPar, " iN ",numINodesRankPar,'i',i

      !deallocate(nodeOwned)
      !deallocate(nodeInBoundary)

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      write(aux_string_mshrank_deb,'(I0)') mshRank
      file_name_deb = 'boundaryNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numBoundNodesRankPar
         iNodeG=boundaryNodes(iNode)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         
         write(1,fmt_csv) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

      file_name_deb = 'boundaryMpiNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numMpiBoundNodesRankPar
         iNodeG=mpiBoundaryNodes(iNode)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         
         write(1,fmt_csv) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

      file_name_deb = 'boundaryMshNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numMshBoundNodesRankPar
         iNodeG=mshBoundaryNodes(iNode)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         
         write(1,fmt_csv) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

#endif
   end subroutine get_rankPartitionBoundaryNodes_in_parallel

   subroutine new_subroutine_without_name(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,face2ijk,listNodesInRank,connecInRank,nodeOwnedCnt,listBoundsInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
      implicit none
      integer,intent(in) :: numElemsInRank,numNodesInRank,numBoundsInRank
      integer,intent(in) :: iElemL,face2ijk(npbou),listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode),nodeOwnedCnt(numNodesInRank),listBoundsInRank(numBoundsInRank),boundFacesInRank(numBoundsInRank,npbou)
      integer,intent(in) :: numBoundFacesToCheck,boundFacesToCheck(maxBoundsPerElem)
      logical,intent(out) :: nodeInBoundary(numNodesInRank),nodeInMshBoundary(numNodesInRank),nodeInMpiBoundary(numNodesInRank)
      integer :: ii,jj,ind,ind_f,iFace,iNodeG,iNodeG_inFace,iNodeG_ToCheck,indexNode,indexNode_inFace
      integer :: iBound,iBoundG,iElemG,indexBound
      logical :: isMshBound
      integer, parameter :: checkFacePos = 6
      
      !# 1.Whatever faceFront ------------------------------------------------------
      ind = face2ijk(checkFacePos)
      iNodeG_inFace = connecInRank(iElemL,gmsh2ijk(ind))
      indexNode_inFace = binarySearch_int_i(listNodesInRank,iNodeG_inFace)
      if(nodeOwnedCnt(indexNode_inFace).eq.1) then !this face is boundary

         !----------------------------------------------------------------
         isMshBound=.false.
         do iBound=1,numBoundFacesToCheck
            iBoundG = boundFacesToCheck(iBound)
            indexBound = binarySearch_int_i(listBoundsInRank,iBoundG)
            checkLoop: do jj=13,16
               iNodeG_ToCheck=boundFacesInRank(indexBound,jj)
               if(iNodeG_inFace.eq.iNodeG_ToCheck) then
                  isMshBound = .true.
                  exit checkLoop
               end if
               !write(*,*) '[',mpi_rank,']iElemL',iElemL,'check iBoundG',iBoundG,'jj',jj,'iNodeG2c',iNodeG_ToCheck,'iNodeGinF',iNodeG_inFace
            end do checkLoop
            !write(*,*) '[',mpi_rank,']iElemL',iElemL,'check iBoundG',iBoundG,'isMshBound',isMshBound
         end do
         !----------------------------------------------------------------
         !----------------------------------------------------------------
         do iFace=1,npbou
            ind_f = gmsh2ijk(face2ijk(iFace))
            iNodeG = connecInRank(iElemL,ind_f)
            indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
            nodeInBoundary(indexNode) = .true.
            if(isMshBound) then
               nodeInMshBoundary(indexNode) = .true.
            else 
               nodeInMpiBoundary(indexNode) = .true.
            end if
         end do
         !----------------------------------------------------------------
      end if

   end subroutine new_subroutine_without_name

   subroutine define_parallelNodePartitioning(numMshRanks2Part,iNodeStartPar)
      implicit none
      integer,intent(in) :: numMshRanks2Part
      integer,dimension(0:numMshRanks2Part-1),intent(out) :: iNodeStartPar
      integer,dimension(0:numMshRanks2Part-1) :: iNodeEndPar,vectorNumNodesRankPar2send,vectorNumNodesRankPar2rcv
      integer :: window_id
      integer :: iMshRank,mshRank,mpiRank,auxNodeCnt

      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      !allocate(vectorNumNodesRankPar(0:mpi_size-1))

      vectorNumNodesRankPar2send(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vectorNumNodesRankPar2send(mshRank) = numNodesRankPar(iMshRank)
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vectorNumNodesRankPar2send(0),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0

      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vectorNumNodesRankPar2rcv(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,&
                     1,MPI_INTEGER,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------

      !write(*,*) 'rank[',mpi_rank,']vNNR2s', vectorNumNodesRankPar2send(:),'vNNR2r',vectorNumNodesRankPar2rcv(:),'nNrp',numNodesRankPar(:)

      auxNodeCnt = 1
      do iMshRank=0,numMshRanks2Part-1
      !do iRank=0,mpi_size-1
         iNodeStartPar(iMshRank) = auxNodeCnt
         iNodeEndPar(iMshRank)   = iNodeStartPar(iMshRank) + vectorNumNodesRankPar2rcv(iMshRank) - 1
         auxNodeCnt = iNodeEndPar(iMshRank) + 1
      end do

      totalNumNodesPar = iNodeEndPar(numMshRanks2Part-1)

      allocate(rankNodeStart(numMshRanksInMpiRank))
      allocate(rankNodeEnd(numMshRanksInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         !vectorNumNodesRankPar(mshRank) = numNodesRankPar(iMshRank)
         rankNodeStart(iMshRank) = iNodeStartPar(mshRank)
         rankNodeEnd(iMshRank)   = iNodeEndPar(mshRank)
         write(*,*) 'mshrank',mshrank,'[',mpi_rank,']nodesInRank',numNodesRankPar(iMshRank),'iNodeS',rankNodeStart(iMshRank),'iNodeE',rankNodeEnd(iMshRank)
         !write(*,*) 'mshrank',mshrank,'[',mpi_rank,']start',iNodeStartPar(:),'end',iNodeEndPar(:),'tNNP',totalNumNodesPar
      end do
   end subroutine define_parallelNodePartitioning

   subroutine reorder_nodes_in_mshRank(mshRank,numMshRanks2Part,numElemsInRank,numNodesInRank,numBoundsInRank,elemGid,listNodesInRank,connecInRank,iNodeStartPar,&
               globalIdSrl,globalIdPar,connecVTK,connecParOrig,boundFacesInRank)
      implicit none
      integer,intent(in) :: mshRank,numMshRanks2Part,numElemsInRank,numNodesInRank,numBoundsInRank
      integer,intent(in) :: elemGid(numElemsInRank),listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode)
      integer,dimension(0:numMshRanks2Part-1),intent(in) :: iNodeStartPar

      integer,intent(out) :: globalIdSrl(numNodesInRank),globalIdPar(numNodesInRank),connecVTK(numElemsInRank*nnode),connecParOrig(numElemsInRank,nnode)
      integer,intent(inout) :: boundFacesInRank(numBoundsInRank,npbou)

      integer :: m,indConn,indexIJK,indexVTK,indexGMSH,indexNew,nodeIndexCnt,indPosListNodes
      integer :: iNodeL,iNodeGpar,iNodeGsrl,iElemL,iElemG,iNode,iBound

      integer,dimension(nnode) :: auxNodeNewOrderInElem,auxNewOrderIndex,auxVTKorder,auxCGNSorder
      integer,dimension(numNodesInRank) :: isNodeAdded
      integer,dimension(npbou) :: aux_ibound_GtoL
      integer,dimension(numNodesInRank,2) :: globalIdSrlOrdered

      !$acc kernels
      isNodeAdded(:)=-1
      globalIdSrlOrdered(:,:)=-1
      !$acc end kernels
      nodeIndexCnt = 0
      indConn = -1

      !@TODO SUPER VERY FUCKING IMPORTANT
      !--------------------------------------------------------------------------------------------------------
      ! now this only works using gmsh2ijk, because the old part of the code is using the atoijk set in 
      ! the subroutine set_hex64_list in mod_elem.f90 file
      ! atoijk is the same than in gmsh2ijk
      ! first do a first migration of the code and once it work, think how to allow different node ordering
      ! BUT FIX IT!!!!!!

      call generate_new_nodeOrder_and_connectivity(gmsh2ijk,auxNewOrderIndex,auxVTKorder,auxCGNSorder)
      !call generate_new_nodeOrder_and_connectivity(cgns2ijk,auxNewOrderIndex,auxCGNSorder)
      !call generate_new_nodeOrder_and_connectivity(dummy2ijk,auxNewOrderIndex,auxCGNSorder)

      !----------------------------------------------------------------------------------------------------------

      do iElemL=1,numElemsInRank
         !iElemG = (iElemL-1) + rankElemStart
         iElemG = elemGid(iElemL)
         auxNodeNewOrderInElem(:)=0

         do indexIJK=1,nnode
            indexGMSH = gmsh2ijk(indexIJK)
            iNodeGsrl = connecInRank(iElemL,indexGMSH)
            !iNodeGsrl = connecGMSH(iElemG,indexGMSH)

            indexNew = auxNewOrderIndex(indexIJK)
            
            auxNodeNewOrderInElem(indexNew) = iNodeGsrl
         end do

         do m=1,nnode
            iNodeGsrl = auxNodeNewOrderInElem(m)
            indPosListNodes = binarySearch_int_i(listNodesInRank,iNodeGsrl)

            if(isNodeAdded(indPosListNodes) < 0) then !node not added put it in the list
               nodeIndexCnt=nodeIndexCnt+1
               iNodeL = nodeIndexCnt

               isNodeAdded(indPosListNodes)  = iNodeL

               iNodeGPar = iNodeL + iNodeStartPar(mpi_rank) - 1

               globalIdSrl(iNodeL) = iNodeGsrl
               globalIdPar(iNodeL) = iNodeGPar

               globalIdSrlOrdered(iNodeL,1) = iNodeGsrl
               globalIdSrlOrdered(iNodeL,2) = iNodeL
            else 
               iNodeL = isNodeAdded(indPosListNodes)
               !iNodeGPar = globalIdPar(iNodeL)
            endif

            connecParOrig(iElemL,m) = iNodeL

            indexVTK = auxVTKorder(m)
            indConn = (iElemL-1)*nnode + indexVTK
            connecVTK(indConn) = iNodeL
         end do
         !write(*,*) '[',mpi_rank,']iElemG ',iElemG,' connecParOrig ',connecParOrig(iElemL,:)
      end do
      !if(mpi_rank.eq.2)write(*,*) 'check globalIdSrl[',mpi_rank,']',globalIdSrl(:)

      call quicksort_matrix_int(globalIdSrlOrdered,1)

      do iBound=1,numBoundsInRank
         aux_ibound_GtoL(:) = 0
         !write(*,*) '1.reorder ibound[',mpi_rank,']npbouG',boundFacesInRank(iBound,:)
         do iNode=1,npbou
            iNodeGsrl = boundFacesInRank(iBound,iNode)
            indexNew = binarySearch_int_i(globalIdSrlOrdered(:,1),iNodeGsrl)
            iNodeL = globalIdSrlOrdered(indexNew,2)
            !if(mpi_rank.eq.2)write(*,*) 'ibound[',mpi_rank,']',iBound,'iNodeGsrl',iNodeGsrl,'iNodeL',iNodeL
            aux_ibound_GtoL(iNode) = iNodeL
         end do

         boundFacesInRank(iBound,:) = aux_ibound_GtoL(:)
         !write(*,*) '2.reorder ibound[',mpi_rank,']npbouG',boundFacesInRank(iBound,:)
      end do

   end subroutine reorder_nodes_in_mshRank

   subroutine generate_new_nodeOrder_and_connectivity(newOrderIJK,auxNewOrderIndex,auxVTKorder,auxCGNSorder)
      integer,intent(in)  :: newOrderIJK(:)
      integer,intent(out) :: auxNewOrderIndex(:),auxVTKorder(:),auxCGNSorder(:)
      integer :: i,j,k,indexIJK,indexNew,indexVTK,indexCGNS

      !!!$acc kernels
      do k = 0,porder
         do i = 0,porder
            do j = 0,porder
               indexIJK = ((porder+1)**2)*k+(porder+1)*i+j+1

               indexNew = newOrderIJK(indexIJK)
               indexVTK = vtk2ijk(indexIJK)
               indexCGNS = cgns2ijk(indexIJK) !posicio requerida en el connec de cgns

               auxNewOrderIndex(indexIJK) = indexNew
               auxVTKorder(indexNew) = indexVTK
               auxCGNSorder(indexNew) = indexCGNS

               !write(*,*) 'test->indexIJK ', indexIJK, ' iNew ', indexNew,' aux ', auxCGNSorder(indexNew)
            end do
         end do
      end do
      !!!$acc end kernels
   end subroutine generate_new_nodeOrder_and_connectivity

   subroutine generate_mpi_comm_scheme_parallel(numMshRanks2Part,numMpiBoundaryNodes,mpiBoundaryNodes_jv)
      !generate a matrix with the comm schemes for shared nodes between procs
      integer,intent(in)  :: numMshRanks2Part,numMpiBoundaryNodes(numMshRanksInMpiRank)
      type(jagged_vector_int),intent(in) :: mpiBoundaryNodes_jv
      type(jagged_vector_int) :: mpiBoundaryNodesAll_jv
      integer, allocatable :: mpiBoundaryNodesAll(:)
      integer :: iMshRank,mshRank,iMshRankTrgt,mpiRank,iAux,memPos,memSize,memDisp
      integer :: totalNumMpiBoundaryNodes
      integer :: mshRankOrig,mshRankTrgt,numBNOrig,numBNTrgt
      integer :: iNode,iPos,iNodeL,iNodeGSrl
      integer, dimension(numMshRanks2Part) :: numMpiBoundaryNodesAll,memDispMpiBN
      integer, dimension(0:numMshRanks2Part-1) :: auxCommSchemeNumNodes
      integer,allocatable :: auxCommSchemeMemPos(:,:),auxCommSchemeMemPosAll(:)
      
      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      !------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) ' # Generating MPI Comm scheme...'

      !---------------------------------------------------------------------
      !1. sharing the number of mpi boundary nodes in all the msh ranks

      numMpiBoundaryNodesAll(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         numMpiBoundaryNodesAll(mshRank+1) = numMpiBoundaryNodes(iMshRank)
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(numMpiBoundaryNodesAll,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         !mshRank = iMshRank
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = iMshRank-1

         call MPI_Get(numMpiBoundaryNodesAll(iMshRank),1,MPI_INTEGER,mpiRank,target_displacement,&
                     1,MPI_INTEGER,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------


      !2. Allocating vector to store all the mpi boundary nodes of the domain 
      totalNumMpiBoundaryNodes=0
      memDispMpiBN(1)=0
      do iMshRank=1,numMshRanks2Part
         totalNumMpiBoundaryNodes = totalNumMpiBoundaryNodes + numMpiBoundaryNodesAll(iMshRank)
         if(iMshRank.gt.1) memDispMpiBN(iMshRank)=memDispMpiBN(iMshRank-1)+numMpiBoundaryNodesAll(iMshRank-1)
      end do
      !write(*,*) 'rank[',mpi_rank,']numMBNA', numMpiBoundaryNodesAll(:),'total',totalNumMpiBoundaryNodes,'memDisp',memDispMpiBN(:)

      ! fill my own mpi boundary nodes
      allocate(mpiBoundaryNodesAll(totalNumMpiBoundaryNodes))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         !write(*,*) 'rank[',mpi_rank,']iMshRank',iMshRank,'mshRank',mshRank,'memDisp',memDispMpiBN(mshRank+1)

         do iAux=1,numMpiBoundaryNodes(iMshRank)
            memPos = memDispMpiBN(mshRank+1) + iAux
            mpiBoundaryNodesAll(memPos) = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*totalNumMpiBoundaryNodes

      call MPI_Win_create(mpiBoundaryNodesAll,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      do iMshRank=1,numMshRanks2Part
         memSize = numMpiBoundaryNodesAll(iMshRank)
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         memPos = memDispMpiBN(iMshRank)+1
         target_displacement = memDispMpiBN(iMshRank)
         !write(*,*) '[',mpi_rank,']memSize',memSize,'mshRank',mshRank,'mpiRank',mpiRank,'td',target_displacement

         if(mpi_rank.ne.mpiRank) then
            call MPI_Get(mpiBoundaryNodesAll(memPos),memSize,MPI_INTEGER,mpiRank,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
         end if
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !---------------------------------------------------------------------------------------------------------

      allocate(mpiBoundaryNodesAll_jv%vector(numMshRanks2Part))
      do iMshRank=1,numMshRanks2Part
         numBNTrgt = numMpiBoundaryNodesAll(iMshRank)
         allocate(mpiBoundaryNodesAll_jv%vector(iMshRank)%elems(numBNTrgt))
         !memPos = memDispMpiBN(iMshRank)

         do iAux=1,numBNTrgt
            memPos = memDispMpiBN(iMshRank) + iAux
            !mpiBoundaryNodesAll(memPos) = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iAux)
            mpiBoundaryNodesAll_jv%vector(iMshRank)%elems(iAux) = mpiBoundaryNodesAll(memPos) 
         end do

      end do
      deallocate(mpiBoundaryNodesAll)

      !---------------------------------------------------------------------------------------------------------

      allocate(numNodesToComm(numMshRanksInMpiRank))
      allocate(numRanksWithComms(numMshRanksInMpiRank))
      allocate(matrixCommScheme_jm%matrix(numMshRanksInMpiRank))
      allocate(ranksToComm_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemSize_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemPosInLoc_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemPosInNgb_jv%vector(numMshRanksInMpiRank))

      allocate(auxCommSchemeMemPos(numMshRanksInMpiRank,numMshRanks2Part))
      allocate(auxCommSchemeMemPosAll(numMshRanks2Part**2))


      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numBNOrig = numMpiBoundaryNodes(iMshRank)

         numNodesToComm(iMshRank)=0
         auxCommSchemeNumNodes(:)=0
         auxCommSchemeMemPos(iMshRank,:)=0

         do iNode=1,numBNOrig
            iNodeGSrl = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)

            do iMshRankTrgt=1,numMshRanks2Part
               mshRankTrgt = iMshRankTrgt-1
               
               if(mshRankTrgt.ne.mshRankOrig) then
                  iPos = binarySearch_int_i(mpiBoundaryNodesAll_jv%vector(iMshRankTrgt)%elems,iNodeGSrl)
                  if(iPos.ne.0) then
                     numNodesToComm(iMshRank) = numNodesToComm(iMshRank) + 1
                     auxCommSchemeNumNodes(mshRankTrgt) = auxCommSchemeNumNodes(mshRankTrgt)+1

                     !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBNO',numBNOrig,'mshRankT',mshRankTrgt,'memPos',memPos,'numBNT',numBNTrgt
                     !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'mshRankT',mshRankTrgt,'iNGO',iNodeGOrig,'iNGT',iNodeGTrgt
                  end if

                  !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBNO',numBNOrig,'mshRankT',mshRankTrgt,'memPos',memPos,'numBNT',numBNTrgt
               end if
            end do
         end do
         
         ! setting numRanksWithsComms
         numRanksWithComms(iMshRank) = 0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               numRanksWithComms(iMshRank) = numRanksWithComms(iMshRank) + 1
            end if
         end do

         allocate(ranksToComm_jv%vector(iMshRank)%elems(numRanksWithComms(iMshRank)))
         allocate(commsMemSize_jv%vector(iMshRank)%elems(numRanksWithComms(iMshRank)))
         allocate(commsMemPosInLoc_jv%vector(iMshRank)%elems(numRanksWithComms(iMshRank)))
         allocate(commsMemPosInNgb_jv%vector(iMshRank)%elems(numRanksWithComms(iMshRank)))

         iAux=0
         memPos=1
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            auxCommSchemeMemPos(iMshRank,iMshRankTrgt)=memPos
            if(auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               iAux=iAux+1
               mshRankTrgt = iMshRankTrgt-1
               ranksToComm_jv%vector(iMshRank)%elems(iAux) = mshRankTrgt

               commsMemSize_jv%vector(iMshRank)%elems(iAux) = auxCommSchemeNumNodes(mshRankTrgt)
               commsMemPosInLoc_jv%vector(iMshRank)%elems(iAux) = memPos
               memPos=memPos+auxCommSchemeNumNodes(mshRankTrgt)
            end if
         end do

         !write(*,*) '[',mpi_rank,']mshRank',mashRankOrig,'csStartEnd->',commSchemeStartEndNodes(:)
         write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'csNumNodes->',auxCommSchemeNumNodes(:)
         write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'numRanksWC->',numRanksWithComms(iMshRank)
         write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'ranksToComm->',ranksToComm_jv%vector(iMshRank)%elems(:)
         write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInLoc->',commsMemPosInLoc_jv%vector(iMshRank)%elems(:)
         write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemSize->',commsMemSize_jv%vector(iMshRank)%elems(:)

         !write(*,*) 'rank[',mpi_rank,']mshRankOrig',mshRankOrig,'numNodesToComm',numNodesToComm(iMshRank),'numBNOrig',numBNOrig,'cSNN',auxCommSchemeNumNodes(:),&
         !          'numRanksWithComms',numRanksWithComms(iMshRank),'r2C',ranksToComm_jv%vector(iMshRank)%elems(:)

         allocate(matrixCommScheme_jm%matrix(iMshRank)%elems(numNodesToComm(iMshRank),3))
         iAux=0
         do iNode=1,numBNOrig
            iNodeGSrl = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)

            do iMshRankTrgt=1,numMshRanks2Part
               mshRankTrgt = iMshRankTrgt-1
               
               if(mshRankTrgt.ne.mshRankOrig) then
                  iPos = binarySearch_int_i(mpiBoundaryNodesAll_jv%vector(iMshRankTrgt)%elems,iNodeGSrl)
                  if(iPos.ne.0) then
                     iAux=iAux+1
                     
                     iNodeL = binarySearch_int_i(globalIdSrl_jv%vector(iMshRank)%elems,iNodeGSrl)
                     matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,1) = iNodeL
                     matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,2) = iNodeGSrl
                     matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,3) = mshRankTrgt
                  end if
               end if
            end do
         end do
      end do

      auxCommSchemeMemPosAll(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numMshRanks2Part
            iPos = mshRank*numMshRanks2Part+iMshRankTrgt
            auxCommSchemeMemPosAll(iPos) = auxCommSchemeMemPos(iMshRank,iMshRankTrgt)
         end do
      end do
      !write(*,*) '1.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*(numMshRanks2Part**2)

      call MPI_Win_create(auxCommSchemeMemPosAll,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         !mshRank = iMshRank
         mpiRank = mapMshRankToMpiRank(iMshRank)
         memSize = numMshRanks2Part
         memPos = (iMshRank-1)*numMshRanks2Part+1
         target_displacement = (iMshRank-1)*numMshRanks2Part

         call MPI_Get(auxCommSchemeMemPosAll(memPos),memSize,MPI_INTEGER,mpiRank,target_displacement,&
                     memSize,MPI_INTEGER,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      !write(*,*) '2.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numRanksWithComms(iMshRank)
            mshRankTrgt=ranksToComm_jv%vector(iMshRank)%elems(iMshRankTrgt)
            iPos = mshRankTrgt*numMshRanks2Part+mshRankOrig+1
            commsMemPosInNgb_jv%vector(iMshRank)%elems(iMshRankTrgt) = auxCommSchemeMemPosAll(iPos)
            !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'iMRT',iMshRankTrgt,'mRT',mshRankTrgt,'iPos',iPos
         end do
         write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInNgb->',commsMemPosInNgb_jv%vector(iMshRank)%elems(:)
      end do

      deallocate(auxCommSchemeMemPos)
      deallocate(auxCommSchemeMemPosAll)

   end subroutine generate_mpi_comm_scheme_parallel

   subroutine set_nodesCoordinates(mshRank,numElemsInRank,numNodesInRank,globalIdSrl,listNodesInRank,connecInRank,coordInRank,coordPar)
      implicit none
      integer, intent(in) :: mshRank,numElemsInRank,numNodesInRank,globalIdSrl(numNodesInRank),listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode)
      real(rp),intent(in) :: coordInRank(numNodesInRank,3)
      real(rp),intent(out) :: coordPar(numNodesInRank,3)
      integer :: iPos,iNodeL,iNodeGsrl

      if(mpi_rank.eq.0) write(*,*) ' # Creating nodes coordinates...'
      !------------------------------------------------------

      !- create the coordinate data for this process
      !!!$acc parallel loop
      do iNodeL=1,numNodesInRank
         iNodeGSrl=globalIdSrl(iNodeL)
         iPos = binarySearch_int_i(listNodesInRank,iNodeGSrl)
         !write (*,*) 'mshRank',mshRank,'iL', iNodeL,'iG', iNodeGSrl,'iPos',iPos,'[',coordInRank(iPos,1),']','[',coordInRank(iPos,2),']'
         coordPar(iNodeL,1) = coordInRank(iPos,1)
         coordPar(iNodeL,2) = coordInRank(iPos,2)
         coordPar(iNodeL,3) = coordInRank(iPos,3)
      end do
      !!!$acc end parallel loop

      call interpolateOriginalCoordinates(numElemsInRank,numNodesInRank,coordPar)

   end subroutine set_nodesCoordinates

   subroutine interpolateOriginalCoordinates(numElemsInRank,numNodesInRank,coordPar)
      integer,intent(in) :: numElemsInRank,numNodesInRank
      real(rp),intent(inout) :: coordPar(numNodesInRank,3)
      real(rp), allocatable :: aux_1(:,:)
      integer(4) :: iElem,iNode,idime
      
      write(*,*) 'REMEMBER TO FULLY IMPLEMENT THIS FUNC! NEED TO DEFINE THE Ngp_l'
      if(mpi_rank.eq.0) write(*,*) "--| Interpolating nodes coordinates..."
      allocate(aux_1(numNodesInRank,ndime))
      aux_1(:,:) = coordPar(:,:)
      do ielem = 1,numElemsInRank
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               !call var_interpolate(aux_1(connecParOrig(ielem,:),idime),Ngp_l(inode,:),coordPar(connecParOrig(ielem,inode),idime))
            end do
         end do
      end do
      deallocate(aux_1)
      !call overwrite_coordinates_hdf5()

   end subroutine interpolateOriginalCoordinates

   subroutine create_working_lists_parallel(isPeriodic,mshRank,numElemsInRank,numNodesInRank,globalIdSrl,connecParOrig,nPerInRank,masSlaInRank,connecParWork,numWorkingNodesInRank,workingNodesInRank)
      implicit none
      logical,intent(in) :: isPeriodic
      integer,intent(in) :: mshRank,numElemsInRank,numNodesInRank,globalIdSrl(numElemsInRank),connecParOrig(numElemsInRank,nnode)
      integer,intent(in) :: nPerInRank,masSlaInRank(nPerInRank,2)
      integer,intent(out) :: connecParWork(numElemsInRank,nnode),numWorkingNodesInRank
      integer,intent(inout),allocatable :: workingNodesInRank(:)

      integer :: iNodeL,iElem,iPer,iAux
      integer :: iNodeL_Per,iNodeL_Per_Pair
      integer, allocatable :: aux_workingNodesInRank(:)

      if(mpi_rank.eq.0) write(*,*) ' # Creating working lists...'

      !$acc kernels
      connecParWork(:,:) = connecParOrig(:,:)
      !$acc end kernels

      if(isPeriodic) then !do all stuff in case mesh is periodic
         write(*,*) 'TO BE CHECK! subroutine create_working_list_parallel() for periodic cases'
#if 1
         !----------------------------------------------------------------
         !-------------  CONNEC   ----------------------------------------

         !$acc kernels
         do iElem = 1,numElemsInRank
            do iAux = 1,nnode
               iNodeL = connecParWork(iElem,iAux)
               !iNodeG = globalIdSrl(iNodeL)
               do iPer = 1,nPerInRank
                  iNodeL_Per = masSlaInRank(iPer,2)
                  if (iNodeL .eq. iNodeL_Per) then
                     iNodeL_Per_Pair = masSlaInRank(iPer,1)
                     connecParWork(iElem,iAux) = iNodeL_Per_Pair
                  end if
               end do
            end do
         end do
         !$acc end kernels

         numWorkingNodesInRank = numNodesInRank - nPerInRank

         allocate(aux_workingNodesInRank(numNodesInRank))
         allocate(workingNodesInRank(numWorkingNodesInRank))

         !$acc parallel loop
         do iNodeL = 1,numNodesInRank
            aux_workingNodesInRank(iNodeL) = iNodeL
         end do
         !$acc end parallel loop

         do iPer = 1,nPerInRank
            iNodeL_Per = masSlaInRank(iPer,2)
            do iNodeL = 1,numNodesInRank
               !iNodeG = globalIdSrl(iNodeL)
               if (iNodeL_Per .eq. iNodeL) then
                  aux_workingNodesInRank(iNodeL) = 0
                  exit
               end if
            end do
         end do

         !TODO: GPU-PARALLELIZE THIS LOOP
         iAux = 0
         do iNodeL = 1,numNodesInRank
            if (aux_workingNodesInRank(iNodeL) .ne. 0) then
               iAux = iAux+1
               workingNodesInRank(iAux) = aux_workingNodesInRank(iNodeL)
               !write(*,*)'wNP(',iAux,')',workingNodesInRank(iAux)
            end if
         end do

         deallocate(aux_workingNodesInRank)
#endif
      else !non-periodic meshes
         numWorkingNodesInRank = numNodesInRank
         allocate(workingNodesInRank(numWorkingNodesInRank))

         !$acc parallel loop
         do iNodeL = 1,numWorkingNodesInRank
            workingNodesInRank(iNodeL) = iNodeL
         end do
         !$acc end parallel loop
      end if

   end subroutine create_working_lists_parallel

end module mod_meshConversorTool