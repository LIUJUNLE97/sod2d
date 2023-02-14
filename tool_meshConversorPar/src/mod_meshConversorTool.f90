module mod_meshConversorTool
   use mod_constants
   use mod_mpi
   use mod_utils
   use mod_gmsh_indices
   use elem_hex
   use quadrature_rules
   use mod_hdf5
   use iso_c_binding
   implicit none

!-----------------------------------   
#define _CHECK_ 0
#define int_size 4
!-----------------------------------   
#ifndef GEMPAINTERFACE
   interface
      subroutine gempa_do_partition(numElemsMpiRankToPart,numMshRanksToPart,x,y,z,weights,part) bind(c)
        import c_int, c_double
        integer(kind=c_int), intent(in), value :: numElemsMpiRankToPart
        integer(kind=c_int), intent(in), value :: numMshRanksToPart
        real(kind=c_double), intent(in) :: x(numElemsMpiRankToPart)
        real(kind=c_double), intent(in) :: y(numElemsMpiRankToPart)
        real(kind=c_double), intent(in) :: z(numElemsMpiRankToPart)
        integer(kind=c_int), intent(in) :: weights(numElemsMpiRankToPart)
        integer(kind=c_int), intent(out) :: part(numElemsMpiRankToPart)
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
integer, parameter :: fileId_geo = 99, fileId_bou = 98, fileId_dims = 97
! ################################################################################################
!     FOR DEBUG
character(128) :: file_name_deb,aux_string_mpirank_deb,aux_string_mshrank_deb
character(*), parameter :: fmt_csv_deb = '(1x,*(g0,","))'
! ################################################################################################

contains

   subroutine read_gmsh_h5_file_and_do_partitioning_in_parallel(gmsh_filePath,gmsh_fileName,mesh_h5_filePath,mesh_h5_fileName,numMshRanks2Part)
      implicit none
      character(500), intent(in)  :: gmsh_filePath,gmsh_fileName,mesh_h5_filePath,mesh_h5_fileName
      logical :: isPeriodic=.false.,isBoundaries=.false.
      integer(4),intent(in) :: numMshRanks2Part

      integer(4) :: numElemsMpiRank,numNodesMpiRank,numBoundsMpiRank,numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4), allocatable :: listElemsMpiRank(:),listNodesMpiRank(:),connecMpiRank(:,:)
      integer(4), allocatable :: linkedPerElemsSrl(:,:),listPerElemsSrl(:),masSlaNodesSrl(:,:)
      real(rp), allocatable :: coordNodesMpiRank(:,:)
      integer(4),allocatable :: numBoundaryNodesMshRank(:),numMpiBoundaryNodes(:),numMshBoundaryNodes(:),numInnerNodes(:)
      integer(4),allocatable :: numDoFMshRank(:)

      type(jagged_vector_int) :: listNodesMshRank_jv,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,boundaryNodes_jv,mshBoundaryNodes_jv,mpiBoundaryNodes_jv,dofNodes_jv
      type(jagged_matrix_int) :: connecMshRank_jm,listElemsBoundsMshRank_jm,connecOrigBoundFacesMshRank_jm,connecBoundFacesMshRank_jm,globalIdSrlOrdered_jm
      type(jagged_matrix_rp) :: coordMshRank_jm
      integer,dimension(0:numMshRanks2Part-1) :: iNodeStartPar
      integer,dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumNodesToCommMshRank,vecNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecBndNumMshRanksWithComms
      integer,dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank
      real(rp) :: Ngp_l(ngaus,nnode)
      integer(4) :: iMshRank,mshRank
      integer(hid_t) :: gmsh_h5_fileId, sod2dmsh_h5_fileId
      real(8),dimension(10) :: start_time,end_time,elapsed_time_r

      ! ################################################################################################
      ! ----------------- VARS for new Par mesh FORMAT -------------------------------------------------
      ! ################################################################################################

      integer(int_size) :: numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numPerFacesGmsh,numPerLinkedNodesGmsh,numNodesParTotal!,numPerNodes
      integer(int_size), allocatable :: numElemsMshRank(:),mshRankElemStart(:),mshRankElemEnd(:)
      integer(int_size), allocatable :: mshRankNodeStart(:), mshRankNodeEnd(:),numNodesMshRank(:)
      integer(int_size), allocatable :: numBoundFacesMshRank(:)
      integer(int_size) :: maxBoundCode

      type(jagged_vector_int) :: elemGidMshRank_jv,globalIdSrl_jv,globalIdPar_jv

      type(jagged_vector_int) :: connecVTK_jv
      type(jagged_matrix_int) :: connecParOrig_jm,connecParWork_jm

      type(jagged_matrix_rp) :: coordPar_jm

      type(jagged_vector_int) :: workingNodesPar_jv
      integer(int_size),allocatable :: numWorkingNodesMshRank(:)

      type(jagged_matrix_int) :: masSlaRankPar_jm
      integer(int_size),allocatable :: numPerNodesMshRank(:)

      integer(int_size) :: mshRankInMpiRankStart,mshRankInMpiRankEnd,numMshRanksInMpiRank,maxNumMshRanks
      integer(int_size), allocatable :: mshRanksInMpiRank(:),mapMshRankToMpiRank(:)

! ################################################################################################
! ------------------------ VARS for MPI COMMS ----------------------------------------------------
! ################################################################################################

      type(jagged_matrix_int) :: matrixCommScheme_jm
      type(jagged_vector_int) :: commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv
      integer(int_size),allocatable :: numNodesToCommMshRank(:),numMshRanksWithComms(:)

      type(jagged_matrix_int) :: bnd_matrixCommScheme_jm
      type(jagged_vector_int) :: bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv
      integer(int_size),allocatable :: bnd_numNodesToCommMshRank(:),bnd_numMshRanksWithComms(:)

! ################################################################################################
! ################################################################################################

      start_time(1) = MPI_Wtime()

      call init_hdf5_interface()

#if _CHECK_
       if(mpi_rank.eq.0) write(*,*) "Be Careful! _CHECK_ flag is ON! Debug files (*.csv) will be created!"     
#endif

      call open_gmsh_h5_file(gmsh_filePath,gmsh_fileName,gmsh_h5_fileId)

      call read_dims_gmsh_h5_file(gmsh_h5_fileId,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numPerFacesGmsh,numPerLinkedNodesGmsh,isPeriodic,isBoundaries)

      start_time(2) = MPI_Wtime()
      call read_elems_nodes_gmsh_h5_file_in_parallel(gmsh_h5_fileId,isPeriodic,numElemsGmsh,numNodesGmsh,numPerFacesGmsh,numPerLinkedNodesGmsh,&
                     numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank,connecMpiRank,coordNodesMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl)
      end_time(2) = MPI_Wtime()
      elapsed_time_r(2) = end_time(2) - start_time(2)
                     
      call distribute_ranks2Part_in_mpiRank(numMshRanks2Part,mshRankInMpiRankStart,mshRankInMpiRankEnd,numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank,mapMshRankToMpiRank)

      start_time(3) = MPI_Wtime()
      call do_element_partitioning_gempa_in_parallel(gmsh_h5_fileId,isPeriodic,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,&
                     numMshRanks2Part,numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank,mapMshRankToMpiRank,numElemsMpiRank,listElemsMpiRank,&
                     numNodesMpiRank,listNodesMpiRank,connecMpiRank,coordNodesMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl,&
                     numElemsMshRank,mshRankElemStart,mshRankElemEnd,numNodesMshRank,numBoundFacesMshRank,numPerNodesMshRank,maxBoundCode,&
                     elemGidMshRank_jv,listNodesMshRank_jv,connecMshRank_jm,listElemsBoundsMshRank_jm,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,&
                     connecOrigBoundFacesMshRank_jm,coordMshRank_jm,masSlaRankPar_jm)
      end_time(3) = MPI_Wtime()
      elapsed_time_r(3) = end_time(3) - start_time(3)

      if(isPeriodic) then
         !since now i have the masSlaRankPar_jm with global nodes ids i can deallocate the serial one
         deallocate(masSlaNodesSrl)
         numMasSlaNodesSrl=0
      end if

      !-----------------------------------------------------------------------
      allocate(boundaryNodes_jv%vector(numMshRanksInMpiRank))
      allocate(dofNodes_jv%vector(numMshRanksInMpiRank))
      allocate(mpiBoundaryNodes_jv%vector(numMshRanksInMpiRank))
      allocate(mshBoundaryNodes_jv%vector(numMshRanksInMpiRank))

      allocate(numMpiBoundaryNodes(numMshRanksInMpiRank))
      allocate(numMshBoundaryNodes(numMshRanksInMpiRank))
      allocate(numInnerNodes(numMshRanksInMpiRank))

      allocate(numBoundaryNodesMshRank(numMshRanksInMpiRank))
      allocate(numDoFMshRank(numMshRanksInMpiRank))

      allocate(mshRankNodeStart(numMshRanksInMpiRank))
      allocate(mshRankNodeEnd(numMshRanksInMpiRank))

      numMpiBoundaryNodes(:) = 0
      numMshBoundaryNodes(:) = 0
      numInnerNodes(:) = 0
      numBoundaryNodesMshRank(:) = 0
      numDoFMshRank(:) = 0
      mshRankNodeStart(:) = 0
      mshRankNodeEnd(:) = 0

      !-----------------------------------------------------------------------

      start_time(4) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Reading Msh Partitions Boundary Nodes"
      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank) 

         call get_rankPartitionBoundaryNodes_in_parallel(mshRank,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),&
                  listNodesMshRank_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems,&
                  listElemsBoundsMshRank_jm%matrix(iMshRank)%elems,listBoundFacesMshRank_jv%vector(iMshRank)%elems,connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,&
                  numPerNodesMshRank(iMshRank),masSlaRankPar_jm%matrix(iMshRank)%elems,&
                  numMshBoundaryNodes(iMshRank),numMpiBoundaryNodes(iMshRank),numInnerNodes(iMshRank),mpiBoundaryNodes_jv%vector(iMshRank)%elems,mshBoundaryNodes_jv%vector(iMshRank)%elems)
      end do
      end_time(4) = MPI_Wtime()
      elapsed_time_r(4) = end_time(4) - start_time(4)

      call define_parallelNodePartitioning(numMshRanks2Part,numMshRanksInMpiRank,numNodesMshRank,mshRanksInMpiRank,mapMshRankToMpiRank,mshRankNodeStart,mshRankNodeEnd,iNodeStartPar,numNodesParTotal)

      !--------------------------------------------------------------------------------------
      call evalShapeFunctions_Ngp_l(Ngp_l)

      allocate(globalIdSrl_jv%vector(numMshRanksInMpiRank))
      allocate(globalIdSrlOrdered_jm%matrix(numMshRanksInMpiRank))
      allocate(globalIdPar_jv%vector(numMshRanksInMpiRank))
      allocate(connecVTK_jv%vector(numMshRanksInMpiRank))
      allocate(connecParOrig_jm%matrix(numMshRanksInMpiRank))
      allocate(connecParWork_jm%matrix(numMshRanksInMpiRank))
      allocate(coordPar_jm%matrix(numMshRanksInMpiRank))
      allocate(numWorkingNodesMshRank(numMshRanksInMpiRank))
      allocate(workingNodesPar_jv%vector(numMshRanksInMpiRank))
      allocate(connecBoundFacesMshRank_jm%matrix(numMshRanksInMpiRank))

      start_time(5) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Reordering nodes, interpolating coordinates and creating working lists..."

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)

         allocate(globalIdSrl_jv%vector(iMshRank)%elems(numNodesMshRank(iMshRank)))
         allocate(globalIdSrlOrdered_jm%matrix(iMshRank)%elems(numNodesMshRank(iMshRank),2))
         allocate(globalIdPar_jv%vector(iMshRank)%elems(numNodesMshRank(iMshRank)))

         allocate(connecVTK_jv%vector(iMshRank)%elems(numElemsMshRank(iMshRank)*nnode))
         allocate(connecParOrig_jm%matrix(iMshRank)%elems(numElemsMshRank(iMshRank),nnode))
         allocate(connecParWork_jm%matrix(iMshRank)%elems(numElemsMshRank(iMshRank),nnode))

         allocate(coordPar_jm%matrix(iMshRank)%elems(numNodesMshRank(iMshRank),3))

         allocate(connecBoundFacesMshRank_jm%matrix(iMshRank)%elems(numBoundFacesMshRank(iMshRank),npbou))

         !$acc kernels
         globalIdSrl_jv%vector(iMshRank)%elems(:) = -1
         globalIdSrlOrdered_jm%matrix(iMshRank)%elems(:,:) = -1
         globalIdPar_jv%vector(iMshRank)%elems(:) = -1
         connecVTK_jv%vector(iMshRank)%elems(:) = 0
         connecParOrig_jm%matrix(iMshRank)%elems(:,:) = 0
         connecParWork_jm%matrix(iMshRank)%elems(:,:) = 0
         coordPar_jm%matrix(iMshRank)%elems(:,:) = 0.0_rp
         connecBoundFacesMshRank_jm%matrix(iMshRank)%elems(:,:) = 0
         !$acc end kernels

         !reordering nodes
         call reorder_nodes_in_mshRank(mshRank,numMshRanks2Part,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),numPerNodesMshRank(iMshRank),&
                           elemGidMshRank_jv%vector(iMshRank)%elems,listNodesMshRank_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,iNodeStartPar,&
                           globalIdSrl_jv%vector(iMshRank)%elems,globalIdSrlOrdered_jm%matrix(iMshRank)%elems,globalIdPar_jv%vector(iMshRank)%elems,&
                           connecVTK_jv%vector(iMshRank)%elems,connecParOrig_jm%matrix(iMshRank)%elems,&
                           connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,masSlaRankPar_jm%matrix(iMshRank)%elems)

         call set_nodesCoordinates(mshRank,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),globalIdSrl_jv%vector(iMshRank)%elems,&
            listNodesMshRank_jv%vector(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems,Ngp_l,connecParOrig_jm%matrix(iMshRank)%elems,&
            coordPar_jm%matrix(iMshRank)%elems)

         call create_working_lists_parallel(isPeriodic,mshRank,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),&
               globalIdSrl_jv%vector(iMshRank)%elems,connecParOrig_jm%matrix(iMshRank)%elems,numPerNodesMshRank(iMshRank),masSlaRankPar_jm%matrix(iMshRank)%elems,&
               connecParWork_jm%matrix(iMshRank)%elems,connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,connecBoundFacesMshRank_jm%matrix(iMshRank)%elems,&
               numWorkingNodesMshRank(iMshRank),workingNodesPar_jv%vector(iMshRank)%elems)

      end do

      end_time(5) = MPI_Wtime()
      elapsed_time_r(5) = end_time(5) - start_time(5)

      start_time(6) = MPI_Wtime()
      call generate_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numMpiBoundaryNodes,mpiBoundaryNodes_jv,globalIdSrl_jv,globalIdSrlOrdered_jm,&
                                    matrixCommScheme_jm,commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv,numNodesToCommMshRank,numMshRanksWithComms)

      end_time(6) = MPI_Wtime()
      elapsed_time_r(6) = end_time(6) - start_time(6)

      start_time(7) = MPI_Wtime()

      call generate_dof_and_boundary_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,&
               numNodesMshRank,numMshBoundaryNodes,mshBoundaryNodes_jv,numMpiBoundaryNodes,mpiBoundaryNodes_jv,globalIdSrl_jv,globalIdSrlOrdered_jm,&
               numBoundaryNodesMshRank,boundaryNodes_jv,numDoFMshRank,dofNodes_jv,&
               bnd_matrixCommScheme_jm,bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv,bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms)

      end_time(7) = mpi_wtime()
      elapsed_time_r(7) = end_time(7) - start_time(7)

      !-------------------------------------------------------------------------------------------

      deallocate(listElemsMpiRank)
      deallocate(listNodesMpiRank)
      !deallocate(listElemsBoundsMpiRank)
      deallocate(connecMpiRank)
      deallocate(coordNodesMpiRank)

      !----------------------------------------------------------------------------------------------
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numWorkingNodesMshRank,vecNumWorkingNodes)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numNodesToCommMshRank,vecNumNodesToCommMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numMshRanksWithComms,vecNumMshRanksWithComms)

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numBoundFacesMshRank,vecNumBoundFacesMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numDoFMshRank,vecNumDoFMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numBoundaryNodesMshRank,vecNumBoundaryNodesMshRank)

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,bnd_numNodesToCommMshRank,vecBndNumNodesToCommMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,bnd_numMshRanksWithComms,vecBndNumMshRanksWithComms)

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numPerNodesMshRank,vecNumPerNodesMshRank)

      !write(*,*) 'l1[',mpi_rank,']vnwn',vecNumWorkingNodes(:),'vnntcmr',vecNumNodesToCommMshRank(:),'vnmrwc',vecNumMshRanksWithComms(:)
      !write(*,*) 'l2[',mpi_rank,']vnbf',vecNumBoundFacesMshRank(:),'vndof',vecNumDoFMshRank(:),'vnbn',vecNumBoundaryNodesMshRank(:)

      !----------------------------------------------------------------------------------------------
      start_time(8) = MPI_Wtime()
      if(mpi_rank.eq.0) write(*,*) "--| Saving Partitioned HDF5 mesh file:",mesh_h5_fileName

      call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName,numMshRanks2Part)

      call create_hdf5_meshFile_from_tool(sod2dmsh_h5_fileId)
      call create_hdf5_groups_datasets_in_meshFile_from_tool(sod2dmsh_h5_fileId,isPeriodic,isBoundaries,numMshRanks2Part,numElemsGmsh,numNodesParTotal,&
               vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         call write_mshRank_data_in_hdf5_meshFile_from_tool(sod2dmsh_h5_fileId,mshRank,numMshRanks2Part,isPeriodic,isBoundaries,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numNodesParTotal,&
            numElemsMshRank(iMshRank),mshRankElemStart(iMshRank),mshRankElemEnd(iMshRank),mshRankNodeStart(iMshRank),&
            mshRankNodeEnd(iMshRank),numNodesMshRank(iMshRank),numWorkingNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),numBoundaryNodesMshRank(iMshRank),numDoFMshRank(iMshRank),maxBoundCode,&
            elemGidMshRank_jv%vector(iMshRank)%elems,globalIdSrl_jv%vector(iMshRank)%elems,globalIdPar_jv%vector(iMshRank)%elems,&
            connecVTK_jv%vector(iMshRank)%elems,connecParOrig_jm%matrix(iMshRank)%elems,connecParWork_jm%matrix(iMshRank)%elems,coordPar_jm%matrix(iMshRank)%elems,workingNodesPar_jv%vector(iMshRank)%elems,&
            boundaryNodes_jv%vector(iMshRank)%elems,dofNodes_jv%vector(iMshRank)%elems,boundFacesCodesMshRank_jv%vector(iMshRank)%elems,connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,connecBoundFacesMshRank_jm%matrix(iMshRank)%elems,&
            numPerNodesMshRank(iMshRank),masSlaRankPar_jm%matrix(iMshRank)%elems,&
            numNodesToCommMshRank(iMshRank),numMshRanksWithComms(iMshRank),matrixCommScheme_jm%matrix(iMshRank)%elems,commsMemPosInLoc_jv%vector(iMshRank)%elems,&
            commsMemSize_jv%vector(iMshRank)%elems,commsMemPosInNgb_jv%vector(iMshRank)%elems,ranksToComm_jv%vector(iMshRank)%elems,&
            bnd_numNodesToCommMshRank(iMshRank),bnd_numMshRanksWithComms(iMshRank),bnd_matrixCommScheme_jm%matrix(iMshRank)%elems,bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems,&
            bnd_commsMemSize_jv%vector(iMshRank)%elems,bnd_commsMemPosInNgb_jv%vector(iMshRank)%elems,bnd_ranksToComm_jv%vector(iMshRank)%elems,&
            vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank) 
      end do

      do iMshRank=(numMshRanksInMpiRank+1),maxNumMshRanks
         !write(*,*) 'FAKE-rank[',mpi_rank,']doing',iMshRank,'max',maxNumMshRanks
         call dummy_write_mshRank_data_in_hdf5_meshFile_from_tool(sod2dmsh_h5_fileId,numMshRanks2Part,isPeriodic,isBoundaries)
      end do


      call close_hdf5_meshFile_from_tool(sod2dmsh_h5_fileId)

      end_time(8) = mpi_wtime()
      elapsed_time_r(8) = end_time(8) - start_time(8)


      end_time(1) = MPI_Wtime()
      elapsed_time_r(1) = end_time(1) - start_time(1)
      if(mpi_rank.eq.0) then
         write(*,*) "--| Everything Done!"
         write(*,*) "------------------------------------------------"
         write(*,*) 'Total Mesh generation time',elapsed_time_r(1)
         write(*,*) ' -ReadElemsNodesAndPer',elapsed_time_r(2)
         write(*,*) ' -ElementPartition',elapsed_time_r(3)
         write(*,*) ' -RankBoundaries',elapsed_time_r(4)
         write(*,*) ' -ReOrder,Coords,Working',elapsed_time_r(5)
         write(*,*) ' -MpiCommScheme',elapsed_time_r(6)
         write(*,*) ' -DoF & BoundCommScheme',elapsed_time_r(7)
         write(*,*) ' -SavingHDF5',elapsed_time_r(8)
      end if

      call end_hdf5_interface()

   end subroutine read_gmsh_h5_file_and_do_partitioning_in_parallel
!---------------------------------------------------------------------------------------------------------------------------
   subroutine open_gmsh_h5_file(gmsh_h5_filePath,gmsh_h5_fileName,gmsh_h5_fileId)
      implicit none
      character(len=*), intent(in) :: gmsh_h5_filePath,gmsh_h5_fileName
      integer(hid_t), intent(out) :: gmsh_h5_fileId
      character(512) :: gmsh_h5_fullFileName
      integer(hid_t) :: file_id,plist_id
      integer :: h5err

      gmsh_h5_fullFileName = trim(adjustl(gmsh_h5_filePath))//trim(adjustl(gmsh_h5_fileName))//'.h5'

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      if(mpi_rank.eq.0) write(*,*) '# Opening gmsh h5 mesh file:', gmsh_h5_fullFileName

      call h5fopen_f(gmsh_h5_fullFileName, H5F_ACC_RDWR_F,gmsh_h5_fileId,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot open h5 meshfile ',trim(adjustl(gmsh_h5_fullFileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)


      !call h5fclose_f(file_id,h5err)
   end subroutine open_gmsh_h5_file

   subroutine read_dims_gmsh_h5_file(gmsh_h5_fileId,numElems,numNodes,numBoundFaces,numPerFaces,numPerLinks,isPeriodic,isBoundaries)
      implicit none     
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      integer(4), intent(out) :: numElems,numNodes,numBoundFaces,numPerFaces,numPerLinks
      logical, intent(out) :: isPeriodic,isBoundaries

      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hid_t) :: dtype
      integer :: ms_rank,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      dtype = H5T_NATIVE_INTEGER
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = 0
      allocate(aux_array(1))

      dsetname = '/dims/numNodes'
      call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numNodes=aux_array(1)

      dsetname = '/dims/numElements'
      call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numElems=aux_array(1)

      dsetname = '/dims/numBoundaryFaces'
      call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numBoundFaces=aux_array(1)

      dsetname = '/dims/numPeriodicFaces'
      call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numPerFaces=aux_array(1)

      dsetname = '/dims/numPeriodicLinks'
      call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numPerLinks=aux_array(1)

      if(numBoundFaces.ne.0) isBoundaries = .true.
      if(numPerFaces.ne.0) isPeriodic = .true.

      if(mpi_rank.eq.0) write(*,*) 'numNodes',numNodes,'numElems',numElems,'numBoundFaces',numBoundFaces,'numPerFaces',numPerFaces,'numPerLinks',numPerLinks

      deallocate(aux_array)
    
   end subroutine read_dims_gmsh_h5_file


   subroutine read_elem_connec_and_nodes_coords_from_gmsh_h5_file(gmsh_h5_fileId,numElemsSrl,numNodesSrl,numElemsInRank,listElemsRank,numNodesInRank,connecRank,listNodesRank,coordNodesRank)
      implicit none
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      integer,intent(in) :: numElemsSrl,numNodesSrl,numElemsInRank,listElemsRank(numElemsInRank)
      integer,intent(out) :: numNodesInRank,connecRank(numElemsInRank,nnode)
      integer,allocatable,intent(inout) :: listNodesRank(:)
      real(rp),allocatable,intent(inout) :: coordNodesRank(:,:)
      integer(4) :: iElem,iNode,iNodeG,prevNodeG,nodeCnt,iAux,jAux,iDim
      integer(4) :: elemToRead,nodeToRead
      !-----------------------------------------------------------
      character(128) :: dsetname
      integer,parameter :: ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer :: h5err
      !------------------------------------------------------------
      integer,allocatable :: rawNodeListRank(:)
      real(8), dimension(5) :: start_time,end_time,elapsed_time,elapsed_time_m

      !if(mpi_rank.eq.0)write(*,*) "--| Reading element table..."
      if(mpi_rank.eq.0) write(*,*) ' -numElems2read',numElemsInRank
      start_time(1) = MPI_Wtime()

      dsetname = '/connec'
      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      !haig de jugar amb el ratio fs_dims i ms_dims

      !call read_elem_connec_by_coords_selection(dset_id,fspace_id,plist_id,numElemsInRank,listElemsRank,connecRank)
      call read_elem_connec_by_chunks(dset_id,fspace_id,plist_id,numElemsSrl,numElemsInRank,listElemsRank,connecRank)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      !----------------------------------------------------------------------------------------------------------
      start_time(2) = MPI_Wtime()
      !----------------------------------------------------------------------------------------------------------
      ! REORDERING NODE LIST

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
      !if(mpi_rank.eq.0) write(*,*) 'listNodesRank',listNodesRank(:)

      deallocate(rawNodeListRank)

      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      start_time(3) = MPI_Wtime()
      ! Nodal coordinates section
      !------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0)write(*,*) "--| Reading coordinates..."

      dsetname = '/coords'

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call read_nodes_coords_by_coords_selection(dset_id,fspace_id,plist_id,numNodesInRank,listNodesRank,coordNodesRank)
      !call read_nodes_coords_by_chunks(dset_id,fspace_id,plist_id,numNodesSrl,numNodesInRank,listNodesRank,coordNodesRank)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(3) = MPI_Wtime()
      elapsed_time(3) = end_time(3) - start_time(3)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),3,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '   1.Read Elems',elapsed_time_m(1)
         write(*,*) '   2.Sort Nodes',elapsed_time_m(2)
         write(*,*) '   3.Read Coord',elapsed_time_m(3)
      end if

   end subroutine read_elem_connec_and_nodes_coords_from_gmsh_h5_file
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_elem_connec_by_chunks(dset_id,fspace_id,plist_id,numElemsSrl,numElemsInRank,listElemsRank,connecRank)
      implicit none
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer,intent(in) :: numElemsSrl,numElemsInRank,listElemsRank(numElemsInRank)
      integer,intent(out) :: connecRank(numElemsInRank,nnode)
      integer(4) :: nextElemToRead,elemCnt,iElem,iElemG
      integer :: maxRows2read,numChunks,iChunk
      !-----------------------------------------------------------
      integer,parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: mspace_id
      integer :: h5err
      integer(hid_t) :: dtype
      integer, allocatable :: vecChunks(:),auxConnec(:,:)
      !----------------------------------------------------------------------------------------------------------
      !TESTING ZONE READ ELEMENTS BY CHUNKS!
      maxRows2read = 20000

      numChunks = ceiling(real(numElemsSrl)/real(maxRows2read))
      if(mpi_rank.eq.0) write(*,*) ' -Reading elems by chunks | numChunks',numChunks,'maxRows2read',maxRows2read

      ms_dims(1) = nnode
      ms_dims(2) = 1
      ms_offset(1) = 0
      ms_offset(2) = 0
      dtype = H5T_NATIVE_INTEGER

      allocate(vecChunks(numChunks))
      call distribution_algorithm(numElemsSrl,numChunks,vecChunks)

      elemCnt=0
      iElemG=0
      nextElemToRead=listElemsRank(elemCnt+1)
      ms_offset(2) = 0
      do iChunk=1,numChunks
         ms_dims(1) = nnode
         ms_dims(2) = vecChunks(iChunk)
         allocate(auxConnec(ms_dims(1),ms_dims(2)))
         !write(*,*) 'iChunk',iChunk,'ms_dims',ms_dims,'ms_offset',ms_offset(2),'+',ms_offset(2) + ms_dims(2)!,'fsdims',fs_dims

         call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         call h5dread_f(dset_id,dtype,auxConnec,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5sclose_f(mspace_id,h5err)

         do iElem=1,ms_dims(2)
            iElemG=iElemG+1
            if(iElemG.eq.nextElemToRead) then
               elemCnt=elemCnt+1
               connecRank(elemCnt,1:nnode) = auxConnec(:,iElem)
               !write(*,*) 'connecRank[',mpi_rank,'](',elemCnt,')',connecRank(elemCnt,:)
               if(elemCnt.lt.numElemsInRank) nextElemToRead=listElemsRank(elemCnt+1)
            end if
         end do

         deallocate(auxConnec)
         ms_offset(2) = ms_offset(2) + ms_dims(2)
      end do
      !write(*,*) 'elemCnt',elemCnt,'numElemsRank',numElemsInRank

      deallocate(vecChunks)
      !----------------------------------------------------------------------------------------------------------

   end subroutine read_elem_connec_by_chunks
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_elem_connec_by_coords_selection(dset_id,fspace_id,plist_id,numElemsInRank,listElemsRank,connecRank)
      implicit none
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer,intent(in) :: numElemsInRank,listElemsRank(numElemsInRank)
      integer,intent(out) :: connecRank(numElemsInRank,nnode)
      !-----------------------------------------------------------------------
      integer, allocatable :: auxConnec(:,:)
      integer :: iAux,iElem,iNode,elemToRead
      !-----------------------------------------------------------------------
      integer,parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hid_t) :: mspace_id
      integer :: h5err
      integer(hid_t) :: dtype
      integer(hsize_t) :: ms_numElems
      integer(hsize_t),allocatable :: ms_coords(:,:)
      !-----------------------------------------------------------------------

      dtype = H5T_NATIVE_INTEGER

      ms_dims(1) = nnode
      ms_dims(2) = numElemsInRank
      ms_numElems = nnode*numElemsInRank
      !write(*,*) 'ms_numElems(',mpi_rank,')',ms_numElems
      allocate(ms_coords(ms_rank,nnode*numElemsInRank))
      allocate(auxConnec(nnode,numElemsInRank))

      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 

      !default
      ms_coords(1,:) = 1
      ms_coords(2,:) = 1

      iAux=0
      do iElem=1,numElemsInRank
         elemToRead=listElemsRank(iElem)
         do iNode=1,nnode
            iAux=iAux+1
            ms_coords(1,iAux) = iNode
            ms_coords(2,iAux) = elemToRead
         end do
      end do

      call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err) 
      call h5dread_f(dset_id,dtype,auxConnec,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      do iElem=1,numElemsInRank
         connecRank(iElem,:) = auxConnec(:,iElem)
      end do

      deallocate(auxConnec)
      deallocate(ms_coords)

      call h5sclose_f(mspace_id,h5err)

   end subroutine read_elem_connec_by_coords_selection
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_nodes_coords_by_chunks(dset_id,fspace_id,plist_id,numNodesSrl,numNodesInRank,listNodesRank,coordNodesRank)
      implicit none
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer,intent(in) :: numNodesSrl,numNodesInRank,listNodesRank(numNodesInRank)
      real(rp),intent(inout) :: coordNodesRank(:,:)
      !-----------------------------------------------------------------------
      real(8),allocatable :: auxCoords(:,:)
      integer(4) :: nextNodeToRead,nodeCnt,iNode,iNodeG
      !-----------------------------------------------------------------------
      integer,parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: mspace_id
      integer :: h5err
      integer(hid_t) :: dtype
      integer, allocatable :: vecChunks(:),auxConnec(:,:)
      !----------------------------------------------------------------------------------------------------------
      integer :: maxRows2read,numChunks,iChunk

      maxRows2read = 1000000

      numChunks = ceiling(real(numNodesSrl)/real(maxRows2read))
      if(mpi_rank.eq.0) write(*,*) ' -Reading nodescoords by chunks | numChunks',numChunks,'maxRows2read',maxRows2read

      ms_dims(1) = ndime
      ms_dims(2) = 1
      ms_offset(1) = 0
      ms_offset(2) = 0
      dtype = H5T_NATIVE_DOUBLE

      allocate(vecChunks(numChunks))
      call distribution_algorithm(numNodesSrl,numChunks,vecChunks)

      nodeCnt=0
      iNodeG=0
      nextNodeToRead=listNodesRank(nodeCnt+1)
      ms_offset(2) = 0
      do iChunk=1,numChunks
         ms_dims(1) = ndime
         ms_dims(2) = vecChunks(iChunk)
         allocate(auxCoords(ms_dims(1),ms_dims(2)))
         !write(*,*) 'iChunk',iChunk,'ms_dims',ms_dims,'ms_offset',ms_offset(2),'+',ms_offset(2) + ms_dims(2),'fsdims',fs_dims

         call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err)
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         call h5dread_f(dset_id,dtype,auxCoords,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5sclose_f(mspace_id,h5err)

         do iNode=1,ms_dims(2)
            iNodeG=iNodeG+1
            if(iNodeG.eq.nextNodeToRead) then
               nodeCnt=nodeCnt+1
               coordNodesRank(nodeCnt,1:ndime) = auxCoords(:,iNode)

               if(nodeCnt.lt.numNodesInRank) nextNodeToRead=listNodesRank(nodeCnt+1)
            end if
         end do

         deallocate(auxCoords)
         ms_offset(2) = ms_offset(2) + ms_dims(2)
      end do

      deallocate(vecChunks)
      !----------------------------------------------------------------------------------------------------------

   end subroutine read_nodes_coords_by_chunks
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_nodes_coords_by_coords_selection(dset_id,fspace_id,plist_id,numNodesInRank,listNodesRank,coordNodesRank)
      implicit none
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer,intent(in) :: numNodesInRank,listNodesRank(numNodesInRank)
      real(rp),intent(inout) :: coordNodesRank(:,:)
      !-----------------------------------------------------------------------
      real(8),allocatable :: auxCoords(:,:)
      integer :: iAux,iElem,iNode,iDim,nodeToRead
      !-----------------------------------------------------------------------
      integer,parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hid_t) :: mspace_id
      integer :: h5err
      integer(hid_t) :: dtype
      integer(hsize_t) :: ms_numElems
      integer(hsize_t),allocatable :: ms_coords(:,:)
      !-----------------------------------------------------------------------
      
      dtype = H5T_NATIVE_DOUBLE

      ms_dims(1) = ndime
      ms_dims(2) = numNodesInRank

      ms_numElems = ms_dims(1)*ms_dims(2)
      allocate(ms_coords(ms_rank,ms_dims(1)*ms_dims(2)))
      allocate(auxCoords(ms_dims(1),ms_dims(2)))

      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 

      !default
      ms_coords(1,:) = 1
      ms_coords(2,:) = 1

      iAux=0
      do iNode=1,numNodesInRank
         nodeToRead=listNodesRank(iNode)
         do iDim=1,ndime
            iAux=iAux+1
            ms_coords(1,iAux) = iDim
            ms_coords(2,iAux) = nodeToRead
         end do
      end do

      call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err) 
      call h5dread_f(dset_id,dtype,auxCoords,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      do iNode=1,numNodesInRank
         coordNodesRank(iNode,:) = auxCoords(:,iNode)
      end do

      deallocate(auxCoords)
      deallocate(ms_coords)

      call h5sclose_f(mspace_id,h5err)

   end subroutine read_nodes_coords_by_coords_selection
!-------------------------------------------------------------------------------------------------------------------------------------------

   subroutine read_periodic_faces_and_links_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,listElemsInRank,connecInRank,&
                        numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl)
      implicit none
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      integer(4),intent(in) :: numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,listElemsInRank(numElemsInRank),connecInRank(numElemsInRank,nnode)
      integer(4),intent(out) :: numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4),allocatable,intent(inout) :: linkedPerElemsSrl(:,:),listPerElemsSrl(:),masSlaNodesSrl(:,:)     
      integer :: numPerFacesInRank,h5err
      integer(4),allocatable :: listElemsPerFacesInRank(:,:),listPerFacesInRank(:),connecPerFacesInRank(:,:)

      character(128) :: dsetname2read

      real(8), dimension(2) :: start_time,end_time,elapsed_time,elapsed_time_m

      if(mpi_rank.eq.0) write(*,*) "--| Reading periodic boundaries..."

      start_time(1) = MPI_Wtime()
      dsetname2read = '/periodicFaces'
      call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,numPerFacesSrl,numElemsInRank,listElemsInRank,connecInRank,&
                                                                        numPerFacesInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank)
      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      if(mpi_rank.eq.0) write(*,*) "--| Reading periodic links..."

      start_time(2) = MPI_Wtime()
      call read_periodic_links_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,numPerFacesInRank,&
                     listElemsInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank,numLinkedPerElemsSrl,linkedPerElemsSrl,&
                     numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl)
      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '   1.Time ReadPerBounds',elapsed_time(1)
         write(*,*) '   2.Time ReadPerLinks',elapsed_time(2)
      end if
      !-----------------------------------------------

   end subroutine read_periodic_faces_and_links_from_gmsh_h5_file_in_parallel

   subroutine read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname,numFacesSrl,numElemsInRank,listElemsInRank,connecInRank,numFacesInRank,listElemsFacesInRank,listFacesInRank,connecFacesInRank)
      implicit none
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      character(*), intent(in) :: dsetname
      integer(4), intent(in) :: numFacesSrl,numElemsInRank,listElemsInRank(numElemsInRank),connecInRank(numElemsInRank,nnode)
      integer, intent(out) :: numFacesInRank
      integer, allocatable, intent(inout) :: listElemsFacesInRank(:,:),listFacesInRank(:),connecFacesInRank(:,:)

      integer(4) :: iFace,iChunk,iFaceG,numFacesToRead,iElem,ind_gmsh,iElemG,iAux,jAux
      integer(4) :: faces2readInChunk,maxFaces2read,iFaceNodes(npbou)
      integer(4),allocatable :: auxFacesInRank(:,:)
      integer(4) :: iNodeG_inFace,iNodeG_inElem,nodeCnt
      integer :: f_iNodeG(4),e_iNodeG(8) !for the vertex of squares of the face and the element
      logical :: vertexFound

      integer,parameter :: ms_rank=2,ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer :: h5err
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dtype

      real(8), dimension(3) :: start_time,end_time,elapsed_time,elapsed_time_m

      !for the moment I will use this auxiliar array, but if we see that is allocating too much memory
      !maybe we can try another strategy, without allocating an array of numBoundFaceSrl
      !but doing a double loop
      !1st loop: determine the bounds that are in the msh rank
      !2nd loop: store only the nodes in the rank
      !OR ANOTHER POSSIBLE STRATEGY:
      ! Determine a maxNumBoundsToRead and then read by chunks
      !

      allocate(listElemsFacesInRank(numElemsInRank,maxBoundsPerElem))
      listElemsFacesInRank(:,:) = 0

      !---------------------------------------------------
      numFacesInRank=0
      !do the reading by chunks-------
      maxFaces2read =1000

      dtype = H5T_NATIVE_INTEGER
      ms_dims(1) = npbou
      ms_dims(2) = 0
      ms_offset(1) = 0
      ms_offset(2) = 0
      !do iChunk=1,numChunks

      faces2readInChunk=numFacesSrl
      allocate(auxFacesInRank(npbou+1,faces2readInChunk))
      auxFacesInRank(:,:) = 0

      ms_dims(2) = faces2readInChunk

      start_time(1) = MPI_Wtime()

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      !call h5dread_f(dset_id,dtype,connecRank(iElem,:),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id)
      call h5dread_f(dset_id,dtype,auxFacesInRank(1:npbou,1:faces2readInChunk),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      start_time(2) = MPI_Wtime()
      !read the chunk
      iChunk=0
      do iFace=1,faces2readInChunk
         iFaceG=iFace !REVISAR AIXO QUAN IMPLEMENTI ELS CHUNKS A FULL
      !do iFace = 1,numFacesSrl
         !ms_offset(2) = iFace-1
         !call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,iFaceNodes)
         !-----------------------------------------------------
         iFaceNodes(:)=auxFacesInRank(1:npbou,iFace)
         !fill the corners of the face to check
         do iAux=1,4
            ind_gmsh = gmsh2ij(posFaceVertices(iAux))
            f_iNodeG(iAux) = iFaceNodes(ind_gmsh)
         end do
         !write(*,*) '[',mpi_rank,']iFace',iFace,' -> f_iNodeG ',f_iNodeG(:)

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
               numFacesInRank=numFacesInRank+1
               !write(*,*) '[',mpi_rank,']iFace',iFace,' -> elem ',iElemG
               !auxFacesInRank(numFacesInRank,1)=iFace
               !auxFacesInRank(numFacesInRank,2:npbou+1)=iFaceNodes(1:npbou)
               auxFacesInRank(npbou+1,iFace)=iFaceG!numFacesInRank! iFaceNodes(1:npbou)

               addLoop : do iAux = 1,maxBoundsPerElem
                  if(listElemsFacesInRank(iElem,iAux).eq.0) then
                     !write(*,*) '[',mpi_rank,'] adding iFace',iFace,'to elem',iElemG
                     listElemsFacesInRank(iElem,iAux) = iFace

                     exit addLoop
                  end if
               end do addLoop
            end if

         end do
      end do
      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)
      !write(*,*) 'numFacesInRank',numFacesInRank

      allocate(listFacesInRank(numFacesInRank))
      allocate(connecFacesInRank(numFacesInRank,npbou))

      !allocate(boundFacesCodesInRank(numBoundsInRank))
      !maxBoundCodeInRank = 0
      !do iFace=1,numFacesInRank
      start_time(3) = MPI_Wtime()

      iAux=0
      do iFace=1,faces2readInChunk
         iFaceG=auxFacesInRank(npbou+1,iFace)
         if(iFaceG.ne.0) then
            iAux=iAux+1
            listFacesInRank(iAux)     = iFaceG!auxFacesInRank(iFace,1)
            connecFacesInRank(iAux,:) = auxFacesInRank(1:npbou,iFace)
         end if
         !boundFacesCodesInRank(iBound) = auxBoundFacesInRank(iBound,2)
         !maxBoundCodeInRank = max(maxBoundCodeInRank,boundFacesCodesInRank(iBound))
         !write(*,*) '[',mpi_rank,']bfmpirank(',iBound,')',boundFacesInRank(iBound,:)
      end do
      end_time(3) = MPI_Wtime()
      elapsed_time(3) = end_time(3) - start_time(3)

      deallocate(auxFacesInRank)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),3,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '     1.Bounds(readHDF5)',elapsed_time_m(1)
         write(*,*) '     2.Bounds(linkElemAndFace)',elapsed_time_m(2)
         write(*,*) '     3.Bounds(genListAndConn)',elapsed_time_m(3)
      end if


   end subroutine read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel

   subroutine read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname,numBoundsInRank,listBoundFacesInRank,boundFacesCodesInRank,maxBoundCodeInRank)
      implicit none
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      character(*),intent(in) :: dsetname
      integer(4),intent(in) :: numBoundsInRank,listBoundFacesInRank(numBoundsInRank)
      integer(4),intent(out) :: maxBoundCodeInRank
      integer,allocatable,intent(inout) :: boundFacesCodesInRank(:)
      
      integer(4) :: iBound     
      !integer(4) :: boundCodeFaceId(1),boundId

      !-----------------------------------------------------------
      integer,parameter :: ms_rank=1,ds_rank=1
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer :: h5err
      integer(hsize_t) :: ms_numElems
      integer(hid_t) :: dtype
      integer(hsize_t),allocatable :: ms_coords(:,:)
      !------------------------------------------------------------

      allocate(boundFacesCodesInRank(numBoundsInRank))
      boundFacesCodesInRank(:) = 0
      !---------------------------------------------------------------
      !TEST ZONE FOR ELEMENT SELECTION!
      dtype = H5T_NATIVE_INTEGER

      if(numBoundsInRank.ne.0) then
         ms_dims(1) = numBoundsInRank
      else
         ms_dims(1) = 1 !dummy acces case numBoundsInRank=0
      end if
      ms_numElems = ms_dims(1)
      allocate(ms_coords(ms_rank,ms_dims(1)))

      !default
      ms_coords(1,:) = 1

      do iBound=1,numBoundsInRank
         ms_coords(1,iBound) = listBoundFacesInRank(iBound)
      end do

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err) 
      call h5dread_f(dset_id,dtype,boundFacesCodesInRank,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      deallocate(ms_coords)
      !deallocate(ms_coords1d)
      !write(*,*) 'fs_dims',fs_dims(:),'max_dims',fs_maxdims(:)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      !END TEST ZONE FOR ELEMENT SELECTION!
      !---------------------------------------------------------------
      
      maxBoundCodeInRank = 0
      do iBound = 1,numBoundsInRank
         maxBoundCodeInRank = max(maxBoundCodeInRank,boundFacesCodesInRank(iBound))
      end do
      !write(*,*) '(',mpi_rank,')maxBC',maxBoundCodeInRank,'numBiR',numBoundsInRank

   end subroutine read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel

   subroutine read_periodic_links_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,numPerFacesInRank,&
               listElemsInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank,numLinkedPerElemsSrl,linkedPerElemsAll,&
               numPerElemsSrl,listPerElemsAll,numMasSlaNodesSrl,masSlaNodesAll)
      implicit none
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      integer(4),intent(in) :: numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,numPerFacesInRank
      integer(4),intent(in) :: listElemsInRank(numElemsInRank),listElemsPerFacesInRank(numElemsInRank,maxBoundsPerElem)
      integer(4),intent(in) :: listPerFacesInRank(numPerFacesInRank),connecPerFacesInRank(numPerFacesInRank,npbou)
      integer(4),intent(out) :: numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4), allocatable,intent(inout) :: linkedPerElemsAll(:,:),listPerElemsAll(:),masSlaNodesAll(:,:)

      integer(4), allocatable ::listPerElemsInRank(:)!linkedPerFacesAll(:,:),

      integer(4), parameter :: numVert = 4
      character(128) :: dsetname
      integer(4), allocatable :: refNodesPerFacesInRank(:),perFaceElemInRank(:),refNodesPerFacesAll(:),perFaceElemAll(:),listPerFacesAll(:)!,mpiRankPerFacesAll(:)
      logical, allocatable :: masterPerFacesInRank(:),masterPerFacesAll(:)
      integer(4), dimension(0:mpi_size-1) :: vecNumPerFacesInRank,vecNumPerElemsInRank

      integer(4) :: iLink,iElem,iPerFace,iFaceG,iBound,ind_gmsh,iElemG,iPos,iPosFace
      !integer(4) :: perFaceId,perFaceIdM,perFaceIdS
      integer(4) :: iFaceMaster,iFaceSlave,iFaceGmstr,iFaceGslv,iNodeGM,iNodeGS,elemGidPerFaceM,elemGidPerFaceS
      integer(4) :: iAux,jAux,iVert,jVert
      integer(4) :: iNodeG,vertNodeCnt,linkedPerElemsCnt,numPerElemsInRank,iNodeGtoReplace
      logical :: vertNodeFound
      integer(4) :: pairPerLinkNodes(2)
      integer(4),allocatable :: matPerLinkNodes(:,:),matPerLinkNodesT(:,:)

      integer,parameter :: ms_rank=2,ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer :: h5err
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dtype
      !integer :: ms_rank,h5err
      !integer(hsize_t), dimension(2) :: ms_dims
      !integer(hssize_t), dimension(1) :: ms_offset

      integer(4) :: memPos,memSize,window_id1,window_id2,mpiRankTrgt,iMpiRank
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size1,win_buffer_size2,trgt_disp

      real(8), dimension(8) :: start_time,end_time,elapsed_time,elapsed_time_m
      !----------------------------------------------------------------------------------------------------------
      !  0. Allocate vectors

      allocate(refNodesPerFacesInRank(numPerFacesInRank*numVert))
      allocate(perFaceElemInRank(numPerFacesInRank))
      allocate(masterPerFacesInRank(numPerFacesInRank))
      allocate(refNodesPerFacesAll(numPerFacesSrl*numVert))
      allocate(perFaceElemAll(numPerFacesSrl))
      allocate(masterPerFacesAll(numPerFacesSrl))
      allocate(listPerFacesAll(numPerFacesSrl))
      !allocate(mpiRankPerFacesAll(numPerFacesSrl))

      !$acc kernels
      refNodesPerFacesInRank(:) = -1
      refNodesPerFacesAll(:) = -1
      perFaceElemInRank(:) = -1
      perFaceElemAll(:) = -1
      masterPerFacesInRank(:) = .true.
      masterPerFacesAll(:) = .true.
      !mpiRankPerFacesAll(:) = -1
      !$acc end kernels

      !----------------------------------------------------------------------------------------------------------

      call MPI_Allgather(numPerFacesInRank,1,MPI_INTEGER,vecNumPerFacesInRank,1,MPI_INTEGER,MPI_COMM_WORLD,mpi_err)
      !write(*,*) '2.vecNumPerFaces',vecNumPerFacesInRank(:)
      !write(*,*) '[',mpi_rank,']',listPerFacesInRank(:)
      !if(mpi_rank.eq.0) write(*,*) '1.refNodesPerFacesAll',refNodesPerFacesAll(:)
      !if(mpi_rank.eq.0) write(*,*) '1.masterPerFacesAll',masterPerFacesAll(:)

      !----------------------------------------------------------------------------------------------------------
      !  2. Second read the periodic links from h5 file
      if(mpi_rank.eq.0) write(*,*) "  |-> Reading periodic links from h5 file..."
      start_time(2) = MPI_Wtime()
      !IMPORTANT: for the moment I will this auxiliar matrix with the full size numPerLinkedNodesSrl for performance reasons
      !if in a future this can be a memory bound problem (which should not be, since the numPerLinkedNodes should not be very big compared with the inner nodes)
      !rethink the method, maybe do the reading line by line of the hdf5 for all the faces or do the method by chunks

      allocate(matPerLinkNodes(numPerLinkedNodesSrl,2))
      allocate(matPerLinkNodesT(2,numPerLinkedNodesSrl))

      dsetname = 'periodicLinks'

      dtype = H5T_NATIVE_INTEGER
      ms_dims(1) = 2
      ms_dims(2) = numPerLinkedNodesSrl
      ms_offset(1) = 0
      ms_offset(2) = 0

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      !call h5dread_f(dset_id,dtype,connecRank(iElem,:),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
      call h5dread_f(dset_id,dtype,matPerLinkNodesT(1:2,1:numPerLinkedNodesSrl),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      do iLink=1,numPerLinkedNodesSrl
         matPerLinkNodes(iLink,:) = matPerLinkNodesT(:,iLink)
      end do

      deallocate(matPerLinkNodesT)

      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      !----------------------------------------------------------------------------------------------------------
      !  3. Fill the refNodesPerFacesAll matrix
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling refNodesPerFacesAll..."
      start_time(3) = MPI_Wtime()
      call quicksort_matrix_int(matPerLinkNodes,1)

      do iPerFace=1,numPerFacesInRank
         !perFaceId = listPerFacesInRank(iPerFace)
         perFaceLoop : do iVert=1,numVert
            ind_gmsh = gmsh2ij(posFaceInnerNodes(iVert))
            iNodeG = connecPerFacesInRank(iPerFace,ind_gmsh)
            iAux = (iPerFace-1)*numVert + iVert
            iPos = binarySearch_int_i(matPerLinkNodes(:,1),iNodeG)
            !write(*,*) 'ipf',iPerFace,'iAux',iAux,'iNodeG',iNodeG,'iPos',iPos
            if(iPos.ne.0) then !found
               refNodesPerFacesInRank(iAux) = iNodeG
            else
               masterPerFacesInRank(iPerFace) = .false.
               exit perFaceLoop
            end if
         end do perFaceLoop
      end do 

      !write(*,*) '[',mpi_rank,']list',listPerFacesInRank(:),'masterPerFace',masterPerFacesInRank(:)
      !then, the one not founds in master, must be in the 'slave column'

      call quicksort_matrix_int(matPerLinkNodes,2)

      do iPerFace=1,numPerFacesInRank
         !perFaceId = listPerFacesInRank(iPerFace)
         if(not(masterPerFacesInRank(iPerFace))) then
            do iVert=1,numVert
               ind_gmsh = gmsh2ij(posFaceInnerNodes(iVert))
               iNodeG = connecPerFacesInRank(iPerFace,ind_gmsh)
               iAux = (iPerFace-1)*numVert + iVert
               iPos = binarySearch_int_i(matPerLinkNodes(:,2),iNodeG)
               refNodesPerFacesInRank(iAux) = matPerLinkNodes(iPos,1)
            end do
         end if
      end do
      end_time(3) = MPI_Wtime()
      elapsed_time(3) = end_time(3) - start_time(3)
      if(mpi_rank.eq.0) write(*,*) "  |-> Looking for and setting the periodic master nodes..."
      start_time(4) = MPI_Wtime()

#if 0
      !----------------------------------------------------------------------------------------------------------
      !  1. Share the mpi rank owning each periodic face
      start_time(1) = MPI_Wtime()
      !--------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesSrl
      call MPI_Win_create(mpiRankPerFacesAll,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      do mpiRankTrgt=0,mpi_size-1
         do iPerFace=1,numPerFacesInRank
            perFaceId = listPerFacesInRank(iPerFace)
            trgt_disp = perFaceId-1

            call MPI_Put(mpi_rank,1,MPI_INTEGER,mpiRankTrgt,trgt_disp,1,MPI_INTEGER,window_id1,mpi_err)
         end do
      end do

      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !write(*,*) '[',mpi_rank,']',mpiRankPerFacesAll(:)
      !--------------------------------------------------------------------------------------------------------------------
      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)
#endif

      !do iPerFace=1,numPerFacesInRank
      !   write(*,*) '[',mpi_rank,']perFace(',listPerFacesInRank(iPerFace),'own',ownedPerFacesInRank(iPerFace,:),'oth',otherPerFacesInRank(iPerFace,:),'M/S',masterPerFacesInRank(iPerFace)
      !end do

      do iLink=1,numPerLinkedNodesSrl
         iNodeG = matPerLinkNodes(iLink,1)
         call find_masterNode_recursive(iNodeG,numPerLinkedNodesSrl,matPerLinkNodes,iNodeGtoReplace)
         
         if(iNodeGtoReplace.ne.0) then
            matPerLinkNodes(iLink,1) = iNodeGtoReplace
         end if
      end do

      !do iLink=1,numPerLinkedNodesSrl
      !   write(*,*) 'n',matPerLinkNodes(iLink,:),'o',matPerLinkNodesC(iLink,:)
      !end do

      numMasSlaNodesSrl = numPerLinkedNodesSrl
      do iLink=1,numPerLinkedNodesSrl-1
         if((matPerLinkNodes(iLink,1).eq.matPerLinkNodes(iLink+1,1)).and.(matPerLinkNodes(iLink,2).eq.matPerLinkNodes(iLink+1,2))) then
            numMasSlaNodesSrl = numMasSlaNodesSrl-1
         end if
      end do
      !write(*,*) 'numMasSlaNodesSrl',numMasSlaNodesSrl,'numPLN',numPerLinkedNodesSrl

      allocate(masSlaNodesAll(numMasSlaNodesSrl,2))

      iAux  = 1 
      iLink = 1
      masSlaNodesAll(iAux,:) = matPerLinkNodes(iLink,:)
      !write(*,*) 'iAux',iAux,'msns',masSlaNodesAll(iAux,:)
      do iLink=2,numPerLinkedNodesSrl
         if(not((masSlaNodesAll(iAux,1).eq.matPerLinkNodes(iLink,1)).and.(masSlaNodesAll(iAux,2).eq.matPerLinkNodes(iLink,2)))) then
            iAux=iAux+1
            masSlaNodesAll(iAux,:) = matPerLinkNodes(iLink,:)
            !write(*,*) 'iAux',iAux,'msns',masSlaNodesAll(iAux,:)
         end if
      end do

      deallocate(matPerLinkNodes)
      !deallocate(matPerLinkNodesC)

      !-----------------------------------------------------------------------------------------------------------
      ! X. fill the vector perFaceElemInRank with the local elems
      end_time(4) = MPI_Wtime()
      elapsed_time(4) = end_time(4) - start_time(4)
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling vector perFaceElemInRank with local Elems..."
      start_time(5) = MPI_Wtime()

      numPerElemsInRank=0
      do iElem=1,numElemsInRank
         do iBound=1,maxBoundsPerElem
            iElemG = listElemsInRank(iElem)
            iFaceG = listElemsPerFacesInRank(iElem,iBound)
            iPosFace = binarySearch_int_i(listPerFacesInRank,iFaceG)
            if(iFaceG.ne.0) then
               perFaceElemInRank(iPosFace) = iElemG
               !perFaceElemAll(perFaceId) = iElemG
            end if
         end do
         if(listElemsPerFacesInRank(iElem,1).ne.0) then
            numPerElemsInRank=numPerElemsInRank+1
         end if
      end do

      call MPI_Allgather(numPerElemsInRank,1,MPI_INTEGER,vecNumPerElemsInRank,1,MPI_INTEGER,MPI_COMM_WORLD,mpi_err)
      
      numPerElemsSrl=0
      do iMpiRank=0,mpi_size-1
         numPerElemsSrl=numPerElemsSrl+vecNumPerElemsInRank(iMpiRank)
      end do

#if 0
      write(*,*) '[',mpi_rank,']nPFiR',numPerFacesInRank,'nPFS',numPerFacesSrl
      write(*,*) '[',mpi_rank,']nPEiR',numPerElemsInRank,'nPES',numPerElemsSrl
      call MPI_Barrier(MPI_COMM_WORLD,h5err)
      if(mpi_rank.eq.0) then
         write(*,*) 'list(:)',listPerFacesInRank
         write(*,*) 'perFaceElemInRank(:)',perFaceElemInRank
      end if
      call MPI_Barrier(MPI_COMM_WORLD,h5err)
      if(mpi_rank.eq.3) then
         write(*,*) 'list(:)',listPerFacesInRank
         write(*,*) 'perFaceElemInRank(:)',perFaceElemInRank
      end if
      call MPI_Barrier(MPI_COMM_WORLD,h5err)
      call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
#endif

      end_time(5) = MPI_Wtime()
      elapsed_time(5) = end_time(5) - start_time(5)
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling vectors listPerElems (Rank & All)..."
      start_time(6) = MPI_Wtime()
      !TODO: NEED TO IMPROVE THESE PART OF COMMS FOR LARGE DATA ANT LARGE NUM OF PROCS
      !NOT EFFICIENT AT ALL --> IS A BOTTLENECK!
      !--------------------------------------------------------------------------------

      allocate(listPerElemsInRank(numPerElemsInRank))
      allocate(listPerElemsAll(numPerElemsSrl))
      iAux=0
      do iElem=1,numElemsInRank
         if(listElemsPerFacesInRank(iElem,1).ne.0) then
            iElemG = listElemsInRank(iElem)
            iAux=iAux+1
            listPerElemsInRank(iAux) = iElemG
         end if
      end do
      !write(*,*) '[',mpi_rank,']lPEiR',listPerElemsInRank(:)
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerElemsInRank
      call MPI_Win_create(listPerElemsInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerElemsInRank(iMpiRank)
         call MPI_Get(listPerElemsAll(memPos),memSize,MPI_INTEGER,iMpiRank,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------

      deallocate(listPerElemsInRank)
      call quicksort_array_int(listPerElemsAll) !i think is not necessary... but just in case
      !write(*,*) '[',mpi_rank,']lPESrl',listPerElemsAll(:)

      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*(numPerFacesInRank*numVert)
      call MPI_Win_create(refNodesPerFacesInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)*numVert
         call MPI_Get(refNodesPerFacesAll(memPos),memSize,MPI_INTEGER,iMpiRank,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
#if 0
      win_buffer_size1 = mpi_integer_size*(numPerFacesSrl*numVert)
      call MPI_Win_create(refNodesPerFacesAll,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      do perFaceId=1,numPerFacesSrl
         !perFaceId = listPerFacesInRank(iPerFace)
         mpiRankTrgt = mpiRankPerFacesAll(perFaceId)

         if(mpi_rank .ne. mpiRankTrgt) then
            iAux = (perFaceId-1)*numVert! + iVert
            memPos = iAux+1
            trgt_disp = iAux
            memSize = numVert

            call MPI_Get(refNodesPerFacesAll(memPos),memSize,MPI_INTEGER,mpiRankTrgt,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)

         end if
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
#endif
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesInRank
      call MPI_Win_create(perFaceElemInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)
         call MPI_Get(perFaceElemAll(memPos),memSize,MPI_INTEGER,iMpiRank,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'perFaceElemAll',perFaceElemAll
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesInRank
      call MPI_Win_create(masterPerFacesInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)
         call MPI_Get(masterPerFacesAll(memPos),memSize,MPI_INTEGER,iMpiRank,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'masterPerFacesAll',masterPerFacesAll

      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesInRank
      call MPI_Win_create(listPerFacesInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)
         call MPI_Get(listPerFacesAll(memPos),memSize,MPI_INTEGER,iMpiRank,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'listPerFacesAll',listPerFacesAll

#if 0
      win_buffer_size1 = mpi_integer_size*(numPerFacesSrl)
      call MPI_Win_create(perFaceElemAll,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)
      if(mpi_rank.eq.0) write(*,*) "  |---> Comm5.b.."
      win_buffer_size2 = mpi_integer_size*(numPerFacesSrl)
      call MPI_Win_create(masterPerFacesAll,win_buffer_size2,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id2,mpi_err)
      call MPI_Win_fence(0,window_id2,mpi_err)

      do perFaceId=1,numPerFacesSrl
         mpiRankTrgt = mpiRankPerFacesAll(perFaceId)
         if(mpi_rank .ne. mpiRankTrgt) then
            memPos = perFaceId
            trgt_disp = perFaceId-1
            memSize = 1

            call MPI_Get(perFaceElemAll(memPos),memSize,MPI_INTEGER,mpiRankTrgt,trgt_disp,memSize,MPI_INTEGER,window_id1,mpi_err)
            call MPI_Get(masterPerFacesAll(memPos),memSize,MPI_INTEGER,mpiRankTrgt,trgt_disp,memSize,MPI_INTEGER,window_id2,mpi_err)
         end if
         
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      call MPI_Win_fence(0,window_id2,mpi_err)
      call MPI_Win_free(window_id2,mpi_err)
#endif
      !TODO: TILL HERE TO IMPROVE THIS PART
      !-------------------------------------------------
      !--------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) '2.refNodesPerFacesAll',refNodesPerFacesAll(:)
      !if(mpi_rank.eq.0) write(*,*) '2.masterPerFacesAll',masterPerFacesAll(:)
      !if(mpi_rank.eq.0) write(*,*) '2.perFaceElemAll',perFaceElemAll(:)

      end_time(6) = MPI_Wtime()
      elapsed_time(6) = end_time(6) - start_time(6)
      if(mpi_rank.eq.0) write(*,*) "  |-> Finding the matches: the linked periodic faces!..."
      start_time(7) = MPI_Wtime()
      !------------------------------------------------------------------------------------------
      ! finding the 'matches', ie the linked periodic faces!
      numLinkedPerElemsSrl = numPerFacesSrl/2
      !allocate(linkedPerFacesAll(numLinkedPerElemsSrl,2))
      allocate(linkedPerElemsAll(numLinkedPerElemsSrl,2))

      linkedPerElemsCnt=0
      do iFaceMaster=1,numPerFacesSrl
         if(masterPerFacesAll(iFaceMaster)) then 
            id2loop: do iFaceSlave=1,numPerFacesSrl
               if(not(masterPerFacesAll(iFaceSlave))) then
                  vertNodeCnt=0
                  iVertLoop: do iVert=1,numVert
                     iAux = (iFaceMaster-1)*numVert + iVert
                     vertNodeFound = .false.
                     iNodeGM = refNodesPerFacesAll(iAux)
                     jVertLoop: do jVert=1,numVert
                        jAux = (iFaceSlave-1)*numVert + jVert
                        iNodeGS = refNodesPerFacesAll(jAux)
                        if(iNodeGM.eq.iNodeGS) then
                           vertNodeFound = .true.
                           vertNodeCnt = vertNodeCnt+1
                           exit jVertLoop
                        end if
                     end do jVertLoop
                     if(not(vertNodeFound)) exit iVertLoop
                  end do iVertLoop
                  if(vertNodeCnt.eq.4) then
                     !if(mpi_rank.eq.0) write(*,*) '(M)pFId',perFaceIdM,'(S)pFId',perFaceIdS
                     linkedPerElemsCnt=linkedPerElemsCnt+1
                     !linkedPerFacesAll(linkedPerElemsCnt,1) = perFaceIdM
                     !linkedPerFacesAll(linkedPerElemsCnt,2) = perFaceIdS
                     iFaceGmstr = listPerFacesAll(iFaceMaster)
                     iFaceGslv  = listPerFacesAll(iFaceSlave)
                     !if(mpi_rank.eq.0) write(*,*) 'matching iFaceGmstr',iFaceGmstr,'iFaceGslv',iFaceGslv

                     linkedPerElemsAll(linkedPerElemsCnt,1) = perFaceElemAll(iFaceMaster)!elemGidPerFaceM
                     linkedPerElemsAll(linkedPerElemsCnt,2) = perFaceElemAll(iFaceSlave)!elemGidPerFaceS

                     exit id2loop
                  end if
               end if
            end do id2loop
         end if
      end do

      if(linkedPerElemsCnt.ne.numLinkedPerElemsSrl) then
         write(*,*) 'ERROR! Crashing in read_periodic_links_from_gmsh_h5_file_in_parallel! EXIT!'
         write(*,*) 'linkedPerElemsCnt',linkedPerElemsCnt,'NOT EQUAL TO numLinkedPerElemsSrl',numLinkedPerElemsSrl
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call quicksort_matrix_int(linkedPerElemsAll,1)

#if 0
      if(mpi_rank.eq.0) then
         do iPos=1,numLinkedPerElemsSrl
            !write(*,*) 'lPFS',linkedPerFacesAll(auxCnt,:)
            write(*,*) '0.lPES',linkedPerElemsAll(iPos,:)
         end do
      end if
      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
      if(mpi_rank.eq.1) then
         do iPos=1,numLinkedPerElemsSrl
            !write(*,*) 'lPFS',linkedPerFacesAll(auxCnt,:)
            write(*,*) '1.lPES',linkedPerElemsAll(iPos,:)
         end do
      end if
#endif

      !----------------------------------------------------------------------------------------------------------
      end_time(7) = MPI_Wtime()
      elapsed_time(7) = end_time(7) - start_time(7)
      if(mpi_rank.eq.0) write(*,*) "  |-> Done periodic links stuff!!..."

      deallocate(refNodesPerFacesAll)
      deallocate(listPerFacesAll)
      deallocate(perFaceElemAll)
      deallocate(masterPerFacesAll)

      elapsed_time(1) = end_time(1) - start_time(1)
      elapsed_time(2) = end_time(2) - start_time(2)
      elapsed_time(3) = end_time(3) - start_time(3)
      elapsed_time(4) = end_time(4) - start_time(4)
      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),7,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '     1.PerLinks(s1)',elapsed_time_m(1)
         write(*,*) '     2.PerLinks(readHDF5)',elapsed_time_m(2)
         write(*,*) '     3.refNodesPerFacesAll',elapsed_time_m(3)
         write(*,*) '     4.Doing PerMaster',elapsed_time_m(4)
         write(*,*) '     5.Vecs.perFaceElemAll',elapsed_time_m(5)
         write(*,*) '     6.Vecs.listPerElems',elapsed_time_m(6)
         write(*,*) '     7.FindingMatches',elapsed_time_m(7)
      end if

   end subroutine read_periodic_links_from_gmsh_h5_file_in_parallel

   subroutine unfold_linked_elems(numElems2PartInRank,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,elemPart)
      implicit none
      integer(4),intent(inout) :: numElems2PartInRank
      integer(4),intent(in) :: numLinkedPerElemsSrl,numPerElemsSrl
      integer(4),intent(in) :: linkedPerElemsSrl(numLinkedPerElemsSrl,2),listPerElemsSrl(numPerElemsSrl)
      integer,allocatable,intent(inout) :: elemPart(:,:)
      integer,allocatable :: elemPartAux(:,:)
      integer(4) :: numAdditionalElemsInRank,unfoldedElems(numPerElemsSrl)
      integer(4) :: iElem,iElemG,mshRank,i,j

      numAdditionalElemsInRank=0
      !allocate(unfoldedElems(numPerElemsSrl))
      unfoldedElems(:)=-1

      do iElem=1,numElems2PartInRank
         iElemG  = elemPart(iElem,1)
         mshRank = elemPart(iElem,2)

         call find_linked_elems_inPar(mshRank,iElemG,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,unfoldedElems,numAdditionalElemsInRank)
      end do

      !do iElemG=1,size(unfoldedElems)
      !   iRank = unfoldedElems(iElemG)
      !   if(iRank.ge.0) then
      !      write(*,*) '[',mpi_rank,'] unfolded Elem ',iElemG,' rank ', iRank
      !   end if
      !end do
      !write(*,*) 'unfolded[',mpi_rank,']',unfoldedElems(:)

      !now we have the 'unfolded elems' and the rank that they belong!
      !we need to add it to the elemPart matrix(:,2)!!!
      !write(*,*) 'numAdditional[',mpi_rank,']',numAdditionalElemsInRank

      allocate(elemPartAux(numElems2PartInRank,2))
      !$acc kernels
      elemPartAux(:,1)=elemPart(:,1)
      elemPartAux(:,2)=elemPart(:,2)
      !$acc end kernels

      deallocate(elemPart)
      numElems2PartInRank = numElems2PartInRank + numAdditionalElemsInRank
      allocate(elemPart(numElems2PartInRank,2))

      i=0
      do j=1,size(elemPartAux(:,1))
         i=i+1
         elemPart(i,:)=elemPartAux(i,:)
      end do

      deallocate(elemPartAux)

      do iElem=1,numPerElemsSrl
         mshRank = unfoldedElems(iElem)
         if(mshRank.ge.0) then
            i=i+1
            iElemG = listPerElemsSrl(iElem)
            elemPart(i,1) = iElemG
            elemPart(i,2) = mshRank
            !write(*,*) '[',mpi_rank,'] unf.Elem ',iElemG,' rank ', mshRank,' i ',i
         end if
      end do

   end subroutine unfold_linked_elems

   recursive subroutine find_linked_elems_inPar(mshRank,iElemG,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,unfoldedElems,numAddElemsInRank)
      implicit none
      integer, intent(in) :: mshRank,iElemG,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl(numLinkedPerElemsSrl,2),listPerElemsSrl(numPerElemsSrl)
      integer, intent(inout) :: unfoldedElems(numPerElemsSrl),numAddElemsInRank
      integer :: iPosEL1ref,iPosEL1ini,iPosEL1end,iPosEL1,iPosEL2,iElemGLink!,iEL2
      logical :: keepLooping         

      iPosEL1ref = binarySearch_int_i(linkedPerElemsSrl(:,1),iElemG)

      if(iPosEL1ref.ne.0) then   !#checking if the elem is in the linkedElems list
         iPosEL1ini=iPosEL1ref
         iPosEL1end=iPosEL1ref
         keepLooping = .true.
         do while(keepLooping)
            if(iPosEL1ini.eq.1) then
               keepLooping = .false.
            else
               if(iElemG.eq.linkedPerElemsSrl(iPosEL1ini-1,1)) then
                  iPosEL1ini = iPosEL1ini - 1
                  keepLooping = .true.
               else
                  keepLooping = .false.
               end if
            end if
         end do

         keepLooping = .true.
         do while(keepLooping)
            if(iPosEL1end.eq.numLinkedPerElemsSrl) then
               keepLooping = .false.
            else
               if(iElemG.eq.linkedPerElemsSrl(iPosEL1end+1,1)) then
                  iPosEL1end = iPosEL1end + 1
                  keepLooping = .true.
               else
                  keepLooping = .false.
               end if
            end if
         end do
         !write(*,*) '[',mpi_rank,']iELini',iPosEL1ini,'iELend',iPosEL1end

         do iPosEL1=iPosEL1ini,iPosEL1end
            iElemGLink = linkedPerElemsSrl(iPosEL1,2)
            !write(*,*) 'iElemG',iElemG,'iPosEl1',iPosEL1,'El',linkedPerElemsSrl(iPosEL1,1),'+1',linkedPerElemsSrl(iPosEL1+1,1)
            !write(*,*) 'iElemG',iElemG,'linked to', iElemGLink
         
            iPosEL2 = binarySearch_int_i(listPerElemsSrl,iElemGLink)
            if(unfoldedElems(iPosEL2).le.0) then
               unfoldedElems(iPosEL2) = mshRank
               numAddElemsInRank=numAddElemsInRank+1
               call find_linked_elems_inPar(mshRank,iElemGLink,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,unfoldedElems,numAddElemsInRank)
            end if        
         end do

      end if

   end subroutine find_linked_elems_inPar

   recursive subroutine find_masterNode_recursive(iNodeG,numPerLinkedNodesSrl,matPerLinkNodes,iNodeGtoReplace)
      implicit none
      integer(4),intent(in) :: iNodeG,numPerLinkedNodesSrl,matPerLinkNodes(numPerLinkedNodesSrl,2) !this matPerLinksNodes must be sorted by slave col
      integer(4),intent(out) :: iNodeGtoReplace
      integer(4) :: iPos,iPosNew,iNodeGnew

      iNodeGtoReplace = 0
      iPos = binarySearch_int_i(matPerLinkNodes(:,2),iNodeG)
      if(iPos.ne.0) then
         !if(mpi_rank.eq.0) write(*,*) 'iNG(M)',iNodeG,'is in slave',matPerLinkNodes(iPos,:)
         iNodeGnew = matPerLinkNodes(iPos,1)
         iPosNew   = binarySearch_int_i(matPerLinkNodes(:,2),iNodeGnew)
         if(iPosNew.ne.0) then
            !if(mpi_rank.eq.0) write(*,*) 'iNGnew(M)',iNodeGnew,'is in slave too!',matPerLinkNodes(iPosNew,:)
            call find_masterNode_recursive(iNodeGnew,numPerLinkedNodesSrl,matPerLinkNodes,iNodeGtoReplace)
         else
            iNodeGtoReplace = iNodeGnew
         end if
      end if

   end subroutine find_masterNode_recursive


   subroutine read_elems_nodes_gmsh_h5_file_in_parallel(gmsh_h5_fileId,isPeriodic,numElemsSrl,numNodesSrl,numPerFacesSrl,numPerLinkedNodesSrl,&
                     numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank,connecMpiRank,coordNodesMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl)
      implicit none
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      logical, intent(in)        :: isPeriodic
      integer(4), intent(in)     :: numElemsSrl,numNodesSrl,numPerFacesSrl,numPerLinkedNodesSrl
      integer(4), intent(out)  :: numElemsMpiRank,numNodesMpiRank,numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4), allocatable,intent(inout) :: listElemsMpiRank(:),listNodesMpiRank(:),connecMpiRank(:,:)
      integer(4), allocatable,intent(inout) :: linkedPerElemsSrl(:,:),listPerElemsSrl(:),masSlaNodesSrl(:,:)
      real(rp), allocatable, intent(inout) :: coordNodesMpiRank(:,:)

      integer(4)                 :: iElemMpiRankStart,iElemMpiRankEnd
      integer(4)                 :: auxCnt,iElem,h5err

      real(8),dimension(5) :: start_time,end_time,elapsed_time,elapsed_time_m


      !-----------------------------------------------------------
      !TESTING STUFF!!!
      real(4),allocatable :: rand_r(:)
      integer(4) :: iniRand,endRand
      character(128) :: dsetname
      integer,parameter :: ms_rank=2,ds_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer(hsize_t) :: ms_numElems
      integer(hid_t) :: dtype
      integer, allocatable :: vecChunks(:),auxConnec(:,:)
      integer :: maxRows2read,numChunks,iChunk
      !------------------------------------------------------------

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

      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
      start_time(1) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Reading element table and coordinates for Mpi Ranks"
      call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(gmsh_h5_fileId,numElemsSrl,numNodesSrl,numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,connecMpiRank,listNodesMpiRank,coordNodesMpiRank)
      end_time(1) = MPI_Wtime()

      !FOR TESTING!!!!
#if 0

      deallocate(listNodesMpiRank)
      deallocate(coordNodesMpiRank)

      allocate(rand_r(numElemsMpiRank))

      call random_number(rand_r)
      iniRand = 1
      endRand = numElemsSrl
      listElemsMpiRank(:) = iniRand + floor((endRand+1-iniRand)*rand_r(:))  ! We want to choose one from m-n+1 integers

      call quicksort_array_int(listElemsMpiRank)

      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
      start_time(3) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Testing 2... "
      call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(gmsh_h5_fileId,numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,connecMpiRank,listNodesMpiRank,coordNodesMpiRank)
      end_time(3) = MPI_Wtime()

      deallocate(rand_r(numElemsMpiRank))
      deallocate(listNodesMpiRank)
      deallocate(coordNodesMpiRank)

      auxCnt=0
      do iElem = iElemMpiRankStart,iElemMpiRankEnd
         auxCnt=auxCnt+1
         listElemsMpiRank(auxCnt) = iElem
      end do

      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
      start_time(4) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Testing 3... "
      call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(gmsh_h5_fileId,numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,connecMpiRank,listNodesMpiRank,coordNodesMpiRank)
      end_time(4) = MPI_Wtime()




      elapsed_time(1) = end_time(1) - start_time(1)
      elapsed_time(3) = end_time(3) - start_time(3)
      elapsed_time(4) = end_time(4) - start_time(4)
      elapsed_time(5) = end_time(5) - start_time(5)
      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),5,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '#OFst ReadElemNode',elapsed_time_m(1)
         write(*,*) '#Rand ReadElemNode',elapsed_time_m(3)
         write(*,*) '#OSnd ReadElemNode',elapsed_time_m(4)
         write(*,*) '#Chnk ReadElemNode',elapsed_time_m(5)
      end if


#endif



      !------------------------------------------------------------------------------------------
      ! Periodic faces

      start_time(2) = MPI_Wtime()
      if(isPeriodic) then
         call read_periodic_faces_and_links_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,numElemsMpiRank,listElemsMpiRank,connecMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl)
      end if
      end_time(2) = MPI_Wtime()

      elapsed_time(1) = end_time(1) - start_time(1)
      elapsed_time(2) = end_time(2) - start_time(2)
      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) ' 1.Time ReadElemNode',elapsed_time_m(1)
         write(*,*) ' 2.Time ReadPeriodic',elapsed_time_m(2)
      end if

   end subroutine read_elems_nodes_gmsh_h5_file_in_parallel


   subroutine do_element_partitioning_gempa_in_parallel(gmsh_h5_fileId,isPeriodic,numElemsSrl,numNodesSrl,numBoundFacesSrl,&
                  numMshRanks2Part,numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank,mapMshRankToMpiRank,numElemsMpiRank,listElemsMpiRank,&
                  numNodesMpiRank,listNodesMpiRank,connecMpiRank,coordNodesMpiRank,&
                  numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl,&
                  numElemsMshRank,mshRankElemStart,mshRankElemEnd,numNodesMshRank,numBoundFacesMshRank,numPerNodesMshRank,maxBoundCode,&
                  elemGid_jv,listNodesMshRank_jv,connecMshRank_jm,listElemsBoundsMshRank_jm,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,&
                  connecBoundFacesMshRank_jm,coordMshRank_jm,masSlaRankPar_jm)
      implicit none
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      logical,intent(in) :: isPeriodic
      integer(4),intent(in) :: numElemsSrl,numNodesSrl,numBoundFacesSrl
      integer(4),intent(in) :: numElemsMpiRank,numNodesMpiRank,numMshRanks2Part,numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4),intent(in) :: numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanksInMpiRank)
      integer(4),intent(in) :: listElemsMpiRank(numElemsMpiRank),listNodesMpiRank(numNodesMpiRank),connecMpiRank(numElemsMpiRank,nnode)
      integer(4),intent(in) :: linkedPerElemsSrl(numLinkedPerElemsSrl,2),listPerElemsSrl(numPerElemsSrl),masSlaNodesSrl(numMasSlaNodesSrl,2)
      real(rp),intent(in) :: coordNodesMpiRank(numNodesMpiRank,ndime)
      integer(4),intent(out) :: maxBoundCode
      integer(4),intent(out),allocatable :: numElemsMshRank(:),mshRankElemStart(:),mshRankElemEnd(:),numNodesMshRank(:),numBoundFacesMshRank(:),numPerNodesMshRank(:)
      type(jagged_vector_int),intent(out) :: elemGid_jv,listNodesMshRank_jv,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv
      type(jagged_matrix_int),intent(out) :: connecMshRank_jm,listElemsBoundsMshRank_jm,connecBoundFacesMshRank_jm,masSlaRankPar_jm
      type(jagged_matrix_rp),intent(out) :: coordMshRank_jm

      integer(4),allocatable :: dummyListElems(:),dummyListNodes(:),dummyListBoundFaces(:),dummyBoundFacesCodes(:)
      integer(4),allocatable :: dummyConnec(:,:),dummyListElemsBounds(:,:),dummyConnecBoundFaces(:,:)
      real(rp),allocatable :: dummyCoordNodes(:,:)
      integer(4) :: dummyNumNodes,dummyNumBoundFaces

      real(8), allocatable, dimension(:) :: x,y,z
      integer, allocatable :: listElems2Part(:),weightElems2Part(:),unfoldedElems(:)
      integer, allocatable :: linkedElems(:,:),elemPart(:,:),elemPartAux(:,:)

      integer :: ii,jj,kk,m,iElem,iElem2Part,iNode,iElemG,iNodeG,iBound,iPos
      integer :: mshRank,mpiRank,iMshRank,iMpiRank
      integer :: numElems2Part,numValues2get,aux_numElemsMshRank,aux_sum_NBP,aux_sum_NB_mpiRank
      integer :: maxBoundCodeMpiRank,maxBoundCodeMshRank
      character(128) :: dsetname2read

      real(8) :: x_a,y_a,z_a
      integer,dimension(0:numMshRanks2Part-1) :: numElems2mshRankInMpiRank,numElemsInMshRank
      integer,dimension(0:numMshRanks2Part-1,0:mpi_size-1) :: matNumElems2mshRankInMpiRank

      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size,target_displacement

      if(mpi_rank.eq.0) write(*,*) '--| Doing element partitioning GEMPA in parallel...'

      !1. obtain the list element 2 part using gempa
      if(isPeriodic) then
         call get_listElems2Par_periodic_in_parallel(numElemsMpiRank,numLinkedPerElemsSrl,numPerElemsSrl,listElemsMpiRank,linkedPerElemsSrl,listElems2Part,weightElems2Part,numElems2Part)
      else
         call get_listElems2Par_in_parallel(numElemsMpiRank,listElemsMpiRank,listElems2Part,weightElems2Part,numElems2Part)
      end if

      allocate(elemPart(numElems2Part,4))
      allocate(x(numElems2Part))
      allocate(y(numElems2Part))
      allocate(z(numElems2Part))

      do iElem2Part=1,numElems2Part
         elemPart(iElem2Part,1) = listElems2Part(iElem2Part)
         elemPart(iElem2Part,3) = weightElems2Part(iElem2Part)
         iElemG = listElems2Part(iElem2Part)
         iElem = binarySearch_int_i(listElemsMpiRank,iElemG)
         elemPart(iElem2Part,4) = iElem

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

         x(iElem2Part) = x_a/8.0d0
         y(iElem2Part) = y_a/8.0d0
         z(iElem2Part) = z_a/8.0d0
         !write(*,*) '[',mpi_rank,']iElemG',elemPart(iElem2Part,1),'w',elemPart(iElem2Part,3),'x',x(iElem2Part),'y',y(iElem2Part),'z',z(iElem2Part)
         !write(*,*) '[',mpi_rank,']iElemG',elemPart(iElem2Part,1),'z',z(iElem2Part)
      end do

      deallocate(listElems2Part)

#ifndef GEMPAINTERFACE
      !---- CALLING GEMPA in PARALLEL-----------------------
      call gempa_do_partition(numElems2Part,numMshRanks2Part,x,y,z,elemPart(:,3),elemPart(:,2))
      !---------------------------------------------------------------------
#else
      write(*,*) 'FATAL ERROR! Trying to call do_element_partitioning_gempa() for program compiled with flag GEMPAINTERFACE!'
      call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
#endif

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'gempaPartitionInfo_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,weight,iElem'
      do ii=1,numElems2Part
         write(1,fmt_csv_deb) x(ii),y(ii),z(ii),elemPart(ii,1),elemPart(ii,2),elemPart(ii,3),elemPart(ii,4)
      end do
      close(1)
#endif

      deallocate(x)
      deallocate(y)
      deallocate(z)

      if(isPeriodic) then
         ! if parallel, now we have to put the rank to the 'slave' elements who were not included in the partitioning process
         ! since they were linked to a 'master' element
         call unfold_linked_elems(numElems2Part,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,elemPart)
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

#if _CHECK_
      if(not(isPeriodic)) then
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'elemPartition_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,iElem'
      do ii=1,numElems2Part
         iElemG   = elemPart(ii,1)
         mshRank  = elemPart(ii,2)-1 
         !iElem    = elemPart(ii,4)

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

         write(1,fmt_csv_deb) x_a,y_a,z_a,iElemG,mshRank,iElem
      end do
      close(1)
      end if
#endif

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(numElems2mshRankInMpiRank,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
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
      if(not(isPeriodic)) then
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

         write(1,fmt_csv_deb) x_a,y_a,z_a,iElemG,mshRank,iElem
      end do
      close(1)
      end if
#endif



      !--------------------------------------------------------------------------------------

      numElemsInMshRank(:)=0

      do mpiRank=0,mpi_size-1
         do mshRank=0,numMshRanks2Part-1
            numElemsInMshRank(mshRank) = numElemsInMshRank(mshRank) + matNumElems2mshRankInMpiRank(mshRank,mpiRank)
         end do
      end do

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
      win_buffer_size = mpi_integer_size*numElems2Part

      call MPI_Win_create(elemPart(:,1),win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
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
            !write(1,fmt_csv_deb) iElem,elemGid(iElem,iMshRank)
            write(1,fmt_csv_deb) iElem,elemGid_jv%vector(iMshRank)%elems(iElem)
         end do
         close(1)
      end do
#endif
!---------------------------------------------------------------------------------
      deallocate(elemPart)

      allocate(numElemsMshRank(numMshRanksInMpiRank))
      allocate(mshRankElemStart(numMshRanksInMpiRank))
      allocate(mshRankElemEnd(numMshRanksInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank= mshRanksInMpiRank(iMshRank)
         numElemsMshRank(iMshRank) = numElemsInMshRank(mshRank)
      end do
      !write(*,*) 'numElemsMshRank[',mpi_rank,']',numElemsMshRank(:)

      mshRankElemStart(:)=1
      do iMshRank=1,numMshRanksInMpiRank
         mshRank= mshRanksInMpiRank(iMshRank)
         do jj=0,mshRank-1
            mshRankElemStart(iMshRank)=mshRankElemStart(iMshRank)+numElemsInMshRank(jj)
         end do
         mshRankElemEnd(iMshRank) = mshRankElemStart(iMshRank) + numElemsInMshRank(jj)-1
      end do
      !write(*,*) 'numElemsMshRank[',mpi_rank,']',numElemsMshRank(:),'eStart',mshRankElemStart(:),'eEnd ',mshRankElemEnd(:)

!------------------------------------------------------------------------------

      allocate(numNodesMshRank(numMshRanksInMpiRank))
      allocate(numBoundFacesMshRank(numMshRanksInMpiRank))
      allocate(listNodesMshRank_jv%vector(numMshRanksInMpiRank))
      allocate(connecMshRank_jm%matrix(numMshRanksInMpiRank))!,boundMshRank_jm,bou_codesMshRank_jm
      allocate(coordMshRank_jm%matrix(numMshRanksInMpiRank))
      
      allocate(numPerNodesMshRank(numMshRanksInMpiRank))
      allocate(masSlaRankPar_jm%matrix(numMshRanksInMpiRank))

      !if(numBoundFacesSrl.ne.0) then
      allocate(listBoundFacesMshRank_jv%vector(numMshRanksInMpiRank))
      allocate(boundFacesCodesMshRank_jv%vector(numMshRanksInMpiRank))
      allocate(listElemsBoundsMshRank_jm%matrix(numMshRanksInMpiRank))
      allocate(connecBoundFacesMshRank_jm%matrix(numMshRanksInMpiRank))
      !end if

      numNodesMshRank(:)  = 0
      numBoundFacesMshRank(:) = 0
      numPerNodesMshRank(:) = 0

      maxBoundCodeMshRank = 0
      maxBoundCodeMpiRank = 0
      maxBoundCode=0

      if(mpi_rank.eq.0) then
         write(*,*) '--| Reading element table and coordinates for Msh Ranks | numMshRanks(0)',numMshRanksInMpiRank,'max',maxNumMshRanks
         write(*,*) '  |-> numMshRanks(rank=0)',numMshRanksInMpiRank,'max',maxNumMshRanks
      end if

      do iMshRank=1,numMshRanksInMpiRank
         !write(*,*) 'REAL-rank[',mpi_rank,']doing',iMshRank,'max',maxNumMshRanks
         aux_numElemsMshRank=numElemsMshRank(iMshRank)
         mshRank= mshRanksInMpiRank(iMshRank)
         !write(*,*) '#Reading stuff for mshRank',mshRank,'in mpi_rank',mpi_rank,'(iMshRank',iMshRank,')','aux_numElemsMshRank',aux_numElemsMshRank

         allocate(connecMshRank_jm%matrix(iMshRank)%elems(aux_numElemsMshRank,nnode))

         call quicksort_array_int(elemGid_jv%vector(iMshRank)%elems)

         call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(gmsh_h5_fileId,numElemsSrl,numNodesSrl,aux_numElemsMshRank,elemGid_jv%vector(iMshRank)%elems,numNodesMshRank(iMshRank),&
                                          connecMshRank_jm%matrix(iMshRank)%elems,listNodesMshRank_jv%vector(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems)



         if(numBoundFacesSrl.ne.0) then

            dsetname2read = '/boundFaces'
            call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,numBoundFacesSrl,aux_numElemsMshRank,&
                     elemGid_jv%vector(iMshRank)%elems,connecMshRank_jm%matrix(iMshRank)%elems,&
                     numBoundFacesMshRank(iMshRank),listElemsBoundsMshRank_jm%matrix(iMshRank)%elems,listBoundFacesMshRank_jv%vector(iMshRank)%elems,connecBoundFacesMshRank_jm%matrix(iMshRank)%elems)

            dsetname2read = '/boundFacesId'
            call read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,numBoundFacesMshRank(iMshRank),&
                     listBoundFacesMshRank_jv%vector(iMshRank)%elems,boundFacesCodesMshRank_jv%vector(iMshRank)%elems,maxBoundCodeMshRank)

            maxBoundCodeMpiRank = max(maxBoundCodeMpiRank,maxBoundCodeMshRank)
            !if(mpi_rank.eq.0) write(*,*) ' maxBoundCode',maxBoundCodeMpiRank
         else
            allocate(listBoundFacesMshRank_jv%vector(iMshRank)%elems(numBoundFacesMshRank(iMshRank)))
            allocate(boundFacesCodesMshRank_jv%vector(iMshRank)%elems(numBoundFacesMshRank(iMshRank)))
            allocate(connecBoundFacesMshRank_jm%matrix(iMshRank)%elems(numBoundFacesMshRank(iMshRank),npbou))
            allocate(listElemsBoundsMshRank_jm%matrix(iMshRank)%elems(numElemsMshRank(iMshRank),maxBoundsPerElem))
            listElemsBoundsMshRank_jm%matrix(iMshRank)%elems(:,:) = 0
         end if

      end do

      aux_numElemsMshRank=1
      allocate(dummyListElems(1))
      dummyListElems(1)=1
      allocate(dummyConnec(1,nnode))

      do iMshRank=(numMshRanksInMpiRank+1),maxNumMshRanks
         !write(*,*) 'FAKE-rank[',mpi_rank,']doing',iMshRank,'max',maxNumMshRanks

         call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(gmsh_h5_fileId,numElemsSrl,numNodesSrl,aux_numElemsMshRank,dummyListElems,dummyNumNodes,&
                                          dummyConnec,dummyListNodes,dummyCoordNodes)

         if(numBoundFacesSrl.ne.0) then
            dsetname2read = '/boundFaces'
            call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,numBoundFacesSrl,aux_numElemsMshRank,&
                     dummyListElems,dummyConnec,dummyNumBoundFaces,dummyListElemsBounds,dummyListBoundFaces,dummyConnecBoundFaces)

            dsetname2read = '/boundFacesId'
            call read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,dummyNumBoundFaces,&
                     dummyListBoundFaces,dummyBoundFacesCodes,maxBoundCodeMshRank)

            deallocate(dummyListElemsBounds)
            deallocate(dummyListBoundFaces)
            deallocate(dummyConnecBoundFaces)
            deallocate(dummyBoundFacesCodes)
         end if

         deallocate(dummyListNodes)
         deallocate(dummyCoordNodes)

      end do
      deallocate(dummyListElems)
      deallocate(dummyConnec)

      call MPI_Allreduce(maxBoundCodeMpiRank,maxBoundCode,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD,mpi_err)

      !fer el check aqui de que s'han llegit correctament totes les boundaries
      if(numBoundFacesSrl.ne.0) then
         aux_sum_NB_mpiRank=0
         do iMshRank=1,numMshRanksInMpiRank
            aux_sum_NB_mpiRank = aux_sum_NB_mpiRank + numBoundFacesMshRank(iMshRank)
         end do

         call MPI_Allreduce(aux_sum_NB_mpiRank,aux_sum_NBP,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,mpi_err)

         if(aux_sum_NBP.ne.numBoundFacesSrl) then
            write(*,*) '[',mpi_rank,']Error in do_element_partitioning_gempa_in_parallel(..) after read_boundaries(..): aux_sum_NBP',aux_sum_NBP,' not equal to  numBoundFacesSrl',numBoundFacesSrl
            call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
         end if
      end if

      !create here masSlaNodesMshRank

      if(isPeriodic) then
         do iMshRank=1,numMshRanksInMpiRank
            mshRank= mshRanksInMpiRank(iMshRank)
            call generate_masSlaRankPar(mshRank,numNodesMshRank(iMshRank),listNodesMshRank_jv%vector(iMshRank)%elems,numMasSlaNodesSrl,masSlaNodesSrl,&
                        numPerNodesMshRank(iMshRank),masSlaRankPar_jm%matrix(iMshRank)%elems)
         end do

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         write(aux_string_mshrank_deb,'(I0)') mshRank
         file_name_deb = 'masSla_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
         open(1,file=file_name_deb)
         write(1,*) 'X,Y,Z,iNodeG,MS'
         do ii=1,numPerNodesMshRank(iMshRank)
            iNodeG=masSlaRankPar_jm%matrix(iMshRank)%elems(ii,1)
            iPos = binarySearch_int_i(listNodesMshRank_jv%vector(iMshRank)%elems,iNodeG)
            x_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,1)
            y_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,2)
            z_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,3)
            write(1,fmt_csv_deb) x_a,y_a,z_a,iNodeG,'1'
         end do
         do ii=1,numPerNodesMshRank(iMshRank)
            iNodeG=masSlaRankPar_jm%matrix(iMshRank)%elems(ii,2)
            iPos = binarySearch_int_i(listNodesMshRank_jv%vector(iMshRank)%elems,iNodeG)
            x_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,1)
            y_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,2)
            z_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,3)
            write(1,fmt_csv_deb) x_a,y_a,z_a,iNodeG,'2'
         end do
         close(1)
      end do
#endif
      end if

   end subroutine do_element_partitioning_gempa_in_parallel
!------------------------------------------------------------------------------------------------

   subroutine generate_masSlaRankPar(mshRank,numNodesInRank,listNodesInRank,numMasSlaNodesSrl,masSlaNodesSrl,numMasSlaNodesInRank,masSlaNodesInRank)
      implicit none
      integer,intent(in) :: mshRank,numNodesInRank,listNodesInRank(numNodesInRank),numMasSlaNodesSrl,masSlaNodesSrl(numMasSlaNodesSrl,2)
      integer,intent(out) :: numMasSlaNodesInRank
      integer,allocatable,intent(out) :: masSlaNodesInRank(:,:)
      integer :: iLink,iNodeGM,iNodeGS,iPos,iNode,iAux
      integer :: numPerNodesInRank

      numMasSlaNodesInRank = 0
      do iLink = 1,numMasSlaNodesSrl
         !iNodeGM = masSlaNodesSrl(iLink,1)
         iNodeGS = masSlaNodesSrl(iLink,2)

         !iPosM = binarySearch_int_i(listNodesInRank,iNodeGM)
         !if(iPosM.ne.0) nodeIsPerInRank(iPosM) = .true.

         iPos = binarySearch_int_i(listNodesInRank,iNodeGS)
         if(iPos.ne.0) numMasSlaNodesInRank = numMasSlaNodesInRank + 1! nodeIsPerInRank(iPosS) = .true.
      end do

      allocate(masSlaNodesInRank(numMasSlaNodesInRank,2))
      iAux = 0
      do iLink = 1,numMasSlaNodesSrl
         iNodeGM = masSlaNodesSrl(iLink,1)
         iNodeGS = masSlaNodesSrl(iLink,2)

         !iPosM = binarySearch_int_i(listNodesInRank,iNodeGM)
         !if(iPosM.ne.0) nodeIsPerInRank(iPosM) = .true.

         iPos = binarySearch_int_i(listNodesInRank,iNodeGS)
         if(iPos.ne.0) then
            iAux=iAux+1
            masSlaNodesInRank(iAux,1) = iNodeGM
            masSlaNodesInRank(iAux,2) = iNodeGS
            !write(*,*) 'mshRank',mshRank,'MS',masSlaNodesInRank(iAux,:)
         end if
      end do

      !write(*,*) 'numMasSlaRank(',mshRank,')',numMasSlaNodesInRank,'numMasSlaSrl',numMasSlaNodesSrl

   end subroutine generate_masSlaRankPar

   subroutine do_element_partitioning_serial(numElems2Par,iElemStart,iElemEnd,iElemsInRank)
      implicit none
      integer, intent(in) :: numElems2Par
      integer, intent(out) :: iElemStart,iElemEnd,iElemsInRank
      integer, dimension(0:mpi_size-1) :: vecElemsInMpiRank
      integer :: iMpiRank

      if(numElems2Par.lt.mpi_size) then
         write(*,*) 'ERROR! The total number of elements to par:',numElems2Par,'is bigger than the number of cpus used to do the partition',mpi_size,&
                     'this is CRAZY BRO! The tool is going to crash! Use a reasonable number of CPUs (smaller than num of elements)'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call distribution_algorithm(numElems2Par,mpi_size,vecElemsInMpiRank)

      iElemsInRank = vecElemsInMpiRank(mpi_rank)
      iElemStart=1
      do iMpiRank=0,(mpi_rank-1) !find the iElemStart
         iElemStart = iElemStart + vecElemsInMpiRank(iMpiRank)
      end do
      iElemEnd = iElemStart + (iElemsInRank - 1)

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

   subroutine get_listElems2Par_periodic_in_parallel(numElemsMpiRank,numPerLinkedPerFacesSrl,numPerElemsSrl,listElemsMpiRank,linkedPerElemsSrl,listElems2Part,weightElems2Part,numElems2Part)
      implicit none
      integer(4), intent(in)  :: numElemsMpiRank,numPerLinkedPerFacesSrl,numPerElemsSrl,listElemsMpiRank(numElemsMpiRank),linkedPerElemsSrl(numPerLinkedPerFacesSrl,2)
      integer, allocatable, intent(out) :: listElems2Part(:),weightElems2Part(:)
      integer, intent(out) :: numElems2Part
      logical :: isElem2PartMpiRank(numElemsMpiRank)

      integer :: iElem,iElemG,iPos,iPerElem,auxCnt

      isElem2PartMpiRank(:) = .true.

      do iElem=1,numPerLinkedPerFacesSrl
         iElemG = linkedPerElemsSrl(iElem,2)

         iPos = binarySearch_int_i(listElemsMpiRank,iElemG)
         !write(*,*) '[',mpi_rank,']iElem',iElem,'iElemG',iElemG,'iPos',iPos
         if(iPos.ne.0) then !the element is in the slave rank, so not to part
            isElem2PartMpiRank(iPos) = .false.
         end if
      end do

      numElems2Part=0
      do iElem=1,numElemsMpiRank
         if(isElem2PartMpiRank(iElem)) then
            numElems2Part=numElems2Part+1
         end if
      end do

      allocate(listElems2Part(numElems2Part))
      allocate(weightElems2Part(numElems2Part))

      auxCnt=0
      do iElem=1,numElemsMpiRank
         if(isElem2PartMpiRank(iElem)) then
            auxCnt=auxCnt+1
            listElems2Part(auxCnt)   = listElemsMpiRank(iElem)
            weightElems2Part(auxCnt) = 1
            do iPerElem=1,numPerLinkedPerFacesSrl
               iElemG = linkedPerElemsSrl(iPerElem,1)
               if(listElems2Part(auxCnt).eq.iElemG) then
                  weightElems2Part(auxCnt) = weightElems2Part(auxCnt) + 1
               end if
            end do
         end if
      end do

      !if(mpi_rank.eq.1) then
      !   write(*,*) 'listElems2Part',listElems2Part(:)
      !   write(*,*) 'weightElems2Part',weightElems2Part(:)
      !end if
      !write(*,*) '[',mpi_rank,']numElems2Part',numElems2Part,'auxCnt',auxCnt

   end subroutine get_listElems2Par_periodic_in_parallel

   subroutine distribute_ranks2Part_in_mpiRank(numMshRanks2Part,iMshRankStart,iMshRankEnd,numRanksMpiRank,maxNumRanksMpiRank,ranksMpiRank,mapRanksToMpiRank)
      implicit none
      integer, intent(in) :: numMshRanks2Part
      integer, intent(out) :: iMshRankStart,iMshRankEnd,numRanksMpiRank,maxNumRanksMpiRank
      integer, allocatable, intent(inout) :: ranksMpiRank(:),mapRanksToMpiRank(:)
      integer, dimension(0:mpi_size-1) :: vecRanksMpiRank
      integer :: iMpiRank,iMshRankCnt,ii

      maxNumRanksMpiRank=1
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
            maxNumRanksMpiRank = max(maxNumRanksMpiRank,vecRanksMpiRank(iMpiRank))
         end do

      else
         !possible first approach...
         !to use ONLY the number of CPUs=numMshRanks2Part and the others do nothing...
         !lets try

         if(mpi_rank.eq.0) then
            write(*,*) 'Tool not optimized for this case numMshRanks2Part<mpi_size -> ',numMshRanks2Part,'<',mpi_size
            write(*,*) 'Maybe it is better to run it with numMshRanks2Part=mpi_size... just suggesting'
         end if

         if(mpi_rank.lt.numMshRanks2Part) then
            numRanksMpiRank = 1
            iMshRankStart   = mpi_rank
            iMshRankEnd     = mpi_rank

            allocate(ranksMpiRank(numRanksMpiRank))
            ranksMpiRank(1) = mpi_rank
         else
            numRanksMpiRank = 0
            iMshRankStart   = 0
            iMshRankEnd     = 0

            allocate(ranksMpiRank(numRanksMpiRank))
         endif

         allocate(mapRanksToMpiRank(numMshRanks2Part))
         do iMpiRank=0,(numMshRanks2Part-1)
            mapRanksToMpiRank(iMpiRank+1)=iMpiRank
         end do

         !write(*,*) 'numMshRanks2Part<mpi_size ... haig de pensar en aquest cas...'
         !call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      !write(*,*) '#rank',mpi_rank,'numRanksMpiRank',numRanksMpiRank,'ranksInP',ranksMpiRank(:),'mapRanksToMpiRank',mapRanksToMpiRank(:)

   end subroutine distribute_ranks2Part_in_mpiRank

   subroutine get_rankPartitionBoundaryNodes_in_parallel(mshRank,numElemsInRank,numNodesInRank,numBoundsInRank,listNodesInRank,connecInRank,coordInRank,&
               listElemsBoundsInRank,listBoundFacesInRank,connecBoundFacesInRank,numMasSlaNodesInRank,masSlaNodesInRank,&
               numMshBoundNodesRankPar,numMpiBoundNodesRankPar,numInnerNodesRankPar,mpiBoundaryNodes,mshBoundaryNodes)
      implicit none
      integer, intent(in) :: mshRank,numElemsInRank,numNodesInRank,numBoundsInRank,listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode)
      integer, intent(in) :: numMasSlaNodesInRank,masSlaNodesInRank(numMasSlaNodesInRank,2)
      real(rp),intent(in) :: coordInRank(numNodesInRank,3)
      integer, intent(in) :: listElemsBoundsInRank(numElemsInRank,maxBoundsPerElem),listBoundFacesInRank(numBoundsInRank),connecBoundFacesInRank(numBoundsInRank,npbou)
      
      integer, intent(out) :: numMshBoundNodesRankPar,numMpiBoundNodesRankPar,numInnerNodesRankPar
      integer, allocatable, intent(out) :: mpiBoundaryNodes(:),mshBoundaryNodes(:)
      
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

      do ind=1,numMasSlaNodesInRank
         iNodeG = masSlaNodesInRank(ind,1)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         nodeOwnedCnt(indexNode) = 2

         iNodeG = masSlaNodesInRank(ind,2)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         nodeOwnedCnt(indexNode) = 0
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

         call find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceFront2ijk, listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceLeft2ijk,  listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceTop2ijk,   listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceBack2ijk,  listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceRight2ijk, listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
         call find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,faceBottom2ijk,listNodesInRank,connecInRank,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)

      end do

      !numBoundNodesRankPar=0
      numMshBoundNodesRankPar=0
      numMpiBoundNodesRankPar=0
      numInnerNodesRankPar=0

      do iNode=1,numNodesInRank
         if(nodeInBoundary(iNode)) then !node is boundary
            !numBoundNodesRankPar=numBoundNodesRankPar+1
            if(nodeInMshBoundary(iNode)) numMshBoundNodesRankPar=numMshBoundNodesRankPar+1
            if(nodeInMpiBoundary(iNode)) numMpiBoundNodesRankPar=numMpiBoundNodesRankPar+1
         else !node is inner
            numInnerNodesRankPar=numInnerNodesRankPar+1
         end if
      end do

      !write(*,*) '1.#rank[',mpi_rank,'] numNodesInRank',numNodesInRank,'bMshN',numMshBoundNodesRankPar,'bMpiN',numMpiBoundNodesRankPar,'iN',numInnerNodesRankPar
      !allocate(boundaryNodes(numBoundNodesRankPar))
      allocate(mpiBoundaryNodes(numMpiBoundNodesRankPar))
      allocate(mshBoundaryNodes(numMshBoundNodesRankPar))

      auxCnt=0
      mpiAuxCnt=0
      mshAuxCnt=0
      do iNode=1,numNodesInRank
         iNodeG=listNodesInRank(iNode)
         if(nodeInBoundary(iNode)) then !node is boundary
            auxCnt=auxCnt+1
            !boundaryNodes(auxCnt) = iNodeG
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

      !write(*,*) '#rank[',mpi_rank,']  nodesInRank ',numNodesMshRank, " bN ", numBNodesRankPar, " iN ",numINodesRankPar,'i',i

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      write(aux_string_mshrank_deb,'(I0)') mshRank
      !file_name_deb = 'boundaryNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      !open(1, file=file_name_deb)
      !write(1,*) 'X,Y,Z,iNodeG'
      !do iNode=1,numBoundNodesRankPar
      !   iNodeG=boundaryNodes(iNode)
      !   indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
      !   
      !   write(1,fmt_csv_deb) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      !end do
      !close(1)

      file_name_deb = 'boundaryMpiNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numMpiBoundNodesRankPar
         iNodeG=mpiBoundaryNodes(iNode)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         
         write(1,fmt_csv_deb) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

      file_name_deb = 'boundaryMshNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numMshBoundNodesRankPar
         iNodeG=mshBoundaryNodes(iNode)
         indexNode = binarySearch_int_i(listNodesInRank,iNodeG)
         
         write(1,fmt_csv_deb) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

#endif
   end subroutine get_rankPartitionBoundaryNodes_in_parallel

   subroutine find_mpi_and_msh_boundaries(numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,face2ijk,listNodesInRank,connecInRank,nodeOwnedCnt,listBoundsInRank,boundFacesInRank,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary)
      implicit none
      integer,intent(in) :: numElemsInRank,numNodesInRank,numBoundsInRank
      integer,intent(in) :: iElemL,face2ijk(npbou),listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode),nodeOwnedCnt(numNodesInRank),listBoundsInRank(numBoundsInRank),boundFacesInRank(numBoundsInRank,npbou)
      integer,intent(in) :: numBoundFacesToCheck,boundFacesToCheck(maxBoundsPerElem)
      logical,intent(out) :: nodeInBoundary(numNodesInRank),nodeInMshBoundary(numNodesInRank),nodeInMpiBoundary(numNodesInRank)
      integer :: ii,jj,ind,ind_f,iFace,iNodeG,iNodeG_inFace,iNodeG_ToCheck,indexNode,indexNode_inFace
      integer :: iBound,iBoundG,iElemG,indexBound
      logical :: isMshBound
      integer, parameter :: checkFacePos = 6
      
      !# 1.Whatever face ------------------------------------------------------
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

   end subroutine find_mpi_and_msh_boundaries

   subroutine define_parallelNodePartitioning(numMshRanks2Part,numMshRanksInMpiRank,numNodesMshRank,mshRanksInMpiRank,mapMshRankToMpiRank,mshRankNodeStart,mshRankNodeEnd,iNodeStartPar,numNodesParTotal)
      implicit none
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in) :: numMshRanksInMpiRank,numNodesMshRank(numMshRanksInMpiRank),mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanksInMpiRank)
      integer,intent(out) :: mshRankNodeStart(numMshRanksInMpiRank),mshRankNodeEnd(numMshRanksInMpiRank)
      integer,dimension(0:numMshRanks2Part-1),intent(out) :: iNodeStartPar
      integer,intent(out) :: numNodesParTotal
      integer,dimension(0:numMshRanks2Part-1) :: iNodeEndPar,vecNumNodesMshRank
      integer :: window_id
      integer :: iMshRank,mshRank,mpiRank,auxNodeCnt

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numNodesMshRank,vecNumNodesMshRank)
      !write(*,*) 'rank[',mpi_rank,']vNNR2s', vectorNumNodesRankPar2send(:),'vNNR2r',vectorNumNodesRankPar2rcv(:),'nNrp',numNodesMshRank(:)

      auxNodeCnt = 1
      do iMshRank=0,numMshRanks2Part-1
      !do iRank=0,mpi_size-1
         iNodeStartPar(iMshRank) = auxNodeCnt
         iNodeEndPar(iMshRank)   = iNodeStartPar(iMshRank) + vecNumNodesMshRank(iMshRank) - 1
         auxNodeCnt = iNodeEndPar(iMshRank) + 1
      end do

      numNodesParTotal = iNodeEndPar(numMshRanks2Part-1)

      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         !vectorNumNodesRankPar(mshRank) = numNodesMshRank(iMshRank)
         mshRankNodeStart(iMshRank) = iNodeStartPar(mshRank)
         mshRankNodeEnd(iMshRank)   = iNodeEndPar(mshRank)
         !write(*,*) 'mshrank',mshrank,'[',mpi_rank,']nodesInRank',numNodesMshRank(iMshRank),'iNodeS',mshRankNodeStart(iMshRank),'iNodeE',mshRankNodeEnd(iMshRank)
         !write(*,*) 'mshrank',mshrank,'[',mpi_rank,']start',iNodeStartPar(:),'end',iNodeEndPar(:),'tNNP',numNodesParTotal
      end do
   end subroutine define_parallelNodePartitioning

   subroutine reorder_nodes_in_mshRank(mshRank,numMshRanks2Part,numElemsInRank,numNodesInRank,numBoundsInRank,numMasSlaNodesInRank,elemGidInRank,listNodesInRank,connecInRank,iNodeStartPar,&
               globalIdSrlInRank,globalIdSrlOrderedInRank,globalIdParInRank,connecVTKInRank,connecParOrigInRank,boundFacesInRank,masSlaNodesInRank)
      implicit none
      integer,intent(in) :: mshRank,numMshRanks2Part,numElemsInRank,numNodesInRank,numBoundsInRank,numMasSlaNodesInRank
      integer,intent(in) :: elemGidInRank(numElemsInRank),listNodesInRank(numNodesInRank),connecInRank(numElemsInRank,nnode)
      integer,dimension(0:numMshRanks2Part-1),intent(in) :: iNodeStartPar

      integer,intent(out) :: globalIdSrlInRank(numNodesInRank),globalIdParInRank(numNodesInRank),connecVTKInRank(numElemsInRank*nnode),connecParOrigInRank(numElemsInRank,nnode)
      integer,intent(inout) :: boundFacesInRank(numBoundsInRank,npbou),masSlaNodesInRank(numMasSlaNodesInRank,2)
      integer,intent(out),dimension(numNodesInRank,2) :: globalIdSrlOrderedInRank

      integer :: m,indConn,indexIJK,indexVTK,indexGMSH,indexNew,nodeIndexCnt,indPosListNodes
      integer :: iNodeL,iNodeGpar,iNodeGsrl,iElemL,iElemG,iNode,iBound,iLink

      integer,dimension(nnode) :: auxNodeNewOrderInElem,auxNewOrderIndex,auxVTKorder,auxCGNSorder
      integer,dimension(numNodesInRank) :: isNodeAdded
      integer :: aux_ibound_GtoL(npbou), aux_iMasSla_GtoL(2)

      !$acc kernels
      isNodeAdded(:)=-1
      globalIdSrlOrderedInRank(:,:)=-1
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
         !iElemG = (iElemL-1) + mshRankElemStart
         iElemG = elemGidInRank(iElemL)
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

               iNodeGPar = iNodeL + iNodeStartPar(mshRank) - 1

               globalIdSrlInRank(iNodeL) = iNodeGsrl
               globalIdParInRank(iNodeL) = iNodeGPar

               globalIdSrlOrderedInRank(iNodeL,1) = iNodeGsrl
               globalIdSrlOrderedInRank(iNodeL,2) = iNodeL
            else 
               iNodeL = isNodeAdded(indPosListNodes)
               !iNodeGPar = globalIdPar(iNodeL)
            endif

            connecParOrigInRank(iElemL,m) = iNodeL

            indexVTK = auxVTKorder(m)
            indConn = (iElemL-1)*nnode + indexVTK
            connecVTKInRank(indConn) = iNodeL
         end do
         !write(*,*) '[',mpi_rank,']iElemG ',iElemG,' connecParOrig ',connecParOrig(iElemL,:)
      end do
      !if(mpi_rank.eq.2)write(*,*) 'check globalIdSrl[',mpi_rank,']',globalIdSrl(:)

      call quicksort_matrix_int(globalIdSrlOrderedInRank,1)

      do iBound=1,numBoundsInRank
         aux_ibound_GtoL(:) = 0
         !write(*,*) '1.reorder ibound[',mpi_rank,']npbouG',boundFacesInRank(iBound,:)
         do iNode=1,npbou
            iNodeGsrl = boundFacesInRank(iBound,iNode)
            indexNew = binarySearch_int_i(globalIdSrlOrderedInRank(:,1),iNodeGsrl)
            iNodeL = globalIdSrlOrderedInRank(indexNew,2)
            !if(mpi_rank.eq.2)write(*,*) 'ibound[',mpi_rank,']',iBound,'iNodeGsrl',iNodeGsrl,'iNodeL',iNodeL
            aux_ibound_GtoL(iNode) = iNodeL
         end do

         boundFacesInRank(iBound,:) = aux_ibound_GtoL(:)
         !write(*,*) '2.reorder ibound[',mpi_rank,']npbouG',boundFacesInRank(iBound,:)
      end do

      !to be done here: modificar el masSlaNodesInRank
      !----------------------------------------------------------------
      do iLink=1,numMasSlaNodesInRank
         !master
         iNodeGsrl = masSlaNodesInRank(iLink,1)
         indexNew = binarySearch_int_i(globalIdSrlOrderedInRank(:,1),iNodeGsrl)
         iNodeL = globalIdSrlOrderedInRank(indexNew,2)
         aux_iMasSla_GtoL(1) = iNodeL
         !slave
         iNodeGsrl = masSlaNodesInRank(iLink,2)
         indexNew = binarySearch_int_i(globalIdSrlOrderedInRank(:,1),iNodeGsrl)
         iNodeL = globalIdSrlOrderedInRank(indexNew,2)
         aux_iMasSla_GtoL(2) = iNodeL

         masSlaNodesInRank(iLink,:) = aux_iMasSla_GtoL(:)
      end do

      !reordeno el vector per la columna slaves, ja que son unics
      call quicksort_matrix_int(masSlaNodesInRank,2)
      !----------------------------------------------------------------

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

#define _NEWMETHOD_ 1

   subroutine generate_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numMpiBoundaryNodes,mpiBoundaryNodes_jv,globalIdSrl_jv,globalIdSrlOrdered_jm,&
               matrixCommScheme_jm,commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv,numNodesToCommMshRank,numMshRanksWithComms)
      !generate a matrix with the comm schemes for shared nodes between procs
      integer,intent(in) :: numMshRanks2Part,numMpiBoundaryNodes(numMshRanksInMpiRank)
      integer,intent(in) :: numMshRanksInMpiRank,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part)
      type(jagged_vector_int),intent(in) :: mpiBoundaryNodes_jv,globalIdSrl_jv
      type(jagged_matrix_int),intent(in) :: globalIdSrlOrdered_jm
      type(jagged_matrix_int),intent(inout) :: matrixCommScheme_jm
      type(jagged_vector_int),intent(inout) :: commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv
      integer(int_size),allocatable,intent(inout) :: numNodesToCommMshRank(:),numMshRanksWithComms(:)

      integer :: iMshRank,mshRank,iMshRankTrgt,mpiRank,iAux,memPos,memSize,memDisp
      integer :: mshRankOrig,mshRankTrgt,numBNOrig,numBNTrgt,mpiRankTrgt
      integer :: iNode,iPos,iNodeL,iNodeGSrl,iNodeLPos

      integer, dimension(0:numMshRanks2Part-1) :: auxCommSchemeNumNodes
      integer,allocatable :: auxCommSchemeMemPos(:,:),auxCommSchemeMemPosAll(:)

      integer, dimension(0:numMshRanks2Part-1) :: vecNumMpiBN,vecMemDispMpiBN
      integer,allocatable :: arrayMpiBNInMpiRank(:),auxMpiBNMshRank(:)
      integer :: numMpiBNInMpiRank,auxMemDispMpiBN
      
      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      !------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) "--| Generating MPI Comm scheme"

      !NEW METHOD
      !---------------------------------------------------------------------

      vecNumMpiBN(:) = 0
      vecMemDispMpiBN(:) = 0
      numMpiBNInMpiRank = 0

      auxMemDispMpiBN = 0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vecNumMpiBN(mshRank) = numMpiBoundaryNodes(iMshRank)
         numMpiBNInMpiRank = numMpiBNInMpiRank + numMpiBoundaryNodes(iMshRank)
         if(iMshRank.ne.1) auxMemDispMpiBN = auxMemDispMpiBN + numMpiBoundaryNodes(iMshRank-1)
         vecMemDispMpiBN(mshRank) = auxMemDispMpiBN
      end do
      
      !if(mpi_rank.eq.0) then
      !   write(*,*) 'numMpiBNInMpiRank',numMpiBNInMpiRank
      !   write(*,*) 'vecNumMpiBN(',mpi_rank,')',vecNumMpiBN(:)
      !   write(*,*) 'vecMemDispMpiBN(',mpi_rank,')',vecMemDispMpiBN(:)
      !endif

      allocate(arrayMpiBNInMpiRank(numMpiBNInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         do iAux=1,numMpiBoundaryNodes(iMshRank)
            memPos = vecMemDispMpiBN(mshRank) + iAux
            arrayMpiBNInMpiRank(memPos) = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      !if(mpi_rank.eq.0) then
      !   write(*,*) 'arrayMpiBNInMpiRank(',mpi_rank,')',arrayMpiBNInMpiRank(:)
      !endif

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vecNumMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecNumMpiBN(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      call MPI_Win_create(vecMemDispMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecMemDispMpiBN(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'vecNumMpiBN',vecNumMpiBN(:)
      !if(mpi_rank.eq.0) write(*,*) 'vecMemDispMpiBN',vecMemDispMpiBN(:)

      allocate(numNodesToCommMshRank(numMshRanksInMpiRank))
      allocate(numMshRanksWithComms(numMshRanksInMpiRank))
      allocate(matrixCommScheme_jm%matrix(numMshRanksInMpiRank))
      allocate(ranksToComm_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemSize_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemPosInLoc_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemPosInNgb_jv%vector(numMshRanksInMpiRank))
      allocate(auxCommSchemeMemPos(numMshRanksInMpiRank,numMshRanks2Part))

      !--------------------------------------------------------------------------------------------
      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

      win_buffer_size = mpi_integer_size*numMpiBNInMpiRank
      call MPI_Win_create(arrayMpiBNInMpiRank,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)

      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numBNOrig = numMpiBoundaryNodes(iMshRank)

         numNodesToCommMshRank(iMshRank)=0
         auxCommSchemeNumNodes(:)=0
         auxCommSchemeMemPos(iMshRank,:)=0

         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then

               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispMpiBN(mshRankTrgt)
               memSize = vecNumMpiBN(mshRankTrgt)
               allocate(auxMpiBNMshRank(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxMpiBNMshRank,memSize,MPI_INTEGER,mpiRankTrgt,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               !if(mpi_rank.eq.0 .and. mshRankTrgt.eq.1) write(*,*) 'memSize',memSize,'td',vecMemDispMpiBN(mshRankTrgt),'auxMpiBN(',mpi_rank,')',auxMpiBNMshRank

               do iNode=1,numBNOrig
                  iNodeGSrl = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)

                  iPos = binarySearch_int_i(auxMpiBNMshRank,iNodeGSrl)
                  if(iPos.ne.0) then
                     numNodesToCommMshRank(iMshRank) = numNodesToCommMshRank(iMshRank) + 1
                     auxCommSchemeNumNodes(mshRankTrgt) = auxCommSchemeNumNodes(mshRankTrgt)+1

                     !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBNO',numBNOrig,'mshRankT',mshRankTrgt,'memPos',memPos,'numBNT',numBNTrgt
                     !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'mshRankT',mshRankTrgt,'iNGO',iNodeGOrig,'iNGT',iNodeGTrgt
                  end if
                  !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBNO',numBNOrig,'mshRankT',mshRankTrgt,'memPos',memPos,'numBNT',numBNTrgt
               end do

               deallocate(auxMpiBNMshRank)
            end if

         end do

         ! setting numRanksWithsComms
         numMshRanksWithComms(iMshRank) = 0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               numMshRanksWithComms(iMshRank) = numMshRanksWithComms(iMshRank) + 1
            end if
         end do

         allocate(ranksToComm_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))
         allocate(commsMemSize_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))
         allocate(commsMemPosInLoc_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))
         allocate(commsMemPosInNgb_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))

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

         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'csNumNodes->',auxCommSchemeNumNodes(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'numRanksWC->',numMshRanksWithComms(iMshRank)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'ranksToComm->',ranksToComm_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInLoc->',commsMemPosInLoc_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemSize->',commsMemSize_jv%vector(iMshRank)%elems(:)

         !write(*,*) 'rank[',mpi_rank,']mshRankOrig',mshRankOrig,'numNodesToCommMshRank',numNodesToCommMshRank(iMshRank),'numBNOrig',numBNOrig,'cSNN',auxCommSchemeNumNodes(:),&
         !          'numMshRanksWithComms',numMshRanksWithComms(iMshRank),'r2C',ranksToComm_jv%vector(iMshRank)%elems(:)
         allocate(matrixCommScheme_jm%matrix(iMshRank)%elems(numNodesToCommMshRank(iMshRank),3))
         iAux=0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispMpiBN(mshRankTrgt)
               memSize = vecNumMpiBN(mshRankTrgt)
               allocate(auxMpiBNMshRank(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxMpiBNMshRank,memSize,MPI_INTEGER,mpiRankTrgt,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numBNOrig
                  iNodeGSrl = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)

                  iPos = binarySearch_int_i(auxMpiBNMshRank,iNodeGSrl)
                  if(iPos.ne.0) then
                     iAux=iAux+1
                     
                     iNodeLPos = binarySearch_int_i(globalIdSrlOrdered_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
                     iNodeL = globalIdSrlOrdered_jm%matrix(iMshRank)%elems(iNodeLPos,2)
                     !write(*,*) 'iNodeL',iNodeL,'iNodeG',iNodeGSrl,'mshRank',mshRankTrgt
                     matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,1) = iNodeL
                     matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,2) = iNodeGSrl
                     matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,3) = mshRankTrgt
                  end if

               end do

               deallocate(auxMpiBNMshRank)
            end if
         end do

      end do

      call MPI_Win_free(window_id,mpi_err)

      !---------------------------------------------------------------------

      allocate(auxCommSchemeMemPosAll(numMshRanks2Part**2))

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
      win_buffer_size = mpi_integer_size*(numMshRanks2Part**2)

      call MPI_Win_create(auxCommSchemeMemPosAll,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mpiRank = mapMshRankToMpiRank(iMshRank)
         memSize = numMshRanks2Part
         memPos = (iMshRank-1)*numMshRanks2Part+1
         target_displacement = (iMshRank-1)*numMshRanks2Part

         call MPI_Get(auxCommSchemeMemPosAll(memPos),memSize,MPI_INTEGER,mpiRank,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------
      !write(*,*) '2.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numMshRanksWithComms(iMshRank)
            mshRankTrgt=ranksToComm_jv%vector(iMshRank)%elems(iMshRankTrgt)
            iPos = mshRankTrgt*numMshRanks2Part+mshRankOrig+1
            commsMemPosInNgb_jv%vector(iMshRank)%elems(iMshRankTrgt) = auxCommSchemeMemPosAll(iPos)
            !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'iMRT',iMshRankTrgt,'mRT',mshRankTrgt,'iPos',iPos
         end do
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInNgb->',commsMemPosInNgb_jv%vector(iMshRank)%elems(:)
      end do

      deallocate(auxCommSchemeMemPos)
      deallocate(auxCommSchemeMemPosAll)

   end subroutine generate_mpi_comm_scheme_parallel

   subroutine generate_dof_and_boundary_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,&
               numNodesMshRank,numMshBoundaryNodes,mshBoundaryNodes_jv,numMpiBoundaryNodes,mpiBoundaryNodes_jv,globalIdSrl_jv,globalIdSrlOrdered_jm,&
               numBoundaryNodesMshRank,boundaryNodes_jv,numDoFMshRank,dofNodes_jv,&
               bnd_matrixCommScheme_jm,bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv,bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms)
      !generate a matrix with the comm schemes for shared nodes between procs
      integer,intent(in) :: numMshRanks2Part,numNodesMshRank(numMshRanksInMpiRank),numMpiBoundaryNodes(numMshRanksInMpiRank),numMshBoundaryNodes(numMshRanksInMpiRank)
      integer,intent(in) :: numMshRanksInMpiRank,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part)
      type(jagged_vector_int),intent(in) :: mshBoundaryNodes_jv,mpiBoundaryNodes_jv,globalIdSrl_jv
      type(jagged_matrix_int),intent(in) :: globalIdSrlOrdered_jm
      integer,dimension(numMshRanksInMpiRank),intent(out) :: numBoundaryNodesMshRank,numDoFMshRank
      type(jagged_vector_int),intent(inout) :: boundaryNodes_jv,dofNodes_jv
      type(jagged_matrix_int),intent(inout) :: bnd_matrixCommScheme_jm
      type(jagged_vector_int),intent(inout) :: bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv
      integer(int_size),allocatable,intent(inout) :: bnd_numNodesToCommMshRank(:),bnd_numMshRanksWithComms(:)

      integer :: iMshRank,mshRank,iMshRankTrgt,mpiRank,memPos,memSize,memDisp
      integer :: iAux,auxCntBoundaryNode,auxCnt1,auxCnt2
      integer :: mshRankOrig,mshRankTrgt,mpiRankTrgt,numMpiBNOrig,numMpiBNTrgt,numMshBNOrig,numMshBNTrgt,numBndMpiNodeOrig,numBndMpiNodeTrgt
      integer :: iNode,iNodeBnd,iNodeMpi,iPos,iNodeL,iNodeGSrl,iNodeLPos,iPosMshBN

      integer :: numBndMpiBoundaryNodes(numMshRanksInMpiRank),numBndMpiBoundaryNodesAll(numMshRanks2Part)
      integer, dimension(0:numMshRanks2Part-1) :: bnd_auxCommSchemeNumNodes
      integer,allocatable :: bnd_auxCommSchemeMemPos(:,:),bnd_auxCommSchemeMemPosAll(:)
      logical,allocatable :: auxListBN(:)

      type(jagged_vector_int) :: bndMpiBoundaryNodes_jv
      integer, dimension(0:numMshRanks2Part-1) :: vecNumMshBN,vecMemDispMshBN,vecNumBndMpiBN,vecMemDispBndMpiBN
      integer,allocatable :: arrayMshBNInMpiRank(:),auxMshBNMshRank(:),arrayBndMpiBNInMpiRank(:),auxBndMpiBNMshRank(:)
      integer :: numMshBNInMpiRank,numBndMpiBNInMpiRank

      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      !------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) "--| Generating DoF & Boundary MPI Comm scheme"

      !NEW APPROACH
      !---------------------------------------------------------------------
      vecNumMshBN(:) = 0
      vecMemDispMshBN(:) = 0
      numMshBNInMpiRank = 0

      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vecNumMshBN(mshRank) = numMshBoundaryNodes(iMshRank)
         numMshBNInMpiRank = numMshBNInMpiRank + numMshBoundaryNodes(iMshRank)
         if(iMshRank.ne.1) vecMemDispMshBN(mshRank) = vecMemDispMshBN(mshRank) + numMshBoundaryNodes(iMshRank-1)
      end do
      !write(*,*) 'numMshBNInMpiRank',numMshBNInMpiRank
      !write(*,*) 'vecNumMshBN(',mpi_rank,')',vecNumMshBN(:)
      !write(*,*) 'vecMemDispMshBN(',mpi_rank,')',vecMemDispMshBN(:)
      allocate(arrayMshBNInMpiRank(numMshBNInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         do iAux=1,numMshBoundaryNodes(iMshRank)
            memPos = vecMemDispMshBN(mshRank) + iAux
            arrayMshBNInMpiRank(memPos) = mshBoundaryNodes_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vecNumMshBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecNumMshBN(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      call MPI_Win_create(vecMemDispMshBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecMemDispMshBN(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------
      !lets determine ALL the boundary nodes in each rank
      !the boundary nodes in a rank will be:
      ! 1. the boundary nodes in the container mshBoundaryNodes_jv of each rank
      ! 2. plus the nodes that are in mpiBoundaryNode_jv in the rank and in the boundary nodes container mshBoundaryNodes_jv in the other ranks 

      !--------------------------------------------------------------------------------------------
      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
      win_buffer_size = mpi_integer_size*numMshBNInMpiRank
      call MPI_Win_create(arrayMshBNInMpiRank,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numMshBNOrig = numMshBoundaryNodes(iMshRank)
         numMpiBNOrig = numMpiBoundaryNodes(iMshRank)

         auxCntBoundaryNode=0

         allocate(auxListBN(numNodesMshRank(iMshRank)))
         auxListBN(:)=.false.
         !1.
         do iNode=1,numMshBNOrig
            iNodeGSrl = mshBoundaryNodes_jv%vector(iMshRank)%elems(iNode)
            iNodeLPos = binarySearch_int_i(globalIdSrlOrdered_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
            iNodeL    = globalIdSrlOrdered_jm%matrix(iMshRank)%elems(iNodeLPos,2)
            
            auxListBN(iNodeL) = .true.
            auxCntBoundaryNode = auxCntBoundaryNode + 1
         end do

         !2.
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispMshBN(mshRankTrgt)
               memSize = vecNumMshBN(mshRankTrgt)
               allocate(auxMshBNMshRank(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxMshBNMshRank,memSize,MPI_INTEGER,mpiRankTrgt,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numMpiBNOrig
                  iNodeGSrl = mpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)
                  iPosMshBN = binarySearch_int_i(mshBoundaryNodes_jv%vector(iMshRank)%elems,iNodeGSrl)
                  if(iPosMshBN.eq.0) then !if is not already in a boundary face
                     iPos = binarySearch_int_i(auxMshBNMshRank,iNodeGSrl)
                     if(iPos.ne.0) then
                        iNodeLPos = binarySearch_int_i(globalIdSrlOrdered_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
                        iNodeL    = globalIdSrlOrdered_jm%matrix(iMshRank)%elems(iNodeLPos,2)
                        if(not(auxListBN(iNodeL))) then
                           auxListBN(iNodeL) = .true.
                           auxCntBoundaryNode = auxCntBoundaryNode + 1
                        end if
                     end if
                  end if
               end do

               deallocate(auxMshBNMshRank)

            end if
         end do

         !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numMshBN',numMshBNOrig,'numMpiBN',numMpiBNOrig,'addCntBN',auxCntBoundaryNode
         numBoundaryNodesMshRank(iMshRank) = auxCntBoundaryNode
         numDoFMshRank(iMshRank) = numNodesMshRank(iMshRank) - numBoundaryNodesMshRank(iMshRank)
         !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBN',numBoundaryNodesMshRank(iMshRank),'numDoF',numDoFMshRank(iMshRank),'numN',numNodesMshRank(iMshRank)

         allocate(boundaryNodes_jv%vector(iMshRank)%elems(numBoundaryNodesMshRank(iMshRank)))
         allocate(dofNodes_jv%vector(iMshRank)%elems(numDoFMshRank(iMshRank)))

         auxCnt1=0
         auxCnt2=0
         do iNodeL=1,numNodesMshRank(iMshRank)
            if(auxListBN(iNodeL)) then
               auxCnt1=auxCnt1+1
               boundaryNodes_jv%vector(iMshRank)%elems(auxCnt1) = iNodeL
            else
               auxCnt2=auxCnt2+1
               dofNodes_jv%vector(iMshRank)%elems(auxCnt2) = iNodeL
            end if
         end do

         deallocate(auxListBN)
      end do

      call MPI_Win_free(window_id,mpi_err)
      !---------------------------------------------------------------------------------------------------------
      ! END DOF SECTION!!!
      !---------------------------------------------------------------------------------------------------------

      allocate(bndMpiBoundaryNodes_jv%vector(numMshRanksInMpiRank))
      !generar el la jagged matrix amb els nou mpi boundary nodes, a partir d'aqui el proces es identic que abans
      !vaja, straightforward! :)

      do iMshRank=1,numMshRanksInMpiRank
         
         allocate(auxListBN(numBoundaryNodesMshRank(iMshRank)))
         auxListBN(:) = .false.
         auxCnt1=0
         do iNodeBnd=1,numBoundaryNodesMshRank(iMshRank)
            iNodeL = boundaryNodes_jv%vector(iMshRank)%elems(iNodeBnd)
            iNodeGSrl = globalIdSrl_jv%vector(iMshRank)%elems(iNodeL)
            iPos = binarySearch_int_i(mpiBoundaryNodes_jv%vector(iMshRank)%elems,iNodeGSrl)
            if(iPos.ne.0) then   
               auxCnt1=auxCnt1+1
               auxListBN(iNodeBnd) = .true.
            end if
         end do
         numBndMpiBoundaryNodes(iMshRank) = auxCnt1
         !write(*,*) 'mshRankO',mshRankOrig,'numBndMpi',numBndMpiBoundaryNodes(iMshRank)
      
         allocate(bndMpiBoundaryNodes_jv%vector(iMshRank)%elems(numBndMpiBoundaryNodes(iMshRank)))

         auxCnt1=0
         do iNodeBnd=1,numBoundaryNodesMshRank(iMshRank)
            if(auxListBN(iNodeBnd)) then
               iNodeL = boundaryNodes_jv%vector(iMshRank)%elems(iNodeBnd)
               iNodeGSrl = globalIdSrl_jv%vector(iMshRank)%elems(iNodeL)
               auxCnt1=auxCnt1+1
               !if(mpi_rank.eq.3)write(*,*) 'rank[',mpi_rank,']aux',auxCnt1,'iNodeL',iNodeL,'iNodeGSrl',iNodeGSrl
               bndMpiBoundaryNodes_jv%vector(iMshRank)%elems(auxCnt1) = iNodeGSrl
            end if
         end do
         deallocate(auxListBN)

         call quicksort_array_int(bndMpiBoundaryNodes_jv%vector(iMshRank)%elems)
      end do

!--------------------------------------------------------------------------------------------

      vecNumBndMpiBN(:) = 0
      vecMemDispBndMpiBN(:) = 0
      numBndMpiBNInMpiRank = 0

      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vecNumBndMpiBN(mshRank) = numBndMpiBoundaryNodes(iMshRank)
         numBndMpiBNInMpiRank = numBndMpiBNInMpiRank + numBndMpiBoundaryNodes(iMshRank)
         if(iMshRank.ne.1) vecMemDispBndMpiBN(mshRank) = vecMemDispBndMpiBN(mshRank) + numBndMpiBoundaryNodes(iMshRank-1)
      end do
      allocate(arrayBndMpiBNInMpiRank(numBndMpiBNInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         do iAux=1,numBndMpiBoundaryNodes(iMshRank)
            memPos = vecMemDispBndMpiBN(mshRank) + iAux
            arrayBndMpiBNInMpiRank(memPos) = bndMpiBoundaryNodes_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vecNumBndMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecNumBndMpiBN(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      call MPI_Win_create(vecMemDispBndMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecMemDispBndMpiBN(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'vecNumBndMpiBN',vecNumBndMpiBN(:)
      !if(mpi_rank.eq.0) write(*,*) 'vecMemDispBndMpiBN',vecMemDispBndMpiBN(:)
      !---------------------------------------------------------------------------------------------------------
      allocate(bnd_numNodesToCommMshRank(numMshRanksInMpiRank))
      allocate(bnd_numMshRanksWithComms(numMshRanksInMpiRank))
      allocate(bnd_matrixCommScheme_jm%matrix(numMshRanksInMpiRank))
      allocate(bnd_ranksToComm_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_commsMemSize_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_commsMemPosInLoc_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_commsMemPosInNgb_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_auxCommSchemeMemPos(numMshRanksInMpiRank,numMshRanks2Part))

      !---------------------------------------------------------------------------------------------------------------
      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
      win_buffer_size = mpi_integer_size*numBndMpiBNInMpiRank
      call MPI_Win_create(arrayBndMpiBNInMpiRank,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
      !---------------------------------------------------------------------------------------------------------------

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numBndMpiNodeOrig = numBndMpiBoundaryNodes(iMshRank)

         bnd_numNodesToCommMshRank(iMshRank)=0
         bnd_auxCommSchemeNumNodes(:)=0
         bnd_auxCommSchemeMemPos(iMshRank,:)=0

         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
               
            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispBndMpiBN(mshRankTrgt)
               memSize = vecNumBndMpiBN(mshRankTrgt)
               allocate(auxBndMpiBNMshRank(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxBndMpiBNMshRank,memSize,MPI_INTEGER,mpiRankTrgt,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numBndMpiNodeOrig
                  iNodeGSrl = bndMpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)
                  iPos = binarySearch_int_i(auxBndMpiBNMshRank,iNodeGSrl)
                  if(iPos.ne.0) then
                     bnd_numNodesToCommMshRank(iMshRank) = bnd_numNodesToCommMshRank(iMshRank) + 1
                     bnd_auxCommSchemeNumNodes(mshRankTrgt) = bnd_auxCommSchemeNumNodes(mshRankTrgt)+1
                  end if
               end do

               deallocate(auxBndMpiBNMshRank)

            end if
         end do

         ! setting numRanksWithsComms
         bnd_numMshRanksWithComms(iMshRank) = 0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(bnd_auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               bnd_numMshRanksWithComms(iMshRank) = bnd_numMshRanksWithComms(iMshRank) + 1
            end if
         end do

         allocate(bnd_ranksToComm_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))
         allocate(bnd_commsMemSize_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))
         allocate(bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))
         allocate(bnd_commsMemPosInNgb_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))

         iAux=0
         memPos=1
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            bnd_auxCommSchemeMemPos(iMshRank,iMshRankTrgt)=memPos
            if(bnd_auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               iAux=iAux+1
               mshRankTrgt = iMshRankTrgt-1
               bnd_ranksToComm_jv%vector(iMshRank)%elems(iAux) = mshRankTrgt

               bnd_commsMemSize_jv%vector(iMshRank)%elems(iAux) = bnd_auxCommSchemeNumNodes(mshRankTrgt)
               bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems(iAux) = memPos
               memPos=memPos+bnd_auxCommSchemeNumNodes(mshRankTrgt)
            end if
         end do


         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_csNumNodes->',bnd_auxCommSchemeNumNodes(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_numNodesToComm',bnd_numNodesToCommMshRank(iMshRank)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_numRanksWC->',bnd_numMshRanksWithComms(iMshRank)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_ranksToComm->',bnd_ranksToComm_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_commsMemPosInLoc->',bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_commsMemSize->',bnd_commsMemSize_jv%vector(iMshRank)%elems(:)

         allocate(bnd_matrixCommScheme_jm%matrix(iMshRank)%elems(bnd_numNodesToCommMshRank(iMshRank),3))
         iAux=0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispBndMpiBN(mshRankTrgt)
               memSize = vecNumBndMpiBN(mshRankTrgt)
               allocate(auxBndMpiBNMshRank(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxBndMpiBNMshRank,memSize,MPI_INTEGER,mpiRankTrgt,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numBndMpiNodeOrig
                  iNodeGSrl = bndMpiBoundaryNodes_jv%vector(iMshRank)%elems(iNode)
                  iPos = binarySearch_int_i(auxBndMpiBNMshRank,iNodeGSrl)
                  if(iPos.ne.0) then
                     iAux=iAux+1

                     iNodeLPos = binarySearch_int_i(globalIdSrlOrdered_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
                     iNodeL = globalIdSrlOrdered_jm%matrix(iMshRank)%elems(iNodeLPos,2)

                     bnd_matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,1) = iNodeL
                     bnd_matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,2) = iNodeGSrl
                     bnd_matrixCommScheme_jm%matrix(iMshRank)%elems(iAux,3) = mshRankTrgt
                  end if
               end do

               deallocate(auxBndMpiBNMshRank)
            end if
         end do
      end do

      call MPI_Win_free(window_id,mpi_err)

      !---------------------------------------------------------------------------------------------------------

      allocate(bnd_auxCommSchemeMemPosAll(numMshRanks2Part**2))

      bnd_auxCommSchemeMemPosAll(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numMshRanks2Part
            iPos = mshRank*numMshRanks2Part+iMshRankTrgt
            bnd_auxCommSchemeMemPosAll(iPos) = bnd_auxCommSchemeMemPos(iMshRank,iMshRankTrgt)
         end do
      end do
      !write(*,*) '1.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      ! Create the window
      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*(numMshRanks2Part**2)

      call MPI_Win_create(bnd_auxCommSchemeMemPosAll,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         !mshRank = iMshRank
         mpiRank = mapMshRankToMpiRank(iMshRank)
         memSize = numMshRanks2Part
         memPos = (iMshRank-1)*numMshRanks2Part+1
         target_displacement = (iMshRank-1)*numMshRanks2Part

         call MPI_Get(bnd_auxCommSchemeMemPosAll(memPos),memSize,MPI_INTEGER,mpiRank,target_displacement,memSize,MPI_INTEGER,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      !write(*,*) '2.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,bnd_numMshRanksWithComms(iMshRank)
            mshRankTrgt=bnd_ranksToComm_jv%vector(iMshRank)%elems(iMshRankTrgt)
            iPos = mshRankTrgt*numMshRanks2Part+mshRankOrig+1
            bnd_commsMemPosInNgb_jv%vector(iMshRank)%elems(iMshRankTrgt) = bnd_auxCommSchemeMemPosAll(iPos)
            !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'iMRT',iMshRankTrgt,'mRT',mshRankTrgt,'iPos',iPos
         end do
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInNgb->',commsMemPosInNgb_jv%vector(iMshRank)%elems(:)
      end do

      deallocate(bnd_auxCommSchemeMemPos)
      deallocate(bnd_auxCommSchemeMemPosAll)

   end subroutine generate_dof_and_boundary_mpi_comm_scheme_parallel

!------------------------------------------------------------------------------------------------------------------------------------

   subroutine set_nodesCoordinates(mshRank,numElemsInRank,numNodesInRank,globalIdSrlInRank,listNodesInRank,coordInRank,Ngp_l,connecParOrigMshRank,coordParMshRank)
      implicit none
      integer, intent(in) :: mshRank,numElemsInRank,numNodesInRank,globalIdSrlInRank(numNodesInRank),listNodesInRank(numNodesInRank),connecParOrigMshRank(numElemsInRank,nnode)
      real(rp),intent(in) :: coordInRank(numNodesInRank,3),Ngp_l(ngaus,nnode)
      real(rp),intent(out) :: coordParMshRank(numNodesInRank,3)
      integer :: iPos,iNodeL,iNodeGsrl

      !if(mpi_rank.eq.0) write(*,*) ' # Creating nodes coordinates...'
      !------------------------------------------------------

      !- create the coordinate data for this process
      !!!$acc parallel loop
      do iNodeL=1,numNodesInRank
         iNodeGSrl=globalIdSrlInRank(iNodeL)
         iPos = binarySearch_int_i(listNodesInRank,iNodeGSrl)
         !write (*,*) 'mshRank',mshRank,'iL', iNodeL,'iG', iNodeGSrl,'iPos',iPos,'[',coordInRank(iPos,1),']','[',coordInRank(iPos,2),']'
         coordParMshRank(iNodeL,1) = coordInRank(iPos,1)
         coordParMshRank(iNodeL,2) = coordInRank(iPos,2)
         coordParMshRank(iNodeL,3) = coordInRank(iPos,3)
      end do
      !!!$acc end parallel loop

      call interpolateOriginalCoordinates(numElemsInRank,numNodesInRank,Ngp_l,connecParOrigMshRank,coordParMshRank)

   end subroutine set_nodesCoordinates

   subroutine interpolateOriginalCoordinates(numElemsInRank,numNodesInRank,Ngp_l,connecParOrigMshRank,coordParMshRank)
      integer,intent(in) :: numElemsInRank,numNodesInRank,connecParOrigMshRank(numElemsInRank,nnode)
      real(rp),intent(in) :: Ngp_l(ngaus,nnode)
      real(rp),intent(inout) :: coordParMshRank(numNodesInRank,3)
      real(rp), allocatable :: aux_1(:,:)
      integer(4) :: iElem,iNode,idime
      
      !if(mpi_rank.eq.0) write(*,*) "--| Interpolating nodes coordinates..."
      allocate(aux_1(numNodesInRank,ndime))
      aux_1(:,:) = coordParMshRank(:,:)
      do ielem = 1,numElemsInRank
         do inode = (2**ndime)+1,nnode
            do idime = 1,ndime
               call var_interpolate(aux_1(connecParOrigMshRank(ielem,:),idime),Ngp_l(inode,:),coordParMshRank(connecParOrigMshRank(ielem,inode),idime))
            end do
         end do
      end do
      deallocate(aux_1)

   end subroutine interpolateOriginalCoordinates

   subroutine create_working_lists_parallel(isPeriodic,mshRank,numElemsInRank,numNodesInRank,numBoundFacesInRank,globalIdSrlInRank,connecParOrigInRank,&
            nPerInRank,masSlaInRank,connecParWorkInRank,connecOrigBoundFacesInRank,connecBoundFacesInRank,numWorkingNodesInRank,workingNodesInRank)
      implicit none
      logical,intent(in) :: isPeriodic
      integer,intent(in) :: mshRank,numElemsInRank,numNodesInRank,numBoundFacesInRank,globalIdSrlInRank(numElemsInRank),connecParOrigInRank(numElemsInRank,nnode),connecOrigBoundFacesInRank(numBoundFacesInRank,npbou)
      integer,intent(in) :: nPerInRank,masSlaInRank(nPerInRank,2)
      integer,intent(out) :: connecParWorkInRank(numElemsInRank,nnode),numWorkingNodesInRank,connecBoundFacesInRank(numBoundFacesInRank,npbou)
      integer,intent(inout),allocatable :: workingNodesInRank(:)
      integer :: iNode,iNodeL,iElem,iBound,iPer,iAux,iPos
      integer :: iNodeL_Per,iNodeL_Per_Pair
      integer, allocatable :: aux_workingNodesInRank(:)

      !if(mpi_rank.eq.0) write(*,*) ' # Creating working lists...'

      !$acc kernels
      connecParWorkInRank(:,:)    = connecParOrigInRank(:,:) 
      connecBoundFacesInRank(:,:) = connecOrigBoundFacesInRank(:,:)   
      !$acc end kernels

      if(isPeriodic) then !do all stuff in case mesh is periodic
         !write(*,*) 'TO BE CHECK! subroutine create_working_list_parallel() for periodic cases'
         !----------------------------------------------------------------
         !-------------  CONNEC   ----------------------------------------

#if 1
         do iElem = 1,numElemsInRank
            do iNode = 1,nnode
               iNodeL = connecParWorkInRank(iElem,iNode)
               iPos = binarySearch_int_i(masSlaInRank(:,2),iNodeL)
               if(iPos.ne.0) then !this node is slave, change!
                  iNodeL_Per_Pair = masSlaInRank(iPos,1)
                  connecParWorkInRank(iElem,iNode) = iNodeL_Per_Pair
               end if
            end do
         end do
#else
         !$acc kernels
         do iElem = 1,numElemsInRank
            do iAux = 1,nnode
               iNodeL = connecParWorkInRank(iElem,iAux)
               !iNodeG = globalIdSrl(iNodeL)
               do iPer = 1,nPerInRank
                  iNodeL_Per = masSlaInRank(iPer,2)
                  if (iNodeL .eq. iNodeL_Per) then
                     iNodeL_Per_Pair = masSlaInRank(iPer,1)
                     connecParWorkInRank(iElem,iAux) = iNodeL_Per_Pair
                  end if
               end do
            end do
         end do
         !$acc end kernels
#endif

         do iBound = 1,numBoundFacesInRank
            do iNode = 1,npbou
               iNodeL = connecBoundFacesInRank(iBound,iNode)
               iPos = binarySearch_int_i(masSlaInRank(:,2),iNodeL)
               if(iPos.ne.0) then !this node is slave, change!
                  iNodeL_Per_Pair = masSlaInRank(iPos,1)
                  connecBoundFacesInRank(iBound,iNode) = iNodeL_Per_Pair
               end if
            end do
         end do

         !----------------------------------------------------------------
         !----------------WORKING NODES-----------------------------------

         numWorkingNodesInRank = numNodesInRank - nPerInRank

         allocate(aux_workingNodesInRank(numNodesInRank))
         allocate(workingNodesInRank(numWorkingNodesInRank))

         !$acc parallel loop
         do iNodeL = 1,numNodesInRank
            aux_workingNodesInRank(iNodeL) = iNodeL
         end do
         !$acc end parallel loop

#if 1
         do iNodeL = 1,numNodesInRank
            iPos = binarySearch_int_i(masSlaInRank(:,2),iNodeL)
            if(iPos.ne.0) then !this node is slave, change!
               aux_workingNodesInRank(iNodeL) = 0
            end if
         end do
#else
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
#endif
         !TODO: GPU-PARALLELIZE THIS LOOP
         iAux = 0
         do iNodeL = 1,numNodesInRank
            if(aux_workingNodesInRank(iNodeL).ne.0) then
               iAux = iAux+1
               workingNodesInRank(iAux) = aux_workingNodesInRank(iNodeL)
               !write(*,*)'wNP(',iAux,')',workingNodesInRank(iAux)
            end if
         end do

         deallocate(aux_workingNodesInRank)

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

   subroutine get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,int2comm,vectorOut)
      implicit none
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in) :: numMshRanksInMpiRank,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanksInMpiRank),int2comm(numMshRanksInMpiRank)
      integer,dimension(0:numMshRanks2Part-1),intent(out) :: vectorOut

      integer,dimension(0:numMshRanks2Part-1) :: auxVector2send
      integer :: window_id
      integer :: iMshRank,mshRank,mpiRank,auxNodeCnt

      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      auxVector2send(:)=0
      vectorOut(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         auxVector2send(mshRank) = int2comm(iMshRank)
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(auxVector2send(0),win_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0

      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vectorOut(mshRank),1,MPI_INTEGER,mpiRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------


   end subroutine get_vector_with_mshRank_values_for_numMshRanks2Part

   subroutine evalShapeFunctions_Ngp_l(Ngp_l)
      real(rp),intent(out) :: Ngp_l(ngaus,nnode)
      real(rp) :: s, t, z
      integer(4) :: igaus
      integer(4) :: atoIJK(nnode)
      integer(4) :: listHEX08((porder**ndime),2**ndime)
      real(rp) :: xgp(ngaus,ndime), wgp(ngaus)
      real(rp) :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1),dNgp_l(ndime,nnode,ngaus)

      !*********************************************************

      call set_hex64_lists(atoIJK,listHEX08)
      call chebyshev_hex(atoIJK,xgp,wgp)

      do igaus = 1,ngaus
         s = xgp(igaus,1)
         t = xgp(igaus,2)
         z = xgp(igaus,3)
         call hex64(s,t,z,atoIJK,Ngp(igaus,:),dNgp(:,:,igaus),Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
      end do

   end subroutine evalShapeFunctions_Ngp_l

!---------------------------------------------------------------------------------------------------------------------------
!-------------------------DEPRECATED CODE SECTION---------------------------------------------------------------------------
#if 0

#if 0
BEGIN WHEN READING LINE BY LINE.. FUCKING SLOW
      !----------------------------------------------------------------------------------------------------------
      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      !call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write

      do iElem=1,numElemsInRank
         elemToRead=listElemsRank(iElem)
         ms_offset(2) = elemToRead-1
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         !call h5dread_f(dset_id,dtype,connecRank(iElem,:),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5dread_f(dset_id,dtype,connecRank(iElem,:),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id)

         !connecRank(iElem,1:nnode) = auxConnec(:)
         !if(mpi_rank.eq.0)write(*,*) 'iElem',iElem,'e2r',elemToRead,'auxConnec',auxConnec(:)
         !if(mpi_rank.eq.0)write(*,*) 'connecRank',connecRank(iElem,:)
      end do

      !call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      !----------------------------------------------------------------------------------------------------------
END WHEN READING LINE BY LINE.. FUCKING SLOW
#endif

#if 0
BEGIN WHEN READING LINE BY LINE.. FUCKING SLOW
      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      !call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write

      do iNode=1,numNodesInRank
         nodeToRead=listNodesRank(iNode)
         ms_offset(2) = nodeToRead-1
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         !call h5dread_f(dset_id,dtype,auxCoords,fs_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5dread_f(dset_id,dtype,auxCoords,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id)
         coordNodesRank(iNode,:) = auxCoords(:)
         !if(mpi_rank.eq.0)write(*,*) 'n2r',nodeToRead,'coords',coordNodesRank(iNode,:)
      end do

      !call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
END WHEN READING LINE BY LINE.. FUCKING SLOW
#endif



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

      ! Boundary nodes section
      !------------------------------------------------------------------------------------------

      if(isPeriodic) then
         if(mpi_rank.eq.0)write(*,*) "--| Shit for periodic...."
         write(*,*) 'ABORTING, SHIT FOR PERIODIC YET TO BE IMPLEMENTED!'
         call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
         !call open_fix_bou_file(file_path,file_name)
         !call read_boundaries(numBoundFacesSrl,numElemsMpiRank,listElemsMpiRank,connecMpiRank,numBoundsMpiRank,listElemsBoundsMpiRank,boundFacesMpiRank)

         !call close_fix_bou_file()
      end if

      call close_geo_dat_file()

      !------------------------------------------------------------------------------------------

      if(mpi_rank.eq.0)write(*,*) "--| End of reading GEO.DAT file in parallel!"
      if(mpi_rank.eq.0)write(*,*) "--------------------------------------------"
    
   end subroutine read_geo_dat_file_in_parallel

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

      deallocate(rawNodeListRank)

   end subroutine read_elem_connec_and_nodes_coords_from_geo_file

   subroutine read_boundaries(numBoundFacesSrl,numElemsInRank,listElemsInRank,connecInRank,numBoundsInRank,maxBoundCodeInRank,listElemsBoundsInRank,listBoundFacesInRank,boundFacesCodesInRank,boundFacesInRank)
      implicit none
      integer,intent(in) :: numBoundFacesSrl,numElemsInRank,listElemsInRank(numElemsInRank),connecInRank(numElemsInRank,nnode)
      integer(4), intent(out) :: numBoundsInRank,maxBoundCodeInRank
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

      !write(*,*) '[',mpi_rank,'] numBoundsInRank',numBoundsInRank
      !do iElem=1,numElemsInRank
      !   if(listElemsBoundsInRank(iElem,1).ne.0) then
      !      write(*,*) 'i',iElem,'iG',listElemsInRank(iElem),'bounds',listElemsBoundsInRank(iElem,:)
      !   end if
      !end do

      allocate(listBoundFacesInRank(numBoundsInRank))
      allocate(boundFacesCodesInRank(numBoundsInRank))
      allocate(boundFacesInRank(numBoundsInRank,npbou))

      maxBoundCodeInRank = 0
      do iBound=1,numBoundsInRank
         listBoundFacesInRank(iBound)  = auxBoundFacesInRank(iBound,1)
         boundFacesCodesInRank(iBound) = auxBoundFacesInRank(iBound,2)
         boundFacesInRank(iBound,:)    = auxBoundFacesInRank(iBound,3:npbou+2)

         maxBoundCodeInRank = max(maxBoundCodeInRank,boundFacesCodesInRank(iBound))
         !write(*,*) '[',mpi_rank,']bfmpirank(',iBound,')',boundFacesInRank(iBound,:)
      end do

      deallocate(auxBoundFacesInRank)

   end subroutine read_boundaries

#endif
!---------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------


end module mod_meshConversorTool
