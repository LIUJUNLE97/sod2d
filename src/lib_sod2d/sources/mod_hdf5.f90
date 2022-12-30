module mod_hdf5
   use hdf5
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms
   use mod_mesh_boundaries
   use mod_comms_boundaries
   implicit none

   character(256) :: meshFile_h5_name,base_resultsFile_h5_name,base_avgResultsFile_h5_name
   integer(hid_t) :: meshFile_h5_id,resultsFile_h5_id

contains

   subroutine init_hdf5_interface()
      implicit none
      integer :: h5err

      !.init h5 interface
      call h5open_f(h5err)

      !call set_hdf5_meshFile_name(mesh_filePath,mesh_fileName)
      !call set_hdf5_baseResultsFile_name(results_filePath,results_fileName,mesh_fileName)
   end subroutine init_hdf5_interface

   subroutine end_hdf5_interface()
      implicit none
      integer :: h5err

      !close h5 interface
      call h5close_f(h5err)
   end subroutine end_hdf5_interface

   subroutine set_hdf5_meshFile_name(file_path,file_name,numRanks)
      implicit none
      character(len=*), intent(in) :: file_path,file_name
      integer,intent(in) :: numRanks
      character(len=12) :: aux_numRanks

      write(aux_numRanks,'(I0)') numRanks
      meshFile_h5_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'-'//trim(aux_numRanks)//'.h5'
   end subroutine set_hdf5_meshFile_name

   subroutine set_hdf5_baseResultsFile_name(res_filePath,res_fileName,mesh_fileName,numRanks)
      implicit none
      character(len=*), intent(in) :: res_filePath,res_fileName,mesh_fileName
      integer,intent(in) :: numRanks
      character(len=12) :: aux_numRanks

      write(aux_numRanks,'(I0)') mpi_size
      base_resultsFile_h5_name = trim(adjustl(res_filePath))//trim(adjustl(res_fileName))//'_'&
         //trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_'
      base_avgResultsFile_h5_name = trim(adjustl(res_filePath))//trim(adjustl(res_fileName))//'_AVG_'&
         //trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_'

   end subroutine set_hdf5_baseResultsFile_name

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

   subroutine read_alyaMesh_part_and_create_hdf5Mesh(gmsh_file_path,gmsh_file_name,periodic)
      implicit none     
      character(len=*), intent(in) :: gmsh_file_path,gmsh_file_name
      logical, intent(in) :: periodic

      !-- read the alya mesh fesh files in GMSH/ALYA FORMAT
      call read_alya_mesh_files(gmsh_file_path,gmsh_file_name,periodic)

      !----- Do mesh partitioning!
      call do_mesh_partitioning()

      !----- for boundaries
      if(totalNumBoundsSrl.ne.0) then
         call splitBoundary_inPar()
         call generate_boundary_mpi_comm_scheme()
      end if

      !----- Deallocate alya/gmsh arrays
      call deallocate_read_alya_mesh_arrays()

      !----- Create HDF5 File
      call create_hdf5_meshFile()

   end subroutine read_alyaMesh_part_and_create_hdf5Mesh

   subroutine create_hdf5_meshFile()
      implicit none
      integer(hid_t) :: file_id,plist_id,dset_id,dspace_id,mspace_id,group_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(256) :: groupname,dsetname

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(meshFile_h5_name,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)
      !call h5pget_driver_f(fapl_id, driver_id, hdferror)
      !call check("h5pget_driver_f", hdferror, nerrors)
      !--------------------------------------------------------------------------------
      groupname = '/Coords'
      call create_group_hdf5(file_id,groupname)
      !call h5gcreate_f(file_id,groupname,group_id,h5err)
      !call h5gclose_f(group_id, h5err)

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar

      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      dtype = H5T_NATIVE_REAL

      dsetname = '/Coords/X'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coord_x)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,1))

      dsetname = '/Coords/Y'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,2))

      dsetname = '/Coords/Z'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,3))

      !--------------------------------------------------------------------------------
      groupname = '/globalIds'
      call create_group_hdf5(file_id,groupname)

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar

      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      dtype = H5T_NATIVE_INTEGER

      dsetname = '/globalIds/globalIdSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrl)

      dsetname = '/globalIds/globalIdPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdPar)

      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      dsetname = '/globalIds/elemGid'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,elemGid)

      !--------------------------------------------------------------------------------
      !save connectivity
      call save_connectivity_hdf5(file_id)

      !--------------------------------------------------------------------------------
      !save parallel data
      if(mpi_size.ge.2) then
         call save_parallel_data_hdf5(file_id)
      end if
      !--------------------------------------------------------------------------------
      !save periodic data
      if(isMeshPeriodic) then
         call save_periodic_data_hdf5(file_id)
      end if
      !--------------------------------------------------------------------------------
      !save boundary data
      if(isMeshBoundaries) then
         call save_boundary_data_hdf5(file_id)
         if(mpi_size.ge.2) then
            call save_parallel_data_boundary_hdf5(file_id)
         end if
      end if
      !--------------------------------------------------------------------------------

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine create_hdf5_meshFile

   subroutine create_hdf5_meshFile_from_tool(hdf5_file_id)
      implicit none
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer(hid_t) :: plist_id
      integer :: h5err
      !-------------------------------------------------------------------------------

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(meshFile_h5_name,H5F_ACC_TRUNC_F,hdf5_file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

   end subroutine

   subroutine close_hdf5_meshFile_from_tool(hdf5_file_id)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      integer :: h5err

      !close the file.
      call h5fclose_f(hdf5_file_id,h5err)
   end subroutine close_hdf5_meshFile_from_tool

   subroutine create_hdf5_groups_datasets_in_meshFile_from_tool(hdf5_file_id,isPeriodic,isBoundaries,numMshRanks2Part,numElemsGmsh,numNodesParTotal,&
               vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer,intent(in) :: numMshRanks2Part,numElemsGmsh,numNodesParTotal
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      integer(hid_t) :: dset_id,dspace_id,group_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims
      integer :: ds_rank,ms_rank,h5err
      character(256) :: groupname,dsetname

      !--------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numNodesParTotal
      
      dtype = H5T_NATIVE_REAL

      groupname = '/Coords'
      call create_group_hdf5(hdf5_file_id,groupname)
      dsetname = '/Coords/X'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/Coords/Y'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/Coords/Z'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
  
      dtype = H5T_NATIVE_INTEGER

      groupname = '/globalIds'
      call create_group_hdf5(hdf5_file_id,groupname)
      dsetname = '/globalIds/globalIdSrl'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/globalIds/globalIdPar'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = numElemsGmsh
      dsetname = '/globalIds/elemGid'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      call create_groups_datasets_connectivity_workingNodes_hdf5(hdf5_file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)

      !if(mpi_size.ge.2) then
         call create_groups_datasets_parallel_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumMshRanksWithComms,vecNumNodesToCommMshRank)
         !if(isBoundaries) then
            !call create_groups_datasets_parallel_data_boundary_firstBlock_hdf5(hdf5_file_id,numMshRanks2Part)
         !end if
      !end if

      !if(isPeriodic) then
         !call create_groups_datasets_periodic_data_firstBlock_hdf5(hdf5_file_id,numMshRanks2Part)
      !end if

      if(isBoundaries) then
         call create_groups_datasets_boundary_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
         !if(mpi_size.ge.2) then
         call create_groups_datasets_parallel_data_boundary_hdf5(hdf5_file_id,numMshRanks2Part,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank)
         !end if
      end if


      !-----------------------------------------------------------------------------------------------------
#if 0
      do iMshRank = 1,numMshRanks

         ms_rank = 1
         ms_dims(1) = numNodesRankPar
         ms_offset(1) = rankNodeStart-1

         dsetname = '/Coords/X'
         call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,1))

         dsetname = '/Coords/Y'
         call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,2))

         dsetname = '/Coords/Z'
         call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,3))

         !--------------------------------------------------------------------------------

         ms_rank = 1
         ms_dims(1) = numNodesRankPar
         ms_offset(1) = rankNodeStart-1

         dsetname = '/globalIds/globalIdSrl'
         call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrl)

         dsetname = '/globalIds/globalIdPar'
         call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdPar)

         ms_dims(1) = numElemsRankPar
         ms_offset(1) = rankElemStart-1

         dsetname = '/globalIds/elemGid'
         call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,elemGid)

         !--------------------------------------------------------------------------------
         !save connectivity
         call save_connectivity_hdf5(file_id)

         !--------------------------------------------------------------------------------
         !save parallel data
         if(mpi_size.ge.2) then
            call save_parallel_data_hdf5(file_id)
         end if
         !--------------------------------------------------------------------------------
         !save periodic data
         if(isMeshPeriodic) then
            call save_periodic_data_hdf5(file_id)
         end if
         !--------------------------------------------------------------------------------
         !save boundary data
         if(isMeshBoundaries) then
            call save_boundary_data_hdf5(file_id)
            if(mpi_size.ge.2) then
               call save_parallel_data_boundary_hdf5(file_id)
            end if
         end if
         !--------------------------------------------------------------------------------

      end do
#endif


   end subroutine create_hdf5_groups_datasets_in_meshFile_from_tool

   subroutine create_groups_datasets_connectivity_workingNodes_hdf5(file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part,numElemsGmsh
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      groupname = trim('/Connectivity')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      !-----------------------------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numElemsGmsh*nnode

      dsetname = '/Connectivity/connecParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Connectivity/connecParWork'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Connectivity/connecVTK'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !-----------------------------------------------------------------------------------------------------
      ds_dims(1) = numMshRanks2Part

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !-----------------------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumWorkingNodes(mshRank)
      end do
      ds_dims(1) = accumVal

      dsetname = '/Connectivity/workingNodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !-----------------------------------------------------------------------------------------------------

   end subroutine create_groups_datasets_connectivity_workingNodes_hdf5
   
   subroutine create_groups_datasets_parallel_data_hdf5(file_id,numMshRanks2Part,vecNumMshRanksWithComms,vecNumNodesToCommMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumMshRanksWithComms,vecNumNodesToCommMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      groupname = trim('/Parallel_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part

      dsetname = '/Parallel_data/rankNodeStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/rankNodeEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/rankElemStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/rankElemEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !----------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumMshRanksWithComms(mshRank)
      end do
      ds_dims(1) = accumVal

      dsetname = '/Parallel_data/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !----------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumNodesToCommMshRank(mshRank)
      end do
      ds_dims(1) = accumVal

      dsetname = '/Parallel_data/matrixCommScheme_iNodeL'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/matrixCommScheme_iNodeGSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data/matrixCommScheme_ngbRank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !----------------------------------------------------------------------------------------

   end subroutine create_groups_datasets_parallel_data_hdf5

   subroutine create_groups_datasets_parallel_data_boundary_hdf5(file_id,numMshRanks2Part,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      groupname = trim('/Parallel_data_boundary')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part

      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !----------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecBndNumMshRanksWithComms(mshRank)
      end do
      ds_dims(1) = accumVal

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !----------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecBndNumNodesToCommMshRank(mshRank)
      end do
      ds_dims(1) = accumVal

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeL'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeGSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Parallel_data_boundary/matrixCommScheme_ngbRank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !----------------------------------------------------------------------------------------

   end subroutine create_groups_datasets_parallel_data_boundary_hdf5

   subroutine create_groups_datasets_periodic_data_firstBlock_hdf5(file_id,mshRank,numMshRanks2Part)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: mshRank,numMshRanks2Part
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Periodic_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mshRank
      allocate(aux_array(1))

      dsetname = '/Periodic_data/nPerRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=nPerRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set nPerRankPar of all ranks
      dsetname = '/Periodic_data/nPerRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=nPerRankPar

      dsetname = '/Periodic_data/masSlaRankPar1'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaRankPar(:,1))

      dsetname = '/Periodic_data/masSlaRankPar2'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaRankPar(:,2))
      !--------------------------------------------------------------------------------------------------------

   end subroutine create_groups_datasets_periodic_data_firstBlock_hdf5

   subroutine create_groups_datasets_boundary_data_hdf5(file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,h5err      
      integer :: mshRank,accumVal

      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = 1

      dsetname = '/Boundary_data/numBoundCodes'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = numMshRanks2Part

      dsetname = '/Boundary_data/numBoundsRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Boundary_data/ndofRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !--------------------------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumBoundFacesMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      !--------------------------------------------------------------------------------------------------------
   
      dsetname = '/Boundary_data/bouCodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = ds_dims(1)*npbou

      dsetname = '/Boundary_data/boundPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Boundary_data/boundParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !--------------------------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumDoFMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      !--------------------------------------------------------------------------------------------------------

      dsetname = '/Boundary_data/ldofPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !--------------------------------------------------------------------------------------------------------
      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumBoundaryNodesMshRank(mshRank)
      end do
      ds_dims(1) = accumVal
      !--------------------------------------------------------------------------------------------------------

      dsetname = '/Boundary_data/lbnodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !--------------------------------------------------------------------------------------------------------
   end subroutine create_groups_datasets_boundary_data_hdf5

#if 0
   subroutine save_boundary_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Boundary_data/numBoundsRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundsRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundCodes'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundCodes
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/ndofRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=ndofRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundaryNodesRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numBoundsRankPar of all ranks
      dsetname = '/Boundary_data/numBoundsRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numBoundsRankPar

      !write(*,*) 'ms_offset ',ms_offset(1),' ms_dims(1) ',ms_dims(1),' ds_dims ',ds_dims(1)

      dsetname = '/Boundary_data/bouCodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bouCodesPar)

      ds_dims(1) = ds_dims(1)*npbou!totalNumBoundsSrl*npbou
      ms_dims(1) = ms_dims(1)*npbou!numBoundsRankPar*npbou
      ms_offset(1) = ms_offset(1)*npbou

      deallocate(aux_array)
      allocate(aux_array(numBoundsRankPar*npbou))

      i=1
      do iBound=1,numBoundsRankPar
         do m=1,npbou
            aux_array(i)=boundPar(iBound,m)
            i=i+1
         end do
      end do

      dsetname = '/Boundary_data/boundPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set ndofRankPar of all ranks
      dsetname = '/Boundary_data/ndofRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=ndofRankPar

      dsetname = '/Boundary_data/ldofPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ldofPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numBoundaryNodesRankPar of all ranks
      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numBoundaryNodesRankPar

      dsetname = '/Boundary_data/lbnodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,lbnodesPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
   end subroutine save_boundary_data_hdf5
#endif


   subroutine write_mshRank_data_in_hdf5_meshFile_from_tool(hdf5_file_id,mshRank,numMshRanks2Part,isPeriodic,isBoundaries,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numNodesParTotal,&
         numElemsMshRank,mshRankElemStart,mshRankElemEnd,mshRankNodeStart,mshRankNodeEnd,numNodesMshRank,numWorkingNodesMshRank,numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank,maxBoundCode,&
         elemGidMshRank,globalIdSrlMshRank,globalIdParMshRank,connecVTKMshRank,connecParOrigMshRank,connecParWorkMshRank,coordParMshRank,workingNodesMshRank,&
         boundaryNodesMshRank,dofNodesMshRank,boundFacesCodesMshRank,boundFacesOrigMshRank,boundFacesMshRank,&
         numNodesToCommMshRank,numMshRanksWithComms,matrixCommSchemeMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank,&
         bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms,bnd_matrixCommSchemeMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank,&
         vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer(4),intent(in) :: mshRank,numMshRanks2Part,numElemsGmsh,numNodesGmsh,numBoundFacesGmsh,numNodesParTotal!,numPerNodes
      integer(4),intent(in) :: numElemsMshRank,mshRankElemStart,mshRankElemEnd
      integer(4),intent(in) :: mshRankNodeStart,mshRankNodeEnd,numNodesMshRank,numWorkingNodesMshRank
      integer(4),intent(in) :: numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank
      integer(4),intent(in) :: maxBoundCode

      integer(4),intent(in) :: elemGidMshRank(numElemsMshRank),globalIdSrlMshRank(numNodesMshRank),globalIdParMshRank(numNodesMshRank),workingNodesMshRank(numWorkingNodesMshRank)
      integer(4),intent(in) :: connecVTKMshRank(numElemsMshRank*nnode)
      integer(4),intent(in) :: connecParOrigMshRank(numElemsMshRank,nnode),connecParWorkMshRank(numElemsMshRank,nnode)
      real(4),intent(in)    :: coordParMshRank(numNodesMshRank,3)

      integer(4),intent(in) :: boundaryNodesMshRank(numBoundaryNodesMshRank),dofNodesMshRank(numDoFMshRank)
      integer(4),intent(in) :: boundFacesCodesMshRank(numBoundFacesMshRank),boundFacesOrigMshRank(numBoundFacesMshRank,npbou),boundFacesMshRank(numBoundFacesMshRank,npbou)
      integer(4),intent(in) :: numNodesToCommMshRank,numMshRanksWithComms
      integer(4),intent(in),dimension(numNodesToCommMshRank,3) :: matrixCommSchemeMshRank
      integer(4),intent(in),dimension(numMshRanksWithComms) :: commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank
      integer(4),intent(in) :: bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms
      integer(4),intent(in),dimension(bnd_numNodesToCommMshRank,3) :: bnd_matrixCommSchemeMshRank
      integer(4),intent(in),dimension(bnd_numMshRanksWithComms) :: bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank

      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      !-------------------------------------------------------------------------------------------------------------------------------
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer :: i,m,iBound,iElemL,ms_rank
      integer(hssize_t), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)
      !-------------------------------------------------------------------------------------------------------------------------------

      ms_rank = 1
      ms_dims(1) = numNodesMshRank
      ms_offset(1) = mshRankNodeStart-1

      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrlMshRank)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdParMshRank)

      ms_dims(1) = numElemsMshRank
      ms_offset(1) = mshRankElemStart-1

      dsetname = '/globalIds/elemGid'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,elemGidMshRank)

      !---------------------------------------------------------------------------------------------------------------------
      ms_rank = 1
      ms_dims(1) = numNodesMshRank
      ms_offset(1) = mshRankNodeStart-1

      dsetname = '/Coords/X'
      call write_dataspace_fp32_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,coordParMshRank(:,1))

      dsetname = '/Coords/Y'
      call write_dataspace_fp32_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,coordParMshRank(:,2))

      dsetname = '/Coords/Z'
      call write_dataspace_fp32_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,coordParMshRank(:,3))

      !---------------------------------------------------------------------------------------------------------------------

      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mshRank
      allocate(aux_array(1))

      dsetname = '/Parallel_data/rankNodeStart'
      aux_array(1)=mshRankNodeStart
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankNodeEnd'
      aux_array(1)=mshRankNodeEnd
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankElemStart'
      aux_array(1)=mshRankElemStart
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankElemEnd'
      aux_array(1)=mshRankElemEnd
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/numRanksWithComms'
      aux_array(1)=numMshRanksWithComms
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/numNodesToComm'
      aux_array(1)=numNodesToCommMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)

      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecNumMshRanksWithComms(i)
      end do
      ms_dims(1)=numMshRanksWithComms
      !write(*,*) '[',mshRank,']ms_offset',ms_offset(1),'ms_dims(1)',ms_dims(1)!,' ds_dims ',ds_dims(1)

      dsetname = '/Parallel_data/ranksToComm'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,ranksToCommMshRank)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInLocMshRank)

      dsetname = '/Parallel_data/commsMemSize'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemSizeMshRank)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInNgbMshRank)

      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecNumNodesToCommMshRank(i)
      end do
      ms_dims(1)=numNodesToCommMshRank
      !write(*,*) '[',mshRank,']ms_offset',ms_offset(1),'ms_dims(1)',ms_dims(1)!,' ds_dims ',ds_dims(1)

      dsetname = '/Parallel_data/matrixCommScheme_iNodeL'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommSchemeMshRank(:,1))

      dsetname = '/Parallel_data/matrixCommScheme_iNodeGSrl'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommSchemeMshRank(:,2))

      dsetname = '/Parallel_data/matrixCommScheme_ngbRank'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommSchemeMshRank(:,3))

      !---------------------------------------------------------------------------------------------------------------------

      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mshRank
      allocate(aux_array(1))

      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      aux_array(1)=bnd_numMshRanksWithComms
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      aux_array(1)=bnd_numNodesToCommMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)

      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecBndNumMshRanksWithComms(i)
      end do
      ms_dims(1)=bnd_numMshRanksWithComms

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_ranksToCommMshRank)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInLocMshRank)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemSizeMshRank)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInNgbMshRank)

      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecBndNumNodesToCommMshRank(i)
      end do
      ms_dims(1)=bnd_numNodesToCommMshRank

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeL'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommSchemeMshRank(:,1))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeGSrl'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommSchemeMshRank(:,2))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_ngbRank'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommSchemeMshRank(:,3))

      !------------------------------------------------------------------------------------------------------------------------

      !------------------------------------------------------------------------------------------------------------------------

      allocate(aux_array(1))
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mshRank

      if(mshRank.eq.0) then
         dsetname = '/Boundary_data/numBoundCodes'
         aux_array(1)=maxBoundCode
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      end if

      dsetname = '/Boundary_data/numBoundsRankPar'
      aux_array(1)=numBoundFacesMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/ndofRankPar'
      aux_array(1)=numDoFMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      aux_array(1)=numBoundaryNodesMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------
      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecNumDoFMshRank(i)
      end do
      ms_dims(1)=numDoFMshRank

      dsetname = '/Boundary_data/ldofPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,dofNodesMshRank)
      !------------------------------------------------------------------------------------------------------------------------
      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecNumBoundaryNodesMshRank(i)
      end do
      ms_dims(1)=numBoundaryNodesMshRank

      dsetname = '/Boundary_data/lbnodesPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,boundaryNodesMshRank)
      !------------------------------------------------------------------------------------------------------------------------
      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecNumBoundFacesMshRank(i)
      end do
      ms_dims(1)=numBoundFacesMshRank

      dsetname = '/Boundary_data/bouCodesPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,boundFacesCodesMshRank)
      !------------------------------------------------------------------------------------------------------------------------
      allocate(aux_array(numBoundFacesMshRank*npbou))

      i=1
      do iBound=1,numBoundFacesMshRank
         do m=1,npbou
            aux_array(i)=boundFacesMshRank(iBound,m)
            i=i+1
         end do
      end do

      ms_dims(1) = ms_dims(1)*npbou
      ms_offset(1) = ms_offset(1)*npbou

      dsetname = '/Boundary_data/boundPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      i=1
      do iBound=1,numBoundFacesMshRank
         do m=1,npbou
            aux_array(i)=boundFacesOrigMshRank(iBound,m)
            i=i+1
         end do
      end do

      dsetname = '/Boundary_data/boundParOrig'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mshRank
      allocate(aux_array(1))

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      aux_array(1)=numWorkingNodesMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------
      !  SAVING connecParOrig(:,:)
      allocate(aux_array(numElemsMshRank*nnode))
      ms_rank = 1
      ms_dims(1) = numElemsMshRank*nnode
      ms_offset(1) = (mshRankElemStart-1)*nnode
      !-----------------------------------------------------------------------------------------------------
      i=1
      do iElemL=1,numElemsMshRank
         do m=1,nnode
            aux_array(i)=connecParOrigMshRank(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParOrig'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParWork(:,:)
      i=1
      do iElemL=1,numElemsMshRank
         do m=1,nnode
            aux_array(i)=connecParWorkMshRank(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParWork'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      deallocate(aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      dsetname = '/Connectivity/connecVTK'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,connecVTKMshRank)
      !-----------------------------------------------------------------------------------------------------
      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+vecNumWorkingNodes(i)
      end do
      ms_dims(1)=numWorkingNodesMshRank

      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,workingNodesMshRank)
      !-----------------------------------------------------------------------------------------------------




   end subroutine write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine load_hdf5_meshFile()
      implicit none     
      character(256) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id,fspace_id
      integer :: h5err
      integer(hsize_t),dimension(1) :: fs_dims,fs_maxdims

      !.init h5 interface
      !call h5open_f(h5err)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      if(mpi_rank.eq.0) write(*,*) '# Loading hdf5 mesh: ',meshFile_h5_name

      call h5fopen_f(meshFile_h5_name, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !-----------------------------------------------------------------------------------------------

      dsetname = '/globalIds/globalIdSrl'
      call h5dopen_f(file_id, dsetname, dset_id, h5err)
      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)
      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      totalNumNodesPar = fs_dims(1)

      dsetname = '/globalIds/elemGid'
      call h5dopen_f(file_id, dsetname, dset_id, h5err)
      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)
      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      totalNumElements = fs_dims(1)

      !write(*,*) 'totalNumNodesPar ',totalNumNodesPar,' totalNumElements ',totalNumElements

#if 0
      ds_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1
      dtype = H5T_NATIVE_INTEGER

      dsetname = '/globalIds/globalIdSrl'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrl)

      dsetname = '/globalIds/globalIdPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdPar)

      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      dsetname = '/globalIds/elemGid'
#endif

      !-----------------------------------------------------------------------------------------------
      !load the parallel data
      if(mpi_size.ge.2) then
         call load_parallel_data_hdf5(file_id)
      else !only 1 rank
         numNodesRankPar = totalNumNodesPar
         numElemsRankPar = totalNumElements
         rankNodeStart = 1
         rankNodeEnd = numNodesRankPar
         rankElemStart = 1
         rankElemEnd = numElemsRankPar
      end if
      !-----------------------------------------------------------------------------------------------
      !load periodic data
      call load_periodic_data_hdf5(file_id)

      !-----------------------------------------------------------------------------------------------
      !load boundary data
      call load_boundary_data_hdf5(file_id)
      if((isMeshBoundaries).and.(mpi_size.ge.2)) then
         call load_parallel_data_boundary_hdf5(file_id)
      end if

      !-----------------------------------------------------------------------------------------------
      !load the coordinates
      call load_coordinates_hdf5(file_id)

      !-----------------------------------------------------------------------------------------------
      !load connectivity
      call load_connectivity_hdf5(file_id)

      !--------------------------------------------------------------------------------
      !load globalIds
      call load_globalIds_hdf5(file_id)

      !close h5 file
      call h5fclose_f(file_id,h5err)

   end subroutine load_hdf5_meshFile

   subroutine close_hdf5_meshFile()
      implicit none
      integer :: h5err

      !close the file.
      call h5fclose_f(meshFile_h5_id,h5err)
   end subroutine close_hdf5_meshFile

   subroutine create_group_hdf5(file_id,groupname)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*), intent(in) :: groupname
      integer(hid_t) :: group_id
      integer :: h5err
      
      call h5gcreate_f(file_id,groupname,group_id,h5err)
      call h5gclose_f(group_id, h5err)
   end subroutine create_group_hdf5

   subroutine create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id
      integer :: h5err

      ! Create the data space for the  dataset. 
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      ! Create the dataset with default properties.
      call h5dcreate_f(file_id,dsetname,dtype,dspace_id,dset_id,h5err)
      
      !write(*,*) 'create dsetname ',dsetname, ' dset_id ',dset_id,' dspace_id ',dspace_id

      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine create_dataspace_hdf5

   subroutine create_chunked_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,chunk_dims,dtype)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank),chunk_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id,plist_id
      integer :: h5err

      ! Create the data space for the  dataset. 
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)

   call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,h5err)
   call h5pset_chunk_f(plist_id,ds_rank,chunk_dims,h5err)
   call h5dcreate_f(file_id, dsetname,dtype,dspace_id,dset_id, h5err,plist_id)

      ! Create the dataset with default properties.
      !call h5dcreate_f(file_id,dsetname,dtype,dspace_id,dset_id,h5err)

      !write(*,*) 'create dsetname ',dsetname, ' dset_id ',dset_id,' dspace_id ',dspace_id
   call h5pclose_f(plist_id,h5err)
      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine create_chunked_dataspace_hdf5

   subroutine write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      real(4),intent(in) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_REAL


      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dwrite_f(dset_id,dtype,data,fs_dims,h5err,&
                      file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine write_dataspace_fp32_hyperslab_parallel

   subroutine write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      real(8),intent(in) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_DOUBLE

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dwrite_f(dset_id,dtype,data,fs_dims,h5err,&
                      file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine write_dataspace_fp64_hyperslab_parallel

   subroutine write_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      integer(1),intent(in) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_STD_U8LE

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dwrite_f(dset_id,dtype,data,fs_dims,h5err,&
                      file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine write_dataspace_int1_hyperslab_parallel

   subroutine read_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      integer(1),intent(out) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_STD_I8LE

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine read_dataspace_int1_hyperslab_parallel

   subroutine write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      integer(4),intent(in) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_INTEGER

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dwrite_f(dset_id,dtype,data,fs_dims,h5err,&
                      file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine write_dataspace_int4_hyperslab_parallel

   subroutine read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      integer(4),intent(out) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_INTEGER

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine read_dataspace_int4_hyperslab_parallel

   subroutine write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      integer(8),intent(in) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_STD_I64LE

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dwrite_f(dset_id,dtype,data,fs_dims,h5err,&
                      file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine write_dataspace_int8_hyperslab_parallel

   subroutine read_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      integer(8),intent(out) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_STD_I64LE

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine read_dataspace_int8_hyperslab_parallel

   subroutine read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      real(4),intent(out) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_REAL

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine read_dataspace_fp32_hyperslab_parallel

   subroutine read_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset 
      real(8),intent(out) :: data(:)
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_DOUBLE

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine read_dataspace_fp64_hyperslab_parallel

   subroutine overwrite_coordinates_hdf5()
      implicit none
      integer(hid_t) :: file_id,plist_id,dtype
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hssize_t), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(128) :: dsetname

      !---------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! open file collectively
      call h5fopen_f(meshFile_h5_name, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot open meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !-------------------------------------------
      dtype = H5T_NATIVE_REAL
      ds_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1
      !-------------------------------------------

      dsetname = '/Coords/X'
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,1))

      dsetname = '/Coords/Y'
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,2))

      dsetname = '/Coords/Z'
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,3))


      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine overwrite_coordinates_hdf5

   subroutine save_connectivity_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank
      integer :: iElemL,i,m,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Connectivity')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = totalNumElements*nnode
      ms_rank = 1
      ms_dims(1) = numElemsRankPar*nnode
      ms_offset(1) = (rankElemStart-1)*nnode

      allocate(aux_array(numElemsRankPar*nnode))
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParOrig(:,:)
      i=1
      do iElemL=1,numElemsRankPar
         do m=1,nnode
            aux_array(i)=connecParOrig(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParWork(:,:)
      i=1
      do iElemL=1,numElemsRankPar
         do m=1,nnode
            aux_array(i)=connecParWork(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParWork'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      dsetname = '/Connectivity/connecVTK'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,connecVTK)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParCGNS(:)
      dsetname = '/Connectivity/connecCGNS'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,connecCGNS)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      !  SAVING numWorkingNodesRankPar
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numWorkingNodesRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !--------------------------------------------------------------------------------------------------------
      !  SAVING workingNodesPar(:)
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numWorkingNodesRankPar of all ranks
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numWorkingNodesRankPar

      dsetname = '/Connectivity/workingNodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,workingNodesPar)
      !--------------------------------------------------------------------------------------------------------

   end subroutine save_connectivity_hdf5

   subroutine load_connectivity_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank
      integer :: iElemL,i,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading connectivity data hdf5...'

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = totalNumElements*nnode
      ms_rank = 1
      ms_dims(1) = numElemsRankPar*nnode
      ms_offset(1) = (rankElemStart-1)*nnode

      !-------------------------------------------------------------------------------------------------------
      allocate( connecParOrig(numElemsRankPar,nnode) )
      allocate( connecParWork(numElemsRankPar,nnode) )
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecParOrig(:,:)
      allocate(aux_array(numElemsRankPar*nnode))

      dsetname = '/Connectivity/connecParOrig'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      i=1
      do iElemL=1,numElemsRankPar
         do m=1,nnode
            connecParOrig(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecParWork(:,:)
      dsetname = '/Connectivity/connecParWork'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      i=1
      do iElemL=1,numElemsRankPar
         do m=1,nnode
            connecParWork(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do

      deallocate(aux_array)
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecVTK(:)
      allocate( connecVTK(numElemsRankPar*nnode) )
      dsetname = '/Connectivity/connecVTK'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,connecVTK)
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecCGNS(:)
      allocate( connecCGNS(numElemsRankPar*nnode) )
      dsetname = '/Connectivity/connecCGNS'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,connecCGNS)

      !-------------------------------------------------------------------------------------------------------
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))
      !-------------------------------------------------------------------------------------------------------
      !LOADING numWorkingNodesRankPar
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numWorkingNodesRankPar=aux_array(1)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      !LOADING workingNodesPar
      allocate(workingNodesPar(numWorkingNodesRankPar))
      allocate(aux_array(mpi_size))
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numWorkingNodesRankPar of all ranks
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      ds_dims(1)=0
      do i=1,mpi_size
         ds_dims(1)=ds_dims(1)+aux_array(i)
      end do

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numWorkingNodesRankPar

      dsetname = '/Connectivity/workingNodesPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,workingNodesPar)
      deallocate(aux_array)
      !-------------------------------------------------------------------------------------------------------

   end subroutine load_connectivity_hdf5

   subroutine save_parallel_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Parallel_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Parallel_data/rankNodeStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankNodeStart
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankNodeEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankNodeEnd
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankElemStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankElemStart
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankElemEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankElemEnd
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numRanksWithComms
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numNodesToComm
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)
      
      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numRanksWithComms

      !write(*,*) 'ms_offset ',ms_offset(1),' ms_dims(1) ',ms_dims(1),' ds_dims ',ds_dims(1)

      dsetname = '/Parallel_data/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ranksToComm)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInLoc)

      dsetname = '/Parallel_data/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemSize)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInNgb)

      ds_dims(1) = mpi_size
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      
      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numNodesToComm

      dsetname = '/Parallel_data/matrixCommScheme_iNodeL'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,1))

      dsetname = '/Parallel_data/matrixCommScheme_iNodeGSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,2))

      dsetname = '/Parallel_data/matrixCommScheme_ngbRank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,3))


      deallocate(aux_array)

   end subroutine save_parallel_data_hdf5

   subroutine load_parallel_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Parallel_data/rankNodeStart'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      rankNodeStart=aux_array(1)

      dsetname = '/Parallel_data/rankNodeEnd'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      rankNodeEnd=aux_array(1)

      dsetname = '/Parallel_data/rankElemStart'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      rankElemStart=aux_array(1)

      dsetname = '/Parallel_data/rankElemEnd'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      rankElemEnd=aux_array(1)

      numNodesRankPar = rankNodeEnd - rankNodeStart + 1
      numElemsRankPar  = rankElemEnd - rankElemStart + 1

      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numRanksWithComms=aux_array(1)

      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      numNodesToComm=aux_array(1)

      !write(*,*) 'rns ',rankNodeStart,' rne ',rankNodeEnd,' res ',rankElemStart,' ree ',rankElemEnd
      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)
      
      ds_dims(1)=0
      do i=1,mpi_size
         ds_dims(1)=ds_dims(1)+aux_array(i)
      end do

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numRanksWithComms

      allocate(ranksToComm(numRanksWithComms))
      allocate(commsMemPosInLoc(numRanksWithComms))
      allocate(commsMemPosInNgb(numRanksWithComms))
      allocate(commsMemSize(numRanksWithComms))

      dsetname = '/Parallel_data/ranksToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ranksToComm)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInLoc)

      dsetname = '/Parallel_data/commsMemSize'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemSize)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInNgb)
      
      ds_dims(1) = mpi_size
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      
      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ds_dims(1)=0
      do i=1,mpi_size
         ds_dims(1)=ds_dims(1)+aux_array(i)
      end do

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numNodesToComm

      allocate(matrixCommScheme(numNodesToComm,3))

      dsetname = '/Parallel_data/matrixCommScheme_iNodeL'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,1))

      dsetname = '/Parallel_data/matrixCommScheme_iNodeGSrl'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,2))

      dsetname = '/Parallel_data/matrixCommScheme_ngbRank'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,3))

      deallocate(aux_array)

   end subroutine load_parallel_data_hdf5

!--------------------------------------------------------------------------------------------------------------------------------

   subroutine save_parallel_data_boundary_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Parallel_data_boundary')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=bnd_numRanksWithComms
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=bnd_numNodesToComm
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      
      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=bnd_numRanksWithComms

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_ranksToComm)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInLoc)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemSize)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInNgb)

      ds_dims(1) = mpi_size
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      
      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=bnd_numNodesToComm

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeL'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,1))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeGSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,2))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_ngbRank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,3))

      deallocate(aux_array)

   end subroutine save_parallel_data_boundary_hdf5

   subroutine load_parallel_data_boundary_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      bnd_numRanksWithComms=aux_array(1)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      bnd_numNodesToComm=aux_array(1)

      !write(*,*) 'rns ',rankNodeStart,' rne ',rankNodeEnd,' res ',rankElemStart,' ree ',rankElemEnd
      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)
      
      ds_dims(1)=0
      do i=1,mpi_size
         ds_dims(1)=ds_dims(1)+aux_array(i)
      end do

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=bnd_numRanksWithComms

      allocate(bnd_ranksToComm(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInLoc(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInNgb(bnd_numRanksWithComms))
      allocate(bnd_commsMemSize(bnd_numRanksWithComms))

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_ranksToComm)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInLoc)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemSize)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInNgb)
      
      ds_dims(1) = mpi_size
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      
      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ds_dims(1)=0
      do i=1,mpi_size
         ds_dims(1)=ds_dims(1)+aux_array(i)
      end do

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=bnd_numNodesToComm

      allocate(bnd_matrixCommScheme(bnd_numNodesToComm,3))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeL'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,1))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeGSrl'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,2))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_ngbRank'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,3))

      deallocate(aux_array)

   end subroutine load_parallel_data_boundary_hdf5

!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------

   subroutine save_periodic_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)
      !PERIODIC SHIT NOT SAVED (NOT NEEDED): masSlaSrl(:,:) nPerSrl

      groupname = trim('/Periodic_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Periodic_data/nPerRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=nPerRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set nPerRankPar of all ranks
      dsetname = '/Periodic_data/nPerRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=nPerRankPar

      dsetname = '/Periodic_data/masSlaRankPar1'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaRankPar(:,1))

      dsetname = '/Periodic_data/masSlaRankPar2'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaRankPar(:,2))
      !--------------------------------------------------------------------------------------------------------

   end subroutine save_periodic_data_hdf5

   subroutine load_periodic_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)
      logical :: isPeriodicFolder

      !ELS WORKING NODES SON SEMPRE, NO NOMES PERIODIC
      !PENSAR SI POSAR AQUI EL CONNECORIG!
      groupname = trim('/Periodic_data')

      call h5lexists_f(file_id,groupname,isPeriodicFolder,h5err)

      if(mpi_rank.eq.0) write(*,*) 'Loading Periodic data hdf5. -> isPeriodic:',isPeriodicFolder

      if(isPeriodicFolder) then
         isMeshPeriodic = .true.

         dtype = H5T_NATIVE_INTEGER
         ds_rank = 1
         ds_dims(1) = mpi_size
         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = mpi_rank
         allocate(aux_array(1))

         !--------------------------------------------------------------------------------------------------------
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         nPerRankPar=aux_array(1)

         allocate(masSlaRankPar(nPerRankPar,2))
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
         allocate(aux_array(mpi_size))
         ms_dims(1) = mpi_size
         ms_offset(1) = 0
         !read data set numBoundsRankPar of all ranks
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+aux_array(i)
         end do
         ms_dims(1)=nPerRankPar

         dsetname = '/Periodic_data/masSlaRankPar1'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaRankPar(:,1))

         dsetname = '/Periodic_data/masSlaRankPar2'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaRankPar(:,2))
         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
      else
         isMeshPeriodic = .false.
      end if

   end subroutine load_periodic_data_hdf5
!--------------------------------------------------------------------------------------------------------------------------------
   subroutine save_boundary_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)

      dtype = H5T_NATIVE_INTEGER
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      allocate(aux_array(1))

      dsetname = '/Boundary_data/numBoundsRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundsRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundCodes'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundCodes
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/ndofRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=ndofRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundaryNodesRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numBoundsRankPar of all ranks
      dsetname = '/Boundary_data/numBoundsRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numBoundsRankPar

      !write(*,*) 'ms_offset ',ms_offset(1),' ms_dims(1) ',ms_dims(1),' ds_dims ',ds_dims(1)

      dsetname = '/Boundary_data/bouCodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bouCodesPar)

      ds_dims(1) = ds_dims(1)*npbou!totalNumBoundsSrl*npbou
      ms_dims(1) = ms_dims(1)*npbou!numBoundsRankPar*npbou
      ms_offset(1) = ms_offset(1)*npbou

      deallocate(aux_array)
      allocate(aux_array(numBoundsRankPar*npbou))

      i=1
      do iBound=1,numBoundsRankPar
         do m=1,npbou
            aux_array(i)=boundPar(iBound,m)
            i=i+1
         end do
      end do

      dsetname = '/Boundary_data/boundPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set ndofRankPar of all ranks
      dsetname = '/Boundary_data/ndofRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=ndofRankPar

      dsetname = '/Boundary_data/ldofPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ldofPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = mpi_size
      ms_offset(1) = 0
      !read data set numBoundaryNodesRankPar of all ranks
      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+aux_array(i)
      end do
      ms_dims(1)=numBoundaryNodesRankPar

      dsetname = '/Boundary_data/lbnodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,lbnodesPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
   end subroutine save_boundary_data_hdf5

   subroutine load_boundary_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank,h5err
      integer :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)
      logical :: isBoundaryFolder

      groupname = trim('/Boundary_data')
      !call create_group_hdf5(file_id,groupname)

      call h5lexists_f(file_id,groupname,isBoundaryFolder,h5err)

      isMeshBoundaries = isBoundaryFolder

      if(mpi_rank.eq.0) write(*,*) 'Loading Boundary data hdf5. -> isBoundaryFolder:',isBoundaryFolder

      if(isMeshBoundaries) then
         dtype = H5T_NATIVE_INTEGER
         ds_rank = 1
         ds_dims(1) = mpi_size
         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = mpi_rank
         allocate(aux_array(1))

         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         numBoundsRankPar=aux_array(1)

         dsetname = '/Boundary_data/numBoundCodes'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         numBoundCodes=aux_array(1)

         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         ndofRankPar=aux_array(1)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         numBoundaryNodesRankPar=aux_array(1)

         !--------------------------------------------------------------------------------------------------------
         allocate(boundPar(numBoundsRankPar,npbou))
         allocate(bouCodesPar(numBoundsRankPar))
         allocate(ldofPar(ndofRankPar))
         allocate(lbnodesPar(numBoundaryNodesRankPar))
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
         allocate( aux_array(mpi_size) )
         ms_dims(1) = mpi_size
         ms_offset(1) = 0
         !read data set numBoundsRankPar of all ranks
         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

         !taking advantadge and setting totalNumBoundsSrl
         totalNumBoundsSrl=0
         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
            totalNumBoundsSrl=totalNumBoundsSrl+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+aux_array(i)
         end do
         ms_dims(1)=numBoundsRankPar

         dsetname = '/Boundary_data/bouCodesPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bouCodesPar)
         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
         ds_dims(1) = ds_dims(1)*npbou!totalNumBoundsSrl*npbou
         ms_dims(1) = ms_dims(1)*npbou!numBoundsRankPar*npbou
         ms_offset(1) = ms_offset(1)*npbou

         allocate(aux_array(numBoundsRankPar*nnode))

         dsetname = '/Boundary_data/boundPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         i=1
         do iBound=1,numBoundsRankPar
            do m=1,npbou
               boundPar(iBound,m)=aux_array(i)
               i=i+1
            end do
         end do

         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
         allocate( aux_array(mpi_size) )
         ms_dims(1) = mpi_size
         ms_offset(1) = 0
         !read data set ndofRankPar of all ranks
         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+aux_array(i)
         end do
         ms_dims(1)=ndofRankPar

         dsetname = '/Boundary_data/ldofPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ldofPar)
         !--------------------------------------------------------------------------------------------------------
         ms_dims(1) = mpi_size
         ms_offset(1) = 0
         !read data set numBoundaryNodesRankPar of all ranks
         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+aux_array(i)
         end do
         ms_dims(1)=numBoundaryNodesRankPar

         dsetname = '/Boundary_data/lbnodesPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ldofPar)
         !--------------------------------------------------------------------------------------------------------

         deallocate(aux_array)
      end if
   end subroutine load_boundary_data_hdf5

!--------------------------------------------------------------------------------------------------------------------------------

   subroutine load_coordinates_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank
      integer(HSSIZE_T), dimension(1) :: ms_offset

      !allocate(coord_x(numNodesRankPar))
      !allocate(coord_y(numNodesRankPar))
      !allocate(coord_z(numNodesRankPar))
      !write(*,*) 'allocating coordPar ->',numNodesRankPar,ndime
      allocate(coordPar(numNodesRankPar,ndime)) !only works for ndime=3

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar

      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      dtype = H5T_NATIVE_REAL

      dsetname = '/Coords/X'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,1))

      dsetname = '/Coords/Y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,2))

      dsetname = '/Coords/Z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,3))

      !write(*,*) 'Loading Coordinates hdf5.... -> (10)',coordPar(10,1),coordPar(10,2),coordPar(10,3)

   end subroutine load_coordinates_hdf5

   subroutine load_globalIds_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,ms_rank
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: iNodeL,iNodeGSrl,max_iNodeGSrl_l,max_iNodeGSrl_g

      allocate(globalIdSrl(numNodesRankPar))
      allocate(globalIdPar(numNodesRankPar))
      allocate(elemGid(numElemsRankPar))

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar

      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      dtype = H5T_NATIVE_INTEGER

      dsetname = '/globalIds/globalIdSrl'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrl)

      dsetname = '/globalIds/globalIdPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdPar)

      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      dsetname = '/globalIds/elemGid'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,elemGid)

      !--------------------------------------------------------------------------------------------------
      !setting totalNumNodesSrl
      max_iNodeGSrl_l=0
      do iNodeL=1,numNodesRankPar
         iNodeGSrl = globalIdSrl(iNodeL)
         max_iNodeGSrl_l = max(iNodeGSrl,max_iNodeGSrl_l)
      end do

      call MPI_Allreduce(max_iNodeGSrl_l,max_iNodeGSrl_g,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_err)

      totalNumNodesSrl = max_iNodeGSrl_g
      !write(*,*) 'setting totalNumNodesSrl ',totalNumNodesSrl
      !--------------------------------------------------------------------------------------------------
   end subroutine load_globalIds_hdf5

!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!           RESULTS FILE
!---------------------------------------------------------------------------------------------------------

   subroutine set_hdf5_resultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_resultsFile_h5_name))//trim(aux_step)//'.h5'
   end subroutine set_hdf5_resultsFile_name

   subroutine save_hdf5_resultsFile(iStep,time,rho,u,pr,E,eta,csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs)
      implicit none
      integer,  intent(in) :: iStep
      real(rp), intent(in) :: time
      real(rp), intent(inout), dimension(numNodesRankPar)       :: rho,pr,E,eta,csound,machno,mu_fluid,divU,Qcrit
      real(rp), intent(inout), dimension(numNodesRankPar,ndime) :: u
      real(rp), intent(inout), dimension(numElemsRankPar,ngaus)  :: mu_e,mu_sgs
      real(rp), intent(inout), dimension(numNodesRankPar,ndime) :: gradRho,curlU
      real(rp)               , dimension(numNodesRankPar)       :: envit,mut
      
      integer(hid_t) :: file_id,plist_id,dset_id,dspace_id,mspace_id,group_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(512) :: full_fileName,groupname,dsetname
      integer :: iPer,i,j
      real(rp) :: aux_array(1)

      !------------------------------------------------------------------------------------
      !-------------------------------------------
      ! 1. Adjust mu_sgs and envit to be nodal
      !$acc kernels
      do i = 1,numElemsRankPar
         do j = 1, nnode
            envit(connecParOrig(i,j)) =  mu_e(i,j)    !@REVISAR SI AQUEST ES EL CONNEC QUE NECESSITO AQUI....
            mut(connecParOrig(i,j))   =  mu_sgs(i,j)
         end do
      end do
      !$acc end kernels

      !-------------------------------------------
      !2. If case is periodic, adjust slave nodes

      if (isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            u(masSlaRankPar(iPer,2),1)       = u(masSlaRankPar(iPer,1),1)
            u(masSlaRankPar(iPer,2),2)       = u(masSlaRankPar(iPer,1),2)
            u(masSlaRankPar(iPer,2),3)       = u(masSlaRankPar(iPer,1),3)
            gradRho(masSlaRankPar(iPer,2),1) = gradRho(masSlaRankPar(iPer,1),1)
            gradRho(masSlaRankPar(iPer,2),2) = gradRho(masSlaRankPar(iPer,1),2)
            gradRho(masSlaRankPar(iPer,2),3) = gradRho(masSlaRankPar(iPer,1),3)
            curlU(masSlaRankPar(iPer,2),1)   = curlU(masSlaRankPar(iPer,1),1)
            curlU(masSlaRankPar(iPer,2),2)   = curlU(masSlaRankPar(iPer,1),2)
            curlU(masSlaRankPar(iPer,2),3)   = curlU(masSlaRankPar(iPer,1),3)
            rho(masSlaRankPar(iPer,2))       = rho(masSlaRankPar(iPer,1))
            pr(masSlaRankPar(iPer,2))        = pr(masSlaRankPar(iPer,1))
            E(masSlaRankPar(iPer,2))         = E(masSlaRankPar(iPer,1))
            eta(masSlaRankPar(iPer,2))       = eta(masSlaRankPar(iPer,1))
            csound(masSlaRankPar(iPer,2))    = csound(masSlaRankPar(iPer,1))
            machno(masSlaRankPar(iPer,2))    = machno(masSlaRankPar(iPer,1))
            mu_fluid(masSlaRankPar(iPer,2))  = mu_fluid(masSlaRankPar(iPer,1))
            mut(masSlaRankPar(iPer,2))       = mut(masSlaRankPar(iPer,1))
            envit(masSlaRankPar(iPer,2))     = envit(masSlaRankPar(iPer,1))
            divU(masSlaRankPar(iPer,2))      = divU(masSlaRankPar(iPer,1))
            Qcrit(masSlaRankPar(iPer,2))     = Qcrit(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if

      !------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_resultsFile_name(iStep,full_fileName)
      !full_fileName = 'test.h5'

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create results file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1
      dtype = H5T_NATIVE_REAL

      dsetname = 'rho'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,3))

      dsetname = 'mu_fluid'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,mu_fluid)

      dsetname = 'mut'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,mut)

      dsetname = 'envit'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,envit)

      dsetname = 'gradRho_x'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,gradRho(:,1))

      dsetname = 'gradRho_y'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,gradRho(:,2))

      dsetname = 'gradRho_z'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,gradRho(:,3))

      dsetname = 'curlU_x'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,curlU(:,1))

      dsetname = 'curlU_y'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,curlU(:,2))

      dsetname = 'curlU_z'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,curlU(:,3))

      dsetname = 'pr'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,E)

      dsetname = 'eta'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,eta)

      dsetname = 'csound'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,csound)

      dsetname = 'machno'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,machno)

      dsetname = 'divU'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,divU)

      dsetname = 'Qcrit'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,Qcrit)

      ! ----  time  -----
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif
      aux_array(1) = time

      dsetname = 'time'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine save_hdf5_resultsFile

   subroutine load_hdf5_resultsFile(load_step,time,rho,u,pr,E,mu_e,mu_sgs)
      implicit none
      integer, intent(in) :: load_step
      real(rp), intent(inout) :: time
      real(rp), intent(inout), dimension(numNodesRankPar)       :: rho,pr,E
      real(rp), intent(inout), dimension(numNodesRankPar,ndime) :: u
      real(rp), intent(inout), dimension(numElemsRankPar,ngaus)  :: mu_e,mu_sgs
      
      character(512) :: full_loadFileName
      integer(hid_t) :: file_id,plist_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(128) :: dsetname
      real(rp) :: aux_array(1)
      real(rp), dimension(numNodesRankPar) :: envit, mut
      integer :: i,j
      
      call set_hdf5_resultsFile_name(load_step,full_loadFileName)
      if(mpi_rank.eq.0) write(*,*) '# Loading results file: ',trim(adjustl(full_loadFileName))

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_loadFileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_loadFileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ms_rank = 1
      ds_rank = 1
      dtype = H5T_NATIVE_REAL

      ! ----  read time  -----
      ds_dims = 1
      ms_offset(1) = 0
      ms_dims(1) = 1

      dsetname = 'time'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      time = aux_array(1)

      if(mpi_rank.eq.0) write(*,*) ' - load_step',load_step,'time',time

      ! ----  read arrays  -----

      ds_dims(1) = totalNumNodesPar
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      dsetname = 'rho'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,3))

      dsetname = 'envit'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,envit)

      dsetname = 'mut'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,mut)

      dsetname = 'pr'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,E)

      !$acc kernels
      do i = 1,numElemsRankPar
         do j = 1, nnode
            mu_e(i,j)   = envit(connecParOrig(i,j))
            mu_sgs(i,j) = mut(connecParOrig(i,j))
         end do
      end do
      !$acc end kernels

      call h5fclose_f(file_id,h5err)

   end subroutine load_hdf5_resultsFile

   subroutine load_hdf5_resultsFile_allArrays(load_step,time,rho,u,pr,E,eta,csound,machno,gradRho,&
                                              curlU,divU,Qcrit,mu_fluid,envit,mut)
      implicit none
      integer, intent(in) :: load_step
      real(rp), intent(inout) :: time
      real(rp), intent(inout), dimension(numNodesRankPar)       :: rho,pr,E,eta,csound,machno,mu_fluid,divU,Qcrit
      real(rp), intent(inout), dimension(numNodesRankPar,ndime) :: u
      real(rp), intent(inout), dimension(numNodesRankPar,ndime) :: gradRho,curlU
      real(rp), intent(inout), dimension(numNodesRankPar)       :: envit,mut
      !real(rp), intent(inout), dimension(numElemsRankPar,ngaus)  :: mu_e,mu_sgs

      character(512) :: full_loadFileName
      integer(hid_t) :: file_id,plist_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(128) :: dsetname
      real(rp) :: aux_array(1)
      integer :: i,j
      
      call set_hdf5_resultsFile_name(load_step,full_loadFileName)
      if(mpi_rank.eq.0) write(*,*) '# Loading results file: ',full_loadFileName

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_loadFileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_loadFileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ms_rank = 1
      ds_rank = 1
      dtype = H5T_NATIVE_REAL

      ! ----  read time  -----
      ds_dims = 1
      ms_offset(1) = 0
      ms_dims(1) = 1

      dsetname = 'time'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      time = aux_array(1)

      if(mpi_rank.eq.0) write(*,*) ' - load_step',load_step,'time',time

      ! ----  read arrays  -----

      ds_dims(1) = totalNumNodesPar
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      dsetname = 'rho'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,u(:,3))

      dsetname = 'envit'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,envit)

      dsetname = 'mut'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,mut)

      dsetname = 'pr'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,E)

      dsetname = 'eta'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,eta)

      dsetname = 'csound'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,csound)

      dsetname = 'machno'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,machno)

      dsetname = 'gradRho_x'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,gradRho(:,1))

      dsetname = 'gradRho_y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,gradRho(:,2))

      dsetname = 'gradRho_z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,gradRho(:,3))

      dsetname = 'curlU_x'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,curlU(:,1))

      dsetname = 'curlU_y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,curlU(:,2))

      dsetname = 'curlU_z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,curlU(:,3))

      dsetname = 'mu_fluid'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,mu_fluid)

      dsetname = 'divU'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,divU)

      dsetname = 'Qcrit'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,Qcrit)

      call h5fclose_f(file_id,h5err)

   end subroutine load_hdf5_resultsFile_allArrays

   subroutine set_hdf5_avgResultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//trim(aux_step)//'.h5'
   end subroutine set_hdf5_avgResultsFile_name

   subroutine save_hdf5_avgResultsFile(iStep,avvel,avve2,avrho,avpre,avmueff)
      implicit none
      integer, intent(in) :: iStep
      real(rp), intent(in), dimension(numNodesRankPar)       :: avrho, avpre, avmueff
      real(rp), intent(in), dimension(numNodesRankPar,ndime) :: avvel, avve2

      integer(hid_t) :: file_id,plist_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(512) :: full_fileName,groupname,dsetname
      real(rp) :: aux_array(1)

      !------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_avgResultsFile_name(iStep,full_fileName)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create avg results file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1
      dtype = H5T_NATIVE_REAL

      dsetname = 'avrho'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avrho)

      dsetname = 'avpre'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avpre)

      dsetname = 'avmueff'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avmueff)

      dsetname = 'avvel_x'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avvel(:,1))

      dsetname = 'avvel_y'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avvel(:,2))

      dsetname = 'avvel_z'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avvel(:,3))

      dsetname = 'avve2_x'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avve2(:,1))

      dsetname = 'avve2_y'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avve2(:,2))

      dsetname = 'avve2_z'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avve2(:,3))

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine save_hdf5_avgResultsFile

   subroutine load_hdf5_avgResultsFile(iStep,avvel,avve2,avrho,avpre,avmueff)
      implicit none
      integer, intent(in) :: iStep
      real(rp), intent(inout), dimension(numNodesRankPar)       :: avrho, avpre, avmueff
      real(rp), intent(inout), dimension(numNodesRankPar,ndime) :: avvel, avve2

      integer(hid_t) :: file_id,plist_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(512) :: full_fileName,groupname,dsetname
      real(rp) :: aux_array(1)

      !------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_avgResultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Loading results file: ',trim(adjustl(full_fileName))

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results avg file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ds_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_rank = 1
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1
      dtype = H5T_NATIVE_REAL

      dsetname = 'avrho'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avrho)

      dsetname = 'avpre'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avpre)

      dsetname = 'avmueff'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avmueff)

      dsetname = 'avvel_x'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avvel(:,1))

      dsetname = 'avvel_y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avvel(:,2))

      dsetname = 'avvel_z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avvel(:,3))

      dsetname = 'avve2_x'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avve2(:,1))

      dsetname = 'avve2_y'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avve2(:,2))

      dsetname = 'avve2_z'
      call read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,avve2(:,3))

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine load_hdf5_avgResultsFile

   subroutine close_hdf5_resultsfile()
      implicit none
      integer :: h5err

      !close the file.
      call h5fclose_f(resultsFile_h5_id,h5err)

   end subroutine close_hdf5_resultsfile

!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------
!        VTKHDF5 FILE
!---------------------------------------------------------------------------------------------------------

   subroutine set_vtkhdf_resultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_resultsFile_h5_name))//trim(aux_step)//'-vtk.hdf'
   end subroutine set_vtkhdf_resultsFile_name

   subroutine set_vtkhdf_avgResultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//trim(aux_step)//'-vtk.hdf'
   end subroutine set_vtkhdf_avgResultsFile_name

   subroutine set_vtkhdf_finalAvgResultsFile_name(full_fileName)
      implicit none
      character(len=*), intent(out) :: full_fileName

      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//'final-vtk.hdf'
   end subroutine set_vtkhdf_finalAvgResultsFile_name

   subroutine create_vtkhdf_resultsFile(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(hid_t) :: plist_id,dset_id,dspace_id,mspace_id,group_id,aspace_id,attr_id
      integer(hid_t) :: dtype,atype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims,a_dims
      integer(HSIZE_T), dimension(2) :: ds_dims2d,ms_dims2d
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(HSSIZE_T), dimension(2) :: ms_offset2d
      integer(SIZE_T) :: attr_length
      integer :: ds_rank,ms_rank,a_rank,h5err
      character(512) :: full_fileName,groupname,dsetname
      character(16) :: attr_value
      integer :: ii,iNodeL,iElemL
      integer(1),allocatable :: aux_array_i8(:)
      integer(4),allocatable :: aux_array_i32(:)
      integer(8),allocatable :: aux_array_i64(:)
      real(8),allocatable :: aux_array_r64(:)

      !------------------------------------------------------------------------------------
      !call set_vtkhdf_resultsFile_name(iStep,full_fileName)

      !! Setup file access property list with parallel I/O access.
      !call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      !call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      !! create file collectively
      !call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      !if(h5err .ne. 0) then
      !   write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
      !   call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      !end if
      !call h5pclose_f(plist_id, h5err)

!--------------------------------------------------------------------------------
      groupname = '/VTKHDF'
      !call create_group_hdf5(file_id,groupname)
      call h5gcreate_f(file_id,groupname,group_id,h5err)
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      ! Attribute 'Type'
      call h5screate_f(H5S_SCALAR_F,aspace_id,h5err)
      
      attr_value = "UnstructuredGrid"
      attr_length = len_trim(attr_value)

      a_dims(1) = 1
      call h5tcopy_f(H5T_C_S1,atype,h5err)
      call h5tset_size_f(atype, attr_length,h5err)
      call h5tset_strpad_f(atype,H5T_STR_NULLPAD_F,h5err)

      call h5acreate_f(group_id,'Type',atype,aspace_id,attr_id,h5err)

      call h5awrite_f(attr_id,atype,attr_value,a_dims, h5err)
      call h5tclose_f(atype,h5err)

      call h5aclose_f(attr_id, h5err)
      
      call h5sclose_f(aspace_id, h5err)
!--------------------------------------------------------------------------------
      ! Attribute 'Version'
      a_rank = 1
      a_dims(1) = 2
      call h5screate_simple_f(a_rank,a_dims,aspace_id,h5err)

      call h5acreate_f(group_id,'Version',H5T_NATIVE_INTEGER,aspace_id,attr_id,h5err)

      allocate(aux_array_i32(2))
      aux_array_i32(1) = 1
      aux_array_i32(2) = 0
      call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,aux_array_i32,a_dims,h5err)
      deallocate(aux_array_i32)   

      call h5aclose_f(attr_id,h5err)

      call h5sclose_f(aspace_id,h5err)

      call h5gclose_f(group_id, h5err)
!--------------------------------------------------------------------------------
      groupname = '/VTKHDF/CellData'
      call create_group_hdf5(file_id,groupname)

      groupname = '/VTKHDF/FieldData'
      call create_group_hdf5(file_id,groupname)

      groupname = '/VTKHDF/PointData'
      call create_group_hdf5(file_id,groupname)
      !--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------
      dtype = H5T_NATIVE_DOUBLE
      ds_rank = 2
      ms_rank = 2
      ds_dims2d(1) = 3
      ms_dims2d(1) = 1
      ds_dims2d(2) = totalNumNodesPar
      ms_dims2d(2) = numNodesRankPar
      ms_offset2d(2) = rankNodeStart-1
!--------------------------------------------------------------------------------
      allocate(aux_array_r64(numNodesRankPar))

      dsetname = '/VTKHDF/Points'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      
!-------------------------------------------------------------------------------- 
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(coordPar(iNodeL,1),8)
      end do
      ms_offset2d(1) = 0
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
!-------------------------------------------------------------------------------- 
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(coordPar(iNodeL,2),8)
      end do
      ms_offset2d(1) = 1
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
!-------------------------------------------------------------------------------- 
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(coordPar(iNodeL,3),8)
      end do
      ms_offset2d(1) = 2
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)

      deallocate(aux_array_r64)
      !-----------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = mpi_rank
      dtype = H5T_STD_I64LE
      allocate(aux_array_i64(1))

      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i64(1) = numNodesRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i64)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i64(1) = numElemsRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i64)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i64(1) = numElemsRankPar*nnode
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i64)

      !-----------------------------------------------------------------------------
      deallocate(aux_array_i64)
      allocate(aux_array_i64(numElemsRankPar+1))

      ds_dims(1) = totalNumElements+mpi_size
      ms_dims(1) = numElemsRankPar+1
      ms_offset(1) = rankElemStart-1+mpi_rank

      dsetname = '/VTKHDF/Offsets'
      aux_array_i64(1) = 0
      do iElemL = 2,(numElemsRankPar+1)
         aux_array_i64(iElemL) = aux_array_i64(iElemL-1)+nnode
      end do

      !call create_chunked_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,ms_dims,dtype)
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i64)
      deallocate(aux_array_i64)
      !-----------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      !dtype = H5T_STD_I64LE
      allocate(aux_array_i64(numElemsRankPar*nnode))

      ds_dims(1) = totalNumElements*nnode
      ms_dims(1) = numElemsRankPar*nnode
      ms_offset(1) = (rankElemStart-1)*nnode

      do ii = 1,numElemsRankPar*nnode
         aux_array_i64(ii) = connecVTK(ii)-1
      end do

      dsetname = '/VTKHDF/Connectivity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i64)
      deallocate(aux_array_i64)
      !-----------------------------------------------------------------------------
      allocate(aux_array_i8(numElemsRankPar))
      dtype = H5T_STD_U8LE
      dsetname = '/VTKHDF/Types'
      !$acc kernels
      aux_array_i8(:) = 72
      !$acc end kernels

      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !call create_chunked_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,ms_dims,dtype)
      call write_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      deallocate(aux_array_i8)

      !close the file.
      !call h5fclose_f(file_id,h5err)

   end subroutine create_vtkhdf_resultsFile
            
   subroutine write_vtkhdf_instResultsFile(full_fileName,rho,pr,E,eta,csound,machno,divU,Qcrit,envit,mut,mu_fluid,u,gradRho,curlU)
      implicit none
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: rho,pr,E,eta,csound,machno,divU,Qcrit
      real(rp),dimension(numNodesRankPar),intent(in) :: envit,mut,mu_fluid
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: u,gradRho,curlU

      character(512) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims,a_dims
      integer(HSIZE_T), dimension(2) :: ds_dims2d,ms_dims2d
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(HSSIZE_T), dimension(2) :: ms_offset2d
      integer :: ds_rank,ms_rank,h5err

      integer :: ii,iNodeL,iElemL
      integer(1),allocatable :: aux_array_i8(:)
      integer(4),allocatable :: aux_array_i32(:)
      integer(8),allocatable :: aux_array_i64(:)
      real(8),allocatable :: aux_array_r64(:)

      !------------------------------------------------------------------------------------

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call create_vtkhdf_resultsFile(file_id)

      allocate(aux_array_r64(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = H5T_NATIVE_DOUBLE
      ds_rank = 1
      ms_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      ! ## RHO ##
      dsetname = '/VTKHDF/PointData/rho'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(rho(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## PR ##
      dsetname = '/VTKHDF/PointData/pr'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(pr(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## Ener ##
      dsetname = '/VTKHDF/PointData/Ener'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(E(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)
      
      ! ## eta ##
      dsetname = '/VTKHDF/PointData/eta'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(eta(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## csound ##
      dsetname = '/VTKHDF/PointData/csound'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(csound(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## machno ##
      dsetname = '/VTKHDF/PointData/machno'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(machno(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## divU ##
      dsetname = '/VTKHDF/PointData/divU'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(divU(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## Qcrit ##
      dsetname = '/VTKHDF/PointData/Qcrit'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(Qcrit(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## envit ##
      dsetname = '/VTKHDF/PointData/envit'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(envit(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## mut ##
      dsetname = '/VTKHDF/PointData/mut'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(mut(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## mu_fluid ##
      dsetname = '/VTKHDF/PointData/mu_fluid'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(mu_fluid(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      !-------------------------------------------------------------------------------------------------------
      ds_rank = 2
      ms_rank = 2
      ds_dims2d(1) = 3
      ms_dims2d(1) = 1
      ds_dims2d(2) = totalNumNodesPar
      ms_dims2d(2) = numNodesRankPar
      ms_offset2d(2) = rankNodeStart-1
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/Velocity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 0
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(u(iNodeL,1),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,u(:,1))
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 1
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(u(iNodeL,2),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 2
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(u(iNodeL,3),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/gradRho'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 0
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(gradRho(iNodeL,1),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 1
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(gradRho(iNodeL,2),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 2
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(gradRho(iNodeL,3),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/curlU'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 0
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(curlU(iNodeL,1),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 1
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(curlU(iNodeL,2),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 2
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(curlU(iNodeL,3),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 

      deallocate(aux_array_r64)

      !-------------------------------------------------------------------------------------------------------

      allocate(aux_array_i8(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = H5T_STD_U8LE
      ds_rank = 1
      ms_rank = 1
      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i8(iElemL) = mpi_rank
      end do
      call write_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      deallocate(aux_array_i8)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_instResultsFile

   subroutine write_vtkhdf_avgResultsFile(full_fileName,avrho,avpre,avmueff,avvel,avve2)
      implicit none
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: avrho,avpre,avmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: avvel,avve2

      character(512) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSIZE_T), dimension(2) :: ds_dims2d,ms_dims2d
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(HSSIZE_T), dimension(2) :: ms_offset2d
      integer :: ds_rank,ms_rank,h5err

      integer :: ii,iNodeL,iElemL
      integer(1),allocatable :: aux_array_i8(:)
      integer(4),allocatable :: aux_array_i32(:)
      integer(8),allocatable :: aux_array_i64(:)
      real(8),allocatable :: aux_array_r64(:)

      !------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call create_vtkhdf_resultsFile(file_id)

      allocate(aux_array_r64(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = H5T_NATIVE_DOUBLE
      ds_rank = 1
      ms_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      ! ## avgRHO ##
      dsetname = '/VTKHDF/PointData/avrho'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avrho(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## avgPRE ##
      dsetname = '/VTKHDF/PointData/avpre'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avpre(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      ! ## avgMUEFF ##
      dsetname = '/VTKHDF/PointData/avmueff' 
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avmueff(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)
      
      !-------------------------------------------------------------------------------------------------------
      ds_rank = 2
      ms_rank = 2
      ds_dims2d(1) = 3
      ms_dims2d(1) = 1
      ds_dims2d(2) = totalNumNodesPar
      ms_dims2d(2) = numNodesRankPar
      ms_offset2d(2) = rankNodeStart-1
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/avvel'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 0
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avvel(iNodeL,1),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 1
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avvel(iNodeL,2),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 2
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avvel(iNodeL,3),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/avve2'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 0
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avve2(iNodeL,1),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 1
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avve2(iNodeL,2),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      ms_offset2d(1) = 2
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(avve2(iNodeL,3),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r64)
      !-------------------------------------------------------------------------------- 
      deallocate(aux_array_r64)

      !-------------------------------------------------------------------------------------------------------
      allocate(aux_array_i8(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = H5T_STD_U8LE
      ds_rank = 1
      ms_rank = 1
      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i8(iElemL) = mpi_rank
      end do
      call write_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      deallocate(aux_array_i8)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_avgResultsFile

   subroutine write_vtkhdf_floatFieldFile(full_fileName,floatField)
      implicit none
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: floatField

      character(512) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err

      integer :: ii,iNodeL,iElemL
      integer(1),allocatable :: aux_array_i8(:)
      real(8),allocatable :: aux_array_r64(:)

      !------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call create_vtkhdf_resultsFile(file_id)

      allocate(aux_array_r64(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = H5T_NATIVE_DOUBLE
      ds_rank = 1
      ms_rank = 1
      ds_dims(1) = totalNumNodesPar
      ms_dims(1) = numNodesRankPar
      ms_offset(1) = rankNodeStart-1

      ! ## floatField ##
      dsetname = '/VTKHDF/PointData/floatField'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iNodeL = 1,numNodesRankPar
         aux_array_r64(iNodeL) = real(floatField(iNodeL),8)
      end do
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_r64)

      !-------------------------------------------------------------------------------------------------------
      allocate(aux_array_i8(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = H5T_STD_U8LE
      ds_rank = 1
      ms_rank = 1
      ds_dims(1) = totalNumElements
      ms_dims(1) = numElemsRankPar
      ms_offset(1) = rankElemStart-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i8(iElemL) = mpi_rank
      end do
      call write_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      deallocate(aux_array_i8)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_floatFieldFile

!-----------------------------------------------------------------------------------------------------------------------------

   subroutine save_vtkhdf_instResultsFile(iStep,rho,pr,E,eta,csound,machno,divU,Qcrit,envit,mut,mu_fluid,u,gradRho,curlU)
      implicit none
      integer, intent(in) :: iStep
      real(rp),dimension(numNodesRankPar),intent(in) :: rho,pr,E,eta,csound,machno,divU,Qcrit
      real(rp),dimension(numNodesRankPar),intent(in) :: envit,mut,mu_fluid
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: u,gradRho,curlU
      character(512) :: full_fileName

      !------------------------------------------------------------------------------------
      call set_vtkhdf_resultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF inst results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_instResultsFile(full_fileName,rho,pr,E,eta,csound,machno,divU,Qcrit,envit,mut,mu_fluid,u,gradRho,curlU)

   end subroutine save_vtkhdf_instResultsFile

   subroutine save_vtkhdf_avgResultsFile(iStep,avrho,avpre,avmueff,avvel,avve2)
      implicit none
      integer, intent(in) :: iStep
      real(rp),dimension(numNodesRankPar),intent(in) :: avrho,avpre,avmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: avvel,avve2
      character(512) :: full_fileName

      call set_vtkhdf_avgResultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF avg results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_avgResultsFile(full_fileName,avrho,avpre,avmueff,avvel,avve2)

   end subroutine save_vtkhdf_avgResultsFile

   subroutine save_vtkhdf_finalAvgResultsFile(favrho,favpre,favmueff,favvel,favve2)
      implicit none
      real(rp),dimension(numNodesRankPar),intent(in) :: favrho,favpre,favmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: favvel,favve2
      character(512) :: full_fileName

      call set_vtkhdf_finalAvgResultsFile_name(full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF final avg results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_avgResultsFile(full_fileName,favrho,favpre,favmueff,favvel,favve2)

   end subroutine save_vtkhdf_finalAvgResultsFile

   subroutine save_vtkhdf_flotFieldFile(floatField)
      implicit none
      !integer, intent(in) :: iStep
      real(rp),dimension(numNodesRankPar),intent(in) :: floatField
      character(512) :: full_fileName

      !------------------------------------------------------------------------------------
      call set_vtkhdf_resultsFile_name(0,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF floatField file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_floatFieldFile(full_fileName,floatField)

   end subroutine save_vtkhdf_flotFieldFile

!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------


end module mod_hdf5
