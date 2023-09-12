module mod_hdf5
   use hdf5
   use mod_constants
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms
   use mod_comms_boundaries
   use mod_custom_types
   implicit none

   character(256) :: meshFile_h5_name,surface_meshFile_h5_name
   character(256) :: base_resultsFile_h5_name,base_avgResultsFile_h5_name,base_restartFile_h5_name

   integer(hid_t) :: h5_datatype_uint1,h5_datatype_int1,h5_datatype_int4,h5_datatype_int8
   integer(hid_t) :: h5_datatype_real4,h5_datatype_real8

   real(rp), allocatable :: auxInterpNodeScalarField(:),auxInterpNodeVectorField(:,:)

contains

   subroutine init_hdf5_interface()
      implicit none
      integer :: h5err

      !.init h5 interface
      call h5open_f(h5err)

      h5_datatype_uint1 = H5T_STD_U8LE
      h5_datatype_int1 = H5T_STD_I8LE
      h5_datatype_int4 = H5T_NATIVE_INTEGER
      h5_datatype_int8 = H5T_STD_I64LE

      h5_datatype_real4 = H5T_NATIVE_REAL
      h5_datatype_real8 = H5T_NATIVE_DOUBLE

   end subroutine init_hdf5_interface

   subroutine end_hdf5_interface()
      implicit none
      integer :: h5err

      !close h5 interface
      call h5close_f(h5err)
   end subroutine end_hdf5_interface

   subroutine init_hdf5_auxiliar_saving_arrays()
      implicit none

      if(mesh_isLoaded .eqv. .false.) then
         write(*,*) 'FATAL ERROR! Mesh not loaded when calling init_hdf5_auxiliar_saving_arrays()! CRASHING!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      allocate(auxInterpNodeScalarField(numNodesRankPar))
      allocate(auxInterpNodeVectorField(numNodesRankPar,ndime))
      !$acc enter data create(auxInterpNodeScalarField(:))
      !$acc enter data create(auxInterpNodeVectorField(:,:))

   end subroutine

   subroutine end_hdf5_auxiliar_saving_arrays()
      implicit none

      !$acc exit data delete(auxInterpNodeScalarField(:))
      !$acc exit data delete(auxInterpNodeVectorField(:,:))
      deallocate(auxInterpNodeScalarField)
      deallocate(auxInterpNodeVectorField)

   end subroutine


   subroutine set_hdf5_meshFile_name(file_path,file_name,numRanks)
      implicit none
      character(len=*), intent(in) :: file_path,file_name
      integer,intent(in) :: numRanks
      character(len=12) :: aux_numRanks

      write(aux_numRanks,'(I0)') numRanks
      meshFile_h5_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'-'//trim(aux_numRanks)//'.hdf'
   end subroutine set_hdf5_meshFile_name

   subroutine set_hdf5_surface_meshFile_name()
      implicit none

      surface_meshFile_h5_name = 'surface_'//trim(adjustl(meshFile_h5_name))
   end subroutine set_hdf5_surface_meshFile_name

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

      base_restartFile_h5_name = trim(adjustl(res_filePath))//'restart_'&
         //trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_'

   end subroutine set_hdf5_baseResultsFile_name

!---------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------

   subroutine create_hdf5_file(full_hdf5FileName,hdf5_file_id)
      implicit none
      character(*), intent(in) :: full_hdf5FileName
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer(hid_t) :: plist_id
      integer :: h5err

      !-----------------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_hdf5FileName,H5F_ACC_TRUNC_F,hdf5_file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create results file ',trim(adjustl(full_hdf5FileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

   end subroutine create_hdf5_file

   subroutine open_hdf5_file(full_hdf5FileName,hdf5_file_id)
      implicit none
      character(*), intent(in) :: full_hdf5FileName
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer(hid_t) :: plist_id
      integer :: h5err

      !-----------------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_hdf5FileName,H5F_ACC_RDWR_F,hdf5_file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load hdf5 file ',trim(adjustl(full_hdf5FileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)


      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

   end subroutine open_hdf5_file

   subroutine close_hdf5_file(hdf5_file_id)
      implicit none
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer :: h5err

      !close the file.
      call h5fclose_f(hdf5_file_id,h5err)
   end subroutine close_hdf5_file
#if 0
   subroutine create_hdf5_meshFile_from_tool(hdf5_file_id)
      implicit none
      integer(hid_t),intent(inout) :: hdf5_file_id
      integer(hid_t) :: plist_id
      integer :: h5err
      !-------------------------------------------------------------------------------

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(meshFile_h5_name,H5F_ACC_TRUNC_F,hdf5_file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(app_comm,-1,mpi_err)
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
#endif
   subroutine create_hdf5_groups_datasets_in_meshFile_from_tool(mnnode,mnpbou,hdf5_file_id,isPeriodic,isBoundaries,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8,&
               vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,&
               vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer(4),intent(in) :: numMshRanks2Part,numElemsGmsh
      integer(8),intent(in) :: numNodesParTotal_i8
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      integer(4),intent(in),dimension(0:numMshRanks2Part-1) :: vecNumPerNodesMshRank
      integer(hid_t) :: dset_id,dspace_id,group_id
      integer(hid_t) :: dtype
      integer(hsize_t) :: ds_dims(1),ds_dims2d(2)
      integer(4) :: ds_rank,h5err
      character(256) :: groupname,dsetname

      !--------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numNodesParTotal_i8

      dtype = h5_datatype_int8

      groupname = '/globalIds'
      call create_group_hdf5(hdf5_file_id,groupname)
      dsetname = '/globalIds/globalIdSrl'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/globalIds/globalIdPar'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      dtype = h5_datatype_int4

      ds_dims(1) = numElemsGmsh
      dsetname = '/globalIds/elemGid'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      !--------------------------------------------------------------------------------
      call select_dtype_rp(dtype)
      ds_rank = 2
      ds_dims2d(1) = ndime
      ds_dims2d(2) = numNodesParTotal_i8

      groupname = '/Coords'
      call create_group_hdf5(hdf5_file_id,groupname)
      dsetname = '/Coords/Points'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims2d,dtype)
      !--------------------------------------------------------------------------------
      call create_groups_datasets_connectivity_workingNodes_hdf5(mnnode,hdf5_file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)

      if(isBoundaries) then
         call create_groups_datasets_boundary_data_hdf5(mnpbou,hdf5_file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
      end if

      if(numMshRanks2Part.ge.2) then
         call create_groups_datasets_parallel_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumMshRanksWithComms,vecNumNodesToCommMshRank)
         if(isBoundaries) then
            call create_groups_datasets_parallel_data_boundary_hdf5(hdf5_file_id,numMshRanks2Part,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank)
         end if
      end if

      if(isPeriodic) then
         call create_groups_datasets_periodic_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumPerNodesMshRank)
      end if

      !---------------------------------------------------------
      ! order

      dtype = h5_datatype_int4

      groupname = '/order'
      call create_group_hdf5(hdf5_file_id,groupname)

      ds_rank = 1
      ds_dims(1) = 1
      dsetname = '/order/porder'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = mnnode
      dsetname = '/order/a2ijk'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/gmsh2ijk'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/vtk2ijk'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

      ds_dims(1) = mnpbou
      dsetname = '/order/a2ij'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/gmsh2ij'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/order/vtk2ij'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)

   end subroutine create_hdf5_groups_datasets_in_meshFile_from_tool

   subroutine create_groups_datasets_connectivity_workingNodes_hdf5(mnnode,file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: mnnode,numMshRanks2Part,numElemsGmsh
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,mshRank,accumVal

      groupname = trim('/Connectivity')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      !-----------------------------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = int(numElemsGmsh,hsize_t)*int(mnnode,hsize_t)
      !if(mpi_rank.eq.0) write(*,*) 'debug ds_dims',ds_dims(1),'numElemsGmsh',numElemsGmsh,'mnnode',mnnode

      dsetname = '/Connectivity/connecParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Connectivity/connecParWork'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !dsetname = '/Connectivity/connecVTK'
      !call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !-----------------------------------------------------------------------------------------------------
      ds_dims(1) = numMshRanks2Part

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !-----------------------------------------------------------------------------------------------------
      ds_dims(1)=0
      do mshRank=0,numMshRanks2Part-1
         ds_dims(1)=ds_dims(1)+int(vecNumWorkingNodes(mshRank),hsize_t)
      end do
      !if(mpi_rank.eq.0) write(*,*) 'debug workingNodesPar ds_dims',ds_dims(1)

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

      dtype = h5_datatype_int4
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

      dsetname = '/Parallel_data/nodesToComm'
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

      dtype = h5_datatype_int4
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

      dsetname = '/Parallel_data_boundary/nodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !----------------------------------------------------------------------------------------

   end subroutine create_groups_datasets_parallel_data_boundary_hdf5

   subroutine create_groups_datasets_periodic_data_hdf5(file_id,numMshRanks2Part,vecNumPerNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumPerNodesMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,mshRank,accumVal

      groupname = trim('/Periodic_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part

      dsetname = '/Periodic_data/nPerRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      accumVal=0
      do mshRank=0,numMshRanks2Part-1
         accumVal=accumVal+vecNumPerNodesMshRank(mshRank)
      end do
      ds_dims(1) = accumVal

      dsetname = '/Periodic_data/masSlaRankPar1'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/Periodic_data/masSlaRankPar2'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !--------------------------------------------------------------------------------------------------------

   end subroutine create_groups_datasets_periodic_data_hdf5

   subroutine create_groups_datasets_boundary_data_hdf5(mnpbou,file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: mnpbou,numMshRanks2Part
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims
      integer(hid_t) :: dtype
      integer :: ds_rank,h5err
      integer :: mshRank,accumVal

      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
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

      ds_dims(1) = ds_dims(1)*mnpbou

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

   subroutine write_mshRank_data_in_hdf5_meshFile_from_tool(mporder,mnnode,mnpbou,hdf5_file_id,mshRank,numMshRanks2Part,isPeriodic,isBoundaries,numElemsGmsh,numBoundFacesGmsh,&
         numElemsMshRank,mshRankElemStart,mshRankElemEnd,mshRankNodeStart_i8,mshRankNodeEnd_i8,numNodesMshRank,numWorkingNodesMshRank,numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank,maxBoundCode,&
         a2ijk,a2ij,gmsh2ijk,gmsh2ij,vtk2ijk,vtk2ij,&
         elemGidMshRank,globalIdSrlMshRank_i8,globalIdParMshRank_i8,connecVTKMshRank,connecParOrigMshRank,connecParWorkMshRank,coordParMshRank,workingNodesMshRank,&
         boundaryNodesMshRank,dofNodesMshRank,boundFacesCodesMshRank,boundFacesOrigMshRank,boundFacesMshRank,numPerNodesMshRank,masSlaNodesMshRank,&
         numNodesToCommMshRank,numMshRanksWithComms,nodesToCommMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank,&
         bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms,bnd_nodesToCommMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank,&
         vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer(4),intent(in) :: mporder,mnnode,mnpbou,mshRank,numMshRanks2Part,numElemsGmsh,numBoundFacesGmsh
      integer(4),intent(in) :: numElemsMshRank,mshRankElemStart,mshRankElemEnd
      integer(8),intent(in) :: mshRankNodeStart_i8,mshRankNodeEnd_i8
      integer(4),intent(in) :: numNodesMshRank,numWorkingNodesMshRank,numPerNodesMshRank
      integer(4),intent(in) :: numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank
      integer(4),intent(in) :: maxBoundCode

      integer(4),intent(in),dimension(mnnode) :: a2ijk,gmsh2ijk,vtk2ijk
      integer(4),intent(in),dimension(mnpbou) :: a2ij,gmsh2ij,vtk2ij
      integer(4),intent(in) :: elemGidMshRank(numElemsMshRank),workingNodesMshRank(numWorkingNodesMshRank)
      integer(8),intent(in) :: globalIdSrlMshRank_i8(numNodesMshRank),globalIdParMshRank_i8(numNodesMshRank)
      integer(4),intent(in) :: connecVTKMshRank(numElemsMshRank*mnnode)
      integer(4),intent(in) :: connecParOrigMshRank(numElemsMshRank,mnnode),connecParWorkMshRank(numElemsMshRank,mnnode)
      real(rp),intent(in)   :: coordParMshRank(numNodesMshRank,3)

      integer(4),intent(in) :: boundaryNodesMshRank(numBoundaryNodesMshRank),dofNodesMshRank(numDoFMshRank)
      integer(4),intent(in) :: boundFacesCodesMshRank(numBoundFacesMshRank),boundFacesOrigMshRank(numBoundFacesMshRank,mnpbou),boundFacesMshRank(numBoundFacesMshRank,mnpbou)
      integer(4),intent(in) :: numNodesToCommMshRank,numMshRanksWithComms
      integer(4),intent(in) :: masSlaNodesMshRank(numPerNodesMshRank,2)
      integer(4),intent(in),dimension(numMshRanksWithComms) :: nodesToCommMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank
      integer(4),intent(in) :: bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms
      integer(4),intent(in),dimension(bnd_numMshRanksWithComms) :: bnd_nodesToCommMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank

      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank
      !-------------------------------------------------------------------------------------------------------------------------------
      character(128) :: dsetname
      integer(hsize_t) :: ms_dims(1),ms_dims2d(2)
      integer(4) :: i,m,iBound,iElemL
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(4),allocatable :: aux_array(:)
      integer(8),allocatable :: aux_array_i8(:)
      !-------------------------------------------------------------------------------------------------------------------------------

      ms_dims(1) = int(numNodesMshRank,hsize_t)
      ms_offset(1) = int(mshRankNodeStart_i8,hssize_t)-1

      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,globalIdSrlMshRank_i8)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,globalIdParMshRank_i8)

      ms_dims(1) = int(numElemsMshRank,hsize_t)
      ms_offset(1) = int(mshRankElemStart,hssize_t)-1

      dsetname = '/globalIds/elemGid'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,elemGidMshRank)

      !---------------------------------------------------------------------------------------------------------------------
      ms_dims2d(1) = ndime
      ms_dims2d(2) = int(numNodesMshRank,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(mshRankNodeStart_i8,hssize_t)-1

      dsetname = '/Coords/Points'
      call write_dataspace_2d_tr_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims2d,ms_offset2d,coordParMshRank)
      !---------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------
      if(numMshRanks2Part.ge.2) then
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array_i8(1))

         dsetname = '/Parallel_data/rankNodeStart'
         aux_array_i8(1)=mshRankNodeStart_i8
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

         dsetname = '/Parallel_data/rankNodeEnd'
         aux_array_i8(1)=mshRankNodeEnd_i8
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

         deallocate(aux_array_i8)
         allocate(aux_array(1))

         dsetname = '/Parallel_data/rankElemStart'
         aux_array(1)=mshRankElemStart
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         dsetname = '/Parallel_data/rankElemEnd'
         aux_array(1)=mshRankElemEnd
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         dsetname = '/Parallel_data/numRanksWithComms'
         aux_array(1)=numMshRanksWithComms
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         dsetname = '/Parallel_data/numNodesToComm'
         aux_array(1)=numNodesToCommMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         deallocate(aux_array)

         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumMshRanksWithComms(i),hssize_t)
         end do
         ms_dims(1)=int(numMshRanksWithComms,hsize_t)
         !write(*,*) '[',mshRank,']ms_offset',ms_offset(1),'ms_dims(1)',ms_dims(1)!,' ds_dims ',ds_dims(1)

         dsetname = '/Parallel_data/ranksToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,ranksToCommMshRank)

         dsetname = '/Parallel_data/commsMemPosInLoc'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,commsMemPosInLocMshRank)

         dsetname = '/Parallel_data/commsMemSize'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,commsMemSizeMshRank)

         dsetname = '/Parallel_data/commsMemPosInNgb'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,commsMemPosInNgbMshRank)

         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumNodesToCommMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numNodesToCommMshRank,hsize_t)
         !write(*,*) '[',mshRank,']ms_offset',ms_offset(1),'ms_dims(1)',ms_dims(1)!,' ds_dims ',ds_dims(1)

         dsetname = '/Parallel_data/nodesToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,nodesToCommMshRank)

      end if
      !------------------------------------------------------------------------------------------------------------------------

      if(isBoundaries) then

         dsetname = '/Boundary_data/numBoundCodes'
         if(mshRank.eq.0) then
            ms_dims(1) = 1
            ms_offset(1) = int(mshRank,hssize_t)
            allocate(aux_array(1))
            aux_array(1)=maxBoundCode

            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
            deallocate(aux_array)
         else
            ms_dims(1) = 0
            ms_offset(1) = 0
            allocate(aux_array(0))

            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
            deallocate(aux_array)
         end if

         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array(1))

         dsetname = '/Boundary_data/numBoundsRankPar'
         aux_array(1)=numBoundFacesMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         dsetname = '/Boundary_data/ndofRankPar'
         aux_array(1)=numDoFMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         aux_array(1)=numBoundaryNodesMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         deallocate(aux_array)
         !------------------------------------------------------------------------------------------------------------------------
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumDoFMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numDoFMshRank,hsize_t)

         dsetname = '/Boundary_data/ldofPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,dofNodesMshRank)
         !------------------------------------------------------------------------------------------------------------------------
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumBoundaryNodesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundaryNodesMshRank,hsize_t)

         dsetname = '/Boundary_data/lbnodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,boundaryNodesMshRank)
         !------------------------------------------------------------------------------------------------------------------------
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumBoundFacesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundFacesMshRank,hsize_t)

         dsetname = '/Boundary_data/bouCodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,boundFacesCodesMshRank)
         !------------------------------------------------------------------------------------------------------------------------
         allocate(aux_array(numBoundFacesMshRank*mnpbou))

         i=1
         do iBound=1,numBoundFacesMshRank
            do m=1,mnpbou
               aux_array(i)=boundFacesMshRank(iBound,m)
               i=i+1
            end do
         end do

         ms_dims(1) = int(ms_dims(1),hsize_t)*int(mnpbou,hsize_t)
         ms_offset(1) = ms_offset(1)*int(mnpbou,hssize_t)

         dsetname = '/Boundary_data/boundPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         i=1
         do iBound=1,numBoundFacesMshRank
            do m=1,mnpbou
               aux_array(i)=boundFacesOrigMshRank(iBound,m)
               i=i+1
            end do
         end do

         dsetname = '/Boundary_data/boundParOrig'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------
         if(numMshRanks2Part.ge.2) then
            ms_dims(1) = 1
            ms_offset(1) = int(mshRank,hssize_t)
            allocate(aux_array(1))

            dsetname = '/Parallel_data_boundary/numRanksWithComms'
            aux_array(1)=bnd_numMshRanksWithComms
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

            dsetname = '/Parallel_data_boundary/numNodesToComm'
            aux_array(1)=bnd_numNodesToCommMshRank
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

            deallocate(aux_array)

            ms_offset(1)=0
            do i=0,mshRank-1 !from rank 0 mpi_rank-1
               ms_offset(1)=ms_offset(1)+int(vecBndNumMshRanksWithComms(i),hssize_t)
            end do
            ms_dims(1)=int(bnd_numMshRanksWithComms,hsize_t)

            dsetname = '/Parallel_data_boundary/ranksToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_ranksToCommMshRank)

            dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInLocMshRank)

            dsetname = '/Parallel_data_boundary/commsMemSize'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_commsMemSizeMshRank)

            dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInNgbMshRank)

            ms_offset(1)=0
            do i=0,mshRank-1 !from rank 0 mpi_rank-1
               ms_offset(1)=ms_offset(1)+int(vecBndNumNodesToCommMshRank(i),hssize_t)
            end do
            ms_dims(1)=int(bnd_numNodesToCommMshRank,hsize_t)

            dsetname = '/Parallel_data_boundary/nodesToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,bnd_nodesToCommMshRank)

            !------------------------------------------------------------------------------------------------------------------------

         end if
      end if

      !---------------------------------------------------------------------------------------------------------------------
      ms_dims(1) = 1
      ms_offset(1) = int(mshRank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      aux_array(1)=numWorkingNodesMshRank
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------
      !  SAVING connecParOrig(:,:)
      allocate(aux_array(numElemsMshRank*mnnode))
      ms_dims(1) = int(numElemsMshRank,hsize_t)*int(mnnode,hsize_t)
      ms_offset(1) = int((mshRankElemStart-1),hssize_t)*int(mnnode,hssize_t)
      !-----------------------------------------------------------------------------------------------------
      i=1
      do iElemL=1,numElemsMshRank
         do m=1,mnnode
            aux_array(i)=connecParOrigMshRank(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParOrig'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParWork(:,:)
      i=1
      do iElemL=1,numElemsMshRank
         do m=1,mnnode
            aux_array(i)=connecParWorkMshRank(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParWork'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      deallocate(aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      !dsetname = '/Connectivity/connecVTK'
      !call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,connecVTKMshRank)
      !-----------------------------------------------------------------------------------------------------
      ms_offset(1)=0
      do i=0,mshRank-1 !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(vecNumWorkingNodes(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesMshRank,hsize_t)

      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,workingNodesMshRank)
      !-----------------------------------------------------------------------------------------------------

      !------------------------------------------------------------------------------------------------------------------------
      if(isPeriodic) then
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array(1))

         dsetname = '/Periodic_data/nPerRankPar'
         aux_array(1)= numPerNodesMshRank
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

         deallocate(aux_array)

         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumPerNodesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numPerNodesMshRank,hsize_t)

         dsetname = '/Periodic_data/masSlaRankPar1'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,masSlaNodesMshRank(:,1))

         dsetname = '/Periodic_data/masSlaRankPar2'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,masSlaNodesMshRank(:,2))
         !--------------------------------------------------------------------------------------------------------
      end if

      !------------------------------------------------------------------------------------------------------------------------
      dsetname = '/order/porder'
      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_array(ms_dims(1)))

      if(mshRank.eq.0) then
         aux_array(1) = mporder
      end if

      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)

      !--------------------------------------


      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = mnnode
      endif

      allocate(aux_array(ms_dims(1)))

      !-----------------------------
      dsetname = '/order/a2ijk'
      if(mshRank.eq.0) then
         aux_array(:) = a2ijk(:)
      end if
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------
      dsetname = '/order/gmsh2ijk'
      if(mshRank.eq.0) then
         aux_array(:) = gmsh2ijk(:)
      end if
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------
      dsetname = '/order/vtk2ijk'
      if(mshRank.eq.0) then
         aux_array(:) = vtk2ijk(:)
      end if
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------

      deallocate(aux_array)
      !--------------------------------------

      ms_dims(1) = 0
      ms_offset(1) = 0
      if(mshRank.eq.0) then
         ms_dims(1) = mnpbou
      endif

      allocate(aux_array(ms_dims(1)))

      !-----------------------------
      dsetname = '/order/a2ij'
      if(mshRank.eq.0) then
         aux_array(:) = a2ij(:)
      end if
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------
      dsetname = '/order/gmsh2ij'
      if(mshRank.eq.0) then
         aux_array(:) = gmsh2ij(:)
      end if
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------
      dsetname = '/order/vtk2ij'
      if(mshRank.eq.0) then
         aux_array(:) = vtk2ij(:)
      end if
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------

      deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------

   end subroutine write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine dummy_write_mshRank_data_in_hdf5_meshFile_from_tool(hdf5_file_id,numMshRanks2Part,isPeriodic,isBoundaries)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer(4),intent(in) :: numMshRanks2Part
      !-------------------------------------------------------------------------------------------------------------------------------
      character(128) :: dsetname
      integer(hsize_t) :: ms_dims(1),ms_dims2d(2)
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(4),allocatable :: empty_array_i4(:)
      integer(8),allocatable :: empty_array_i8(:)
      real(rp),allocatable :: empty_array_rp(:),empty_array2d_rp(:,:)

      !-------------------------------------------------------------------------------------------------------------------------------

      allocate(empty_array_i4(0))
      allocate(empty_array_i8(0))
      allocate(empty_array_rp(0))
      allocate(empty_array2d_rp(0,0))

      ms_dims(1) = 0
      ms_offset(1) = 0

      ms_dims2d(1) = 0
      ms_dims2d(2) = 0
      ms_offset2d(1) = 0
      ms_offset2d(2) = 0

      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/globalIds/elemGid'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      !---------------------------------------------------------------------------------------------------------------------

      dsetname = '/Coords/Points'
      call write_dataspace_2d_tr_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims2d,ms_offset2d,empty_array2d_rp)

      !---------------------------------------------------------------------------------------------------------------------
      if(numMshRanks2Part.ge.2) then
         dsetname = '/Parallel_data/rankNodeStart'
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

         dsetname = '/Parallel_data/rankNodeEnd'
         call write_dataspace_1d_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

         dsetname = '/Parallel_data/rankElemStart'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/rankElemEnd'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/numRanksWithComms'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/numNodesToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/ranksToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemPosInLoc'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemSize'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemPosInNgb'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/nodesToComm'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      end if
      !------------------------------------------------------------------------------------------------------------------------

      if(isBoundaries) then

         dsetname = '/Boundary_data/numBoundCodes'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/numBoundsRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/ndofRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         !----------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/ldofPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         !----------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/lbnodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         !----------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/bouCodesPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         !----------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/boundPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/boundParOrig'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      !------------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------
         if(numMshRanks2Part.ge.2) then
            dsetname = '/Parallel_data_boundary/numRanksWithComms'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/numNodesToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/ranksToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemSize'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/nodesToComm'
            call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

            !------------------------------------------------------------------------------------------------------------------------

         end if
      end if

      !---------------------------------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      !----------------------------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/connecParOrig'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      !---------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/connecParWork'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      !---------------------------------------------------------------------------------------------
      !dsetname = '/Connectivity/connecVTK'
      !call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetnamk,ms_dims,ms_offset,empty_array_i4)
      !---------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      !-----------------------------------------------------------------------------------------------------

      !------------------------------------------------------------------------------------------------------------------------
      if(isPeriodic) then

         dsetname = '/Periodic_data/nPerRankPar'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Periodic_data/masSlaRankPar1'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Periodic_data/masSlaRankPar2'
         call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
         !--------------------------------------------------------------------------------------------------------
      end if

      !------------------------------------------------------------------------------------------------------------------------
      dsetname = '/order/porder'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/a2ijk'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/a2ij'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/gmsh2ijk'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/gmsh2ij'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/vtk2ijk'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)

      dsetname = '/order/vtk2ij'
      call write_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,empty_array_i4)
      !------------------------------------------------------------------------------------------------------------------------

      deallocate(empty_array_i4)
      deallocate(empty_array_i8)
      deallocate(empty_array_rp)
      deallocate(empty_array2d_rp)

   end subroutine dummy_write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine load_hdf5_meshFile(mnnode,mnpbou)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou
      character(256) :: groupname,dsetname
      integer(hid_t) :: file_id,dset_id,fspace_id
      integer(4) :: h5err
      integer(hsize_t),dimension(1) :: fs_dims,fs_maxdims

      if(mpi_rank.eq.0) write(*,*) '# Loading hdf5 mesh: ',meshFile_h5_name

      call open_hdf5_file(meshFile_h5_name,file_id)

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

      !--------------------------------------------------------------------------------
      !load porder
      call load_porder_hdf5(file_id,mnnode,mnpbou)

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

         numNodesToComm=0
         numRanksWithComms=0
         allocate(nodesToComm(numNodesToComm))
         allocate(ranksToComm(numRanksWithComms))
         allocate(commsMemPosInLoc(numRanksWithComms))
         allocate(commsMemPosInNgb(numRanksWithComms))
         allocate(commsMemSize(numRanksWithComms))

      end if
      !-----------------------------------------------------------------------------------------------
      !load periodic data
      call load_periodic_data_hdf5(file_id)

      !-----------------------------------------------------------------------------------------------
      !load boundary data
      call load_boundary_data_hdf5(file_id,mnnode,mnpbou)
      if((isMeshBoundaries).and.(mpi_size.ge.2)) then
         call load_parallel_data_boundary_hdf5(file_id)
      else
         bnd_numNodesToComm=0
         bnd_numRanksWithComms=0
         allocate(bnd_nodesToComm(bnd_numNodesToComm))
         allocate(bnd_ranksToComm(bnd_numRanksWithComms))
         allocate(bnd_commsMemPosInLoc(bnd_numRanksWithComms))
         allocate(bnd_commsMemPosInNgb(bnd_numRanksWithComms))
         allocate(bnd_commsMemSize(bnd_numRanksWithComms))
      end if

      !-----------------------------------------------------------------------------------------------
      !load the coordinates
      call load_coordinates_hdf5(file_id)

      !-----------------------------------------------------------------------------------------------
      !load connectivity
      call load_connectivity_hdf5(file_id,mnnode)

      !--------------------------------------------------------------------------------
      !load globalIds
      call load_globalIds_hdf5(file_id)

      !close h5 file
      call close_hdf5_file(file_id)

      if(mpi_rank.eq.0) write(*,*) '# Mesh ',meshFile_h5_name,'succesfully loaded!'
      mesh_isLoaded = .true.

   end subroutine load_hdf5_meshFile

   subroutine create_group_hdf5(file_id,groupname)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*), intent(in) :: groupname
      integer(hid_t) :: group_id
      integer(4) :: h5err

      call h5gcreate_f(file_id,groupname,group_id,h5err)
      call h5gclose_f(group_id, h5err)
   end subroutine create_group_hdf5

   subroutine create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id
      integer(4) :: h5err

      ! Create the data space for the  dataset.
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      ! Create the dataset with default properties.
      call h5dcreate_f(file_id,dsetname,dtype,dspace_id,dset_id,h5err)

      !write(*,*) 'create dsetname ',dsetname, ' dset_id ',dset_id,' dspace_id ',dspace_id

      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine create_dataspace_hdf5

   subroutine create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims, chunk_dims, dtype)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank),max_dims(ds_rank), chunk_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id,plist_id
      integer :: h5err

      ! Create the data space for the  dataset.
      call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err,max_dims)

      ! Create the dataset with default properties.
      call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,h5err)
      call h5pset_chunk_f(plist_id,ds_rank,chunk_dims,h5err)
      call h5dcreate_f(file_id, dsetname,dtype,dspace_id,dset_id, h5err,plist_id)

      !write(*,*) 'create dsetname ',dsetname, ' dset_id ',dset_id,' dspace_id ',dspace_id
      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(dspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine create_dataspace_maxdims_hdf5

   subroutine extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank)
      integer(hid_t) :: dset_id
      integer :: h5err

      ! Open dataset
      call h5dopen_f(file_id, dsetname,dset_id,h5err)
      ! Extend the dataset to ds_dims.
      call h5dextend_f(dset_id,ds_dims,h5err)

      call h5dclose_f(dset_id,h5err)
   end subroutine extend_dataset_hdf5

   subroutine create_chunked_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,chunk_dims,dtype)
      !BE CAREFUL: THIS SUBROUTINE MUST BE DOUBLECHECKED
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ds_rank
      integer(hsize_t),intent(in) :: ds_dims(ds_rank),chunk_dims(ds_rank)
      integer(hid_t),intent(in) :: dtype
      integer(hid_t) :: dset_id,dspace_id,plist_id
      integer(4) :: h5err

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

!-------------------------------------------------------------------------------------------------------------------
   subroutine open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                             dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset
      integer(hid_t),intent(out) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank),intent(out) :: fs_dims,fs_maxdims
      integer(4) :: h5err

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
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

   end subroutine open_create_dataspace_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
   subroutine close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)
      implicit none
      integer(hid_t),intent(inout) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine close_dataspace_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!  FP(rp 4/8)

   subroutine select_dtype_rp(dtype)
      implicit none
      integer(hid_t),intent(inout) :: dtype

      if(rp.eq.4) then
         dtype = h5_datatype_real4
      else if(rp.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in select_dtype_rp! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
   end subroutine select_dtype_rp

   subroutine select_dtype_rp_vtk(dtype)
      implicit none
      integer(hid_t),intent(inout) :: dtype

      if(rp_vtk.eq.4) then
         dtype = h5_datatype_real4
      else if(rp_vtk.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in select_dtype_rp_vtk! rp_vtk is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
   end subroutine select_dtype_rp_vtk

   subroutine write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call select_dtype_rp(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_real_rp_hyperslab_parallel

   subroutine read_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call select_dtype_rp(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_real_rp_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(in) :: array2d(ms_dims(2),ms_dims(1)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      !!!!$acc kernels
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d_tr(ii,jj)=array2d(jj,ii)
         end do
      end do
      !!!!$acc end kernels
      !!!!$acc update device(host(:,:)) !THINK ABOUT IT!!Update the value in the host before saving in hdf5

      call select_dtype_rp(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_2d_tr_real_rp_hyperslab_parallel

   subroutine read_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(out) :: array2d(ms_dims(2),ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      call select_dtype_rp(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

      !!!!$acc update device(array2d(:,:)) !put in device memory the value
      !!!!$acc kernels
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d(jj,ii) = array2d_tr(ii,jj)
         end do
      end do
      !!!!$acc end kernels
      !!!!$acc update host(array2d(:,:)) ! put in device memory the value

   end subroutine read_dataspace_2d_tr_real_rp_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call select_dtype_rp_vtk(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_real_rp_vtk_hyperslab_parallel

   subroutine read_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter:: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call select_dtype_rp_vtk(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_real_rp_vtk_hyperslab_parallel
!--------------------------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(in) :: array2d(ms_dims(2),ms_dims(1)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp_vtk) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      !!!!!$acc kernels !!BE CAREFUL HtoD/DtoH requrired if to be used!!!
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d_tr(ii,jj)=array2d(jj,ii)
         end do
      end do
      !!!!!$acc end kernels

      call select_dtype_rp_vtk(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel

   subroutine read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array2d)
      implicit none
      integer(4),parameter :: ms_rank = 2
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(out) :: array2d(ms_dims(2),ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: ii,jj,h5err
      real(rp_vtk) :: array2d_tr(ms_dims(1),ms_dims(2)) !FORTRAN is COLUMN-MAJOR & HDF5 is ROW-MAJOR

      call select_dtype_rp_vtk(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,array2d_tr,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

      !!!!$acc kernels !!BE CAREFUL HtoD/DtoH requrired if to be used!!!
      do ii=1,ms_dims(1)
         do jj=1,ms_dims(2)
            array2d(jj,ii) = array2d_tr(ii,jj)
         end do
      end do
      !!!$acc end kernels

   end subroutine read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel

!-------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------
!  FP32
   subroutine write_dataspace_1d_fp32_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(4),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_real4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_fp32_hyperslab_parallel

   subroutine read_dataspace_1d_fp32_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(4),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_real4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_fp32_hyperslab_parallel

   subroutine read_dataspace_2d_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer,intent(in) :: ms_rank !assuming ms_rank=fs_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset
      real(4),intent(out) :: data(ms_dims(1), ms_dims(2))
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
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine read_dataspace_2d_fp32_hyperslab_parallel

!-------------------------------------------------------------------------------------------------------------------
!  FP64
   subroutine write_dataspace_1d_fp64_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(8),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_real8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_fp64_hyperslab_parallel

   subroutine read_dataspace_1d_fp64_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(8),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_real8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_fp64_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!  INT1
   subroutine write_dataspace_1d_int1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter:: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_int1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_int1_hyperslab_parallel



   subroutine read_dataspace_1d_int1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_int1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_int1_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  UINT1
   subroutine write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_uint1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_uint1_hyperslab_parallel

   subroutine read_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_uint1,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_uint1_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  INT4
   subroutine write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(4),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_int4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_int4_hyperslab_parallel

   subroutine read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(4),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_int4,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_int4_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  INT8
   subroutine write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(8),intent(in) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_int8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_1d_int8_hyperslab_parallel


   subroutine read_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,array1d)
      implicit none
      integer(4),parameter :: ms_rank = 1
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(8),intent(out) :: array1d(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_int8,array1d,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_1d_int8_hyperslab_parallel

!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------
#if 0
   !DEPRECATED
   subroutine overwrite_coordinates_hdf5()
      implicit none
      integer(hid_t) :: file_id,plist_id
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4) :: ms_rank,h5err
      character(128) :: dsetname

      !---------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! open file collectively
      call h5fopen_f(meshFile_h5_name, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot open meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !-------------------------------------------
      ms_rank = 1
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-------------------------------------------

      dsetname = '/Coords/X'
      call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,coordPar(:,1))

      dsetname = '/Coords/Y'
      call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,coordPar(:,2))

      dsetname = '/Coords/Z'
      call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,coordPar(:,3))


      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine overwrite_coordinates_hdf5
#endif
   subroutine save_connectivity_hdf5(file_id,mnnode)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank
      integer(4) :: iElemL,i,m,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Connectivity')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = totalNumElements*mnnode
      ms_dims(1) = int(numElemsRankPar,hsize_t)*int(mnnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)*int(mnnode,hssize_t)

      allocate(aux_array(numElemsRankPar*mnnode))
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParOrig(:,:)
      i=1
      do iElemL=1,numElemsRankPar
         do m=1,mnnode
            aux_array(i)=connecParOrig(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParWork(:,:)
      i=1
      do iElemL=1,numElemsRankPar
         do m=1,mnnode
            aux_array(i)=connecParWork(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParWork'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      !dsetname = '/Connectivity/connecVTK'
      !call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,connecVTK)
      !-----------------------------------------------------------------------------------------------------

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      !  SAVING numWorkingNodesRankPar
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numWorkingNodesRankPar
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !--------------------------------------------------------------------------------------------------------
      !  SAVING workingNodesPar(:)
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numWorkingNodesRankPar of all ranks
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesRankPar,hsize_t)

      dsetname = '/Connectivity/workingNodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,workingNodesPar)
      !--------------------------------------------------------------------------------------------------------

   end subroutine save_connectivity_hdf5

   subroutine load_connectivity_hdf5(file_id,mnnode)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: iElemL,i,m
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading connectivity data hdf5...'

      ms_dims(1) = int(numElemsRankPar,hsize_t)*int(mnnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)*int(mnnode,hssize_t)

      !-------------------------------------------------------------------------------------------------------
      allocate( connecParOrig(numElemsRankPar,mnnode) )
      allocate( connecParWork(numElemsRankPar,mnnode) )
      !$acc enter data create(connecParOrig(:,:))
      !$acc enter data create(connecParWork(:,:))
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecParOrig(:,:)
      allocate(aux_array(numElemsRankPar*mnnode))

      dsetname = '/Connectivity/connecParOrig'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      i=1
      do iElemL=1,numElemsRankPar
         do m=1,mnnode
            connecParOrig(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do
      !$acc update device(connecParOrig(:,:))
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecParWork(:,:)
      dsetname = '/Connectivity/connecParWork'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      i=1
      do iElemL=1,numElemsRankPar
         do m=1,mnnode
            connecParWork(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do

      deallocate(aux_array)
      !$acc update device(connecParWork(:,:))
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecVTK(:)
      !allocate( connecVTK(numElemsRankPar*mnnode) )
      !dsetname = '/Connectivity/connecVTK'
      !call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,connecVTK)
      !-------------------------------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------------------------------
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))
      !-------------------------------------------------------------------------------------------------------
      !LOADING numWorkingNodesRankPar
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numWorkingNodesRankPar=aux_array(1)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      !LOADING workingNodesPar
      allocate(workingNodesPar(numWorkingNodesRankPar))
      !$acc enter data create(workingNodesPar(:))

      allocate(aux_array(mpi_size))
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numWorkingNodesRankPar of all ranks
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesRankPar,hsize_t)

      dsetname = '/Connectivity/workingNodesPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,workingNodesPar)
      deallocate(aux_array)
      !-------------------------------------------------------------------------------------------------------
      !$acc update device(workingNodesPar(:))

   end subroutine load_connectivity_hdf5
#if 0
   subroutine save_parallel_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,ms_rank,h5err
      integer(4) :: i,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Parallel_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Parallel_data/rankNodeStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankNodeStart
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankNodeEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankNodeEnd
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankElemStart'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankElemStart
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/rankElemEnd'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=rankElemEnd
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numRanksWithComms
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numNodesToComm
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numRanksWithComms,hsize_t)

      !write(*,*) 'ms_offset ',ms_offset(1),' ms_dims(1) ',ms_dims(1),' ds_dims ',ds_dims(1)

      dsetname = '/Parallel_data/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ranksToComm)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInLoc)

      dsetname = '/Parallel_data/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemSize)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInNgb)

      ds_dims(1) = mpi_size
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0

      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numNodesToComm,hsize_t)

      dsetname = '/Parallel_data/matrixCommScheme_iNodeL'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,1))

      dsetname = '/Parallel_data/matrixCommScheme_iNodeGSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,2))

      dsetname = '/Parallel_data/matrixCommScheme_ngbRank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,matrixCommScheme(:,3))


      deallocate(aux_array)

   end subroutine save_parallel_data_hdf5
#endif
   subroutine load_parallel_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Parallel_data/rankNodeStart'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankNodeStart=aux_array(1)

      dsetname = '/Parallel_data/rankNodeEnd'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankNodeEnd=aux_array(1)

      dsetname = '/Parallel_data/rankElemStart'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankElemStart=aux_array(1)

      dsetname = '/Parallel_data/rankElemEnd'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      rankElemEnd=aux_array(1)

      numNodesRankPar = rankNodeEnd - rankNodeStart + 1
      numElemsRankPar  = rankElemEnd - rankElemStart + 1

      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numRanksWithComms=aux_array(1)

      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      numNodesToComm=aux_array(1)

      !write(*,*) 'rns ',rankNodeStart,' rne ',rankNodeEnd,' res ',rankElemStart,' ree ',rankElemEnd
      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numRanksWithComms,hsize_t)

      allocate(ranksToComm(numRanksWithComms))
      allocate(commsMemPosInLoc(numRanksWithComms))
      allocate(commsMemPosInNgb(numRanksWithComms))
      allocate(commsMemSize(numRanksWithComms))
      !$acc enter data create(ranksToComm(:))
      !$acc enter data create(commsMemPosInLoc(:))
      !$acc enter data create(commsMemSize(:))

      dsetname = '/Parallel_data/ranksToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,ranksToComm)
      !$acc update device(ranksToComm(:))

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,commsMemPosInLoc)
      !$acc update device(commsMemPosInLoc(:))

      dsetname = '/Parallel_data/commsMemSize'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,commsMemSize)
      !$acc update device(commsMemSize(:))

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,commsMemPosInNgb)

      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0

      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numNodesToComm,hsize_t)

      allocate(nodesToComm(numNodesToComm))
      !$acc enter data create(nodesToComm(:))

      dsetname = '/Parallel_data/nodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,nodesToComm)
      !$acc update device(nodesToComm(:))

      deallocate(aux_array)

   end subroutine load_parallel_data_hdf5

!--------------------------------------------------------------------------------------------------------------------------------
#if 0
   subroutine save_parallel_data_boundary_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,ms_rank,h5err
      integer(4) :: i,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Parallel_data_boundary')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=bnd_numRanksWithComms
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=bnd_numNodesToComm
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numRanksWithComms,hsize_t)

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_ranksToComm)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInLoc)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemSize)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInNgb)

      ds_dims(1) = int(mpi_size,hsize_t)
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numNodesToComm,hsize_t)

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeL'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,1))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_iNodeGSrl'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,2))

      dsetname = '/Parallel_data_boundary/matrixCommScheme_ngbRank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_matrixCommScheme(:,3))

      deallocate(aux_array)

   end subroutine save_parallel_data_boundary_hdf5
#endif
   subroutine load_parallel_data_boundary_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      bnd_numRanksWithComms=aux_array(1)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      bnd_numNodesToComm=aux_array(1)

      !write(*,*) 'rns ',rankNodeStart,' rne ',rankNodeEnd,' res ',rankElemStart,' ree ',rankElemEnd
      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numRanksWithComms,hsize_t)

      allocate(bnd_ranksToComm(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInLoc(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInNgb(bnd_numRanksWithComms))
      allocate(bnd_commsMemSize(bnd_numRanksWithComms))
      !$acc enter data create(bnd_ranksToComm(:))
      !$acc enter data create(bnd_commsMemPosInLoc(:))
      !$acc enter data create(bnd_commsMemSize(:))

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_ranksToComm)
      !$acc update device(bnd_ranksToComm(:))

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInLoc)
      !$acc update device(bnd_commsMemPosInLoc(:))

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_commsMemSize)
      !$acc update device(bnd_commsMemSize(:))

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_commsMemPosInNgb)

      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numNodesToComm,hsize_t)

      allocate(bnd_nodesToComm(bnd_numNodesToComm))
      !$acc enter data create(bnd_nodesToComm(:))

      dsetname = '/Parallel_data_boundary/nodesToComm'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bnd_nodesToComm)
      !$acc update device(bnd_nodesToComm(:))

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
      integer(4) :: ds_rank,h5err
      integer(4) :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Periodic_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Periodic_data/nPerRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=nPerRankPar
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set nPerRankPar of all ranks
      dsetname = '/Periodic_data/nPerRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(nPerRankPar,hsize_t)

      dsetname = '/Periodic_data/masSlaRankPar1'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,masSlaRankPar(:,1))

      dsetname = '/Periodic_data/masSlaRankPar2'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,masSlaRankPar(:,2))
      !--------------------------------------------------------------------------------------------------------

   end subroutine save_periodic_data_hdf5

   subroutine load_periodic_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,h5err
      integer(4) :: i,accumVal,iBound,m
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

         dtype = h5_datatype_int4
         ds_rank = 1
         ds_dims(1) = mpi_size
         ms_dims(1) = 1
         ms_offset(1) = int(mpi_rank,hssize_t)
         allocate(aux_array(1))

         !--------------------------------------------------------------------------------------------------------
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         nPerRankPar=aux_array(1)

         allocate(masSlaRankPar(nPerRankPar,2))
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
         allocate(aux_array(mpi_size))
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         !read data set numBoundsRankPar of all ranks
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(nPerRankPar,hsize_t)

         dsetname = '/Periodic_data/masSlaRankPar1'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,masSlaRankPar(:,1))

         dsetname = '/Periodic_data/masSlaRankPar2'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,masSlaRankPar(:,2))
         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
      else
         isMeshPeriodic = .false.
      end if

   end subroutine load_periodic_data_hdf5
!--------------------------------------------------------------------------------------------------------------------------------
   subroutine save_boundary_data_hdf5(file_id,mnpbou)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnpbou
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,h5err
      integer(4) :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Boundary_data/numBoundsRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundsRankPar
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundCodes'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundCodes
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/ndofRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=ndofRankPar
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=numBoundaryNodesRankPar
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numBoundsRankPar of all ranks
      dsetname = '/Boundary_data/numBoundsRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numBoundsRankPar,hsize_t)

      !write(*,*) 'ms_offset ',ms_offset(1),' ms_dims(1) ',ms_dims(1),' ds_dims ',ds_dims(1)

      dsetname = '/Boundary_data/bouCodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bouCodesPar)

      ds_dims(1) = int(ds_dims(1),hsize_t)*int(mnpbou,hsize_t)!totalNumBoundsSrl*mnpbou
      ms_dims(1) = int(ms_dims(1),hsize_t)*int(mnpbou,hsize_t)!numBoundsRankPar*mnpbou
      ms_offset(1) = ms_offset(1)*int(mnpbou,hssize_t)

      deallocate(aux_array)
      allocate(aux_array(numBoundsRankPar*mnpbou))

      !boundPar
      i=1
      do iBound=1,numBoundsRankPar
         do m=1,mnpbou
            aux_array(i)=boundPar(iBound,m)
            i=i+1
         end do
      end do

      dsetname = '/Boundary_data/boundPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      !boundParOrig
      i=1
      do iBound=1,numBoundsRankPar
         do m=1,mnpbou
            aux_array(i)=boundParOrig(iBound,m)
            i=i+1
         end do
      end do

      dsetname = '/Boundary_data/boundParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set ndofRankPar of all ranks
      dsetname = '/Boundary_data/ndofRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(ndofRankPar,hsize_t)

      dsetname = '/Boundary_data/ldofPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,ldofPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numBoundaryNodesRankPar of all ranks
      dsetname = '/Boundary_data/numBoundaryNodesRankPar'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      ds_dims(1) = accumVal

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numBoundaryNodesRankPar,hsize_t)

      dsetname = '/Boundary_data/lbnodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,lbnodesPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
   end subroutine save_boundary_data_hdf5

   subroutine load_boundary_data_hdf5(file_id,mnnode,mnpbou)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode,mnpbou
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,h5err
      integer(4) :: i,accumVal,iBound,m,iNodeL
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)
      logical :: isBoundaryFolder

      groupname = trim('/Boundary_data')
      !call create_group_hdf5(file_id,groupname)

      call h5lexists_f(file_id,groupname,isBoundaryFolder,h5err)

      isMeshBoundaries = isBoundaryFolder

      if(mpi_rank.eq.0) write(*,*) 'Loading Boundary data hdf5. -> isBoundaryFolder:',isBoundaryFolder

      if(isMeshBoundaries) then
         dtype = h5_datatype_int4
         ds_rank = 1
         ds_dims(1) = mpi_size
         ms_dims(1) = 1
         ms_offset(1) = 0
         allocate(aux_array(1))

         dsetname = '/Boundary_data/numBoundCodes'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         numBoundCodes=aux_array(1)

         ms_offset(1) = int(mpi_rank,hssize_t)

         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         numBoundsRankPar=aux_array(1)

         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         ndofRankPar=aux_array(1)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
         numBoundaryNodesRankPar=aux_array(1)

         !--------------------------------------------------------------------------------------------------------
         allocate(boundPar(numBoundsRankPar,mnpbou))
         allocate(boundParOrig(numBoundsRankPar,mnpbou))
         allocate(bouCodesPar(numBoundsRankPar))
         allocate(ldofPar(ndofRankPar))
         allocate(lbnodesPar(numBoundaryNodesRankPar))
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
         allocate( aux_array(mpi_size) )
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         !read data set numBoundsRankPar of all ranks
         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)
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
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundsRankPar,hsize_t)

         dsetname = '/Boundary_data/bouCodesPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,bouCodesPar)
         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
         ds_dims(1) = int(ds_dims(1),hsize_t)*int(mnpbou,hsize_t)!totalNumBoundsSrl*mnpbou
         ms_dims(1) = int(ms_dims(1),hsize_t)*int(mnpbou,hsize_t)!numBoundsRankPar*mnpbou
         ms_offset(1) = ms_offset(1)*int(mnpbou,hssize_t)

         allocate(aux_array(numBoundsRankPar*mnnode))

         !boundPar
         dsetname = '/Boundary_data/boundPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

         i=1
         do iBound=1,numBoundsRankPar
            do m=1,mnpbou
               boundPar(iBound,m)=aux_array(i)
               i=i+1
            end do
         end do

         !boundParOrig
         dsetname = '/Boundary_data/boundParOrig'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

         i=1
         do iBound=1,numBoundsRankPar
            do m=1,mnpbou
               boundParOrig(iBound,m)=aux_array(i)
               i=i+1
            end do
         end do

         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
         allocate( aux_array(mpi_size) )
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         !read data set ndofRankPar of all ranks
         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(ndofRankPar,hsize_t)

         dsetname = '/Boundary_data/ldofPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,ldofPar)
         !--------------------------------------------------------------------------------------------------------
         ms_dims(1) = int(mpi_size,hsize_t)
         ms_offset(1) = 0
         !read data set numBoundaryNodesRankPar of all ranks
         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

         ds_dims(1)=0
         do i=1,mpi_size
            ds_dims(1)=ds_dims(1)+aux_array(i)
         end do

         ms_offset(1)=0
         do i=1,(mpi_rank) !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundaryNodesRankPar,hsize_t)

         dsetname = '/Boundary_data/lbnodesPar'
         call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,lbnodesPar)
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
      else
         numBoundsRankPar=0
         allocate(boundPar(numBoundsRankPar,mnpbou))
         allocate(boundParOrig(numBoundsRankPar,mnpbou))
         allocate(bouCodesPar(numBoundsRankPar))

         numBoundaryNodesRankPar=0
         ndofRankPar = numNodesRankPar
         allocate(ldofPar(ndofRankPar))
         allocate(lbnodesPar(numBoundaryNodesRankPar))
         do iNodeL = 1,ndofRankPar
            ldofPar(iNodeL) = iNodeL
         end do
      end if
   end subroutine load_boundary_data_hdf5

!--------------------------------------------------------------------------------------------------------------------------------

   subroutine load_coordinates_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(2) :: ms_dims2d
      integer(hssize_t), dimension(2) :: ms_offset2d

      allocate(coordPar(numNodesRankPar,ndime))
      !$acc enter data create(coordPar(:,:))

      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1

      !dsetname = '/VTKHDF/Points'
      dsetname = '/Coords/Points'
      call read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,coordPar)

      !$acc update device(coordPar(:,:))

   end subroutine load_coordinates_hdf5

   subroutine load_globalIds_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4) :: iNodeL
      integer(8) :: iNodeGSrl,max_iNodeGSrl_l,max_iNodeGSrl_g

      allocate(globalIdSrl(numNodesRankPar))
      allocate(globalIdPar(numNodesRankPar))
      allocate(elemGid(numElemsRankPar))

      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = '/globalIds/globalIdSrl'
      call read_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,globalIdSrl)

      dsetname = '/globalIds/globalIdPar'
      call read_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,globalIdPar)

      ms_dims(1) = int(numElemsRankPar,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      dsetname = '/globalIds/elemGid'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,elemGid)

      !--------------------------------------------------------------------------------------------------
      !setting totalNumNodesSrl
      max_iNodeGSrl_l=0
      do iNodeL=1,numNodesRankPar
         iNodeGSrl = globalIdSrl(iNodeL)
         max_iNodeGSrl_l = max(iNodeGSrl,max_iNodeGSrl_l)
      end do

      call MPI_Allreduce(max_iNodeGSrl_l,max_iNodeGSrl_g,1,mpi_datatype_int8,MPI_MAX,app_comm,mpi_err)

      totalNumNodesSrl = max_iNodeGSrl_g
      !write(*,*) 'setting totalNumNodesSrl ',totalNumNodesSrl
      !--------------------------------------------------------------------------------------------------
   end subroutine load_globalIds_hdf5

   subroutine load_porder_hdf5(file_id,mnnode,mnpbou)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode,mnpbou
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array(:)

      !------------------------------------------------------------------------------------------------
      allocate(aux_array(1))
      ms_dims(1) = 1
      ms_offset(1) = 0

      dsetname = '/order/porder'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array)

      mesh_porder = aux_array(1)
      deallocate(aux_array)

      if(mesh_porder .ne. porder) then
         write(*,*) 'FATAL ERROR! mesh_porder',mesh_porder,' different to porder',porder
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      !------------------------------------------------------------------------------------------------
      allocate(mesh_a2ijk(mnnode))
      allocate(mesh_gmsh2ijk(mnnode))
      allocate(mesh_vtk2ijk(mnnode))
      ms_dims(1) = mnnode
      ms_offset(1) = 0

      dsetname = '/order/a2ijk'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_a2ijk)

      dsetname = '/order/gmsh2ijk'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_gmsh2ijk)

      dsetname = '/order/vtk2ijk'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_vtk2ijk)

      !------------------------------------------------------------------------------------------------
      allocate(mesh_a2ij(mnpbou))
      allocate(mesh_gmsh2ij(mnpbou))
      allocate(mesh_vtk2ij(mnpbou))
      ms_dims(1) = mnpbou
      ms_offset(1) = 0

      dsetname = '/order/a2ij'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_a2ij)

      dsetname = '/order/gmsh2ij'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_gmsh2ij)

      dsetname = '/order/vtk2ij'
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,mesh_vtk2ij)

      !--------------------------------------------------------------------------------------------------
   end subroutine load_porder_hdf5

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
      full_fileName = trim(adjustl(base_resultsFile_h5_name))//trim(aux_step)//'.hdf'
   end subroutine set_hdf5_resultsFile_name

   subroutine set_hdf5_restartFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_restartFile_h5_name))//trim(aux_step)//'.h5'
   end subroutine set_hdf5_restartFile_name

   subroutine set_hdf5_avgResultsFile_name(iStep,full_fileName)
      implicit none
      integer, intent(in) :: iStep
      character(len=*), intent(out) :: full_fileName
      character(len=12) :: aux_step

      write(aux_step,'(I0)') iStep
      full_fileName = trim(adjustl(base_avgResultsFile_h5_name))//trim(aux_step)//'.hdf'
   end subroutine set_hdf5_avgResultsFile_name

   subroutine set_hdf5_surface_resultsFile_name(surf_res_fileName,res_fileName)
      implicit none
      !integer(4),intent(in) :: surfCode
      character(len=*), intent(out) :: surf_res_fileName
      character(len=*), intent(in)  :: res_fileName
      !character(len=12) :: aux_surfCode

      !write(aux_surfCode,'(I0)') surfCode
      !surf_res_fileName = 'surface_'//trim(aux_surfCode)//'_'//trim(adjustl(res_fileName))
      surf_res_fileName = 'surface_'//trim(adjustl(res_fileName))
   end subroutine set_hdf5_surface_resultsFile_name

!----------------------------------------------------------------------------------------------------------------------------------

   subroutine create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(4),intent(in) :: ds_rank
      integer(hsize_t),dimension(ds_rank),intent(in) :: ds_dims
      integer(hid_t) :: dtype

      if(rp_vtk.eq.4) then
         dtype = h5_datatype_real4
      else if(rp_vtk.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in create_dataspace_for_rp_vtk_hdf5! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

   end subroutine create_dataspace_for_rp_vtk_hdf5

   subroutine save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,data_array_rp,isCreateDataspaceOpt)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(1),intent(in) :: ds_dims
      integer(hsize_t),dimension(1),intent(in) :: ms_dims
      integer(hssize_t),dimension(1),intent(in) :: ms_offset
      real(rp),intent(in) :: data_array_rp(ms_dims(1))
      logical, intent(in), optional :: isCreateDataspaceOpt
      integer(4) :: ds_rank = 1 !it is forced
      logical :: isCreateDataspace
      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)
   !---------------------------------------------------------------------------------------------------
      isCreateDataspace = .true.
      if (present(isCreateDataspaceOpt)) then
         isCreateDataspace = isCreateDataspaceOpt
      end if
      if (isCreateDataspace) then
         call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)
      end if

      if(rp .eq. rp_vtk) then
         call write_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
      else
         !copying in the aux array
         allocate(aux_data_array_rp_vtk(ms_dims(1)))
         aux_data_array_rp_vtk(:) = real(data_array_rp(:),rp_vtk)
         call write_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_data_array_rp_vtk)
         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine save_array1D_rp_in_dataset_hdf5_file

   subroutine save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,data_array_rp,isCreateDataspaceOpt)
      implicit none
      integer(4),parameter :: ds_rank = 2 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(ds_rank),intent(in) :: ds_dims2d,ms_dims2d
      integer(hssize_t),dimension(ds_rank),intent(inout) :: ms_offset2d
      logical, intent(in), optional :: isCreateDataspaceOpt
      real(rp),intent(in) :: data_array_rp(ms_dims2d(2),ms_dims2d(1)) !fortran is column-major & hdf5 writes in row-major
      logical :: isCreateDataspace
      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:,:)
   !---------------------------------------------------------------------------------------------------
      isCreateDataspace = .true.
      if (present(isCreateDataspaceOpt)) then
         isCreateDataspace = isCreateDataspaceOpt
      end if
      if (isCreateDataspace) then
         call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
      end if

      if(rp .eq. rp_vtk) then
         call write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp)
      else
         allocate(aux_data_array_rp_vtk(ms_dims2d(2),ms_dims2d(1)))
         aux_data_array_rp_vtk(:,:) = real(data_array_rp(:,:),rp_vtk)
         call write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,aux_data_array_rp_vtk)
         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine save_array2D_tr_rp_in_dataset_hdf5_file
!------------------------------------------------------------------------------------------------------------------------------

   subroutine read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(HSIZE_T),dimension(1),intent(in) :: ms_dims
      integer(HSSIZE_T),dimension(1),intent(in) :: ms_offset
      real(rp),intent(out) :: data_array_rp(ms_dims(1))

      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)
   !---------------------------------------------------------------------------------------------------

      if(rp .eq. rp_vtk) then
         call read_dataspace_1d_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
      else
         !copying in the aux array
         allocate(aux_data_array_rp_vtk(ms_dims(1)))

         call read_dataspace_1d_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_data_array_rp_vtk)
         data_array_rp(:) = real(aux_data_array_rp_vtk(:),rp)

         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine read_array1D_rp_in_dataset_hdf5_file

   subroutine read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(2),intent(in) :: ms_dims2d
      integer(hssize_t),dimension(2),intent(inout) :: ms_offset2d
      real(rp),intent(out) :: data_array_rp(ms_dims2d(2),ms_dims2d(1)) !fortran is column-major & hdf5 writes in row-major

      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:,:)
   !------------------------------------------------------------------v---------------------------------

      if(rp .eq. rp_vtk) then
         call read_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,data_array_rp)
      else
         allocate(aux_data_array_rp_vtk(ms_dims2d(2),ms_dims2d(1)))
         call read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,aux_data_array_rp_vtk)
         data_array_rp(:,:) = real(aux_data_array_rp_vtk(:,:),rp)
         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine read_array2D_tr_rp_in_dataset_hdf5_file

   subroutine save_int4_in_dataset_hdf5_file(file_id,dsetname,int2save)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(4),intent(in) :: int2save
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err
      integer(4),allocatable :: aux_data_array_int4(:)

      !----------------------------------------------------------------------------------------------------------
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_data_array_int4(ms_dims(1)))

      if(mpi_rank.eq.0) then
         aux_data_array_int4(1) = int2save
      end if

      !---------------------------------------------------------------------------------------------------

      dtype = h5_datatype_int4
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,aux_data_array_int4,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

      deallocate(aux_data_array_int4)

   end subroutine save_int4_in_dataset_hdf5_file

   subroutine read_int4_in_dataset_hdf5_file(file_id,dsetname,int2read)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(4),intent(out) :: int2read
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err
      integer(4),allocatable :: aux_data_array_int4(:)

      !----------------------------------------------------------------------------------------------------------
      ms_dims(1) = 1
      ds_dims = 1
      ms_offset(1) = 0

      allocate(aux_data_array_int4(ms_dims(1)))

      !---------------------------------------------------------------------------------------------------

      dtype = h5_datatype_int4

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,aux_data_array_int4,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

      int2read = aux_data_array_int4(1)

      deallocate(aux_data_array_int4)

   end subroutine read_int4_in_dataset_hdf5_file

   subroutine save_real_rp_in_dataset_hdf5_file(file_id,dsetname,real2save)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      real(rp),intent(in) :: real2save
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err

      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)

      !----------------------------------------------------------------------------------------------------------
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif

      allocate(aux_data_array_rp_vtk(ms_dims(1)))

      if(mpi_rank.eq.0) then
         aux_data_array_rp_vtk(1) = real(real2save,rp_vtk)
      end if

      !---------------------------------------------------------------------------------------------------

      call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)

      call select_dtype_rp_vtk(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,aux_data_array_rp_vtk,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

      deallocate(aux_data_array_rp_vtk)

   end subroutine save_real_rp_in_dataset_hdf5_file

   subroutine read_real_rp_in_dataset_hdf5_file(file_id,dsetname,real2read)
      implicit none
      integer(4),parameter :: ds_rank = 1, ms_rank = 1 !it is forced
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      real(rp),intent(out) :: real2read
      integer(hsize_t),dimension(ms_rank) :: ds_dims,ms_dims,fs_dims,fs_maxdims
      integer(hssize_t),dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(4) :: h5err

      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)

      !----------------------------------------------------------------------------------------------------------
      ms_dims(1) = 1
      ds_dims = 1
      ms_offset(1) = 0

      allocate(aux_data_array_rp_vtk(ms_dims(1)))

      !---------------------------------------------------------------------------------------------------

      call select_dtype_rp_vtk(dtype)

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,aux_data_array_rp_vtk,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

      real2read = real(aux_data_array_rp_vtk(1),rp)

      deallocate(aux_data_array_rp_vtk)

   end subroutine read_real_rp_in_dataset_hdf5_file

!----------------------------------------------------------------------------------------------------------------------------------

   subroutine save_hdf5_restartFile(mnnode,mngaus,restartCnt,iStep,flag_walave,flag_buffer,time,rho,u,pr,E,mu_e,mu_t,walave_u,u_buffer)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus
      integer(4),intent(in) :: restartCnt,iStep,flag_walave
      logical,intent(in) :: flag_buffer
      real(rp),intent(in) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u,walave_u,u_buffer
      real(rp),intent(inout),dimension(numElemsRankPar,mngaus) :: mu_e,mu_t

      integer(hid_t) :: file_id,plist_id,dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4) :: ds_rank,ms_rank,h5err
      character(512) :: full_fileName,dsetname

      integer(4) :: aux_array_i4(1)

      real(rp),dimension(numNodesRankPar) :: aux_mu_e,aux_mu_t
      integer(4) :: iElem,iNode

      !$acc kernels
      do iElem = 1,numElemsRankPar
         do iNode = 1,mnnode
            aux_mu_e(connecParOrig(iElem,iNode)) = mu_e(iElem,iNode)
            aux_mu_t(connecParOrig(iElem,iNode)) = mu_t(iElem,iNode)
         end do
      end do
      !$acc end kernels

      !$acc update host(rho(:))
      !$acc update host(u(:,:))
      !$acc update host(pr(:))
      !$acc update host(E(:))
      if(flag_walave == 1) then
         !$acc update host(walave_u(:,:))
      end if

      if(flag_buffer) then
         !$acc update host(u_buffer(:,:))
      end if

      !-----------------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_restartFile_name(restartCnt,full_fileName)

      call create_hdf5_file(full_fileName,file_id)

      !-----------------------------------------------------------------------------------------------
      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------

      dsetname = 'rho'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,3))

      dsetname = 'pr'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,E)

      dsetname = 'mue'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_mu_e)

      dsetname = 'mut'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_mu_t)

      if(flag_walave == 1) then
         dsetname = 'walave_u_x'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,walave_u(:,1))

         dsetname = 'walave_u_y'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,walave_u(:,2))

         dsetname = 'walave_u_z'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,walave_u(:,3))
      end if

      if(flag_buffer) then
         dsetname = 'u_buffer_x'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u_buffer(:,1))

         dsetname = 'u_buffer_y'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u_buffer(:,2))

         dsetname = 'u_buffer_z'
         call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u_buffer(:,3))
      end if

      !-----------------------------------------------------------------------------------------------
      ! ----  time  -----
      dsetname = 'time'
      call save_real_rp_in_dataset_hdf5_file(file_id,dsetname,time)

      !ms_dims(1) = 0
      !ds_dims = 1
      !ms_offset(1) = 0
      !if(mpi_rank.eq.0) then
      !   ms_dims(1) = 1
      !endif
      !aux_array_rp(1) = time
      !call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_array_rp)

      !-----------------------------------------------------------------------------------------------
      ! ---- istep -----
      dsetname = 'istep'
      call save_int4_in_dataset_hdf5_file(file_id,dsetname,iStep)
      !dtype = h5_datatype_int4
      !ds_rank = 1
      !ms_rank = 1
      !aux_array_i4(1) = iStep

      !dsetname = 'istep'
      !call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      !call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i4)

      !close the file.
      call close_hdf5_file(file_id)

   end subroutine save_hdf5_restartFile

   subroutine load_hdf5_restartFile(mnnode,mngaus,restartCnt,load_step,flag_walave,flag_buffer,time,rho,u,pr,E,mu_e,mu_t,walave_u,u_buffer)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,restartCnt,flag_walave
      logical,intent(in) :: flag_buffer
      integer(4),intent(inout) :: load_step
      real(rp),intent(inout) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u,walave_u,u_buffer
      real(rp),intent(inout),dimension(numElemsRankPar,mngaus) :: mu_e,mu_t

      character(512) :: full_restartFileName
      integer(hid_t) :: file_id,plist_id
      integer(hsize_t),dimension(1) :: ms_dims
      integer(hssize_t),dimension(1) :: ms_offset
      integer(4) :: ms_rank,iPer,h5err
      character(128) :: dsetname

      real(rp) :: aux_array_rp(1)
      integer(4) :: aux_array_i4(1)

      real(rp),dimension(numNodesRankPar) :: aux_mu_e,aux_mu_t
      integer(4) :: iElem,iNode
      logical    :: link_exists

      !----------------------------------------------------------------------------------------------------------------------

      if((restartCnt.ne.1).and.(restartCnt.ne.2)) then
         write(*,*) 'FATAL ERROR in load_hdf5_restartFile! restartFile to load must be 1 or 2 and is',restartCnt,'CRASHING!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      call set_hdf5_restartFile_name(restartCnt,full_restartFileName)
      if(mpi_rank.eq.0) write(*,*) '# Loading restart file: ',trim(adjustl(full_restartFileName))

      call open_hdf5_file(full_restartFileName,file_id)

      !ms_rank = 1
      !ms_offset(1) = 0
      !ms_dims(1) = 1
      ! ----  read time  --------------------------------------------------------------------------

      dsetname = 'time'
      call read_real_rp_in_dataset_hdf5_file(file_id,dsetname,time)
      !call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_array_rp)
      !time = aux_array_rp(1)

      ! ----  read istep --------------------------------------------------------------------------
      dsetname = 'istep'
      call read_int4_in_dataset_hdf5_file(file_id,dsetname,load_step)
      !call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i4)
      !load_step = aux_array_i4(1)

      ! ----  read arrays  --------------------------------------------------------------------------
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'rho'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,3))

      dsetname = 'pr'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,E)

      dsetname = 'mue'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_mu_e)

      dsetname = 'mut'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_mu_t)

      if(flag_walave==1) then
         dsetname = 'walave_u_x'
         call h5lexists_f(file_id,dsetname, link_exists, h5err)
         if(h5err /= 0) then
            write(*,*) ' error checking if walave_u_x exists in restart file'
            call MPI_Abort(app_comm,-1,mpi_err)
         end if

         if(link_exists) then
            if(mpi_rank.eq.0) write(111,*) ' walave_u exists'
            call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,walave_u(:,1))
            ! I suposse that if walave_u_x exists then y and z also exist - no need to check again
            dsetname = 'walave_u_y'
            call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,walave_u(:,2))
            dsetname = 'walave_u_z'
            call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,walave_u(:,3))
         else
            if(mpi_rank.eq.0) write(111,*) ' walave_u does not exist - using u'
            walave_u(:,:) = u(:,:)
         end if
         !$acc update device(walave_u(:,:))
      end if

      if(flag_buffer) then
         dsetname = 'u_buffer_x'
         call h5lexists_f(file_id,dsetname, link_exists, h5err)
         if(h5err /= 0) then
            write(*,*) ' error checking if u_buffer_x exists in restart file'
            call MPI_Abort(app_comm,-1,mpi_err)
         end if

         if(link_exists) then
         dsetname = 'u_buffer_x'
         call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u_buffer(:,1))

         dsetname = 'u_buffer_y'
         call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u_buffer(:,2))

         dsetname = 'u_buffer_z'
         call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u_buffer(:,3))
         else
            if(mpi_rank.eq.0) write(111,*) ' u_buffer does not exist in restart file, setting to 0.'
            walave_u(:,:) = 0.0_rp
         end if
         !$acc update device(u_buffer(:,:))
      end if

      !$acc update device(rho(:))
      !$acc update device(pr(:))
      !$acc update device(E(:))
      !$acc update device(u(:,:))

      !$acc kernels
      do iElem = 1,numElemsRankPar
         do iNode = 1, mnnode
            mu_e(iElem,iNode) = aux_mu_e(connecParOrig(iElem,iNode))
            mu_t(iElem,iNode) = aux_mu_t(connecParOrig(iElem,iNode))
         end do
      end do
      !$acc end kernels

      !-------------------------------------------
      !2. If case is periodic, adjust slave nodes

      if (isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            rho(masSlaRankPar(iPer,2)) = rho(masSlaRankPar(iPer,1))
            u(masSlaRankPar(iPer,2),1) = u(masSlaRankPar(iPer,1),1)
            u(masSlaRankPar(iPer,2),2) = u(masSlaRankPar(iPer,1),2)
            u(masSlaRankPar(iPer,2),3) = u(masSlaRankPar(iPer,1),3)
            pr(masSlaRankPar(iPer,2))  = pr(masSlaRankPar(iPer,1))
            E(masSlaRankPar(iPer,2))   = E(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
         if(flag_walave == 1) then
            !$acc parallel loop
            do iPer = 1,nPerRankPar
               walave_u(masSlaRankPar(iPer,2),1) = walave_u(masSlaRankPar(iPer,1),1)
               walave_u(masSlaRankPar(iPer,2),2) = walave_u(masSlaRankPar(iPer,1),2)
               walave_u(masSlaRankPar(iPer,2),3) = walave_u(masSlaRankPar(iPer,1),3)
            end do
            !$acc end parallel loop
         end if
      end if

      call close_hdf5_file(file_id)

   end subroutine load_hdf5_restartFile

!----------------------------------------------------------------------------------------------------------------------------------
#if 0
   subroutine save_hdf5_resultsFile_old(iStep,time,rho,u,pr,E,eta,csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs)
      implicit none
      integer(4), intent(in) :: iStep
      real(rp),intent(in) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E,eta,csound,machno,mu_fluid,divU,Qcrit
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u
      real(rp),intent(inout),dimension(numElemsRankPar,mngaus) :: mu_e,mu_sgs
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: gradRho,curlU
      real(rp),dimension(numNodesRankPar) :: envit,mut

      integer(hid_t) :: file_id,plist_id
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4) :: h5err
      character(512) :: full_fileName,dsetname
      integer(4) :: iPer,i,j
      real(rp) :: aux_array_time(1)

      !------------------------------------------------------------------------------------
      !-------------------------------------------
      ! 1. Adjust mu_sgs and envit to be nodal
      !$acc kernels
      do i = 1,numElemsRankPar
         do j = 1, mnnode
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

      !-----------------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_resultsFile_name(iStep,full_fileName)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)
      !-----------------------------------------------------------------------------------------------
      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------
      dsetname = 'rho'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,3))

      dsetname = 'mu_fluid'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mu_fluid)

      dsetname = 'mut'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mut)

      dsetname = 'envit'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,envit)

      dsetname = 'gradRho_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,gradRho(:,1))

      dsetname = 'gradRho_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,gradRho(:,2))

      dsetname = 'gradRho_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,gradRho(:,3))

      dsetname = 'curlU_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,curlU(:,1))

      dsetname = 'curlU_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,curlU(:,2))

      dsetname = 'curlU_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,curlU(:,3))

      dsetname = 'pr'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,E)

      dsetname = 'eta'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,eta)

      dsetname = 'csound'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,csound)

      dsetname = 'machno'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,machno)

      dsetname = 'divU'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,divU)

      dsetname = 'Qcrit'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,Qcrit)

      ! ----  time  -----
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif
      aux_array_time = time

      dsetname = 'time'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_array_time)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine save_hdf5_resultsFile_old
#endif

   subroutine interpolate_scalarField_in_nodes(mnnode,mngaus,Ngp,connecParW,connecParO,origNodeScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode),origNodeScalarField(numNodesRankPar)
      real(rp),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,igp,inode
      real(rp) :: var_a

      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp(igp,inode)*origNodeScalarField(connecParW(iElem,inode))
            end do
            !$acc atomic write
            interpNodeScalarField(connecParO(iElem,igp)) = var_a
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop

   end subroutine interpolate_scalarField_in_nodes

   subroutine interpolate_vectorField_in_nodes(mnnode,mngaus,Ngp,connecParW,connecParO,origNodeVectorField,interpNodeVectorField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode),origNodeVectorField(numNodesRankPar,ndime)
      real(rp),intent(out) :: interpNodeVectorField(numNodesRankPar,ndime)
      integer(4) :: iElem,igp,inode,idime
      real(rp) :: var_a

      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector collapse(2)
         do igp = 1,mngaus
            do idime = 1,ndime
               var_a = 0.0_rp
               !$acc loop seq
               do inode = 1,mnnode
                  var_a = var_a+Ngp(igp,inode)*origNodeVectorField(connecParW(iElem,inode),idime)
               end do
               !$acc atomic write
               interpNodeVectorField(connecParO(iElem,igp),idime) = var_a
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop

   end subroutine interpolate_vectorField_in_nodes

   subroutine interpolate_elemGpScalarField_in_nodes(mnnode,mngaus,Ngp,connecParW,connecParO,origElemGpScalarField,interpNodeScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode),origElemGpScalarField(numElemsRankPar,mngaus)
      real(rp),intent(out) :: interpNodeScalarField(numNodesRankPar)
      integer(4) :: iElem,igp,inode,iPer
      real(rp) :: var_a

      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar

         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp(igp,inode)*origElemGpScalarField(iElem,inode)
            end do
            !$acc atomic write
            interpNodeScalarField(connecParO(iElem,igp)) = var_a
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop

      if(isMeshPeriodic) then
         !$acc parallel loop
         do iPer = 1,nPerRankPar
            interpNodeScalarField(masSlaRankPar(iPer,2)) = interpNodeScalarField(masSlaRankPar(iPer,1))
         end do
         !$acc end parallel loop
      end if

   end subroutine interpolate_elemGpScalarField_in_nodes

   subroutine interpolate_scalarField_in_elemGp(mnnode,mngaus,Ngp,connecParW,connecParO,origNodeScalarField,interpElemGpScalarField)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,connecParW(numElemsRankPar,mnnode),connecParO(numElemsRankPar,mnnode)
      real(rp),intent(in) :: Ngp(mngaus,mnnode),origNodeScalarField(numNodesRankPar)
      real(rp),intent(out) :: interpElemGpScalarField(numElemsRankPar,mngaus)
      integer(4) :: iElem,igp,inode,iPer
      real(rp) :: var_a

      !$acc parallel loop gang
      do iElem = 1,numElemsRankPar
         !$acc loop vector
         do igp = 1,mngaus
            var_a = 0.0_rp
            !$acc loop seq
            do inode = 1,mnnode
               var_a = var_a+Ngp(igp,inode)*origNodeScalarField(connecParW(iElem,inode))
            end do
            !$acc atomic write
            interpElemGpScalarField(iElem,igp) = var_a
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop

   end subroutine interpolate_scalarField_in_elemGp

   subroutine save_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp_equi,hdf5_fileId,numNodeScalarFields2save,nodeScalarFields2save,nameNodeScalarFields2save,&
                                             numNodeVectorFields2save,nodeVectorFields2save,nameNodeVectorFields2save,&
                                             numElemGpScalarFields2save,elemGpScalarFields2save,nameElemGpScalarFields2save)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus
      real(rp),intent(in) :: Ngp_equi(mngaus,mnnode)
      integer(hid_t),intent(in) :: hdf5_fileId
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      type(ptr_array1d_rp),intent(in) :: nodeScalarFields2save(:)
      type(ptr_array2d_rp),intent(in) :: nodeVectorFields2save(:),elemGpScalarFields2save(:)

      character(128),intent(in)   :: nameNodeScalarFields2save(numNodeScalarFields2save),nameNodeVectorFields2save(numNodeVectorFields2save),nameElemGpScalarFields2save(numElemGpScalarFields2save)

      integer(hsize_t) :: ds_dims(1),ms_dims(1),ds_dims2d(2),ms_dims2d(2)
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(4) :: h5err
      character(512) :: groupname,dsetname
      integer(4) :: iElem,iGp,iPer,iField

      !-----------------------------------------------------------------------------------------------
      !   Creating the VTK-HDF structure
      call create_vtkhdf_unstructuredGrid_struct_for_resultsFile(mnnode,hdf5_fileId)

      !-----------------------------------------------------------------------------------------------
      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------

      !if(mpi_rank.eq.0) write(*,*) 'saving hdf5 resultsfile istep',iStep,'time',time

      groupname = '/VTKHDF/PointData/'

      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numNodeScalarFields2save
         !if(mpi_rank.eq.0) write(*,*) 'saving field',iField,'name',dsetname

         call interpolate_scalarField_in_nodes(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,nodeScalarFields2save(iField)%ptr,auxInterpNodeScalarField)

         !$acc update host(auxInterpNodeScalarField(:))

         dsetname = trim(adjustl(groupname))//trim(nameNodeScalarFields2save(iField))
         call save_array1D_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims,ms_dims,ms_offset,auxInterpNodeScalarField)

      end do

      !--------------------------------------------------------------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------------------
      ds_dims2d(1) = int(ndime,hsize_t)
      ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------

      do iField=1,numNodeVectorFields2save

         call interpolate_vectorField_in_nodes(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,nodeVectorFields2save(iField)%ptr,auxInterpNodeVectorField)

         !$acc update host(auxInterpNodeVectorField(:,:))

         dsetname = trim(adjustl(groupname))//trim(nameNodeVectorFields2save(iField))
         call save_array2D_tr_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxInterpNodeVectorField)
      end do

      do iField=1,numElemGpScalarFields2save
         !if(mpi_rank.eq.0) write(*,*) 'saving field',iField,'name',dsetname

         call interpolate_elemGpScalarField_in_nodes(mnnode,mngaus,Ngp_equi,connecParWork,connecParOrig,elemGpScalarFields2save(iField)%ptr,auxInterpNodeScalarField)

         !$acc update host(auxInterpNodeScalarField(:))

         dsetname = trim(adjustl(groupname))//trim(nameElemGpScalarFields2save(iField))
         call save_array1D_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims,ms_dims,ms_offset,auxInterpNodeScalarField)

      end do

   end subroutine save_hdf5_resultsFile_baseFunc

   subroutine load_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,numNodeScalarFields2load,nodeScalarFields2load,nameNodeScalarFields2load,&
                                             numNodeVectorFields2load,nodeVectorFields2load,nameNodeVectorFields2load,&
                                             numElemGpScalarFields2load,elemGpScalarFields2load,nameElemGpScalarFields2load)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      integer(hid_t),intent(in) :: hdf5_fileId
      integer(4),intent(in) :: numNodeScalarFields2load,numNodeVectorFields2load,numElemGpScalarFields2load
      type(ptr_array1d_rp),intent(inout) :: nodeScalarFields2load(:)
      type(ptr_array2d_rp),intent(inout) :: nodeVectorFields2load(:),elemGpScalarFields2load(:)

      character(128),intent(in)   :: nameNodeScalarFields2load(numNodeScalarFields2load),nameNodeVectorFields2load(numNodeVectorFields2load),nameElemGpScalarFields2load(numElemGpScalarFields2load)

      integer(hsize_t) :: ds_dims(1),ms_dims(1),ds_dims2d(2),ms_dims2d(2)
      integer(hssize_t) :: ms_offset(1),ms_offset2d(2)
      integer(4) :: h5err
      character(512) :: groupname,dsetname
      integer(4) :: iElem,iGp,iPer,iField
      real(rp) :: aux_nodeScalarField(numNodesRankPar)
      real(rp) :: aux_array_time(1)

      !-----------------------------------------------------------------------------------------------
      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------

      groupname = '/VTKHDF/PointData/'

      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numNodeScalarFields2load
         dsetname = trim(adjustl(groupname))//trim(nameNodeScalarFields2load(iField))

         call read_array1D_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ms_dims,ms_offset,auxInterpNodeScalarField)
         !$acc update device(auxInterpNodeScalarField(:))

         call interpolate_scalarField_in_nodes(mnnode,mngaus,Ngp,connecParWork,connecParOrig,auxInterpNodeScalarField,nodeScalarFields2load(iField)%ptr)

      end do

      !--------------------------------------------------------------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------------------
      ds_dims2d(1) = int(ndime,hsize_t)
      ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------

      do iField=1,numNodeVectorFields2load
         dsetname = trim(adjustl(groupname))//trim(nameNodeVectorFields2load(iField))

         call read_array2D_tr_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ms_dims2d,ms_offset2d,auxInterpNodeVectorField)
         !$acc update device(auxInterpNodeVectorField(:,:))

         call interpolate_vectorField_in_nodes(mnnode,mngaus,Ngp,connecParWork,connecParOrig,auxInterpNodeVectorField,nodeVectorFields2load(iField)%ptr)

      end do

      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numElemGpScalarFields2load
         dsetname = trim(adjustl(groupname))//trim(nameElemGpScalarFields2load(iField))
         !if(mpi_rank.eq.0) write(*,*) 'saving field',iField,'name',dsetname

         call read_array1D_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ms_dims,ms_offset,auxInterpNodeScalarField)
         !$acc update device(auxInterpNodeScalarField(:))

         call interpolate_scalarField_in_elemGp(mnnode,mngaus,Ngp,connecParWork,connecParOrig,auxInterpNodeScalarField,elemGpScalarFields2load(iField)%ptr)

      end do

   end subroutine load_hdf5_resultsFile_baseFunc

   subroutine save_instResults_hdf5_file(mnnode,mngaus,Ngp,iStep,time,numNodeScalarFields2save,nodeScalarFields2save,nameNodeScalarFields2save,&
                                          numNodeVectorFields2save,nodeVectorFields2save,nameNodeVectorFields2save,&
                                          numElemGpScalarFields2save,elemGpScalarFields2save,nameElemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: mnnode,mngaus,iStep
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(in) :: time
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      type(ptr_array1d_rp),intent(in) :: nodeScalarFields2save(:)
      type(ptr_array2d_rp),intent(in) :: nodeVectorFields2save(:),elemGpScalarFields2save(:)
      character(128),intent(in)   :: nameNodeScalarFields2save(numNodeScalarFields2save),nameNodeVectorFields2save(numNodeVectorFields2save),nameElemGpScalarFields2save(numElemGpScalarFields2save)

      integer(hid_t) :: hdf5_fileId
      character(512) :: full_hdf5_fileName,dsetname

      !-----------------------------------------------------------------------------------------------

      call set_hdf5_resultsFile_name(iStep,full_hdf5_fileName)

      call create_hdf5_file(full_hdf5_fileName,hdf5_fileId)

      call save_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,numNodeScalarFields2save,nodeScalarFields2save,nameNodeScalarFields2save,&
                                          numNodeVectorFields2save,nodeVectorFields2save,nameNodeVectorFields2save,&
                                          numElemGpScalarFields2save,elemGpScalarFields2save,nameElemGpScalarFields2save)

      dsetname = 'time'
      call save_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,time)

      !----------------------------------------------------------------------------------------------------------

      call close_hdf5_file(hdf5_fileId)

   end subroutine save_instResults_hdf5_file

   subroutine save_avgResults_hdf5_file(mnnode,mngaus,Ngp,restartCnt,initial_avgTime,elapsed_avgTime,numAvgNodeScalarFields2save,avgNodeScalarFields2save,nameAvgNodeScalarFields2save,&
                                       numAvgNodeVectorFields2save,avgNodeVectorFields2save,nameAvgNodeVectorFields2save,&
                                       numAvgElemGpScalarFields2save,avgElemGpScalarFields2save,nameAvgElemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: mnnode,mngaus,restartCnt
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(in) :: initial_avgTime,elapsed_avgTime
      integer(4),intent(in) :: numAvgNodeScalarFields2save,numAvgNodeVectorFields2save,numAvgElemGpScalarFields2save
      type(ptr_array1d_rp),intent(in) :: avgNodeScalarFields2save(:)
      type(ptr_array2d_rp),intent(in) :: avgNodeVectorFields2save(:),avgElemGpScalarFields2save(:)
      character(128),intent(in)   :: nameAvgNodeScalarFields2save(numAvgNodeScalarFields2save),nameAvgNodeVectorFields2save(numAvgNodeVectorFields2save),nameAvgElemGpScalarFields2save(numAvgElemGpScalarFields2save)

      integer(hid_t) :: hdf5_fileId
      character(512) :: full_hdf5_fileName,dsetname
      integer(hsize_t),dimension(1) :: ds_dims,ms_dims
      integer(hssize_t),dimension(1) :: ms_offset
      real(rp) :: aux_array_rp(1)
      !-----------------------------------------------------------------------------------------------
      call set_hdf5_avgResultsFile_name(restartCnt,full_hdf5_fileName)

      call create_hdf5_file(full_hdf5_fileName,hdf5_fileId)

      call save_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,numAvgNodeScalarFields2save,avgNodeScalarFields2save,nameAvgNodeScalarFields2save,&
                                          numAvgNodeVectorFields2save,avgNodeVectorFields2save,nameAvgNodeVectorFields2save,&
                                          numAvgElemGpScalarFields2save,avgElemGpScalarFields2save,nameAvgElemGpScalarFields2save)
      dsetname = 'elapsed_avgTime'
      call save_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,elapsed_avgTime)

      dsetname = 'initial_avgTime'
      call save_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,initial_avgTime)

      !-----------------------------------------------------------------------------------------------

      call close_hdf5_file(hdf5_fileId)

   end subroutine save_avgResults_hdf5_file

   subroutine load_avgResults_hdf5_file(mnnode,mngaus,Ngp,restartCnt,initial_avgTime,elapsed_avgTime,numAvgNodeScalarFields2load,avgNodeScalarFields2load,nameAvgNodeScalarFields2load,&
                                       numAvgNodeVectorFields2load,avgNodeVectorFields2load,nameAvgNodeVectorFields2load,&
                                       numAvgElemGpScalarFields2load,avgElemGpScalarFields2load,nameAvgElemGpScalarFields2load)
      implicit none
      integer(4),intent(in) :: mnnode,mngaus,restartCnt
      real(rp),intent(in) :: Ngp(mngaus,mnnode)
      real(rp),intent(inout) :: initial_avgTime,elapsed_avgTime
      integer(4),intent(in) :: numAvgNodeScalarFields2load,numAvgNodeVectorFields2load,numAvgElemGpScalarFields2load
      type(ptr_array1d_rp),intent(inout) :: avgNodeScalarFields2load(:)
      type(ptr_array2d_rp),intent(inout) :: avgNodeVectorFields2load(:),avgElemGpScalarFields2load(:)
      character(128),intent(in)   :: nameAvgNodeScalarFields2load(numAvgNodeScalarFields2load),nameAvgNodeVectorFields2load(numAvgNodeVectorFields2load),nameAvgElemGpScalarFields2load(numAvgElemGpScalarFields2load)

      integer(hid_t) :: hdf5_fileId
      character(512) :: full_hdf5_fileName,dsetname

      call set_hdf5_avgResultsFile_name(restartCnt,full_hdf5_fileName)

      call open_hdf5_file(full_hdf5_fileName,hdf5_fileId)

      call load_hdf5_resultsFile_baseFunc(mnnode,mngaus,Ngp,hdf5_fileId,numAvgNodeScalarFields2load,avgNodeScalarFields2load,nameAvgNodeScalarFields2load,&
                                          numAvgNodeVectorFields2load,avgNodeVectorFields2load,nameAvgNodeVectorFields2load,&
                                          numAvgElemGpScalarFields2load,avgElemGpScalarFields2load,nameAvgElemGpScalarFields2load)

      !-----------------------------------------------------------------------------------------------

      dsetname = 'elapsed_avgTime'
      call read_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,elapsed_avgTime)

      dsetname = 'initial_avgTime'
      call read_real_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,initial_avgTime)

      !-----------------------------------------------------------------------------------------------

      call close_hdf5_file(hdf5_fileId)

   end subroutine load_avgResults_hdf5_file

   !-------------------------------------------------------------------------------------------------------------------------------
   !de moment la faig generica, ho enfoc a result instantani, despres veure si em funciona i ho puc estendre al avg
   !cool

   subroutine save_surface_instResults_hdf5_file(iStep,numNodeScalarFields2save,nameNodeScalarFields2save,&
                                                  numNodeVectorFields2save,nameNodeVectorFields2save,&
                                                  numElemGpScalarFields2save,nameElemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: iStep
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      character(128),intent(in) :: nameNodeScalarFields2save(numNodeScalarFields2save),nameNodeVectorFields2save(numNodeVectorFields2save),nameElemGpScalarFields2save(numElemGpScalarFields2save)
      character(512) :: res_hdf5_fileName

      call set_hdf5_resultsFile_name(iStep,res_hdf5_fileName)

      call save_surface_results_hdf5_file(res_hdf5_fileName,numNodeScalarFields2save,nameNodeScalarFields2save,&
                                             numNodeVectorFields2save,nameNodeVectorFields2save,&
                                             numElemGpScalarFields2save,nameElemGpScalarFields2save)

   end subroutine save_surface_instResults_hdf5_file

   subroutine save_surface_avgResults_hdf5_file(restartCnt,numAvgNodeScalarFields2save,nameAvgNodeScalarFields2save,&
                                                   numAvgNodeVectorFields2save,nameAvgNodeVectorFields2save,&
                                                   numAvgElemGpScalarFields2save,nameAvgElemGpScalarFields2save)
      implicit none
      integer(4), intent(in) :: restartCnt
      integer(4),intent(in) :: numAvgNodeScalarFields2save,numAvgNodeVectorFields2save,numAvgElemGpScalarFields2save
      character(128),intent(in) :: nameAvgNodeScalarFields2save(numAvgNodeScalarFields2save),nameAvgNodeVectorFields2save(numAvgNodeVectorFields2save),nameAvgElemGpScalarFields2save(numAvgElemGpScalarFields2save)
      character(512) :: res_hdf5_fileName

      call set_hdf5_avgResultsFile_name(restartCnt,res_hdf5_fileName)

      call save_surface_results_hdf5_file(res_hdf5_fileName,numAvgNodeScalarFields2save,nameAvgNodeScalarFields2save,&
                                             numAvgNodeVectorFields2save,nameAvgNodeVectorFields2save,&
                                             numAvgElemGpScalarFields2save,nameAvgElemGpScalarFields2save)

   end subroutine save_surface_avgResults_hdf5_file

   subroutine save_surface_results_hdf5_file(res_hdf5_fileName,numNodeScalarFields2save,nameNodeScalarFields2save,&
                                             numNodeVectorFields2save,nameNodeVectorFields2save,&
                                             numElemGpScalarFields2save,nameElemGpScalarFields2save)
      implicit none
      character(512),intent(in) :: res_hdf5_fileName
      integer(4),intent(in) :: numNodeScalarFields2save,numNodeVectorFields2save,numElemGpScalarFields2save
      character(128),intent(in) :: nameNodeScalarFields2save(numNodeScalarFields2save),nameNodeVectorFields2save(numNodeVectorFields2save),nameElemGpScalarFields2save(numElemGpScalarFields2save)
      integer(4) :: h5err

      integer(hid_t) :: hdf5_fileId
      character(512) :: surf_res_hdf5_fileName,groupname,dsetname
      integer(4) :: iField

      !--------------------------------------------------------------------------------------------------------------------------------------
      call set_hdf5_surface_resultsFile_name(surf_res_hdf5_fileName,res_hdf5_fileName)
      call create_hdf5_file(surf_res_hdf5_fileName,hdf5_fileId)
      call set_vtkhdf_attributes_and_basic_groups(hdf5_fileId)
      !--------------------------------------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------------------------------------
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Points',hdf5_fileId,'/VTKHDF/Points',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/NumberOfPoints',hdf5_fileId,'/VTKHDF/NumberOfPoints',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/NumberOfCells',hdf5_fileId,'/VTKHDF/NumberOfCells',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/NumberOfConnectivityIds',hdf5_fileId,'/VTKHDF/NumberOfConnectivityIds',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/Offsets',hdf5_fileId,'/VTKHDF/Offsets',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/Connectivity',hdf5_fileId,'/VTKHDF/Connectivity',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/Types',hdf5_fileId,'/VTKHDF/Types',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/CellData/mpi_rank',hdf5_fileId,'/VTKHDF/CellData/mpi_rank',h5err)
      call h5lcreate_external_f(surface_meshFile_h5_name,'/VTKHDF/CellData/boundCode',hdf5_fileId,'/VTKHDF/CellData/boundCode',h5err)
      !--------------------------------------------------------------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------------------------------------------------------
      groupname = '/VTKHDF/PointData/'
      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numNodeScalarFields2save
         dsetname = trim(adjustl(groupname))//trim(nameNodeScalarFields2save(iField))
         call h5lcreate_external_f(res_hdf5_fileName,dsetname,hdf5_fileId,dsetname,h5err)
      end do
      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numNodeVectorFields2save
         dsetname = trim(adjustl(groupname))//trim(nameNodeVectorFields2save(iField))
         call h5lcreate_external_f(res_hdf5_fileName,dsetname,hdf5_fileId,dsetname,h5err)
      end do
      !--------------------------------------------------------------------------------------------------------------------------------------
      do iField=1,numElemGpScalarFields2save
         dsetname = trim(adjustl(groupname))//trim(nameElemGpScalarFields2save(iField))
         call h5lcreate_external_f(res_hdf5_fileName,dsetname,hdf5_fileId,dsetname,h5err)
      end do
      !--------------------------------------------------------------------------------------------------------------------------------------

      call close_hdf5_file(hdf5_fileId)

   end subroutine save_surface_results_hdf5_file
   !-------------------------------------------------------------------------------------------------------------------------------

   subroutine save_surface_mesh_hdf5_file(mnpbou,gmsh2ij,vtk2ij)
      implicit none
      integer(4),intent(in) :: mnpbou,gmsh2ij(mnpbou),vtk2ij(mnpbou)
      integer(4) :: ds_rank,h5err
      integer(hsize_t),dimension(1) :: ds_dims,ms_dims
      integer(hssize_t),dimension(1) :: ms_offset

      integer(hid_t) :: hdf5_fileId
      character(512) :: groupname,dsetname
      integer(hid_t) :: dtype
      integer(1),allocatable :: aux_array_i1(:)
      integer(8),allocatable :: aux_array_i8(:)
      integer(4), dimension(0:mpi_size-1) :: vecNumBoundsRankPar
      integer(4) :: mpiRankBoundStart
      integer(4) :: ii,jj,iBound,iRank,iNodeL

      !---------------------------------------------------------------------------------
      call set_hdf5_surface_meshFile_name()
      call create_hdf5_file(surface_meshFile_h5_name,hdf5_fileId)
      call set_vtkhdf_attributes_and_basic_groups(hdf5_fileId)
      !---------------------------------------------------------------------------------

      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Points',hdf5_fileId,'/VTKHDF/Points',h5err)

      !-----------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = int(mpi_size,hsize_t)
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8
      allocate(aux_array_i8(1))

      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesRankPar
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numBoundsRankPar
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numBoundsRankPar*mnpbou
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

      deallocate(aux_array_i8)

      !-----------------------------------------------------------------------------
      call MPI_Allgather(numBoundsRankPar,1,mpi_datatype_int4,vecNumBoundsRankPar,1,mpi_datatype_int4,app_comm,mpi_err)

      mpiRankBoundStart = 1
      do iRank=1,mpi_rank
         mpiRankBoundStart = mpiRankBoundStart +vecNumBoundsRankPar(iRank-1)
      end do
      !------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------

      allocate(aux_array_i8(numBoundsRankPar+1))

      ms_dims(1)   = int(numBoundsRankPar,hsize_t)  + 1
      ms_offset(1) = int(mpiRankBoundStart,hssize_t) - 1 + int(mpi_rank,hssize_t)

      dsetname = '/VTKHDF/Offsets'

      ds_rank    = 1
      ds_dims(1) = int(totalNumBoundsSrl,hsize_t) + int(mpi_size,hsize_t)
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)

      aux_array_i8(1) = 0
      do iBound = 2,(numBoundsRankPar+1)
         aux_array_i8(iBound) = aux_array_i8(iBound-1)+mnpbou
      end do

      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)

      !-----------------------------------------------------------------------------

      allocate(aux_array_i8(numBoundsRankPar*mnpbou))

      dsetname   = '/VTKHDF/Connectivity'

      ds_rank    = 1
      ds_dims(1) = int(totalNumBoundsSrl,hsize_t)  * int(mnpbou,hsize_t)
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)

      ms_dims(1)   = int(numBoundsRankPar,hsize_t)   * int(mnpbou,hsize_t)
      ms_offset(1) = int((mpiRankBoundStart-1),hssize_t)* int(mnpbou,hssize_t)

      !write(*,*) 'EI! THIS FUNC save_surface_mesh_hdf5_file() in mod_hdf5.f90 NOW DOES NOT WORK! FIX IT!'
      !call MPI_Abort(app_comm,-1,mpi_err)

      do iBound=1,numBoundsRankPar
         do ii=1,mnpbou

            iNodeL = boundParOrig(iBound,gmsh2ij(ii))
            jj = (iBound-1)*mnpbou + vtk2ij(ii)
            aux_array_i8(jj) = iNodeL - 1
         end do
      end do

      call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)

      !-----------------------------------------------------------------------------
      allocate(aux_array_i1(numBoundsRankPar))
      dsetname = '/VTKHDF/Types'
      ds_rank    = 1
      ds_dims(1) = int(totalNumBoundsSrl,hsize_t)
      dtype      = h5_datatype_uint1
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)

      aux_array_i1(:) = 70

      ms_dims(1)   = int(numBoundsRankPar ,hsize_t)
      ms_offset(1) = int(mpiRankBoundStart,hssize_t) - 1

      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)

      !------------------------------------------------------------------------------------------------------
      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)

      do iBound = 1,numBoundsRankPar
         aux_array_i1(iBound) = mpi_rank
      end do
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)

      ! ## bou_codes ##
      dsetname = '/VTKHDF/CellData/boundCode'
      call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)

      do iBound = 1,numBoundsRankPar
         aux_array_i1(iBound) = bouCodesPar(iBound)
      end do
      call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      call close_hdf5_file(hdf5_fileId)

   end subroutine save_surface_mesh_hdf5_file

#if 0
   subroutine save_hdf5_avgResultsFile_old(iStep,avvel,avve2,avvex,avrho,avpre,avmueff,avtw)
      implicit none
      integer,intent(in) :: iStep
      real(rp),intent(in),dimension(numNodesRankPar)       :: avrho,avpre,avmueff
      real(rp),intent(in),dimension(numNodesRankPar,ndime) :: avvel,avve2,avvex,avtw

      integer(hid_t) :: file_id,plist_id,dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(4) :: h5err
      character(512) :: full_fileName,dsetname
      real(rp) :: aux_array(1)

      !------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_avgResultsFile_name(iStep,full_fileName)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create avg results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'avrho'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avrho)

      dsetname = 'avpre'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avpre)

      dsetname = 'avmueff'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avmueff)

      dsetname = 'avvel_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvel(:,1))

      dsetname = 'avvel_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvel(:,2))

      dsetname = 'avvel_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvel(:,3))

      dsetname = 'avve2_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avve2(:,1))

      dsetname = 'avve2_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avve2(:,2))

      dsetname = 'avve2_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avve2(:,3))

      dsetname = 'avvex_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvex(:,1))

      dsetname = 'avvex_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvex(:,2))

      dsetname = 'avvex_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvex(:,3))

      dsetname = 'avtw_x'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avtw(:,1))

      dsetname = 'avtw_y'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avtw(:,2))

      dsetname = 'avtw_z'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avtw(:,3))

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine save_hdf5_avgResultsFile_old

   subroutine load_hdf5_avgResultsFile(iStep,avvel,avve2,avvex,avrho,avpre,avmueff,avtw)
      implicit none
      integer, intent(in) :: iStep
      real(rp),intent(inout),dimension(numNodesRankPar)       :: avrho,avpre,avmueff
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: avvel,avve2,avvex,avtw

      integer(hid_t) :: file_id,plist_id,dtype
      integer(HSIZE_T), dimension(1) :: ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer :: ms_rank,h5err
      character(512) :: full_fileName,dsetname
      real(rp) :: aux_array(1)

      !------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_avgResultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Loading results file: ',trim(adjustl(full_fileName))

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results avg file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'avrho'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avrho)

      dsetname = 'avpre'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avpre)

      dsetname = 'avmueff'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avmueff)

      dsetname = 'avvel_x'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvel(:,1))

      dsetname = 'avvel_y'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvel(:,2))

      dsetname = 'avvel_z'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvel(:,3))

      dsetname = 'avve2_x'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avve2(:,1))

      dsetname = 'avve2_y'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avve2(:,2))

      dsetname = 'avve2_z'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avve2(:,3))

      dsetname = 'avvex_x'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvex(:,1))

      dsetname = 'avvex_y'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvex(:,2))

      dsetname = 'avvex_z'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvex(:,3))

      dsetname = 'avtw_x'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avtw(:,1))

      dsetname = 'avtw_y'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avtw(:,2))

      dsetname = 'avtw_z'
      call read_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avtw(:,3))

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine load_hdf5_avgResultsFile

   subroutine close_hdf5_resultsfile()
      implicit none
      integer :: h5err

      !close the file.
      call h5fclose_f(resultsFile_h5_id,h5err)

   end subroutine close_hdf5_resultsfile
#endif
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

   subroutine set_vtkhdf_attributes_and_basic_groups(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(512) :: groupname
      integer(hid_t) :: group_id,aspace_id,attr_id,atype
      integer(size_t) :: attr_length
      integer(hsize_t) :: a_dims(1)
      integer(4) :: a_rank,h5err
      character(16) :: attr_value
      integer(4),allocatable :: aux_array_i4(:)

      !--------------------------------------------------------------------------------
      groupname = '/VTKHDF'
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

      call h5awrite_f(attr_id,atype,attr_value,a_dims,h5err)
      call h5tclose_f(atype,h5err)

      call h5aclose_f(attr_id, h5err)

      call h5sclose_f(aspace_id, h5err)
      !--------------------------------------------------------------------------------
      ! Attribute 'Version'
      a_rank = 1
      a_dims(1) = 2
      call h5screate_simple_f(a_rank,a_dims,aspace_id,h5err)

      call h5acreate_f(group_id,'Version',h5_datatype_int4,aspace_id,attr_id,h5err)

      allocate(aux_array_i4(2))
      aux_array_i4(1) = 1
      aux_array_i4(2) = 0
      call h5awrite_f(attr_id,h5_datatype_int4,aux_array_i4,a_dims,h5err)
      deallocate(aux_array_i4)

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
   end subroutine

   subroutine create_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t),intent(in) :: file_id
      integer(hid_t) :: dtype
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hsize_t), dimension(2) :: ds_dims2d,ms_dims2d
      integer(hssize_t), dimension(1) :: ms_offset
      integer(hssize_t), dimension(2) :: ms_offset2d
      integer(4) :: ds_rank,h5err
      character(512) :: dsetname
      integer(4) :: ii,iElemL
      integer(1),allocatable :: aux_array_i1(:)
      integer(8),allocatable :: aux_array_i8(:)

      !--------------------------------------------------------------------------------

      call set_vtkhdf_attributes_and_basic_groups(file_id)

      !--------------------------------------------------------------------------------
      ds_dims2d(1) = int(ndime,hsize_t)
      ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      ms_dims2d(1) = int(ndime,hsize_t)
      ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      ms_offset2d(1) = 0
      ms_offset2d(2) = int(rankNodeStart,hssize_t)-1

!--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/Points'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,coordPar)

      !-----------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8
      allocate(aux_array_i8(1))

      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numElemsRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numElemsRankPar*mnnode
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      !-----------------------------------------------------------------------------
      deallocate(aux_array_i8)
      allocate(aux_array_i8(numElemsRankPar+1))

      ds_dims(1)   = int(totalNumElements,hsize_t) + int(mpi_size,hsize_t)
      ms_dims(1)   = int(numElemsRankPar,hsize_t)  + 1
      ms_offset(1) = int(rankElemStart,hssize_t)   - 1 + int(mpi_rank,hssize_t)

      dsetname = '/VTKHDF/Offsets'
      aux_array_i8(1) = 0
      do iElemL = 2,(numElemsRankPar+1)
         aux_array_i8(iElemL) = aux_array_i8(iElemL-1)+mnnode
      end do

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)
      !-----------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      allocate(aux_array_i8(numElemsRankPar*mnnode))

      ds_dims(1)   = int(totalNumElements,hsize_t)  * int(mnnode,hsize_t)
      ms_dims(1)   = int(numElemsRankPar,hsize_t)   * int(mnnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)* int(mnnode,hssize_t)

      do ii = 1,numElemsRankPar*mnnode
         aux_array_i8(ii) = connecVTK(ii)-1
      end do

      dsetname = '/VTKHDF/Connectivity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)
      !-----------------------------------------------------------------------------
      allocate(aux_array_i1(numElemsRankPar))
      dtype = h5_datatype_uint1
      dsetname = '/VTKHDF/Types'
      !$acc kernels
      aux_array_i1(:) = 72
      !$acc end kernels

      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t) - 1

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

   end subroutine create_vtkhdf_unstructuredGrid_meshFile

!-------------------------------------------------------------------------------------------------------------------------------------------

   subroutine create_groups_datasets_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode,numMshRanks2Part,numElemsGmsh
      integer(8),intent(in) :: numNodesParTotal_i8
      integer(hid_t) :: dtype
      integer(hsize_t) :: ds_dims(1),ds_dims2d(2)
      integer(4) :: ds_rank,h5err
      character(512) :: dsetname
      !--------------------------------------------------------------------------------

      call set_vtkhdf_attributes_and_basic_groups(file_id)

      !--------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------
      call select_dtype_rp(dtype)
      ds_rank = 2
      ds_dims2d(1) = ndime
      ds_dims2d(2) = numNodesParTotal_i8

      dsetname = '/VTKHDF/Points'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)

      !-----------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numMshRanks2Part
      dtype = h5_datatype_int8

      dsetname = '/VTKHDF/NumberOfPoints'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/VTKHDF/NumberOfCells'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !-----------------------------------------------------------------------------
      ds_dims(1)   = int(numElemsGmsh,hsize_t) + int(numMshRanks2Part,hsize_t)

      dsetname = '/VTKHDF/Offsets'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !-----------------------------------------------------------------------------
      ds_dims(1)   = int(numElemsGmsh,hsize_t)  * int(mnnode,hsize_t)

      dsetname = '/VTKHDF/Connectivity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !-----------------------------------------------------------------------------
      ds_dims(1)   = int(numElemsGmsh,hsize_t)
      dtype = h5_datatype_uint1

      dsetname = '/VTKHDF/Types'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

      !-----------------------------------------------------------------------------
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

   end subroutine create_groups_datasets_vtkhdf_unstructuredGrid_meshFile

   subroutine write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id,mshRank,numMshRanks2Part,numElemsMshRank,mshRankElemStart,mshRankElemEnd,mshRankNodeStart_i8,mshRankNodeEnd_i8,numNodesMshRank,coordVTKMshRank,connecVTKMshRank)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: mnnode,mshRank,numMshRanks2Part
      integer(4),intent(in) :: numElemsMshRank,mshRankElemStart,mshRankElemEnd
      integer(8),intent(in) :: mshRankNodeStart_i8,mshRankNodeEnd_i8
      integer(4),intent(in) :: numNodesMshRank
      real(rp),intent(in)   :: coordVTKMshRank(numNodesMshRank,3)
      integer(4),intent(in) :: connecVTKMshRank(numElemsMshRank*mnnode)

      integer(hsize_t), dimension(1) :: ms_dims
      integer(hsize_t), dimension(2) :: ms_dims2d
      integer(hssize_t), dimension(1) :: ms_offset
      integer(hssize_t), dimension(2) :: ms_offset2d
      character(512) :: dsetname
      integer(4) :: ii,iElemL
      integer(1),allocatable :: aux_array_i1(:)
      integer(8),allocatable :: aux_array_i8(:)

      !--------------------------------------------------------------------------------
      ms_dims2d(1) = ndime
      ms_dims2d(2) = numNodesMshRank
      ms_offset2d(1) = 0
      ms_offset2d(2) = mshRankNodeStart_i8-1
      !--------------------------------------------------------------------------------

      dsetname = '/VTKHDF/Points'
      call write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,coordVTKMshRank)
      !--------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------
      ms_dims(1) = 1
      ms_offset(1) = int(mshRank,hssize_t)

      allocate(aux_array_i8(1))

      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesMshRank
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numElemsMshRank
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numElemsMshRank*mnnode
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      !-----------------------------------------------------------------------------
      deallocate(aux_array_i8)
      allocate(aux_array_i8(numElemsMshRank+1))

      ms_dims(1)   = int(numElemsMshRank,hsize_t)  + 1
      ms_offset(1) = int(mshRankElemStart,hssize_t)   - 1 + int(mshRank,hssize_t)

      dsetname = '/VTKHDF/Offsets'
      aux_array_i8(1) = 0
      do iElemL = 2,(numElemsMshRank+1)
         aux_array_i8(iElemL) = aux_array_i8(iElemL-1)+mnnode
      end do

      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)
      !-----------------------------------------------------------------------------
      allocate(aux_array_i8(numElemsMshRank*mnnode))

      ms_dims(1)   = int(numElemsMshRank,hsize_t)   * int(mnnode,hsize_t)
      ms_offset(1) = int((mshRankElemStart-1),hssize_t)* int(mnnode,hssize_t)

      do ii = 1,numElemsMshRank*mnnode
         aux_array_i8(ii) = connecVTKMshRank(ii)-1
      end do

      dsetname = '/VTKHDF/Connectivity'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)

      !-----------------------------------------------------------------------------
      allocate(aux_array_i1(numElemsMshRank))
      dsetname = '/VTKHDF/Types'
      !$acc kernels
      aux_array_i1(:) = 72
      !$acc end kernels

      ms_dims(1)   = int(numElemsMshRank ,hsize_t)
      ms_offset(1) = int(mshRankElemStart,hssize_t) - 1

      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)

      !------------------------------------------------------------------------------------------------------
      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      do iElemL = 1,numElemsMshRank
         aux_array_i1(iElemL) = mshRank
      end do
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

   end subroutine write_mshRank_data_vtkhdf_unstructuredGrid_meshFile

   subroutine dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(file_id,numMshRanks2Part)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer(4),intent(in) :: numMshRanks2Part

      integer(hsize_t), dimension(1) :: ms_dims
      integer(hsize_t), dimension(2) :: ms_dims2d
      integer(hssize_t), dimension(1) :: ms_offset
      integer(hssize_t), dimension(2) :: ms_offset2d
      character(512) :: dsetname
      integer(1),allocatable :: empty_array_i1(:)
      integer(8),allocatable :: empty_array_i8(:)
      real(rp),allocatable :: empty_array2d_rp(:,:)

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
      allocate(empty_array_i1(0))
      allocate(empty_array_i8(0))
      allocate(empty_array2d_rp(0,0))

      !--------------------------------------------------------------------------------
      ms_dims2d(1) = 0
      ms_dims2d(2) = 0
      ms_offset2d(1) = 0
      ms_offset2d(2) = 0
      !--------------------------------------------------------------------------------

      dsetname = '/VTKHDF/Points'
      call write_dataspace_2d_tr_real_rp_hyperslab_parallel(file_id,dsetname,ms_dims2d,ms_offset2d,empty_array2d_rp)

      !-----------------------------------------------------------------------------
      ms_dims(1) = 0
      ms_offset(1) = 0

      dsetname = '/VTKHDF/NumberOfPoints'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      !---------------------------------------------------------------------
      dsetname = '/VTKHDF/Offsets'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      !-----------------------------------------------------------------------------
      dsetname = '/VTKHDF/Connectivity'
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i8)

      !-----------------------------------------------------------------------------
      dsetname = '/VTKHDF/Types'
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i1)

      !------------------------------------------------------------------------------------------------------
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,empty_array_i1)

      deallocate(empty_array_i1)
      deallocate(empty_array_i8)
      deallocate(empty_array2d_rp)

   end subroutine dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile

   subroutine create_vtkhdf_unstructuredGrid_struct_for_resultsFile(mnnode,file_id)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t),intent(in) :: file_id
      integer(hid_t) :: dtype
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hssize_t), dimension(1) :: ms_offset
      integer(4) :: ds_rank,h5err
      character(512) :: dsetname
      integer(8),allocatable :: aux_array_i8(:)

      call set_vtkhdf_attributes_and_basic_groups(file_id)

      !--------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------
      !dtype = h5_datatype_real8
      !ds_dims2d(1) = 3
      !ms_dims2d(1) = 1
      !ds_dims2d(2) = totalNumNodesPar
      !ms_dims2d(2) = numNodesRankPar
      !ms_offset2d(1) = 0
      !ms_offset2d(2) = rankNodeStart-1
!--------------------------------------------------------------------------------
      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Points',file_id,'/VTKHDF/Points',h5err)

      !-----------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = int(mpi_size,hsize_t)
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8
      allocate(aux_array_i8(1))

      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numElemsRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numElemsRankPar*mnnode
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int8_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

      deallocate(aux_array_i8)

      !-----------------------------------------------------------------------------

      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Offsets',file_id,'/VTKHDF/Offsets',h5err)

      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Connectivity',file_id,'/VTKHDF/Connectivity',h5err)

      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/Types',file_id,'/VTKHDF/Types',h5err)

      call h5lcreate_external_f(meshFile_h5_name,'/VTKHDF/CellData/mpi_rank',file_id,'/VTKHDF/CellData/mpi_rank',h5err)

   end subroutine create_vtkhdf_unstructuredGrid_struct_for_resultsFile

!-------------------------------------------------------------------------------------------------------------------------------------------

   subroutine write_vtkhdf_instResultsFile(mnnode,full_fileName,rho,pr,E,eta,csound,machno,divU,Qcrit,envit,mut,mu_fluid,u,gradRho,curlU)
      implicit none
      integer(4),intent(in) :: mnnode
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: rho,pr,E,eta,csound,machno,divU,Qcrit
      real(rp),dimension(numNodesRankPar),intent(in) :: envit,mut,mu_fluid
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: u,gradRho,curlU

      character(512) :: dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSIZE_T), dimension(2) :: ds_dims2d,ms_dims2d
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(HSSIZE_T), dimension(2) :: ms_offset2d
      integer(4) :: ds_rank,h5err

      integer(4) :: ii,iElemL
      integer(1),allocatable :: aux_array_i1(:)

      !------------------------------------------------------------------------------------

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call create_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id)

      !------------------------------------------------------------------------------------------------------
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      ! ## RHO ##
      dsetname = '/VTKHDF/PointData/rho'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,rho)

      ! ## PR ##
      dsetname = '/VTKHDF/PointData/pr'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,pr)

      ! ## Ener ##
      dsetname = '/VTKHDF/PointData/Ener'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,E)

      ! ## eta ##
      dsetname = '/VTKHDF/PointData/eta'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,eta)

      ! ## csound ##
      dsetname = '/VTKHDF/PointData/csound'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,csound)

      ! ## machno ##
      dsetname = '/VTKHDF/PointData/machno'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,machno)

      ! ## divU ##
      dsetname = '/VTKHDF/PointData/divU'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,divU)

      ! ## Qcrit ##
      dsetname = '/VTKHDF/PointData/Qcrit'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,Qcrit)

      ! ## envit ##
      dsetname = '/VTKHDF/PointData/envit'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,envit)

      ! ## mut ##
      dsetname = '/VTKHDF/PointData/mut'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mut)

      ! ## mu_fluid ##
      dsetname = '/VTKHDF/PointData/mu_fluid'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mu_fluid)

      !-------------------------------------------------------------------------------------------------------
      ds_rank = 2
      ds_dims2d(1) = 3
      ms_dims2d(1) = 1
      ds_dims2d(2) = totalNumNodesPar
      ms_dims2d(2) = numNodesRankPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = rankNodeStart-1
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/Velocity'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,u)
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/gradRho'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,gradRho)
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/curlU'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,curlU)
      !-------------------------------------------------------------------------------------------------------

      allocate(aux_array_i1(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = h5_datatype_uint1
      ds_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_instResultsFile

   subroutine write_vtkhdf_avgResultsFile(mnnode,full_fileName,avrho,avpre,avmueff,avvel,avve2,avvex,avtw)
      implicit none
      integer(4),intent(in) :: mnnode
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: avrho,avpre,avmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: avvel,avve2,avvex,avtw

      character(512) :: dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSIZE_T), dimension(2) :: ds_dims2d,ms_dims2d
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(HSSIZE_T), dimension(2) :: ms_offset2d
      integer(4) :: ds_rank,h5err

      integer(4) :: iElemL
      integer(1),allocatable :: aux_array_i1(:)

      !------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call create_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id)

      !------------------------------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      ! ## avgRHO ##
      dsetname = '/VTKHDF/PointData/avrho'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avrho)

      ! ## avgPRE ##
      dsetname = '/VTKHDF/PointData/avpre'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avpre)

      ! ## avgMUEFF ##
      dsetname = '/VTKHDF/PointData/avmueff'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avmueff)

      !-------------------------------------------------------------------------------------------------------
      ds_rank = 2
      ds_dims2d(1) = 3
      ms_dims2d(1) = 1
      ds_dims2d(2) = totalNumNodesPar
      ms_dims2d(2) = numNodesRankPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = rankNodeStart-1
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/avvel'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avvel)
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/avve2'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avve2)
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/avvex'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avvex)
      !--------------------------------------------------------------------------------
      dsetname = '/VTKHDF/PointData/avtw'
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avtw)
      !--------------------------------------------------------------------------------

      !-------------------------------------------------------------------------------------------------------
      allocate(aux_array_i1(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = h5_datatype_uint1
      ds_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_avgResultsFile

   subroutine write_vtkhdf_realFieldFile(mnnode,full_fileName,realField)
      implicit none
      integer(4),intent(in) :: mnnode
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: realField

      character(512) :: dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer :: ds_rank,h5err

      integer(4) :: iElemL
      integer(1),allocatable :: aux_array_i1(:)

      !------------------------------------------------------------------------------------
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create VTKHDF ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      call create_vtkhdf_unstructuredGrid_meshFile(mnnode,file_id)

      !------------------------------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      ! ## realField ##
      dsetname = '/VTKHDF/PointData/realField'
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,realField)

      !-------------------------------------------------------------------------------------------------------
      allocate(aux_array_i1(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype =  h5_datatype_uint1
      ds_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call write_dataspace_1d_uint1_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_realFieldFile

!-----------------------------------------------------------------------------------------------------------------------------

   subroutine save_vtkhdf_instResultsFile(mnnode,iStep,rho,pr,E,eta,csound,machno,divU,Qcrit,envit,mut,mu_fluid,u,gradRho,curlU)
      implicit none
      integer, intent(in) :: mnnode,iStep
      real(rp),dimension(numNodesRankPar),intent(in) :: rho,pr,E,eta,csound,machno,divU,Qcrit
      real(rp),dimension(numNodesRankPar),intent(in) :: envit,mut,mu_fluid
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: u,gradRho,curlU
      character(512) :: full_fileName

      !------------------------------------------------------------------------------------
      call set_vtkhdf_resultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF inst results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_instResultsFile(mnnode,full_fileName,rho,pr,E,eta,csound,machno,divU,Qcrit,envit,mut,mu_fluid,u,gradRho,curlU)

   end subroutine save_vtkhdf_instResultsFile

   subroutine save_vtkhdf_avgResultsFile(mnnode,iStep,avrho,avpre,avmueff,avvel,avve2,avvex,avtw)
      implicit none
      integer, intent(in) :: mnnode,iStep
      real(rp),dimension(numNodesRankPar),intent(in) :: avrho,avpre,avmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: avvel,avve2,avvex,avtw
      character(512) :: full_fileName

      call set_vtkhdf_avgResultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF avg results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_avgResultsFile(mnnode,full_fileName,avrho,avpre,avmueff,avvel,avve2,avvex,avtw)

   end subroutine save_vtkhdf_avgResultsFile

   subroutine save_vtkhdf_finalAvgResultsFile(mnnode,favrho,favpre,favmueff,favvel,favve2,favvex,favtw)
      implicit none
      integer(4),intent(in) :: mnnode
      real(rp),dimension(numNodesRankPar),intent(in) :: favrho,favpre,favmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: favvel,favve2,favvex,favtw
      character(512) :: full_fileName

      call set_vtkhdf_finalAvgResultsFile_name(full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF final avg results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_avgResultsFile(mnnode,full_fileName,favrho,favpre,favmueff,favvel,favve2,favvex,favtw)

   end subroutine save_vtkhdf_finalAvgResultsFile

   subroutine save_vtkhdf_realFieldFile(mnnode,realField)
      implicit none
      integer(4),intent(in) :: mnnode
      real(rp),dimension(numNodesRankPar),intent(in) :: realField
      character(512) :: full_fileName

      !------------------------------------------------------------------------------------
      call set_vtkhdf_resultsFile_name(0,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF realField file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_realFieldFile(mnnode,full_fileName,realField)

   end subroutine save_vtkhdf_realFieldFile

!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------

   subroutine get_dims(file_id,dsetname, ms_rank, fs_dims, fs_maxdims)
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer, intent(in) :: ms_rank
      integer(hid_t) :: dset_id,fspace_id
      integer :: h5err
      integer(hsize_t), intent(out) :: fs_dims(ms_rank),fs_maxdims(ms_rank)
      integer(hid_t) :: dtype
      dtype = H5T_NATIVE_INTEGER

      call h5dopen_f(file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine get_dims

!-------------------------------WITNESS POINTS-------------------------------!
   subroutine create_witness_hdf5(full_fileName, mnnode, xyz, witel, witxi, shapewit, nwit, nwitPar, witGlob, save_u_i, save_pr, save_rho)
      implicit none
      character(512),intent(in)  :: full_fileName
      integer(4),intent(in)      :: mnnode, nwit, nwitPar
      integer(4),intent(in)      :: witel(nwit), witGlob(nwit)
      real(rp),intent(in)        :: witxi(nwit, ndime), shapewit(nwit,mnnode)
      real(rp),intent(in)        :: xyz(nwit,ndime)
      logical,intent(in)         :: save_u_i, save_pr, save_rho
      integer(4)                 :: aux(1), nwitParAllRanks(mpi_size), nwitOffset=0, inode
      integer(hid_t)             :: file_id,plist_id,dset_id,dspace_id,group_id, dtype
      integer(hsize_t)           :: ds_dims(1),ms_dims(1),max_dims(1),chunk_dims(1)
      integer(hssize_t)          :: ms_offset(1)
      integer(hsize_t)           :: ds_dims2d(2),ms_dims2d(2),max_dims2d(2),chunk_dims2d(2)
      integer(hssize_t)          :: ms_offset2d(2)
      integer(4)                 :: ds_rank,h5err,irank,iwit
      character(256)             :: dsetname
      real(rp)                   :: auxwitxyz(nwitPar, ndime), auxwitxi(nwitPar,ndime), auxshapefunc(nwitPar,mnnode)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !Create dataspece for nwitPar and save it!
      dtype        = h5_datatype_int4
      ds_rank      = 1
      dsetname     = 'nwitPar'
      ds_dims(1)   = mpi_size
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      aux(1)       = nwitPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux)
      call MPI_Barrier(app_comm, mpi_err)

      !Compute sum of nwitPar until that rank and save it on its dataspace!
      dtype        = h5_datatype_int4
      ms_dims(1)   = mpi_size
      ms_offset(1) = 0
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,nwitParAllRanks)
      if (mpi_rank > 0) then
         do irank = 1, mpi_rank
            nwitOffset = nwitOffset + nwitParAllRanks(irank)
         end do
      end if
      ds_rank      = 1
      dsetname     = 'nwitOffset'
      ds_dims(1)   = mpi_size
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      aux(1)       = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,aux)
      call MPI_Barrier(app_comm, mpi_err)

      !Create dataspece for global numeration and save it!
      dtype        = h5_datatype_int4
      dsetname     = 'global'
      ds_rank      = 1
      ds_dims(1)   = nwit
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,witGlob)

      !Create dataspece for element containing the witness and save it!
      dtype        = h5_datatype_int4
      dsetname     = 'element'
      ds_rank      = 1
      ds_dims(1)   = nwit
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,witel)

      !Create dataspace for witness coordinates and save them!
      do iwit = 1, nwitPar
         auxwitxyz(iwit,:) = xyz(iwit,:)
         auxwitxi(iwit,:)  = witxi(iwit,:)
         auxshapefunc(iwit,:) = shapewit(iwit,:)
      end do
      dsetname       = 'xyz'
      ds_dims2d(1)   = ndime
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = ndime
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwitxyz)

      !Create dataspace for witness isoparametric coordinates and save them!
      dsetname       = 'witxi'
      ds_dims2d(1)   = ndime
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = ndime
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwitxi)

      !Create dataspace for the shape functions evaluated on the witness points and save them!
      dsetname       = 'shape_functions'
      ds_dims2d(1)   = mnnode
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = mnnode
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxshapefunc)

      !Create time dataset!
      call select_dtype_rp_vtk(dtype)
      ds_rank       = 1
      dsetname      = 'time'
      ds_dims(1)    = 1
      max_dims(1)   = H5S_UNLIMITED_F
      chunk_dims(1) = 1
      call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)

      !Create istep dataset!
      dtype         = h5_datatype_int4
      ds_rank       = 1
      dsetname      = 'istep'
      ds_dims(1)    = 1
      max_dims(1)   = H5S_UNLIMITED_F
      chunk_dims(1) = 1
      call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)

      !Create dataspaces for the magnitudes to save!
      call select_dtype_rp_vtk(dtype)
      ds_rank         = 2
      ds_dims2d(1)    = 1
      ds_dims2d(2)    = nwit
      max_dims2d(1)   = H5S_UNLIMITED_F
      max_dims2d(2)   = nwit
      chunk_dims2d(1) = 1
      chunk_dims2d(2) = nwit
      if (save_u_i) then
         dsetname = 'u_x'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)

         dsetname = 'u_y'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)

         dsetname = 'u_z'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)
      end if

      if (save_pr) then
         dsetname = 'pr'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)
      end if

      if (save_rho) then
         dsetname = 'rho'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims2d,max_dims2d,chunk_dims2d,dtype)
      end if

      call h5fclose_f(file_id,h5err)

   end subroutine create_witness_hdf5

   subroutine load_witness_hdf5(full_fileName, mnnode, nwit, loadstep, load_stepwit, nwitPar, witel, witxi, shapefunc)
      implicit none
      character(512), intent(in)  :: full_fileName
      integer(4),     intent(in)  :: mnnode, nwit, loadstep
      integer(4),     intent(out) :: witel(nwit)
      real(rp),       intent(out) :: witxi(nwit,ndime), shapefunc(nwit,mnnode)
      integer(4),     intent(out) :: nwitPar, load_stepwit
      integer(hid_t)              :: file_id,plist_id,dset_id,dspace_id,group_id
      integer(hsize_t)            :: ds_dims(1),ms_dims(1),max_dims(1),chunk_dims(1)
      integer(hssize_t)           :: ms_offset(1)
      integer(hsize_t)            :: ds_dims2d(2),ms_dims2d(2),max_dims2d(2),chunk_dims2d(2)
      integer(hssize_t)           :: ms_offset2d(2)
      integer                     :: ms_rank, h5err, iwit, istep
      integer(hsize_t)            :: nsteps(2), maxnsteps(2)
      character(256)              :: dsetname
      real(rp), allocatable       :: auxwitxi(:,:), auxshapefunc(:,:)
      integer(4)                  :: nwitOffset, auxread(1)
      integer(4), allocatable     :: steps(:)

      witel(:)   = 0
      witxi(:,:) = 0.0_rp
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !Read step to continue!
      dsetname     = 'istep'
      ms_rank      = 1
      call get_dims(file_id, dsetname, ms_rank, nsteps, maxnsteps)
      ms_dims(1)   = nsteps(1)
      ms_offset(1) = 0
      allocate(steps(nsteps(1)))
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,steps)
      do istep = nsteps(1), 1, -1
         if (steps(istep) < loadstep) then
            load_stepwit = istep+1
            exit
         end if
      end do
      deallocate(steps)

      !Read nwitPar!
      dsetname     = 'nwitPar'
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,auxread)
      nwitPar = auxread(1)
      allocate(auxwitxi(nwitPar,ndime))
      allocate(auxshapefunc(nwitPar,mnnode))

      !Read nwitOffset!
      dsetname     = 'nwitOffset'
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,auxread)
      nwitOffset = auxread(1)

      !Read elements containing the witness points
      dsetname     = 'element'
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,witel)

      !Read witness isoparametric coordinates!
      dsetname       = 'witxi'
      ms_dims2d(1)   = ndime
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,auxwitxi)

      !Read the shape functions coordinates!
      dsetname       = 'shape_functions'
      ms_dims2d(1)   = mnnode
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = 0
      ms_offset2d(2) = nwitOffset
      call read_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ms_dims2d,ms_offset2d,auxshapefunc)

      do iwit = 1,nwitPar
         witxi(iwit,:) = auxwitxi(iwit,:)
         shapefunc(iwit,:) = auxshapefunc(iwit,:)
      end do

      deallocate(auxwitxi)
      deallocate(auxshapefunc)

      call h5fclose_f(file_id,h5err)

   end subroutine load_witness_hdf5

   subroutine update_witness_hdf5(itewit, leapwitsave, witval, nwit, nwitPar, nvarwit, full_fileName, t, steps, save_u_i, save_pr, save_rho)
      integer(4),     intent(in) :: itewit, nwit, nwitPar, nvarwit, leapwitsave, steps(leapwitsave)
      real(rp),       intent(in) :: witval(leapwitsave, nwitPar, nvarwit), t(leapwitsave)
      logical,        intent(in) :: save_u_i, save_pr, save_rho
      character(512), intent(in) :: full_fileName
      character(256)             :: dsetname
      real(rp)                   :: auxwrite(leapwitsave,nwitPar)
      integer(4)                 :: h5err, iwit, ds_rank, ileap
      integer(4)                 :: nwitOffset, auxread(1)
      integer(hsize_t)           :: ds_dims(1),ms_dims(1),max_dims(1),chunk_dims(1)
      integer(hssize_t)          :: ms_offset(1)
      integer(hsize_t)           :: ds_dims2d(2),ms_dims2d(2),max_dims2d(2),chunk_dims2d(2)
      integer(hssize_t)          :: ms_offset2d(2)
      integer(hid_t)             :: file_id,plist_id, dtype

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_fileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !Read nwitOffset!
      dsetname     = 'nwitOffset'
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,auxread)
      nwitOffset = auxread(1)

      !Save variables!
      ms_dims        = nwitPar
      ds_rank        = 2
      ds_dims2d(1)   = itewit
      ds_dims2d(2)   = nwit
      ms_dims2d(1)   = leapwitsave
      ms_dims2d(2)   = nwitPar
      ms_offset2d(1) = itewit-leapwitsave
      ms_offset2d(2) = nwitOffset
      if (save_u_i) then
         dsetname = 'u_x'
         !$acc kernels
         auxwrite(:,:) = witval(:,:,1)
         !$acc end kernels
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         dsetname = 'u_y'
         !$acc kernels
         auxwrite(:,:) = witval(:,:,2)
         !$acc end kernels
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
         dsetname = 'u_z'
         !$acc kernels
         auxwrite(:,:) = witval(:,:,3)
         !$acc end kernels
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
      end if
      if (save_pr) then
         dsetname = 'pr'
         !$acc kernels
         auxwrite(:,:) = witval(:,:,4)
         !$acc end kernels
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
      end if
      if (save_rho) then
         dsetname = 'rho'
         !$acc kernels
         auxwrite(:,:) = witval(:,:,5)
         !$acc end kernels
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims2d)
         call save_array2D_tr_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,auxwrite,isCreateDataspaceOpt=.False.)
      end if

      !Save time!
      dsetname     = 'time'
      ms_dims(1)   = leapwitsave
      ms_offset(1) = itewit - leapwitsave
      ds_rank      = 1
      ds_dims(1)   = itewit
      call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call save_array1D_rp_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,t,isCreateDataspaceOpt=.False.)

      !Save istep!
      dsetname     = 'istep'
      ms_dims(1)   = leapwitsave
      ms_offset(1) = itewit - leapwitsave
      ds_rank      = 1
      ds_dims(1)   = itewit
      call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call write_dataspace_1d_int4_hyperslab_parallel(file_id,dsetname,ms_dims,ms_offset,steps)
      call h5fclose_f(file_id,h5err)

   end subroutine update_witness_hdf5

end module mod_hdf5
