module mod_hdf5
   use hdf5
   use mod_constants
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms
   use mod_comms_boundaries
   implicit none

   character(256) :: meshFile_h5_name,base_resultsFile_h5_name,base_avgResultsFile_h5_name
   integer(hid_t) :: meshFile_h5_id,resultsFile_h5_id

   integer(hid_t) :: h5_datatype_uint1,h5_datatype_int1,h5_datatype_int4,h5_datatype_int8
   integer(hid_t) :: h5_datatype_real4,h5_datatype_real8

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

   subroutine create_hdf5_groups_datasets_in_meshFile_from_tool(hdf5_file_id,isPeriodic,isBoundaries,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8,&
               vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,&
               vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)
      implicit none
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
      integer(HSIZE_T), dimension(1) :: ds_dims
      integer(4) :: ds_rank,ms_rank,h5err
      character(256) :: groupname,dsetname

      !--------------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numNodesParTotal_i8

      if(rp.eq.4) then
         dtype = h5_datatype_real4
      else if(rp.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in create_hdf5_groups_datasets_in_meshFile_from_tool! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      !dtype = h5_datatype_real

      groupname = '/Coords'
      call create_group_hdf5(hdf5_file_id,groupname)
      dsetname = '/Coords/X'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/Coords/Y'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      dsetname = '/Coords/Z'
      call create_dataspace_hdf5(hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
  
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

      call create_groups_datasets_connectivity_workingNodes_hdf5(hdf5_file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)

      if(isBoundaries) then
         call create_groups_datasets_boundary_data_hdf5(hdf5_file_id,numMshRanks2Part,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank)
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

   end subroutine create_hdf5_groups_datasets_in_meshFile_from_tool

   subroutine create_groups_datasets_connectivity_workingNodes_hdf5(file_id,numMshRanks2Part,numElemsGmsh,vecNumWorkingNodes)
      implicit none
      integer(hid_t),intent(in) :: file_id
      integer,intent(in) :: numMshRanks2Part,numElemsGmsh
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
      ds_dims(1) = int(numElemsGmsh,hsize_t)*int(nnode,hsize_t)
      if(mpi_rank.eq.0) write(*,*) 'debug ds_dims',ds_dims(1),'numElemsGmsh',numElemsGmsh,'nnode',nnode

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
      ds_dims(1)=0
      do mshRank=0,numMshRanks2Part-1
         ds_dims(1)=ds_dims(1)+int(vecNumWorkingNodes(mshRank),hsize_t)
      end do
      if(mpi_rank.eq.0) write(*,*) 'debug workingNodesPar ds_dims',ds_dims(1)
      !ds_dims(1) = accumVal

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

   subroutine write_mshRank_data_in_hdf5_meshFile_from_tool(hdf5_file_id,mshRank,numMshRanks2Part,isPeriodic,isBoundaries,numElemsGmsh,numBoundFacesGmsh,&
         numElemsMshRank,mshRankElemStart,mshRankElemEnd,mshRankNodeStart_i8,mshRankNodeEnd_i8,numNodesMshRank,numWorkingNodesMshRank,numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank,maxBoundCode,&
         elemGidMshRank,globalIdSrlMshRank_i8,globalIdParMshRank_i8,connecVTKMshRank,connecParOrigMshRank,connecParWorkMshRank,coordParMshRank,workingNodesMshRank,&
         boundaryNodesMshRank,dofNodesMshRank,boundFacesCodesMshRank,boundFacesOrigMshRank,boundFacesMshRank,numPerNodesMshRank,masSlaNodesMshRank,&
         numNodesToCommMshRank,numMshRanksWithComms,nodesToCommMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank,&
         bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms,bnd_nodesToCommMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank,&
         vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer(4),intent(in) :: mshRank,numMshRanks2Part,numElemsGmsh,numBoundFacesGmsh
      integer(4),intent(in) :: numElemsMshRank,mshRankElemStart,mshRankElemEnd
      integer(8),intent(in) :: mshRankNodeStart_i8,mshRankNodeEnd_i8
      integer(4),intent(in) :: numNodesMshRank,numWorkingNodesMshRank,numPerNodesMshRank
      integer(4),intent(in) :: numBoundFacesMshRank,numBoundaryNodesMshRank,numDoFMshRank
      integer(4),intent(in) :: maxBoundCode

      integer(4),intent(in) :: elemGidMshRank(numElemsMshRank),workingNodesMshRank(numWorkingNodesMshRank)
      integer(8),intent(in) :: globalIdSrlMshRank_i8(numNodesMshRank),globalIdParMshRank_i8(numNodesMshRank)
      integer(4),intent(in) :: connecVTKMshRank(numElemsMshRank*nnode)
      integer(4),intent(in) :: connecParOrigMshRank(numElemsMshRank,nnode),connecParWorkMshRank(numElemsMshRank,nnode)
      real(rp),intent(in)   :: coordParMshRank(numNodesMshRank,3)

      integer(4),intent(in) :: boundaryNodesMshRank(numBoundaryNodesMshRank),dofNodesMshRank(numDoFMshRank)
      integer(4),intent(in) :: boundFacesCodesMshRank(numBoundFacesMshRank),boundFacesOrigMshRank(numBoundFacesMshRank,npbou),boundFacesMshRank(numBoundFacesMshRank,npbou)
      integer(4),intent(in) :: numNodesToCommMshRank,numMshRanksWithComms
      integer(4),intent(in) :: masSlaNodesMshRank(numPerNodesMshRank,2)
      !integer(4),intent(in),dimension(numNodesToCommMshRank,3) :: matrixCommSchemeMshRank
      integer(4),intent(in),dimension(numMshRanksWithComms) :: nodesToCommMshRank,commsMemPosInLocMshRank,commsMemSizeMshRank,commsMemPosInNgbMshRank,ranksToCommMshRank
      integer(4),intent(in) :: bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms
      !integer(4),intent(in),dimension(bnd_numNodesToCommMshRank,3) :: bnd_matrixCommSchemeMshRank
      integer(4),intent(in),dimension(bnd_numMshRanksWithComms) :: bnd_nodesToCommMshRank,bnd_commsMemPosInLocMshRank,bnd_commsMemSizeMshRank,bnd_commsMemPosInNgbMshRank,bnd_ranksToCommMshRank

      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank
      integer,intent(in),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank
      !-------------------------------------------------------------------------------------------------------------------------------
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,m,iBound,iElemL,ms_rank
      integer(hssize_t), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)
      integer(8),allocatable :: aux_array_i8(:)
      !-------------------------------------------------------------------------------------------------------------------------------

      ms_rank = 1
      ms_dims(1) = int(numNodesMshRank,hsize_t)
      ms_offset(1) = mshRankNodeStart_i8-1

      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrlMshRank_i8)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdParMshRank_i8)

      ms_dims(1) = int(numElemsMshRank,hsize_t)
      ms_offset(1) = int(mshRankElemStart,hssize_t)-1

      dsetname = '/globalIds/elemGid'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,elemGidMshRank)

      !---------------------------------------------------------------------------------------------------------------------
      ms_rank = 1
      ms_dims(1) = int(numNodesMshRank,hsize_t)
      ms_offset(1) = mshRankNodeStart_i8-1

      dsetname = '/Coords/X'
      call write_dataspace_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,coordParMshRank(:,1))

      dsetname = '/Coords/Y'
      call write_dataspace_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,coordParMshRank(:,2))

      dsetname = '/Coords/Z'
      call write_dataspace_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,coordParMshRank(:,3))

      !---------------------------------------------------------------------------------------------------------------------
      if(numMshRanks2Part.ge.2) then
         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array_i8(1))

         dsetname = '/Parallel_data/rankNodeStart'
         aux_array_i8(1)=mshRankNodeStart_i8
         call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

         dsetname = '/Parallel_data/rankNodeEnd'
         aux_array_i8(1)=mshRankNodeEnd_i8
         call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

         deallocate(aux_array_i8)
         allocate(aux_array(1))

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
            ms_offset(1)=ms_offset(1)+int(vecNumMshRanksWithComms(i),hssize_t)
         end do
         ms_dims(1)=int(numMshRanksWithComms,hsize_t)
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
            ms_offset(1)=ms_offset(1)+int(vecNumNodesToCommMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numNodesToCommMshRank,hsize_t)
         !write(*,*) '[',mshRank,']ms_offset',ms_offset(1),'ms_dims(1)',ms_dims(1)!,' ds_dims ',ds_dims(1)

         dsetname = '/Parallel_data/nodesToComm'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,nodesToCommMshRank)

      end if
      !------------------------------------------------------------------------------------------------------------------------

      if(isBoundaries) then

         dsetname = '/Boundary_data/numBoundCodes'
         if(mshRank.eq.0) then
            ms_rank = 1
            ms_dims(1) = 1
            ms_offset(1) = int(mshRank,hssize_t)
            allocate(aux_array(1))
            aux_array(1)=maxBoundCode

            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
            deallocate(aux_array)
         else
            ms_rank = 1
            ms_dims(1) = 0
            ms_offset(1) = 0
            allocate(aux_array(0))

            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
            deallocate(aux_array)
         end if

         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array(1))

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
            ms_offset(1)=ms_offset(1)+int(vecNumDoFMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numDoFMshRank,hsize_t)

         dsetname = '/Boundary_data/ldofPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,dofNodesMshRank)
         !------------------------------------------------------------------------------------------------------------------------
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumBoundaryNodesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundaryNodesMshRank,hsize_t)

         dsetname = '/Boundary_data/lbnodesPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,boundaryNodesMshRank)
         !------------------------------------------------------------------------------------------------------------------------
         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumBoundFacesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundFacesMshRank,hsize_t)

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

         ms_dims(1) = int(ms_dims(1),hsize_t)*int(npbou,hsize_t)
         ms_offset(1) = ms_offset(1)*int(npbou,hssize_t)

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
         if(numMshRanks2Part.ge.2) then
            ms_rank = 1
            ms_dims(1) = 1
            ms_offset(1) = int(mshRank,hssize_t)
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
               ms_offset(1)=ms_offset(1)+int(vecBndNumMshRanksWithComms(i),hssize_t)
            end do
            ms_dims(1)=int(bnd_numMshRanksWithComms,hsize_t)

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
               ms_offset(1)=ms_offset(1)+int(vecBndNumNodesToCommMshRank(i),hssize_t)
            end do
            ms_dims(1)=int(bnd_numNodesToCommMshRank,hsize_t)

            dsetname = '/Parallel_data_boundary/nodesToComm'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_nodesToCommMshRank)

            !------------------------------------------------------------------------------------------------------------------------

         end if
      end if

      !---------------------------------------------------------------------------------------------------------------------
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mshRank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      aux_array(1)=numWorkingNodesMshRank
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      deallocate(aux_array)
      !------------------------------------------------------------------------------------------------------------------------
      !  SAVING connecParOrig(:,:)
      allocate(aux_array(numElemsMshRank*nnode))
      ms_rank = 1
      ms_dims(1) = int(numElemsMshRank,hsize_t)*int(nnode,hsize_t)
      ms_offset(1) = int((mshRankElemStart-1),hssize_t)*int(nnode,hssize_t)
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
         ms_offset(1)=ms_offset(1)+int(vecNumWorkingNodes(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesMshRank,hsize_t)

      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,workingNodesMshRank)
      !-----------------------------------------------------------------------------------------------------

      !------------------------------------------------------------------------------------------------------------------------
      if(isPeriodic) then
         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = int(mshRank,hssize_t)
         allocate(aux_array(1))

         dsetname = '/Periodic_data/nPerRankPar'
         aux_array(1)= numPerNodesMshRank
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         ms_offset(1)=0
         do i=0,mshRank-1 !from rank 0 mpi_rank-1
            ms_offset(1)=ms_offset(1)+int(vecNumPerNodesMshRank(i),hssize_t)
         end do
         ms_dims(1)=int(numPerNodesMshRank,hsize_t)

         dsetname = '/Periodic_data/masSlaRankPar1'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaNodesMshRank(:,1))

         dsetname = '/Periodic_data/masSlaRankPar2'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,masSlaNodesMshRank(:,2))
         !--------------------------------------------------------------------------------------------------------
      end if

   end subroutine write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine dummy_write_mshRank_data_in_hdf5_meshFile_from_tool(hdf5_file_id,numMshRanks2Part,isPeriodic,isBoundaries)
      implicit none
      integer(hid_t),intent(in) :: hdf5_file_id
      logical,intent(in) :: isPeriodic,isBoundaries
      integer(4),intent(in) :: numMshRanks2Part
      !-------------------------------------------------------------------------------------------------------------------------------
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer :: i,m,iBound,iElemL,ms_rank
      integer(hssize_t), dimension(1) :: ms_offset 
      integer(4),allocatable :: empty_array_i4(:)
      integer(8),allocatable :: empty_array_i8(:)
      real(rp),allocatable :: empty_array_rp(:)

      !-------------------------------------------------------------------------------------------------------------------------------

      allocate(empty_array_i4(0))
      allocate(empty_array_i8(0))
      allocate(empty_array_rp(0))
      ms_rank = 1
      ms_dims(1) = 0
      ms_offset(1) = 0

      dsetname = '/globalIds/globalIdSrl'
      call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/globalIds/globalIdPar'
      call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i8)

      dsetname = '/globalIds/elemGid'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

      !---------------------------------------------------------------------------------------------------------------------
      dsetname = '/Coords/X'
      call write_dataspace_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_rp)

      dsetname = '/Coords/Y'
      call write_dataspace_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_rp)

      dsetname = '/Coords/Z'
      call write_dataspace_real_rp_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_rp)

      !---------------------------------------------------------------------------------------------------------------------
      if(numMshRanks2Part.ge.2) then
         dsetname = '/Parallel_data/rankNodeStart'
         call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i8)

         dsetname = '/Parallel_data/rankNodeEnd'
         call write_dataspace_int8_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i8)

         dsetname = '/Parallel_data/rankElemStart'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/rankElemEnd'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/numRanksWithComms'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/numNodesToComm'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/ranksToComm'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemPosInLoc'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemSize'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/commsMemPosInNgb'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Parallel_data/nodesToComm'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

      end if
      !------------------------------------------------------------------------------------------------------------------------

      if(isBoundaries) then

         dsetname = '/Boundary_data/numBoundCodes'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/numBoundsRankPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/ndofRankPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         !------------------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/ldofPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
         !------------------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/lbnodesPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
         !------------------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/bouCodesPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
         !------------------------------------------------------------------------------------------------------------------------
         dsetname = '/Boundary_data/boundPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Boundary_data/boundParOrig'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

      !------------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------
         if(numMshRanks2Part.ge.2) then
            dsetname = '/Parallel_data_boundary/numRanksWithComms'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/numNodesToComm'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/ranksToComm'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemSize'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            dsetname = '/Parallel_data_boundary/nodesToComm'
            call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

            !------------------------------------------------------------------------------------------------------------------------

         end if
      end if

      !---------------------------------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
      !------------------------------------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/connecParOrig'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
      !-----------------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/connecParWork'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      dsetname = '/Connectivity/connecVTK'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
      !-----------------------------------------------------------------------------------------------------

      dsetname = '/Connectivity/workingNodesPar'
      call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
      !-----------------------------------------------------------------------------------------------------

      !------------------------------------------------------------------------------------------------------------------------
      if(isPeriodic) then

         dsetname = '/Periodic_data/nPerRankPar'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Periodic_data/masSlaRankPar1'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)

         dsetname = '/Periodic_data/masSlaRankPar2'
         call write_dataspace_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_rank,ms_dims,ms_offset,empty_array_i4)
         !--------------------------------------------------------------------------------------------------------
      end if

      deallocate(empty_array_i4)
      deallocate(empty_array_i8)
      deallocate(empty_array_rp)

   end subroutine dummy_write_mshRank_data_in_hdf5_meshFile_from_tool

   subroutine load_hdf5_meshFile()
      implicit none     
      character(256) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id,fspace_id
      integer(4) :: h5err
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
      call load_boundary_data_hdf5(file_id)
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
      call load_connectivity_hdf5(file_id)

      !--------------------------------------------------------------------------------
      !load globalIds
      call load_globalIds_hdf5(file_id)

      !close h5 file
      call h5fclose_f(file_id,h5err)

   end subroutine load_hdf5_meshFile

   subroutine close_hdf5_meshFile()
      implicit none
      integer(4) :: h5err

      !close the file.
      call h5fclose_f(meshFile_h5_id,h5err)
   end subroutine close_hdf5_meshFile

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
   subroutine write_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      if(rp.eq.4) then
         dtype = h5_datatype_real4
      else if(rp.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in write_dataspace_real_rp_hyperslab_parallel! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_real_rp_hyperslab_parallel

   subroutine read_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      if(rp.eq.4) then
         dtype = h5_datatype_real4
      else if(rp.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in read_dataspace_real_rp_hyperslab_parallel! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_real_rp_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
   subroutine write_dataspace_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      if(rp_vtk.eq.4) then
         dtype = h5_datatype_real4
      else if(rp_vtk.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in write_dataspace_real_rp_hyperslab_parallel! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,dtype,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_real_rp_vtk_hyperslab_parallel

   subroutine read_dataspace_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(rp_vtk),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      if(rp_vtk.eq.4) then
         dtype = h5_datatype_real4
      else if(rp_vtk.eq.8) then
         dtype = h5_datatype_real8
      else
         write(*,*) 'Fatal error in read_dataspace_real_rp_hyperslab_parallel! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,dtype,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_real_rp_vtk_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------------------------
!  FP32
   subroutine write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      real(4),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_real4,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_fp32_hyperslab_parallel

   subroutine read_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset(ms_rank)
      real(4),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_real4,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_fp32_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  FP64
   subroutine write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset(ms_rank)
      real(8),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_real8,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_fp64_hyperslab_parallel

   subroutine read_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),dimension(ms_rank),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),dimension(ms_rank),intent(in) :: ms_offset(ms_rank)
      real(8),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(4) :: h5err
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_real8,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_fp64_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!  INT1
   subroutine write_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank) 
      integer(1),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_int1,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_int1_hyperslab_parallel

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
   subroutine read_dataspace_int1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_int1,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_int1_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  UINT1
   subroutine write_dataspace_uint1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      
      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_uint1,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_uint1_hyperslab_parallel

   subroutine read_dataspace_uint1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(1),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_uint1,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_uint1_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  INT4
   subroutine write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(4),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_int4,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_int4_hyperslab_parallel

   subroutine read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(4),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_int4,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_int4_hyperslab_parallel
!-------------------------------------------------------------------------------------------------------------------
!  INT8
   subroutine write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2write)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank) 
      integer(8),intent(in) :: data2write(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dwrite_f(dset_id,h5_datatype_int8,data2write,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine write_dataspace_int8_hyperslab_parallel

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
   subroutine read_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data2read)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(len=*),intent(in) :: dsetname
      integer(4),intent(in) :: ms_rank
      integer(hsize_t),intent(in) :: ms_dims(ms_rank)
      integer(hssize_t),intent(in) :: ms_offset(ms_rank)
      integer(8),intent(out) :: data2read(ms_dims(1))
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ms_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err

      call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
                                         dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)

      call h5dread_f(dset_id,h5_datatype_int8,data2read,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)

   end subroutine read_dataspace_int8_hyperslab_parallel

!--------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------

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
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! open file collectively
      call h5fopen_f(meshFile_h5_name, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot open meshfile ',trim(adjustl(meshFile_h5_name))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !-------------------------------------------
      ms_rank = 1
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-------------------------------------------

      dsetname = '/Coords/X'
      call write_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,1))

      dsetname = '/Coords/Y'
      call write_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,2))

      dsetname = '/Coords/Z'
      call write_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,3))


      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine overwrite_coordinates_hdf5

   subroutine save_connectivity_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: groupname,dsetname
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims
      integer(hid_t) :: dtype
      integer(4) :: ds_rank,ms_rank
      integer(4) :: iElemL,i,m,accumVal
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Connectivity')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = totalNumElements*nnode
      ms_rank = 1
      ms_dims(1) = int(numElemsRankPar,hsize_t)*int(nnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)*int(nnode,hssize_t)

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
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !--------------------------------------------------------------------------------------------------------
      !  SAVING workingNodesPar(:)
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesRankPar,hsize_t)

      dsetname = '/Connectivity/workingNodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,workingNodesPar)
      !--------------------------------------------------------------------------------------------------------

   end subroutine save_connectivity_hdf5

   subroutine load_connectivity_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: ms_rank
      integer(4) :: iElemL,i,m
      integer(hssize_t), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading connectivity data hdf5...'

      ms_rank = 1
      ms_dims(1) = int(numElemsRankPar,hsize_t)*int(nnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)*int(nnode,hssize_t)

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

      !-------------------------------------------------------------------------------------------------------
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
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
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numWorkingNodesRankPar of all ranks
      dsetname = '/Connectivity/numWorkingNodesRankPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numWorkingNodesRankPar,hsize_t)

      dsetname = '/Connectivity/workingNodesPar'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,workingNodesPar)
      deallocate(aux_array)
      !-------------------------------------------------------------------------------------------------------

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
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numRanksWithComms,hsize_t)

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
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numNodesToComm,hsize_t)

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
#endif
   subroutine load_parallel_data_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,ms_rank,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
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
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
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

      dsetname = '/Parallel_data/ranksToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ranksToComm)

      dsetname = '/Parallel_data/commsMemPosInLoc'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInLoc)

      dsetname = '/Parallel_data/commsMemSize'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemSize)

      dsetname = '/Parallel_data/commsMemPosInNgb'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,commsMemPosInNgb)
      
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      
      dsetname = '/Parallel_data/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numNodesToComm,hsize_t)

      allocate(nodesToComm(numNodesToComm))

      dsetname = '/Parallel_data/nodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,nodesToComm)

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
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=bnd_numNodesToComm
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      !--------------------------------------------------------------------------------------------------------
      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numRanksWithComms,hsize_t)

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

      ds_dims(1) = int(mpi_size,hsize_t)
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numNodesToComm,hsize_t)

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
#endif
   subroutine load_parallel_data_boundary_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: i,ms_rank,h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      !write(*,*) 'Loading parallel data hdf5...'

      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
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
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      !read data set numRanksWithComms of all ranks
      dsetname = '/Parallel_data_boundary/numRanksWithComms'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
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

      dsetname = '/Parallel_data_boundary/ranksToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_ranksToComm)

      dsetname = '/Parallel_data_boundary/commsMemPosInLoc'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInLoc)

      dsetname = '/Parallel_data_boundary/commsMemSize'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemSize)

      dsetname = '/Parallel_data_boundary/commsMemPosInNgb'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_commsMemPosInNgb)
      
      ms_dims(1) = int(mpi_size,hsize_t)
      ms_offset(1) = 0
      
      dsetname = '/Parallel_data_boundary/numNodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
      !write(*,*) 'rank[',mpi_rank,'] ',aux_array(:)

      ms_offset(1)=0
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(bnd_numNodesToComm,hsize_t)

      allocate(bnd_nodesToComm(bnd_numNodesToComm))

      dsetname = '/Parallel_data_boundary/nodesToComm'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bnd_nodesToComm)

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
      integer(4) :: ds_rank,ms_rank,h5err
      integer(4) :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Periodic_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      allocate(aux_array(1))

      dsetname = '/Periodic_data/nPerRankPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      aux_array(1)=nPerRankPar
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(nPerRankPar,hsize_t)

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
      integer(4) :: ds_rank,ms_rank,h5err
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
         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = int(mpi_rank,hssize_t)
         allocate(aux_array(1))

         !--------------------------------------------------------------------------------------------------------
         dsetname = '/Periodic_data/nPerRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         nPerRankPar=aux_array(1)

         allocate(masSlaRankPar(nPerRankPar,2))
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
         allocate(aux_array(mpi_size))
         ms_dims(1) = int(mpi_size,hsize_t)
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
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(nPerRankPar,hsize_t)

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
      integer(4) :: ds_rank,ms_rank,h5err
      integer(4) :: i,accumVal,iBound,m
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4),allocatable :: aux_array(:)

      groupname = trim('/Boundary_data')
      call create_group_hdf5(file_id,groupname)

      dtype = h5_datatype_int4
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
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
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numBoundsRankPar,hsize_t)

      !write(*,*) 'ms_offset ',ms_offset(1),' ms_dims(1) ',ms_dims(1),' ds_dims ',ds_dims(1)

      dsetname = '/Boundary_data/bouCodesPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bouCodesPar)

      ds_dims(1) = int(ds_dims(1),hsize_t)*int(npbou,hsize_t)!totalNumBoundsSrl*npbou
      ms_dims(1) = int(ms_dims(1),hsize_t)*int(npbou,hsize_t)!numBoundsRankPar*npbou
      ms_offset(1) = ms_offset(1)*int(npbou,hssize_t)

      deallocate(aux_array)
      allocate(aux_array(numBoundsRankPar*npbou))

      !boundPar
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

      !boundParOrig
      i=1
      do iBound=1,numBoundsRankPar
         do m=1,npbou
            aux_array(i)=boundParOrig(iBound,m)
            i=i+1
         end do
      end do

      dsetname = '/Boundary_data/boundParOrig'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(ndofRankPar,hsize_t)

      dsetname = '/Boundary_data/ldofPar'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ldofPar)

      deallocate(aux_array)
      !--------------------------------------------------------------------------------------------------------
      allocate( aux_array(mpi_size) )
      ms_dims(1) = int(mpi_size,hsize_t)
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
         ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
      end do
      ms_dims(1)=int(numBoundaryNodesRankPar,hsize_t)

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
      integer(4) :: ds_rank,ms_rank,h5err
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
         ms_rank = 1
         ms_dims(1) = 1
         ms_offset(1) = 0
         allocate(aux_array(1))

         dsetname = '/Boundary_data/numBoundCodes'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         numBoundCodes=aux_array(1)

         ms_offset(1) = int(mpi_rank,hssize_t)

         dsetname = '/Boundary_data/numBoundsRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         numBoundsRankPar=aux_array(1)

         dsetname = '/Boundary_data/ndofRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         ndofRankPar=aux_array(1)

         dsetname = '/Boundary_data/numBoundaryNodesRankPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
         numBoundaryNodesRankPar=aux_array(1)

         !--------------------------------------------------------------------------------------------------------
         allocate(boundPar(numBoundsRankPar,npbou))
         allocate(boundParOrig(numBoundsRankPar,npbou))
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
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundsRankPar,hsize_t)

         dsetname = '/Boundary_data/bouCodesPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,bouCodesPar)
         deallocate(aux_array)
         !-------------------------------------------------------------------------------------------------------
         ds_dims(1) = int(ds_dims(1),hsize_t)*int(npbou,hsize_t)!totalNumBoundsSrl*npbou
         ms_dims(1) = int(ms_dims(1),hsize_t)*int(npbou,hsize_t)!numBoundsRankPar*npbou
         ms_offset(1) = ms_offset(1)*int(npbou,hssize_t)

         allocate(aux_array(numBoundsRankPar*nnode))

         !boundPar
         dsetname = '/Boundary_data/boundPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         i=1
         do iBound=1,numBoundsRankPar
            do m=1,npbou
               boundPar(iBound,m)=aux_array(i)
               i=i+1
            end do
         end do

         !boundParOrig
         dsetname = '/Boundary_data/boundParOrig'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

         i=1
         do iBound=1,numBoundsRankPar
            do m=1,npbou
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
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

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
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,ldofPar)
         !--------------------------------------------------------------------------------------------------------
         ms_dims(1) = int(mpi_size,hsize_t)
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
            ms_offset(1)=ms_offset(1)+int(aux_array(i),hssize_t)
         end do
         ms_dims(1)=int(numBoundaryNodesRankPar,hsize_t)

         dsetname = '/Boundary_data/lbnodesPar'
         call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,lbnodesPar)
         !--------------------------------------------------------------------------------------------------------
         deallocate(aux_array)
      else 
         numBoundsRankPar=0
         allocate(boundPar(numBoundsRankPar,npbou))
         allocate(boundParOrig(numBoundsRankPar,npbou))
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
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: ms_rank
      integer(HSSIZE_T), dimension(1) :: ms_offset

      allocate(coordPar(numNodesRankPar,ndime)) !only works for ndime=3

      ms_rank = 1
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = '/Coords/X'
      call read_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,1))

      dsetname = '/Coords/Y'
      call read_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,2))

      dsetname = '/Coords/Z'
      call read_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,coordPar(:,3))

   end subroutine load_coordinates_hdf5

   subroutine load_globalIds_hdf5(file_id)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(4) :: ms_rank
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4) :: iNodeL
      integer(8) :: iNodeGSrl,max_iNodeGSrl_l,max_iNodeGSrl_g

      allocate(globalIdSrl(numNodesRankPar))
      allocate(globalIdPar(numNodesRankPar))
      allocate(elemGid(numElemsRankPar))

      ms_rank = 1
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = '/globalIds/globalIdSrl'
      call read_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdSrl)

      dsetname = '/globalIds/globalIdPar'
      call read_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,globalIdPar)

      ms_dims(1) = int(numElemsRankPar,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      dsetname = '/globalIds/elemGid'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,elemGid)

      !--------------------------------------------------------------------------------------------------
      !setting totalNumNodesSrl
      max_iNodeGSrl_l=0
      do iNodeL=1,numNodesRankPar
         iNodeGSrl = globalIdSrl(iNodeL)
         max_iNodeGSrl_l = max(iNodeGSrl,max_iNodeGSrl_l)
      end do

      call MPI_Allreduce(max_iNodeGSrl_l,max_iNodeGSrl_g,1,mpi_datatype_int8,MPI_MAX,MPI_COMM_WORLD,mpi_err)

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
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)

   end subroutine create_dataspace_for_rp_vtk_hdf5

   subroutine save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(1),intent(in) :: ds_dims
      integer(hsize_t),dimension(1),intent(in) :: ms_dims
      integer(hssize_t),dimension(1),intent(in) :: ms_offset 
      real(rp),intent(in) :: data_array_rp(ms_dims(1))

      integer(4) :: ds_rank = 1,ms_rank = 1 !it is forced
      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)
   !---------------------------------------------------------------------------------------------------

      call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims)

      if(rp .eq. rp_vtk) then
         call write_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data_array_rp)
      else
         !copying in the aux array
         allocate(aux_data_array_rp_vtk(ms_dims(1)))
         aux_data_array_rp_vtk(:) = real(data_array_rp(:),rp_vtk)
         call write_dataspace_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_data_array_rp_vtk)
         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine save_array1D_in_dataset_hdf5_file

   subroutine save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(hsize_t),dimension(2),intent(in) :: ds_dims2d
      integer(hsize_t),dimension(2),intent(in) :: ms_dims2d
      integer(hssize_t),dimension(2),intent(inout) :: ms_offset2d 
      real(rp),intent(in) :: data_array_rp(ms_dims2d(2),ms_dims2d(1)) !fortran is column-major & hdf5 writes in row-major

      integer(4) :: ds_rank = 2,ms_rank = 2 !it is forced
      integer(4) :: iCol,h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)
   !---------------------------------------------------------------------------------------------------

      call create_dataspace_for_rp_vtk_hdf5(file_id,dsetname,ds_rank,ds_dims2d)

      if(rp .eq. rp_vtk) then
         do iCol= 1,ds_dims2d(1)
            ms_offset2d(1) = iCol-1
            call write_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,data_array_rp(:,iCol))
         end do
      else
         allocate(aux_data_array_rp_vtk(ms_dims2d(2)))
         do iCol= 1,ds_dims2d(1)
            !copying in the aux array
            aux_data_array_rp_vtk(:) = real(data_array_rp(:,iCol),rp_vtk)
            ms_offset2d(1) = iCol-1
            call write_dataspace_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_data_array_rp_vtk)
         end do
         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine save_array2D_in_dataset_hdf5_file
!------------------------------------------------------------------------------------------------------------------------------

   subroutine read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,data_array_rp)
      implicit none
      integer(hid_t),intent(in) :: file_id
      character(*),intent(in) :: dsetname
      integer(HSIZE_T),dimension(1),intent(in) :: ms_dims
      integer(HSSIZE_T),dimension(1),intent(in) :: ms_offset 
      real(rp),intent(out) :: data_array_rp(ms_dims(1))

      integer(4) :: ms_rank=1 !forced
      integer(4) :: h5err
      real(rp_vtk),allocatable :: aux_data_array_rp_vtk(:)
   !---------------------------------------------------------------------------------------------------

      if(rp .eq. rp_vtk) then
         call read_dataspace_real_rp_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,data_array_rp)
      else
         !copying in the aux array
         allocate(aux_data_array_rp_vtk(ms_dims(1)))

         call read_dataspace_real_rp_vtk_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_data_array_rp_vtk)
         data_array_rp(:) = real(aux_data_array_rp_vtk(:),rp)

         deallocate(aux_data_array_rp_vtk)
      end if

   end subroutine read_array1D_in_dataset_hdf5_file

!----------------------------------------------------------------------------------------------------------------------------------

   subroutine save_hdf5_resultsFile(iStep,time,rho,u,pr,E,eta,csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs)
      implicit none
      integer(4), intent(in) :: iStep
      real(rp),intent(in) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E,eta,csound,machno,mu_fluid,divU,Qcrit
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u
      real(rp),intent(inout),dimension(numElemsRankPar,ngaus) :: mu_e,mu_sgs
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

      !-----------------------------------------------------------------------------------------------
      ! Writing HDF5 Files

      call set_hdf5_resultsFile_name(iStep,full_fileName)

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
      !-----------------------------------------------------------------------------------------------
      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------
      
      dsetname = 'rho'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,u(:,3))

      dsetname = 'mu_fluid'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mu_fluid)

      dsetname = 'mut'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mut)

      dsetname = 'envit'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,envit)

      dsetname = 'gradRho_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,gradRho(:,1))

      dsetname = 'gradRho_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,gradRho(:,2))

      dsetname = 'gradRho_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,gradRho(:,3))

      dsetname = 'curlU_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,curlU(:,1))

      dsetname = 'curlU_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,curlU(:,2))

      dsetname = 'curlU_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,curlU(:,3))

      dsetname = 'pr'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,E)

      dsetname = 'eta'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,eta)

      dsetname = 'csound'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,csound)

      dsetname = 'machno'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,machno)

      dsetname = 'divU'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,divU)

      dsetname = 'Qcrit'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,Qcrit)

      ! ----  time  -----
      ms_dims(1) = 0
      ds_dims = 1
      ms_offset(1) = 0
      if(mpi_rank.eq.0) then
         ms_dims(1) = 1
      endif
      aux_array_time(1) = time

      dsetname = 'time'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,aux_array_time)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine save_hdf5_resultsFile

   subroutine load_hdf5_resultsFile(load_step,time,rho,u,pr,E,mu_e,mu_sgs)
      implicit none
      integer,intent(in) :: load_step
      real(rp),intent(inout) :: time
      real(rp),intent(inout),dimension(numNodesRankPar)       :: rho,pr,E
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u
      real(rp),intent(inout),dimension(numElemsRankPar,ngaus) :: mu_e,mu_sgs
      
      character(512) :: full_loadFileName
      integer(hid_t) :: file_id,plist_id
      integer(HSIZE_T), dimension(1) :: ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4) :: i,j,h5err
      character(128) :: dsetname
      real(rp) :: aux_array_time(1)
      real(rp), dimension(numNodesRankPar) :: envit, mut
      
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

      ! ----  read time  --------------------------------------------------------------------------
      ms_offset(1) = 0
      ms_dims(1) = 1

      dsetname = 'time'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_array_time)
      time = aux_array_time(1)

      if(mpi_rank.eq.0) write(*,*) ' - load_step',load_step,'time',time

      ! ----  read arrays  --------------------------------------------------------------------------

      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'rho'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,3))

      dsetname = 'envit'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,envit)

      dsetname = 'mut'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,mut)

      dsetname = 'pr'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,E)

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

      character(512) :: full_loadFileName
      integer(hid_t) :: file_id,plist_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(4) :: i,j,h5err
      character(128) :: dsetname
      real(rp) :: aux_array_time(1)
      
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

      ! ----  read time  --------------------------------------------------------------------------
      ms_offset(1) = 0
      ms_dims(1) = 1

      dsetname = 'time'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,aux_array_time)
      time = aux_array_time(1)

      if(mpi_rank.eq.0) write(*,*) ' - load_step',load_step,'time',time

      ! ----  read arrays  --------------------------------------------------------------------------

      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'rho'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,rho)

      dsetname = 'u_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,1))

      dsetname = 'u_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,2))

      dsetname = 'u_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,u(:,3))

      dsetname = 'envit'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,envit)

      dsetname = 'mut'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,mut)

      dsetname = 'pr'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,pr)

      dsetname = 'E'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,E)

      dsetname = 'eta'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,eta)

      dsetname = 'csound'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,csound)

      dsetname = 'machno'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,machno)

      dsetname = 'gradRho_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,gradRho(:,1))

      dsetname = 'gradRho_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,gradRho(:,2))

      dsetname = 'gradRho_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,gradRho(:,3))

      dsetname = 'curlU_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,curlU(:,1))

      dsetname = 'curlU_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,curlU(:,2))

      dsetname = 'curlU_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,curlU(:,3))

      dsetname = 'mu_fluid'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,mu_fluid)

      dsetname = 'divU'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,divU)

      dsetname = 'Qcrit'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,Qcrit)

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

   subroutine save_hdf5_avgResultsFile(iStep,avvel,avve2,avvex,avrho,avpre,avmueff,avtw)
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
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(full_fileName,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot create avg results file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'avrho'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avrho)

      dsetname = 'avpre'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avpre)

      dsetname = 'avmueff'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avmueff)

      dsetname = 'avvel_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvel(:,1))

      dsetname = 'avvel_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvel(:,2))

      dsetname = 'avvel_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvel(:,3))

      dsetname = 'avve2_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avve2(:,1))

      dsetname = 'avve2_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avve2(:,2))

      dsetname = 'avve2_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avve2(:,3))

      dsetname = 'avvex_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvex(:,1))

      dsetname = 'avvex_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvex(:,2))

      dsetname = 'avvex_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avvex(:,3))

      dsetname = 'avtw_x'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avtw(:,1))

      dsetname = 'avtw_y'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avtw(:,2))

      dsetname = 'avtw_z'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avtw(:,3))

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine save_hdf5_avgResultsFile

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
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results avg file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      dsetname = 'avrho'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avrho)

      dsetname = 'avpre'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avpre)

      dsetname = 'avmueff'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avmueff)

      dsetname = 'avvel_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvel(:,1))

      dsetname = 'avvel_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvel(:,2))

      dsetname = 'avvel_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvel(:,3))

      dsetname = 'avve2_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avve2(:,1))

      dsetname = 'avve2_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avve2(:,2))

      dsetname = 'avve2_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avve2(:,3))

      dsetname = 'avvex_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvex(:,1))

      dsetname = 'avvex_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvex(:,2))

      dsetname = 'avvex_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avvex(:,3))

      dsetname = 'avtw_x'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avtw(:,1))

      dsetname = 'avtw_y'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avtw(:,2))

      dsetname = 'avtw_z'
      call read_array1D_in_dataset_hdf5_file(file_id,dsetname,ms_dims,ms_offset,avtw(:,3))

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
      integer(hsize_t), dimension(1) :: ds_dims,ms_dims,a_dims
      integer(hsize_t), dimension(2) :: ds_dims2d,ms_dims2d
      integer(hssize_t), dimension(1) :: ms_offset 
      integer(hssize_t), dimension(2) :: ms_offset2d
      integer(size_t) :: attr_length
      integer(4) :: ds_rank,ms_rank,a_rank,h5err
      character(512) :: full_fileName,groupname,dsetname
      character(16) :: attr_value
      integer(4) :: ii,iNodeL,iElemL
      integer(1),allocatable :: aux_array_i1(:)
      integer(4),allocatable :: aux_array_i4(:)
      integer(8),allocatable :: aux_array_i8(:)
      real(8),allocatable :: aux_array_r8(:)

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
      !--------------------------------------------------------------------------------
      dtype = h5_datatype_real8
      ds_rank = 2
      ms_rank = 2
      ds_dims2d(1) = 3
      ms_dims2d(1) = 1
      ds_dims2d(2) = totalNumNodesPar
      ms_dims2d(2) = numNodesRankPar
      ms_offset2d(2) = rankNodeStart-1
!--------------------------------------------------------------------------------
      allocate(aux_array_r8(numNodesRankPar))

      dsetname = '/VTKHDF/Points'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims2d,dtype)
      
!-------------------------------------------------------------------------------- 
      do iNodeL = 1,numNodesRankPar
         aux_array_r8(iNodeL) = real(coordPar(iNodeL,1),8)
      end do
      ms_offset2d(1) = 0
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r8)
!-------------------------------------------------------------------------------- 
      do iNodeL = 1,numNodesRankPar
         aux_array_r8(iNodeL) = real(coordPar(iNodeL,2),8)
      end do
      ms_offset2d(1) = 1
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r8)
!-------------------------------------------------------------------------------- 
      do iNodeL = 1,numNodesRankPar
         aux_array_r8(iNodeL) = real(coordPar(iNodeL,3),8)
      end do
      ms_offset2d(1) = 2
      call write_dataspace_fp64_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims2d,ms_offset2d,aux_array_r8)

      deallocate(aux_array_r8)
      !-----------------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = mpi_size
      ms_rank = 1
      ms_dims(1) = 1
      ms_offset(1) = int(mpi_rank,hssize_t)
      dtype = h5_datatype_int8!H5T_STD_I64LE
      allocate(aux_array_i8(1))

      dsetname = '/VTKHDF/NumberOfPoints'
      aux_array_i8(1) = numNodesRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfCells'
      aux_array_i8(1) = numElemsRankPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      dsetname = '/VTKHDF/NumberOfConnectivityIds'
      aux_array_i8(1) = numElemsRankPar*nnode
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)

      !-----------------------------------------------------------------------------
      deallocate(aux_array_i8)
      allocate(aux_array_i8(numElemsRankPar+1))

      ds_dims(1)   = int(totalNumElements,hsize_t) + int(mpi_size,hsize_t)
      ms_dims(1)   = int(numElemsRankPar,hsize_t)  + 1
      ms_offset(1) = int(rankElemStart,hssize_t)   - 1 + int(mpi_rank,hssize_t)

      dsetname = '/VTKHDF/Offsets'
      aux_array_i8(1) = 0
      do iElemL = 2,(numElemsRankPar+1)
         aux_array_i8(iElemL) = aux_array_i8(iElemL-1)+nnode
      end do

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)
      deallocate(aux_array_i8)
      !-----------------------------------------------------------------------------
      !  SAVING connecParVTK(:)
      allocate(aux_array_i8(numElemsRankPar*nnode))

      ds_dims(1)   = int(totalNumElements,hsize_t)  * int(nnode,hsize_t)
      ms_dims(1)   = int(numElemsRankPar,hsize_t)   * int(nnode,hsize_t)
      ms_offset(1) = int((rankElemStart-1),hssize_t)* int(nnode,hssize_t)

      do ii = 1,numElemsRankPar*nnode
         aux_array_i8(ii) = connecVTK(ii)-1
      end do

      dsetname = '/VTKHDF/Connectivity'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int8_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i8)
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
      
      call write_dataspace_uint1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

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
      integer(4) :: ds_rank,ms_rank,h5err

      integer(4) :: ii,iElemL
      integer(1),allocatable :: aux_array_i1(:)

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

      !------------------------------------------------------------------------------------------------------
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      ! ## RHO ##
      dsetname = '/VTKHDF/PointData/rho'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,rho)

      ! ## PR ##
      dsetname = '/VTKHDF/PointData/pr'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,pr)

      ! ## Ener ##
      dsetname = '/VTKHDF/PointData/Ener'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,E)
     
      ! ## eta ##
      dsetname = '/VTKHDF/PointData/eta'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,eta)

      ! ## csound ##
      dsetname = '/VTKHDF/PointData/csound'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,csound)

      ! ## machno ##
      dsetname = '/VTKHDF/PointData/machno'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,machno)

      ! ## divU ##
      dsetname = '/VTKHDF/PointData/divU'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,divU)

      ! ## Qcrit ##
      dsetname = '/VTKHDF/PointData/Qcrit'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,Qcrit)

      ! ## envit ##
      dsetname = '/VTKHDF/PointData/envit'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,envit)

      ! ## mut ##
      dsetname = '/VTKHDF/PointData/mut'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mut)

      ! ## mu_fluid ##
      dsetname = '/VTKHDF/PointData/mu_fluid'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,mu_fluid)

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
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,u)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/gradRho'
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,gradRho)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/curlU'
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,curlU)
      !-------------------------------------------------------------------------------- 

      !-------------------------------------------------------------------------------------------------------

      allocate(aux_array_i1(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = h5_datatype_uint1
      ds_rank = 1
      ms_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call write_dataspace_uint1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_instResultsFile

   subroutine write_vtkhdf_avgResultsFile(full_fileName,avrho,avpre,avmueff,avvel,avve2,avvex,avtw)
      implicit none
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: avrho,avpre,avmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: avvel,avve2,avvex,avtw

      character(512) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSIZE_T), dimension(2) :: ds_dims2d,ms_dims2d
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer(HSSIZE_T), dimension(2) :: ms_offset2d
      integer(4) :: ds_rank,ms_rank,h5err

      integer(4) :: iElemL
      integer(1),allocatable :: aux_array_i1(:)

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

      !------------------------------------------------------------------------------------------------------
      ds_rank = 1
      ms_rank = 1
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      ! ## avgRHO ##
      dsetname = '/VTKHDF/PointData/avrho'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avrho)

      ! ## avgPRE ##
      dsetname = '/VTKHDF/PointData/avpre'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avpre)

      ! ## avgMUEFF ##
      dsetname = '/VTKHDF/PointData/avmueff' 
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,avmueff)
      
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
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avvel)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/avve2'
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avve2)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/avvex'
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avvex)
      !-------------------------------------------------------------------------------- 
      dsetname = '/VTKHDF/PointData/avtw'
      call save_array2D_in_dataset_hdf5_file(file_id,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,avtw)
      !-------------------------------------------------------------------------------- 

      !-------------------------------------------------------------------------------------------------------
      allocate(aux_array_i1(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype = h5_datatype_uint1
      ds_rank = 1
      ms_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call write_dataspace_uint1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_avgResultsFile

   subroutine write_vtkhdf_realFieldFile(full_fileName,realField)
      implicit none
      character(512),intent(in) :: full_fileName
      real(rp),dimension(numNodesRankPar),intent(in) :: realField

      character(512) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err

      integer(4) :: iElemL
      integer(1),allocatable :: aux_array_i1(:)

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

      !------------------------------------------------------------------------------------------------------
      ds_rank = 1
      ms_rank = 1
      ds_dims(1)   = int(totalNumNodesPar,hsize_t)
      ms_dims(1)   = int(numNodesRankPar ,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1

      ! ## realField ##      
      dsetname = '/VTKHDF/PointData/realField'
      call save_array1D_in_dataset_hdf5_file(file_id,dsetname,ds_dims,ms_dims,ms_offset,realField)

      !-------------------------------------------------------------------------------------------------------
      allocate(aux_array_i1(numNodesRankPar))
      !------------------------------------------------------------------------------------------------------
      dtype =  h5_datatype_uint1
      ds_rank = 1
      ms_rank = 1
      ds_dims(1)   = int(totalNumElements,hsize_t)
      ms_dims(1)   = int(numElemsRankPar ,hsize_t)
      ms_offset(1) = int(rankElemStart,hssize_t)-1

      ! ## mpi_rank ##
      dsetname = '/VTKHDF/CellData/mpi_rank'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do iElemL = 1,numElemsRankPar
         aux_array_i1(iElemL) = mpi_rank
      end do
      call write_dataspace_uint1_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array_i1)

      deallocate(aux_array_i1)

      !close the file.
      call h5fclose_f(file_id,h5err)

   end subroutine write_vtkhdf_realFieldFile

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

   subroutine save_vtkhdf_avgResultsFile(iStep,avrho,avpre,avmueff,avvel,avve2,avvex,avtw)
      implicit none
      integer, intent(in) :: iStep
      real(rp),dimension(numNodesRankPar),intent(in) :: avrho,avpre,avmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: avvel,avve2,avvex,avtw
      character(512) :: full_fileName

      call set_vtkhdf_avgResultsFile_name(iStep,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF avg results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_avgResultsFile(full_fileName,avrho,avpre,avmueff,avvel,avve2,avvex,avtw)

   end subroutine save_vtkhdf_avgResultsFile

   subroutine save_vtkhdf_finalAvgResultsFile(favrho,favpre,favmueff,favvel,favve2,favvex,favtw)
      implicit none
      real(rp),dimension(numNodesRankPar),intent(in) :: favrho,favpre,favmueff
      real(rp),dimension(numNodesRankPar,ndime),intent(in) :: favvel,favve2,favvex,favtw
      character(512) :: full_fileName

      call set_vtkhdf_finalAvgResultsFile_name(full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF final avg results file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_avgResultsFile(full_fileName,favrho,favpre,favmueff,favvel,favve2,favvex,favtw)

   end subroutine save_vtkhdf_finalAvgResultsFile

   subroutine save_vtkhdf_realFieldFile(realField)
      implicit none
      real(rp),dimension(numNodesRankPar),intent(in) :: realField
      character(512) :: full_fileName

      !------------------------------------------------------------------------------------
      call set_vtkhdf_resultsFile_name(0,full_fileName)
      if(mpi_rank.eq.0) write(*,*) '# Saving VTKHDF realField file: ',trim(adjustl(full_fileName))

      call write_vtkhdf_realFieldFile(full_fileName,realField)

   end subroutine save_vtkhdf_realFieldFile

!---------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------

!-------------------------------WITNESS POINTS-------------------------------!
   subroutine create_witness_hdf5(full_fileName, xyz, witel, witxi, shapewit, nwit, nwitPar, witGlob, save_u_i, save_pr, save_rho)
      implicit none
      character(512), intent(in) :: full_fileName
      integer(rp),    intent(in) :: nwit, nwitPar
      integer(rp),    intent(in) :: witel(nwit), witGlob(nwit)
      real(rp),       intent(in) :: witxi(nwit, ndime), shapewit(nwit,nnode)
      real(rp),       intent(in) :: xyz(nwit,ndime)
      logical,        intent(in) :: save_u_i, save_pr, save_rho
      integer(rp)                :: aux(1), nwitParAllRanks(mpi_size), nwitOffset=0, inode
      integer(hid_t)             :: file_id,plist_id,dset_id,dspace_id,group_id, dtype
      integer(HSIZE_T)           :: ds_dims(2), ms_dims(2), max_dims(2), chunk_dims(2)
      integer(HSSIZE_T)          :: ms_offset(2)
      integer                    :: ds_rank, ms_rank, h5err, irank, iwit
      character(256)             :: groupname,dsetname
      real(rp)                   :: auxwitx(nwitPar), auxwity(nwitPar), auxwitz(nwitPar), auxwitxi1(nwitPar), auxwitxi2(nwitPar), auxwitxi3(nwitPar), auxshapefunc(nwitPar) 

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
      dtype = H5T_NATIVE_REAL

      !Create dataspece for nwitPar and save it!
      ds_rank      = 1
      dsetname     = 'nwitPar'
      ds_dims(1)   = mpi_size
      ms_rank      = 1
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      aux(1)       = nwitPar
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux)
      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

      !Compute sum of nwitPar until that rank and save it on its dataspace!
      ms_rank      = 1
      ms_dims(1)   = mpi_size
      ms_offset(1) = 0
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,nwitParAllRanks)
      if (mpi_rank > 0) then
         do irank = 1, mpi_rank
            nwitOffset = nwitOffset + nwitParAllRanks(irank)
         end do 
      end if
      ds_rank      = 1
      dsetname     = 'nwitOffset'
      ds_dims(1)   = mpi_size
      ms_rank      = 1
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      aux(1)       = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux)
      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

      !Create dataspece for global numeration and save it!
      dsetname     = 'global'
      ds_rank      = 1
      ds_dims(1)   = nwit
      ms_rank      = 1
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      !$acc kernels
      auxwitxi1(:) = witxi(:,1)
      auxwitxi2(:) = witxi(:,2)
      auxwitxi3(:) = witxi(:,3)
      auxwitx(:)   = xyz(:,1)
      auxwity(:)   = xyz(:,2)
      auxwitz(:)   = xyz(:,3)
      !$acc end kernels

      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,witGlob)

      !Create dataspece for element containing the witness and save it!
      dsetname     = 'element'
      ds_rank      = 1
      ds_dims(1)   = nwit
      ms_rank      = 1
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,witel)
      
      !Create dataspace for witness coordinates and save them!
      ds_rank      = 2
      dsetname     = 'xyz'
      ds_dims(1)   = ndime
      ds_dims(2)   = nwit
      ms_rank      = 2
      ms_dims(1)   = ndime
      ms_dims(2)   = nwitPar
      ms_offset(2) = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      ms_offset(1) = 0
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwitx)
      ms_offset(1) = 1
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwity)
      ms_offset(1) = 2
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwitz)


      !Create dataspace for witness isoparametric coordinates and save them!
      ds_rank      = 2
      dsetname     = 'witxi'
      ds_dims(1)   = ndime
      ds_dims(2)   = nwit
      ms_rank      = 2
      ms_dims(1)   = ndime
      ms_dims(2)   = nwitPar
      ms_offset(2) = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      ms_offset(1) = 0
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwitxi1)
      ms_offset(1) = 1
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwitxi2)
      ms_offset(1) = 2
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwitxi3)

      !Create dataspace for the shape functions evaluated on the witness points and save them!
      ds_rank      = 2
      dsetname     = 'shape_functions'
      ds_dims(1)   = nnode
      ds_dims(2)   = nwit
      ms_rank      = 2
      ms_dims(1)   = nnode
      ms_dims(2)   = nwitPar
      ms_offset(2) = nwitOffset
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      do inode = 1, nnode
         ms_offset(1) = inode-1
         auxshapefunc(:) = shapewit(:,inode)
         call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxshapefunc)
      end do

      !Create time dataset!
      ds_rank       = 1
      dsetname      = 'time'
      ds_dims(1)    = 1
      max_dims(1)   = H5S_UNLIMITED_F
      chunk_dims(1) = 1
      call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)

      !Create istep dataset!
      ds_rank       = 1
      dsetname      = 'istep'
      ds_dims(1)    = 1
      max_dims(1)   = H5S_UNLIMITED_F
      chunk_dims(1) = 1
      call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)

      !Create dataspaces for the magnitudes to save!
      ds_rank       = 2
      ds_dims(1)    = 1
      ds_dims(2)    = nwit
      max_dims(1)   = H5S_UNLIMITED_F
      max_dims(2)   = nwit
      chunk_dims(1) = 1
      chunk_dims(2) = nwit
      if (save_u_i) then
         dsetname = 'u_x'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)
         
         dsetname = 'u_y'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)

         dsetname = 'u_z'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)
      end if

      if (save_pr) then
         dsetname = 'pr'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)
      end if

      if (save_rho) then
         dsetname = 'rho'
         call create_dataspace_maxdims_hdf5(file_id,dsetname,ds_rank,ds_dims,max_dims,chunk_dims,dtype)
      end if

      call h5fclose_f(file_id,h5err)

   end subroutine create_witness_hdf5

   subroutine load_witness_hdf5(full_fileName, nwit, loadstep, load_stepwit, nwitPar, witel, witxi, shapefunc) 
      implicit none
      character(512), intent(in)  :: full_fileName
      integer(rp),    intent(in)  :: nwit, loadstep
      integer(rp),    intent(out) :: witel(nwit)
      real(rp),       intent(out) :: witxi(nwit,ndime), shapefunc(nwit,nnode)!, t
      integer(rp),    intent(out) :: nwitPar, load_stepwit
      integer(hid_t)              :: file_id,plist_id,dset_id,dspace_id,group_id, dtype
      integer(HSIZE_T)            :: ms_dims(2), max_dims(2)
      integer(HSSIZE_T)           :: ms_offset(2)
      integer                     :: ms_rank, h5err, iwit, istep
      integer(hsize_t)            :: nsteps(2), maxnsteps(2)
      character(256)              :: groupname,dsetname
      integer(rp)                 :: nwitOffset, auxread(1)
      real(rp)                    :: auxwitxi(ndime,nwit), auxshapefunc(nnode, nwit)!, auxt(1)
      integer(rp), allocatable    :: steps(:)

      witel(:)   = 0
      witxi(:,:) = 0.0_rp
      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      dtype = H5T_NATIVE_REAL

      !Read step to continue!
      dsetname     = 'istep'
      ms_rank      = 1
      call get_dims(file_id, dsetname, ms_rank, nsteps, maxnsteps)
      ms_dims(1)   = nsteps(1)
      ms_offset(1) = 0
      allocate(steps(nsteps(1)))
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,steps)
      do istep = nsteps(1), 1, -1
         if (steps(istep) < loadstep) then
            load_stepwit = istep+1
            exit
         end if
      end do
      write(*,*) load_stepwit
      deallocate(steps)
      
      !Read nwitPar!
      dsetname     = 'nwitPar'
      ms_rank      = 1
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxread)
      nwitPar = auxread(1)

      !Read nwitOffset!
      dsetname     = 'nwitOffset'
      ms_rank      = 1
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxread)
      nwitOffset = auxread(1)

      !Read elements containing the witness points
      dsetname     = 'element'
      ms_rank      = 1
      ms_dims(1)   = nwitPar
      ms_offset(1) = nwitOffset
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,witel)
      
      !Read witness isoparametric coordinates!
      dsetname     = 'witxi'
      ms_rank      = 2
      ms_dims(1)   = ndime
      ms_dims(2)   = nwitPar
      ms_offset(1) = 0
      ms_offset(2) = nwitOffset
      call read_dataspace_2d_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwitxi)

      !Read the shape functions coordinates!
      dsetname     = 'shape_functions'
      ms_rank      = 2
      ms_dims(1)   = nnode
      ms_dims(2)   = nwitPar
      ms_offset(1) = 0
      ms_offset(2) = nwitOffset
      call read_dataspace_2d_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxshapefunc)
      do iwit =1,nwitPar
         witxi(iwit,:) = auxwitxi(:,iwit)
         shapefunc(iwit,:) = auxshapefunc(:,iwit)
      end do

      call h5fclose_f(file_id,h5err)

   end subroutine load_witness_hdf5

   subroutine update_witness_hdf5(itewit, leapwitsave, witval, nwit, nwitPar, nvarwit, full_fileName, t, steps, save_u_i, save_pr, save_rho)
      integer(4), intent(in)     :: itewit, nwit, nwitPar, nvarwit, leapwitsave, steps(leapwitsave)
      real(rp), intent(in)       :: witval(leapwitsave, nwitPar, nvarwit), t(leapwitsave) 
      logical, intent(in)        :: save_u_i, save_pr, save_rho
      character(512), intent(in) :: full_fileName
      character(256)             :: dsetname
      real(rp)                   :: auxwrite(nwitPar)
      integer(HSSIZE_T)          :: ms_offset(2)
      integer                    :: ms_rank,h5err, iwit, ds_rank, ileap
      integer(4)                 :: nwitOffset, auxread(1)
      integer(HSIZE_T)           :: ms_dims(2), ds_dims(2)
      integer(hid_t)             :: file_id,plist_id, dtype

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      call h5fopen_f(full_fileName, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot load results file ',trim(adjustl(full_fileName))
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

      !Read nwitOffset!
      dsetname     = 'nwitOffset'
      ms_rank      = 1
      ms_dims(1)   = 1
      ms_offset(1) = mpi_rank
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxread)
      nwitOffset = auxread(1)

      !Save variables!
      ms_rank      = 2
      ms_dims(1)   = 1
      ms_dims(2)   = nwitPar
      ds_rank      = 2
      ds_dims(1)   = itewit
      ds_dims(2)   = nwit
      ms_offset(2) = nwitOffset    
      if (save_u_i) then
         dsetname = 'u_x'
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
         dsetname = 'u_y'
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
         dsetname = 'u_z'
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      end if
      if (save_pr) then
         dsetname = 'pr'
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      end if
      if (save_rho) then
         dsetname = 'rho'
         call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      end if
      do ileap = 1, leapwitsave
         ms_offset(1) = itewit - leapwitsave + ileap - 1
         if (save_u_i) then
            dsetname = 'u_x'
            !$acc kernels
            auxwrite(:) = witval(ileap,:,1)
            !$acc end kernels
            call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwrite)
            dsetname = 'u_y'  
            !$acc kernels
            auxwrite(:) = witval(ileap,:,2)
            !$acc end kernels
            call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwrite)
            dsetname = 'u_z'  
            !$acc kernels
            auxwrite(:) = witval(ileap,:,3)
            !$acc end kernels
            call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwrite)
         end if
         if (save_pr) then
            dsetname = 'pr'  
            !$acc kernels
            auxwrite(:) = witval(ileap,:,4)
            !$acc end kernels
            call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwrite)
         end if
         if (save_rho) then
            dsetname = 'rho'  
            !$acc kernels
            auxwrite(:) = witval(ileap,:,5)
            !$acc end kernels
            call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,auxwrite)
         end if
      end do   

      !Save time!
      dsetname     = 'time'
      ms_rank      = 1
      ms_dims(1)   = leapwitsave
      ms_offset(1) = itewit - leapwitsave
      ds_rank      = 1
      ds_dims(1)   = itewit
      call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call write_dataspace_fp32_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,t)

      !Save istep!
      dsetname     = 'istep'
      ms_rank      = 1
      ms_dims(1)   = leapwitsave
      ms_offset(1) = itewit - leapwitsave
      ds_rank      = 1
      ds_dims(1)   = itewit
      call extend_dataset_hdf5(file_id,dsetname,ds_rank,ds_dims)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,steps)
      call h5fclose_f(file_id,h5err)
   
   end subroutine update_witness_hdf5

end module mod_hdf5