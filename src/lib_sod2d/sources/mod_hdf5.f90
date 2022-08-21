module mod_hdf5
   use hdf5
   use mod_mpi
   use mod_mpi_mesh
   implicit none

   character(256) :: hdf5_mesh_name
   integer(hid_t) :: mesh_file_id

contains

   subroutine set_hdf5_mesh_file_name(file_path,file_name)
      implicit none
      character(len=*), intent(in) :: file_path,file_name
      character(len=12) :: aux_mpiSize

      write(aux_mpiSize,'(I0)') mpi_size
      hdf5_mesh_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'-'//trim(aux_mpiSize)//'.h5'
   end subroutine

   subroutine read_alyaMesh_part_and_create_hdf5Mesh(gmsh_file_path,gmsh_file_name,periodic,&
         mesh_h5_file_path,mesh_h5_file_name)
      implicit none     
      character(len=*), intent(in) :: gmsh_file_path,gmsh_file_name
      logical :: periodic
      character(len=*), intent(in) :: mesh_h5_file_path,mesh_h5_file_name

      !-- read the alya mesh fesh files in GMSH/ALYA FORMAT
      call read_alya_mesh_files(gmsh_file_path,gmsh_file_name,periodic)

      !----- Do mesh partitioning!
      call do_mesh_partitioning()

      !----- Create HDF5 File
      call create_hdf5_meshfile(mesh_h5_file_path,mesh_h5_file_name)

   end subroutine read_alyaMesh_part_and_create_hdf5Mesh

   subroutine create_hdf5_meshfile(file_path,file_name)
      implicit none     
      character(len=*), intent(in) :: file_path,file_name
      integer(hid_t) :: file_id,plist_id,dset_id,dspace_id,mspace_id,group_id
      integer(hid_t) :: dtype
      integer(HSIZE_T), dimension(1) :: ds_dims,ms_dims
      integer(HSSIZE_T), dimension(1) :: ms_offset 
      integer :: ds_rank,ms_rank,h5err
      character(256) :: groupname,dsetname

      call set_hdf5_mesh_file_name(file_path,file_name)

      !.init h5 interface
      call h5open_f(h5err)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      ! create file collectively
      call h5fcreate_f(hdf5_mesh_name,H5F_ACC_TRUNC_F,file_id,h5err,access_prp=plist_id)
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
      ms_dims(1) = numElemsInRank
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
      if(totalNumBoundsSrl.ne.0) then
         call save_boundary_data_hdf5(file_id)
      end if
      !--------------------------------------------------------------------------------

      !close the file.
      call h5fclose_f(file_id,h5err)

      !close h5 interface
      call h5close_f(h5err)

      !--------------------------------------------------------------------------------


      ! Create the data space for the  dataset. 
      !call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      ! Create the dataset with default properties.

      !call h5dcreate_f(file_id,dsetname,H5T_NATIVE_REAL,dspace_id,dset_id,h5err)

      !call h5sclose_f(dspace_id,h5err)
      !call h5dclose_f(dset_id,h5err)

      ! Create the data space for the  dataset. 
      !call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      ! Create the dataset with default properties.

      !call h5dcreate_f(file_id,dsetname,H5T_NATIVE_REAL,dspace_id,dset_id,h5err)

      !call h5sclose_f(dspace_id,h5err)
      !call h5dclose_f(dset_id,h5err)

      ! Create the data space for the  dataset. 
      !call h5screate_simple_f(ds_rank,ds_dims,dspace_id,h5err)
      ! Create the dataset with default properties.

      !call h5dcreate_f(file_id,dsetname,H5T_NATIVE_REAL,dspace_id,dset_id,h5err)

      !call h5sclose_f(dspace_id,h5err)
      !call h5dclose_f(dset_id,h5err)

      !--------------------------------------------------------------------------------
#if 0
      ms_dims(1) = numNodesRankPar! ds_dims(1)/mpi_size
      ms_offset(1) = rankNodeStart-1! ms_dims(1)*mpi_rank
      write(*,*) 'writing -> ',totalNumNodesPar,ms_dims(1),ms_offset(1)
      CALL h5screate_simple_f(ds_rank,ms_dims,mspace_id,h5err) 

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id,dspace_id,h5err)
      CALL h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      CALL h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      ! Write the dataset collectively. 
      CALL h5dwrite_f(dset_id,H5T_NATIVE_REAL,coord_x,ds_dims,h5err, &
                      file_space_id = dspace_id, mem_space_id = mspace_id, xfer_prp = plist_id)

     ! Write the dataset collectively. 
     !CALL h5dwrite_f(dset_id,H5T_NATIVE_REAL,dset,ds_dims,h5err,xfer_prp=plist_id)
     !
     CALL h5sclose_f(dspace_id,h5err)
     call h5sclose_f(mspace_id,h5err)
     CALL h5dclose_f(dset_id,h5err)
      call h5pclose_f(plist_id,h5err)
#endif

   end subroutine create_hdf5_meshfile

   subroutine load_hdf5_meshfile(file_path,file_name)
      implicit none     
      character(len=*), intent(in) :: file_path,file_name
      character(256) :: groupname,dsetname
      integer(hid_t) :: file_id,plist_id,dset_id,fspace_id
      integer :: h5err
      integer(hsize_t),dimension(1) :: fs_dims,fs_maxdims

      call set_hdf5_mesh_file_name(file_path,file_name)

      !.init h5 interface
      call h5open_f(h5err)

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,h5err)

      if(mpi_rank.eq.0) write(*,*) '# Loading hdf5 mesh: ',hdf5_mesh_name

      call h5fopen_f(hdf5_mesh_name, H5F_ACC_RDWR_F,file_id,h5err,access_prp=plist_id)
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
      ms_dims(1) = numElemsInRank
      ms_offset(1) = rankElemStart-1

      dsetname = '/globalIds/elemGid'
#endif

      !-----------------------------------------------------------------------------------------------
      !load the parallel data
      if(mpi_size.ge.2) then
         call load_parallel_data_hdf5(file_id)
      else !only 1 rank
         numNodesRankPar = totalNumNodesPar
         numElemsInRank = totalNumElements
         rankNodeStart = 1
         rankNodeEnd = numNodesRankPar
         rankElemStart = 1
         rankElemEnd = numElemsInRank
      end if
      !-----------------------------------------------------------------------------------------------
      !load periodic data
      call load_periodic_data_hdf5(file_id)

      !-----------------------------------------------------------------------------------------------
      !load boundary data
      call load_boundary_data_hdf5(file_id)

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
      call h5close_f(h5err)

   end subroutine load_hdf5_meshfile


   subroutine close_hdf5_meshfile()
      implicit none
      integer :: h5err

      !close the file.
      call h5fclose_f(mesh_file_id,h5err)

      !close h5 interface
      call h5close_f(h5err)
   end subroutine close_hdf5_meshfile

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
      !write(*,*) 'dsetname ',dsetname, ' dset_id ',dset_id,' fspace_id ',fspace_id,' t_dims ',t_dims

      ! Each process defines dataset in memory and writes it to the hyperslab in the file. 
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) 

      !write(*,*) 'dsetname ',dsetname, ' dset_id ',dset_id,' dspace_id ',dspace_id
      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dwrite_f(dset_id,dtype,data,fs_dims,h5err,&
                      file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)


      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
   end subroutine write_dataspace_fp32_hyperslab_parallel

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
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

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
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine read_dataspace_int4_hyperslab_parallel

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
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5dread_f(dset_id,dtype,data,fs_dims,h5err,&
                     file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

   end subroutine read_dataspace_fp32_hyperslab_parallel

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
      ms_dims(1) = numElemsInRank*nnode
      ms_offset(1) = (rankElemStart-1)*nnode

      allocate(aux_array(numElemsInRank*nnode))
      !-----------------------------------------------------------------------------------------------------
      !  SAVING connecParOrig(:,:)
      i=1
      do iElemL=1,numElemsInRank
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
      do iElemL=1,numElemsInRank
         do m=1,nnode
            aux_array(i)=connecParWork(iElemL,m)
            i=i+1
         end do
      end do

      dsetname = '/Connectivity/connecParWork'
      call create_dataspace_hdf5(file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)
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
      ms_dims(1) = numElemsInRank*nnode
      ms_offset(1) = (rankElemStart-1)*nnode

      !-------------------------------------------------------------------------------------------------------
      allocate( connecParOrig(numElemsInRank,nnode) )
      allocate( connecParWork(numElemsInRank,nnode) )
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecParOrig(:,:)
      allocate(aux_array(numElemsInRank*nnode))

      dsetname = '/Connectivity/connecParOrig'
      call read_dataspace_int4_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,aux_array)

      i=1
      do iElemL=1,numElemsInRank
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
      do iElemL=1,numElemsInRank
         do m=1,nnode
            connecParWork(iElemL,m)=aux_array(i) !it contains the iNodeL
            i=i+1
         end do
      end do

      deallocate(aux_array)
      !-------------------------------------------------------------------------------------------------------
      !LOADING connecCGNS(:)
      allocate( connecCGNS(numElemsInRank*nnode) )

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
      numElemsInRank  = rankElemEnd - rankElemStart + 1

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

      if(isBoundaryFolder) then
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
      allocate(elemGid(numElemsInRank))

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
      ms_dims(1) = numElemsInRank
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

!--------------------------------------------------------------------------------------
!------------- deprecated shit --------------------------------------------------------
!--------------------------------------------------------------------------------------

    subroutine create_hdf5_file()
        CHARACTER(LEN=8), PARAMETER :: filename = "dsetf.h5" ! File name
        CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name

        INTEGER(HID_T) :: file_id       ! File identifier
        INTEGER(HID_T) :: dset_id       ! Dataset identifier
        INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


        INTEGER(HSIZE_T), DIMENSION(2) :: dims = (/4,6/) ! Dataset dimensions
        INTEGER     ::   rank = 2                        ! Dataset rank

        INTEGER     ::   error ! Error flag

        !
        ! Initialize FORTRAN interface.
        !
        CALL h5open_f(error)

        !
        ! Create a new file using default properties.
        !
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        !
        ! Create the dataspace.
        !
        CALL h5screate_simple_f(rank, dims, dspace_id, error)

        !
        ! Create the dataset with default properties.
        !
        CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
             dset_id, error)

        !
        ! End access to the dataset and release resources used by it.
        !
        CALL h5dclose_f(dset_id, error)

        !
        ! Terminate access to the data space.
        !
        CALL h5sclose_f(dspace_id, error)

        !
        ! Close the file.
        !
        CALL h5fclose_f(file_id, error)

        !
        ! Close FORTRAN interface.
        !
        CALL h5close_f(error)


    end subroutine


end module mod_hdf5
