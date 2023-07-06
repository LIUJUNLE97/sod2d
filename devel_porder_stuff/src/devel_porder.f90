program devel_porder
   use mod_mpi
   use mod_read_inputFile
   use mod_mpi_mesh
   use mod_comms
   use mod_hdf5
   use mod_nvtx
   use mod_ijk_indices

   implicit none

   character(512) :: mesh_h5_file_path,mesh_h5_file_name,results_h5_file_path,results_h5_file_name
   logical :: useIntInComms,useRealInComms
   integer(4) :: i,j,k,pIndex,auxCnt
   real(rp),allocatable :: coordsHex(:,:),coordsQuad(:,:)
   real(rp) :: x,y,z

   !------------------------------------------------------------
   integer(4) :: ds_rank,h5err
   integer(hsize_t),dimension(1) :: ds_dims,ms_dims
   integer(hssize_t),dimension(1) :: ms_offset 
   integer(hsize_t),dimension(2) :: ds_dims2d,ms_dims2d
   integer(hssize_t),dimension(2) :: ms_offset2d

   integer(hid_t) :: hdf5_fileId
   character(512) :: groupname,dsetname,outputfilename
   character(12) :: aux_porder
   integer(hid_t) :: dtype
   integer(1),allocatable :: aux_array_i1(:)
   integer(8),allocatable :: aux_array_i8(:)
   integer(4) :: ii,jj,iBound,iRank,iNodeL
   integer(4) :: mporder,mnnode,mngaus,mnpbou
   integer(4),dimension(:),allocatable :: gmsh2ijk,vtk2ijk,gmsh2ij,vtk2ij

   !------------------------------------------------------------

   call init_mpi()
   call init_hdf5_interface()

   !-----------------------------------------------------------------------------

   mporder = 3
   if(mpi_rank.eq.0) write(*,*) 'testing stuff for porder! mporder',mporder

   call initGmshIJKTables(mporder)

   call get_porder_values(mporder,mnnode,mngaus,mnpbou)
   call set_allocate_hexahedronHO_ijk_indices(mporder,gmsh2ijk,vtk2ijk)
   call set_allocate_quadrilateralHO_ij_indices(mporder,gmsh2ij,vtk2ij)



   allocate(coordsHex(mnnode,ndime))
   allocate(coordsQuad(mnpbou,ndime))

   auxCnt = 1
   do i=0,(mporder)
      x = real(i)/real(mporder)*1.0_rp
      do j=0,(mporder)
         y = real(j)/real(mporder)*1.0_rp
         do k=0,(mporder)
            z = real(k)/real(mporder)*1.0_rp
            call vtkHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,pIndex)
            coordsHex(pIndex,1) = x
            coordsHex(pIndex,2) = y
            coordsHex(pIndex,3) = z
            vtk2ijk(auxCnt) = auxCnt!pIndex
            auxCnt = auxCnt + 1
            !write(*,*) i,j,k,pIndex
         enddo
      end do
   end do

   auxCnt = 1
   z = 0.0_rp
   do i=0,(mporder)
      x = real(i)/real(mporder)*1.0_rp
      do j=0,(mporder)
         y = real(j)/real(mporder)*1.0_rp

         call vtkHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,pIndex)
         coordsQuad(pIndex,1) = x
         coordsQuad(pIndex,2) = y
         coordsQuad(pIndex,3) = z
         vtk2ij(auxCnt) = auxCnt!pIndex
         auxCnt = auxCnt + 1
      end do
   end do



   !-----------------------------------------------------------------------------

   write(aux_porder,'(I0)') mporder
   outputfilename = 'hexaedra_vtkhdf_p'//trim(aux_porder)//'.hdf'
   call create_hdf5_file(outputfilename,hdf5_fileId)
   call set_vtkhdf_attributes_and_basic_groups(hdf5_fileId)

   !-----------------------------------------------------------------------------
   ds_rank = 1
   ds_dims(1) = int(mpi_size,hsize_t)
   ms_dims(1) = 1
   ms_offset(1) = int(mpi_rank,hssize_t)
   dtype = h5_datatype_int8
   allocate(aux_array_i8(1))

   dsetname = '/VTKHDF/NumberOfPoints'
   aux_array_i8(1) = mnnode
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   dsetname = '/VTKHDF/NumberOfCells'
   aux_array_i8(1) = 1
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   dsetname = '/VTKHDF/NumberOfConnectivityIds'
   aux_array_i8(1) = 1*mnnode
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   deallocate(aux_array_i8)
   !-----------------------------------------------------------------------------------
   allocate(aux_array_i8(2))
   ds_dims(1)   = 2
   ms_dims(1)   = 2
   ms_offset(1) = 0

   dsetname = '/VTKHDF/Offsets'
   aux_array_i8(1) = 0
   aux_array_i8(2) = aux_array_i8(1)+mnnode
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   deallocate(aux_array_i8)
   !-----------------------------------------------------------------------------
   allocate(aux_array_i8(mnnode))

   ds_dims(1)   = mnnode
   ms_dims(1)   = mnnode
   ms_offset(1) = 0

   do ii = 1,mnnode
      aux_array_i8(ii) = vtk2ijk(ii)-1
   end do

   dsetname = '/VTKHDF/Connectivity'
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
   deallocate(aux_array_i8)

   !-----------------------------------------------------------------------------
   allocate(aux_array_i1(1))
   dtype = h5_datatype_uint1
   dsetname = '/VTKHDF/Types'
   aux_array_i1(:) = 72

   ds_dims(1)   = 1
   ms_dims(1)   = 1
   ms_offset(1) = 0
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)

   deallocate(aux_array_i1)
   !--------------------------------------------------------------------------------
   ds_dims2d(1) = int(ndime,hsize_t)
   ds_dims2d(2) = int(mnnode,hsize_t)
   ms_dims2d(1) = int(ndime,hsize_t)
   ms_dims2d(2) = int(mnnode,hsize_t)
   ms_offset2d(1) = 0
   ms_offset2d(2) = 0!int(0,hssize_t)-1

!-----------------------------------------------------------------------------
   dsetname = '/VTKHDF/Points'
   call save_array2D_tr_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,coordsHex)

   call close_hdf5_file(hdf5_fileId)

!--------------------------------------------------------------------------------------------------------------------------
   !-----------------------------------------------------------------------------

   outputfilename = 'quadrilateral_vtkhdf_p'//trim(aux_porder)//'.hdf'
   call create_hdf5_file(outputfilename,hdf5_fileId)
   call set_vtkhdf_attributes_and_basic_groups(hdf5_fileId)

   !-----------------------------------------------------------------------------
   ds_rank = 1
   ds_dims(1) = int(mpi_size,hsize_t)
   ms_dims(1) = 1
   ms_offset(1) = int(mpi_rank,hssize_t)
   dtype = h5_datatype_int8
   allocate(aux_array_i8(1))

   dsetname = '/VTKHDF/NumberOfPoints'
   aux_array_i8(1) = mnpbou
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   dsetname = '/VTKHDF/NumberOfCells'
   aux_array_i8(1) = 1
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   dsetname = '/VTKHDF/NumberOfConnectivityIds'
   aux_array_i8(1) = 1*mnpbou
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   deallocate(aux_array_i8)
   !-----------------------------------------------------------------------------------
   allocate(aux_array_i8(2))
   ds_dims(1)   = 2
   ms_dims(1)   = 2
   ms_offset(1) = 0

   dsetname = '/VTKHDF/Offsets'
   aux_array_i8(1) = 0
   aux_array_i8(2) = aux_array_i8(1)+mnpbou
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)

   deallocate(aux_array_i8)
   !-----------------------------------------------------------------------------
   allocate(aux_array_i8(mnpbou))

   ds_dims(1)   = mnpbou
   ms_dims(1)   = mnpbou
   ms_offset(1) = 0

   do ii = 1,mnpbou
      aux_array_i8(ii) = vtk2ij(ii)-1
   end do

   dsetname = '/VTKHDF/Connectivity'
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_int8_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
   deallocate(aux_array_i8)

   !-----------------------------------------------------------------------------
   allocate(aux_array_i1(1))
   dtype = h5_datatype_uint1
   dsetname = '/VTKHDF/Types'
   aux_array_i1(:) = 70

   ds_dims(1)   = 1
   ms_dims(1)   = 1
   ms_offset(1) = 0
   call create_dataspace_hdf5(hdf5_fileId,dsetname,ds_rank,ds_dims,dtype)
   call write_dataspace_1d_uint1_hyperslab_parallel(hdf5_fileId,dsetname,ms_dims,ms_offset,aux_array_i1)
   !--------------------------------------------------------------------------------
   ds_dims2d(1) = int(ndime,hsize_t)
   ds_dims2d(2) = int(mnpbou,hsize_t)
   ms_dims2d(1) = int(ndime,hsize_t)
   ms_dims2d(2) = int(mnpbou,hsize_t)
   ms_offset2d(1) = 0
   ms_offset2d(2) = 0!int(0,hssize_t)-1

   deallocate(aux_array_i1)
!-----------------------------------------------------------------------------
   dsetname = '/VTKHDF/Points'
   call save_array2D_tr_rp_in_dataset_hdf5_file(hdf5_fileId,dsetname,ds_dims2d,ms_dims2d,ms_offset2d,coordsQuad)







   call close_hdf5_file(hdf5_fileId)








!-----------------------------------------------------------------------------------------------------------------------------

   ! End hdf5 interface
   call end_hdf5_interface()

   call end_mpi()

end program devel_porder
