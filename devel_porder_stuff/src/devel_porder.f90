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
   real(rp), allocatable :: test_array(:)
   integer(4) :: i,j,k,p_order,pIndex
   !-----------------------------------

   call init_mpi()
   if(mpi_rank.eq.0) write(*,*) 'testing stuff for porder!'

   call init_hdf5_interface()

   p_order = 4

   do i=0,(p_order)
      do j=0,(p_order)
         do k=0,(p_order)
            call vtk_pointIndexFromIJK(i,j,k,p_order,pIndex)
            write(*,*) i,j,k,pIndex
         enddo
      end do
   end do



#if 0
    write(mesh_h5_file_path,*) ""
    write(mesh_h5_file_name,*) "cube_per10"

    write(results_h5_file_path,*) ""
    write(results_h5_file_name,*) "results"

    numIters = 5

    call set_hdf5_meshFile_name(mesh_h5_file_path,mesh_h5_file_name,mpi_size)
    call set_hdf5_baseResultsFile_name(results_h5_file_path,results_h5_file_name,mesh_h5_file_name,mpi_size)

    useIntInComms=.true.
    useRealInComms=.true.

    call load_hdf5_meshfile()
    ! init comms
    call init_comms(useIntInComms,useRealInComms)
    call init_comms_bnd(useIntInComms,useRealInComms)

    allocate(test_array(numNodesRankPar))
    !$acc enter data create(test_array(:))

    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

    if(mpi_rank.eq.0) write(*,*) 'testing new way'
    call nvtxStartRange("new_way")
    do iter=1,numIters

        call nvtxStartRange("new_way_op")
        !$acc parallel loop present(test_array(:))
        do iNode=1,numNodesRankPar
           test_array(iNode) = mpi_rank + 0.5
        end do
        !$acc end parallel loop
        call nvtxEndRange

        call nvtxStartRange("new_way_comm")
        call mpi_halo_atomic_update_real(test_array)
        !call mpi_halo_atomic_update_real_iSendiRcv_devel(test_array)
        call nvtxEndRange
    end do
    call nvtxEndRange

    if(mpi_rank.eq.0) write(*,*) 'test_array(1)',test_array(:)

    !$acc update host(test_array(:))
    if(mpi_rank.eq.0) write(*,*) 'test_array(2)',test_array(:)
    !do i=1,numNodesToComm
    !    iNodeL = ndevel_porderodesToComm(i)
    !    write(*,*) 'test_array(',iNodeL,')'!realField(iNodeL) = realField(iNodeL) + aux_realField_r(i)
    !end do

    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

    if(mpi_rank.eq.0) write(*,*) 'testing old way'
    call nvtxStartRange("old_way")
    do iter=1,numIters

        call nvtxStartRange("old_way_op")
        !$acc parallel loop present(test_array(:))
        do iNode=1,numNodesRankPar
           test_array(iNode) = mpi_rank + 0.5
        end do
        !$acc end parallel loop
        call nvtxEndRange

        call nvtxStartRange("old_way_comm")
        call mpi_halo_atomic_update_real_iSendiRcv(test_array)
        call nvtxEndRange

    end do
    call nvtxEndRange

    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

    !$acc exit data delete(test_array(:))
    deallocate(test_array)

    !call test_mpi_cudaware(100,10)

    ! End comms
    call end_comms()
    call end_comms_bnd()
#endif

    ! End hdf5 interface
   call end_hdf5_interface()

   call end_mpi()

end program devel_porder
