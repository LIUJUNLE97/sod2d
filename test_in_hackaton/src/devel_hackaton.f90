
program devel_hackaton
    use mod_mpi
    use mod_read_inputFile
    use mod_mpi_mesh
    use mod_comms
    use mod_hdf5
    use mod_nvtx
    use mod_test_funcs
    !use openacc
    !use cudafor

    implicit none

    character(512) :: mesh_h5_file_path,mesh_h5_file_name,results_h5_file_path,results_h5_file_name
    logical :: useIntInComms,useRealInComms
    real(rp), allocatable :: test_array(:)
    integer(4) :: iter,numIters,iNode
    !-----------------------------------


    call init_mpi()
    if(mpi_rank.eq.0) write(*,*) 'testing stuff for the hackaton!'

    call init_hdf5_interface()

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
    !    iNodeL = nodesToComm(i)
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


    ! End hdf5 interface
    call end_hdf5_interface()

    ! End comms
    call end_comms()
    call end_comms_bnd()

   call end_mpi()

#if 0
    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, input_file)
        if(mpi_rank.eq.0) write(*,*) '# Input file: ',trim(adjustl(input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call Tool COMMS PERFROMANCE with an input file!'
        call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
    ! Reading the parameters
    call open_inputFile(input_file)
    lineCnt = 1

    !1. numIters--------------------------------------------------------------------------
    parameter2read = 'numIters'
    call read_inputFile_integer(lineCnt,parameter2read,numIters)

    !2. Use INT in COMMS--------------------------------------------------------------------------
    parameter2read = 'useIntInComms'
    call read_inputFile_logical(lineCnt,parameter2read,useIntInComms)

    !3. Use FLOAT in COMMS--------------------------------------------------------------------------
    parameter2read = 'useFloatInComms'
    call read_inputFile_logical(lineCnt,parameter2read,useFloatInComms)

    !4. Use DOUBLE in COMMS--------------------------------------------------------------------------
    parameter2read = 'useDoubleInComms'
    call read_inputFile_logical(lineCnt,parameter2read,useDoubleInComms)

    !5. Use mesh--------------------------------------------------------------------------
    parameter2read = 'useMesh'
    call read_inputFile_logical(lineCnt,parameter2read,useMesh)

    if(useMesh) then

        !6. mesh_h5_file_path--------------------------------------------------------------
        parameter2read = 'mesh_h5_file_path'
        call read_inputFile_string(lineCnt,parameter2read,mesh_h5_file_path)

        !7. mesh_h5_file_name--------------------------------------------------------------
        parameter2read = 'mesh_h5_file_name'
        call read_inputFile_string(lineCnt,parameter2read,mesh_h5_file_name)

    else

        !6. numNodesSrl--------------------------------------------------------------------------
        parameter2read = 'numNodesSrl'
        call read_inputFile_integer(lineCnt,parameter2read,numNodesSrl)

        !7. numNodesB_1r --------------------------------------------------------------------------
        parameter2read = 'numNodesB_1r'
        call read_inputFile_integer(lineCnt,parameter2read,numNodesB_1r)

       if(mpi_rank.eq.0) write(*,*) 'numNodesSrl ',numNodesSrl,' numNodesB_1r ',numNodesB_1r,' numIters ',numIters

    end if
     
    call close_inputFile() 

    if(useMesh) then

        !write(results_h5_file_path,*) ""
        !write(results_h5_file_name,*) "dummy"

        call init_hdf5_interface()
        call set_hdf5_meshFile_name(mesh_h5_file_path,mesh_h5_file_name,mpi_size)
        call load_hdf5_meshfile()
    else    
        call create_dummy_1Dmesh(numNodesSrl,numNodesB_1r)
    end if

   call init_comms_performance(useIntInComms,useFloatInComms,useDoubleInComms)

    !call debug_comms_float()

   call test_comms_performance_float(numIters)


   !call nvtxStartRange("saxpy loop")
   !call do_saxpy_loop(numIters)
   !call nvtxEndRange
   !call test_mpi_cudaware(numNodesSrl,numIters)
   !call do_crazy_mpi_test(numNodesSrl,numIters)

   call end_comms()
#endif


end program devel_hackaton