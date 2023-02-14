
program tool_commsPerfomance
    use mod_mpi
    use mod_read_inputFile
    use mod_mpi_mesh
    use mod_comms
    use mod_comms_performance
    use mod_hdf5
    implicit none

    character(999) :: input_file
    integer :: lineCnt
    integer :: numNodesSrl,numNodesB_1r,numIters
    character(256) :: parameter2read
    character(512) :: mesh_h5_file_path,mesh_h5_file_name
    character(512) :: results_h5_file_path,results_h5_file_name

    logical :: useMesh=.false.
    logical :: useIntInComms=.false.,useFloatInComms=.false.,useDoubleInComms=.false.

    call init_mpi()

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
   call end_mpi()

end program tool_commsPerfomance
