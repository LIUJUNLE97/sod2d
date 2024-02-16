
program tool_copy_results
    use mod_mpi
    use mod_read_inputFile
    use mod_mpi_mesh
    use mod_comms
    use mod_hdf5
    use mod_copy_results
    implicit none

    character(999) :: input_file
    integer(4) :: lineCnt
    character(256) :: parameter2read
    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    character(512) :: results_h5_filePath,results_h5_fileName
    integer(4) :: results_first,results_last,results_step,target_Nprocs

    logical :: useMesh=.false.
    logical :: useIntInComms=.false.,useRealInComms=.false.

    call init_mpi()

    if(mpi_rank.eq.0) then
        write(*,*) '|-- WELCOME TO SOD2D COPY RESULTS TOOL!'
        write(*,*) '|-- Copy results in same mesh and different number of partitions'
    end if

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, input_file)
        if(mpi_rank.eq.0) write(*,*) '# Input file: ',trim(adjustl(input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call tool Copy Results with an input file!'
        call MPI_Abort(app_comm,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
    ! Reading the parameters
    call open_inputFile(input_file)
    lineCnt = 1

    !1. mesh_h5_filePath--------------------------------------------------------------
    parameter2read = 'mesh_h5_filePath'
    call read_inputFile_string(lineCnt,parameter2read,mesh_h5_filePath)

    !2. mesh_h5_fileName--------------------------------------------------------------
    parameter2read = 'mesh_h5_fileName'
    call read_inputFile_string(lineCnt,parameter2read,mesh_h5_fileName)

    !3. results_h5_filePath-----------------------------------------------------------
    parameter2read = 'results_h5_filePath'
    call read_inputFile_string(lineCnt,parameter2read,results_h5_filePath)

    !4. results_h5_fileName-----------------------------------------------------------
    parameter2read = 'results_h5_fileName'
    call read_inputFile_string(lineCnt,parameter2read,results_h5_fileName)

    !5. target_Nprocs-----------------------------------------------------------------
    parameter2read = 'target_Nprocs'
    call read_inputFile_integer(lineCnt,parameter2read,target_Nprocs)

    !6. results_first----------------------------------------------------------------
    parameter2read = 'results_first'
    call read_inputFile_integer(lineCnt,parameter2read,results_first)

    !7. results_last-----------------------------------------------------------------
    parameter2read = 'results_last'
    call read_inputFile_integer(lineCnt,parameter2read,results_last)

    !8. results_step-----------------------------------------------------------------
    parameter2read = 'results_step'
    call read_inputFile_integer(lineCnt,parameter2read,results_step)

    call close_inputFile()
    if(mpi_rank.eq.0) write(*,*) '## End of Reading input file: ',trim(adjustl(input_file))

!---------------------------------------------------------------------------------------------------------
    call copy_results_same_mesh_Npartitions(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName,target_Nprocs,&
                                            results_first,results_last,results_step)
!---------------------------------------------------------------------------------------------------------

    call end_mpi()

end program tool_copy_results