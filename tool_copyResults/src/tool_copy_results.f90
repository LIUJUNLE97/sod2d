
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
    integer(4) :: type_resultsFile,results_first,results_last,results_step,target_Nprocs
    logical :: generateNewMesh=.false.

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

    !3. generate_mesh
    parameter2read = 'generate_mesh'
    call read_inputFile_logical(lineCnt,parameter2read,generateNewMesh)

    !4. results_h5_filePath-----------------------------------------------------------
    parameter2read = 'results_h5_filePath'
    call read_inputFile_string(lineCnt,parameter2read,results_h5_filePath)

    !5. results_h5_fileName-----------------------------------------------------------
    parameter2read = 'results_h5_fileName'
    call read_inputFile_string(lineCnt,parameter2read,results_h5_fileName)

    !6. target_Nprocs-----------------------------------------------------------------
    parameter2read = 'target_Nprocs'
    call read_inputFile_integer(lineCnt,parameter2read,target_Nprocs)

    !7. type_resultsFile-------------------------------------------------------------
    parameter2read = 'type_resultsFile'
    call read_inputFile_integer(lineCnt,parameter2read,type_resultsFile)

    if(type_resultsFile .eq. 1) then
        !8. results_first----------------------------------------------------------------
        parameter2read = 'results_first'
        call read_inputFile_integer(lineCnt,parameter2read,results_first)

        !9. results_last-----------------------------------------------------------------
        parameter2read = 'results_last'
        call read_inputFile_integer(lineCnt,parameter2read,results_last)

        !10. results_step-----------------------------------------------------------------
        parameter2read = 'results_step'
        call read_inputFile_integer(lineCnt,parameter2read,results_step)
    else if(type_resultsFile .eq. 2) then
        !8. resAvg_file----------------------------------------------------------------
        parameter2read = 'resAvg_file'
        call read_inputFile_integer(lineCnt,parameter2read,results_first)

        results_last = results_first
        results_step = 1
    else if(type_resultsFile .eq. 3) then
        !8. restart_file----------------------------------------------------------------
        parameter2read = 'restart_file'
        call read_inputFile_integer(lineCnt,parameter2read,results_first)

        results_last = results_first
        results_step = 1
    else if(type_resultsFile .le. 4) then
        !8. results_file-----------------------------------------------------------------
        parameter2read = 'results_file'
        call read_inputFile_integer(lineCnt,parameter2read,results_first)

        !9. restart_file-----------------------------------------------------------------
        parameter2read = 'restart_file'
        call read_inputFile_integer(lineCnt,parameter2read,results_step)
        
        results_last = results_first
    else
        write(*,*) "Wrong type_resultsFile! Must be 1,2 or 3 (1:inst, 2:avg, 3:restart, 4:inst_to_restart[todo])! Aborting!"
       	call MPI_Abort(app_comm,-1,mpi_err)
    end if

    call close_inputFile()
    if(mpi_rank.eq.0) write(*,*) '## End of Reading input file: ',trim(adjustl(input_file))

!---------------------------------------------------------------------------------------------------------
    call copy_results_same_mesh_Npartitions(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName,target_Nprocs,&
                                            type_resultsFile,generateNewMesh,results_first,results_last,results_step)
!---------------------------------------------------------------------------------------------------------

    call end_mpi()

end program tool_copy_results