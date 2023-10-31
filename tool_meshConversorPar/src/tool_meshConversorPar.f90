
program tool_meshConversorPar
    use mod_mpi
    use mod_read_inputFile
    use mod_meshConversorTool

    implicit none

    character(999) :: input_file
    character(512) :: gmsh_filePath,gmsh_fileName
    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    character(256) :: parameter2read
    integer(4) :: lineCnt,num_partitions
    logical :: lineal_output
    logical :: eval_mesh_quality

!------------------------------------------------------------------------------------------------------

    call init_mpi()

    if(mpi_rank.eq.0) then
        write(*,*) '|-- WELCOME TO THE AWESOME MESH CONVERSION PARALLEL TOOL ! ;)'
        write(*,*) '|-- From GMSH to Sod2D format'
    end if

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, input_file)
        if(mpi_rank.eq.0) write(*,*) 'Input file: ',trim(adjustl(input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call this amazing tool with an input file!!!'
        call MPI_Abort(app_comm,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ! Reading the parameters
    call open_inputFile(input_file)
    lineCnt = 1

    !1. gmsh_filePath--------------------------------------------------------------
    parameter2read = 'gmsh_filePath'
    call read_inputFile_string(lineCnt,parameter2read,gmsh_filePath)

    !2. gmsh_fileName--------------------------------------------------------------
    parameter2read = 'gmsh_fileName'
    call read_inputFile_string(lineCnt,parameter2read,gmsh_fileName)

    !3. mesh_h5_filePath--------------------------------------------------------------
    parameter2read = 'mesh_h5_filePath'
    call read_inputFile_string(lineCnt,parameter2read,mesh_h5_filePath)

    !4. mesh_h5_fileName--------------------------------------------------------------
    parameter2read = 'mesh_h5_fileName'
    call read_inputFile_string(lineCnt,parameter2read,mesh_h5_fileName)

    !5. num_partitions--------------------------------------------------------------------------
    parameter2read = 'num_partitions'
    call read_inputFile_integer(lineCnt,parameter2read,num_partitions)

    !6. lineal_output--------------------------------------------------------------------------
    parameter2read = 'lineal_output'
    call read_inputFile_logical(lineCnt,parameter2read,lineal_output)

    !7. eval_mesh_quality--------------------------------------------------------------------------
    parameter2read = 'eval_mesh_quality'
    call read_inputFile_logical(lineCnt,parameter2read,eval_mesh_quality)

    close(99)
    if(mpi_rank.eq.0) write(*,*) '## End of Reading input file: ',trim(adjustl(input_file))

!---------------------------------------------------------------------------------------------------------

    call read_gmsh_h5_file_and_do_partitioning_in_parallel(gmsh_filePath,gmsh_fileName,mesh_h5_filePath,mesh_h5_fileName,num_partitions,lineal_output,eval_mesh_quality)

!---------------------------------------------------------------------------------------------------------

    call end_mpi()

end program tool_meshConversorPar
