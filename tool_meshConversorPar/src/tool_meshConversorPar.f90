
program tool_meshConversorPar
    use mod_mpi
    use mod_ioutils
    use mod_meshConversorTool
    use json_module

    implicit none

    character(512) :: json_input_file
    character(512) :: gmsh_filePath,gmsh_fileName
    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    integer(4) :: num_partitions,eval_mesh_quality
    logical :: lineal_output
    type(json_file) :: json_f
    logical :: isFound
    character(len=:), allocatable :: str_value
!------------------------------------------------------------------------------------------------------

    call init_mpi()

    if(mpi_rank.eq.0) then
        write(*,*) '|-- WELCOME TO THE AWESOME MESH CONVERSION PARALLEL TOOL ! ;)'
        write(*,*) '|-- From GMSH to Sod2D format'
    end if

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, json_input_file)
        if(mpi_rank.eq.0) write(*,*) 'Input file: ',trim(adjustl(json_input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call this amazing tool with an input file!!!'
        call MPI_Abort(app_comm,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    ! Opening json file
    call open_json_file(json_input_file,json_f)

    !1. gmsh_filePath--------------------------------------------------------------
    call json_f%get("gmsh_filePath",str_value,isFound,"")
    write(gmsh_filePath,*) str_value

    !2. gmsh_fileName--------------------------------------------------------------
    call json_f%get("gmsh_fileName",str_value,isFound,"")
    write(gmsh_fileName,*) str_value

    !3. mesh_h5_filePath--------------------------------------------------------------
    call json_f%get("mesh_h5_filePath",str_value,isFound,"")
    write(mesh_h5_filePath,*) str_value

    !4. mesh_h5_fileName--------------------------------------------------------------
    call json_f%get("mesh_h5_fileName",str_value,isFound,"")
    write(mesh_h5_fileName,*) str_value

    !5. num_partitions--------------------------------------------------------------------------
    call json_f%get("num_partitions",num_partitions,isFound,1);

    !6. lineal_output--------------------------------------------------------------------------
    call json_f%get("lineal_output",lineal_output,isFound,.true.);

    !7. eval_mesh_quality--------------------------------------------------------------------------
    call json_f%get("eval_mesh_quality",eval_mesh_quality,isFound,0);

    ! Closing json file
    call close_json_file(json_f)

!---------------------------------------------------------------------------------------------------------

    call read_gmsh_h5_file_and_do_partitioning_in_parallel(gmsh_filePath,gmsh_fileName,mesh_h5_filePath,mesh_h5_fileName,num_partitions,lineal_output,eval_mesh_quality)

!---------------------------------------------------------------------------------------------------------

    call end_mpi()

end program tool_meshConversorPar
