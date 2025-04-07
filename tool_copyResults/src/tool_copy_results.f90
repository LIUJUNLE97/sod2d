
program tool_copy_results
    use mod_mpi
    use json_module
    use mod_ioutils
    use mod_mpi_mesh
    use mod_comms
    use mod_hdf5
    use mod_copy_results
    implicit none

    character(512) :: json_input_file
    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    character(512) :: results_h5_filePath,results_h5_fileName
    integer(4) :: type_resultsFile,results_first,results_last,results_step,target_Nprocs
    logical :: generateNewMesh=.false.
    type(json_file) :: json_f
    logical :: isFound
    character(len=:), allocatable :: str_value
!------------------------------------------------------------------------------------------------------

    call init_mpi()

    if(mpi_rank.eq.0) then
        write(*,*) '|-- WELCOME TO SOD2D COPY RESULTS TOOL!'
        write(*,*) '|-- Copy results in same mesh and different number of partitions'
    end if

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, json_input_file)
        if(mpi_rank.eq.0) write(*,*) '# Input file: ',trim(adjustl(json_input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call tool Copy Results with an input file!'
        call MPI_Abort(app_comm,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
    ! Opening json file
    call open_json_file(json_input_file,json_f)

    !1. mesh_h5_filePath--------------------------------------------------------------
    call json_f%get("mesh_h5_filePath",str_value,isFound,"")
    write(mesh_h5_filePath,*) str_value

    !2. mesh_h5_fileName--------------------------------------------------------------
    call json_f%get("mesh_h5_fileName",str_value,isFound,"")
    write(mesh_h5_fileName,*) str_value

    !3. generate_mesh
    call json_f%get("generate_mesh",generateNewMesh,isFound,.false.);

    !4. results_h5_filePath-----------------------------------------------------------
    call json_f%get("results_h5_filePath",str_value,isFound,"")
    write(results_h5_filePath,*) str_value

    !5. results_h5_fileName-----------------------------------------------------------
    call json_f%get("results_h5_fileName",str_value,isFound,"")
    write(results_h5_fileName,*) str_value

    !6. target_Nprocs-----------------------------------------------------------------
    call json_f%get("target_Nprocs",target_Nprocs,isFound,1);

    !7. type_resultsFile-------------------------------------------------------------
    call json_f%get("type_resultsFile",type_resultsFile,isFound,1);

    if(type_resultsFile .eq. 1) then
        !8. results_first----------------------------------------------------------------
        call json_f%get("results_first",results_first,isFound,1);

        !9. results_last-----------------------------------------------------------------
        call json_f%get("results_last",results_last,isFound,1);

        !10. results_step-----------------------------------------------------------------
        call json_f%get("results_step",results_step,isFound,1);
    else if(type_resultsFile .eq. 2) then
        !8. resAvg_file----------------------------------------------------------------
        call json_f%get("resAvg_file",results_first,isFound,1);

        results_last = results_first
        results_step = 1
    else if(type_resultsFile .eq. 3) then
        !8. restart_file----------------------------------------------------------------
        call json_f%get("restart_file",results_first,isFound,1);

        results_last = results_first
        results_step = 1
    else if(type_resultsFile .eq. 4) then
        !8. results_file-----------------------------------------------------------------
        call json_f%get("results_file",results_first,isFound,1);

        !9. restart_file-----------------------------------------------------------------
        call json_f%get("restart_file",results_step,isFound,1);
        
        results_last = results_first
    else if(type_resultsFile .le. 5) then
        results_first = 0
        results_step = 0
        results_last = 0
    else
        write(*,*) "Wrong type_resultsFile! Must be 1,2,3,4 or 5 (1:inst, 2:avg, 3:restart, 4:inst_to_restart, 5:mapping)! Aborting!"
       	call MPI_Abort(app_comm,-1,mpi_err)
    end if

    ! Closing json file
    call close_json_file(json_f)

!---------------------------------------------------------------------------------------------------------
    call copy_results_same_mesh_Npartitions(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName,target_Nprocs,&
                                            type_resultsFile,generateNewMesh,results_first,results_last,results_step)
!---------------------------------------------------------------------------------------------------------

    call end_mpi()

end program tool_copy_results