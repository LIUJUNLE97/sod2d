
program tool_commsPerfomance
    use mod_mpi
    use mod_ioutils
    use mod_mpi_mesh
    use mod_comms
    use mod_hdf5
    use mod_comms_performance
    implicit none

    character(512) :: json_input_file
    integer :: numNodesSrl,numNodesB_1r,numIters
    character(512) :: mesh_h5_file_path,mesh_h5_file_name,meshFile_h5_full_name
    character(1024) :: log_file_comms,log_file_append
    character(128) :: aux_string_mpisize,aux_string_numNodes,aux_string_nodesBound
    logical :: useMesh=.false.
    type(json_file) :: json_f
    logical :: isFound
    character(len=:), allocatable :: str_value
!------------------------------------------------------------------------------------------------------

    call init_mpi()

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, json_input_file)
        if(mpi_rank.eq.0) write(*,*) '# Input file: ',trim(adjustl(json_input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call Tool COMMS PERFROMANCE with an input file!'
        call MPI_Abort(app_comm,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
    ! Opening json file
    call open_json_file(json_input_file,json_f)

    !1. numIters--------------------------------------------------------------------------
    call json_f%get("numIters",numIters,isFound,1000)

    !2. Use mesh--------------------------------------------------------------------------
    call json_f%get("useMesh",useMesh,isFound,.true.);

    if(useMesh) then
        !3. mesh_h5_file_path--------------------------------------------------------------
        call json_f%get("mesh_h5_file_path",str_value,isFound,"")
        write(mesh_h5_file_path,*) str_value

        !4. mesh_h5_file_name--------------------------------------------------------------
        call json_f%get("mesh_h5_file_name",str_value,isFound,"")
        write(mesh_h5_file_name,*) str_value
    else
        !3. numNodesSrl--------------------------------------------------------------------------
        call json_f%get("numNodesSrl",numNodesSrl,isFound,100000)

        !4. numNodesB_1r --------------------------------------------------------------------------
        call json_f%get("numNodesB_1r",numNodesB_1r,isFound,100)

       if(mpi_rank.eq.0) write(*,*) 'numNodesSrl ',numNodesSrl,' numNodesB_1r ',numNodesB_1r,' numIters ',numIters
    end if

    ! Closing json file
    call close_json_file(json_f)

    !-----------------------------------------------------------------------------------
    if(useMesh) then
        call init_hdf5_interface()
        call set_hdf5_meshFile_name(mesh_h5_file_path,mesh_h5_file_name,mpi_size,meshFile_h5_full_name)
        call load_hdf5_meshfile(meshFile_h5_full_name)
    else
        call create_dummy_1Dmesh(numNodesSrl,numNodesB_1r)
    end if

    if(mpi_rank.eq.0) then
        write(aux_string_mpisize,'(I0)') mpi_size
        if(useMesh) then
            log_file_append = trim(adjustl(mesh_h5_file_name))
        else
            write(aux_string_numNodes,'(I0)') numNodesSrl
            write(aux_string_nodesBound,'(I0)') numNodesB_1r
            log_file_append = 'mesh1D_'//trim(aux_string_numNodes)//'_bound_'//trim(aux_string_nodesBound)
        end if
        log_file_comms = 'commPerf_'//trim(adjustl(log_file_append))//'-'//trim(aux_string_mpisize)//'.dat'
    end if

   call test_comms_performance_real(numIters,log_file_comms)

   call end_comms()
   call end_mpi()

end program tool_commsPerfomance
