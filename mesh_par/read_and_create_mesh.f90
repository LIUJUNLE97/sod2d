program test_alya_mesh
    use mod_mpi
    use mod_cgns_mesh
    use mod_comms
    implicit none

    real(8) :: x, y, z

    character(500) :: gmsh_file_path, gmsh_file_name
    character(500) :: cgns_file_path, cgns_file_name

    real(8), dimension(:), allocatable :: res_field
    integer :: index_res


!------------------------------------------------------------------------------------------------------
    call init_mpi()

    if(mpi_rank.eq.0) write(*,*) ' #Implementing cgns format...'

    !write(gmsh_file_path,*) "./mesh/"
    write(gmsh_file_path,*) "./mesh_cube50/"
    !write(gmsh_file_path,*) "./mesh2/"
    write(gmsh_file_name,*) "cube" ! Nsys

    write(cgns_file_path,*) ""
    !write(cgns_file_name,*) "cube"
    write(cgns_file_name,*) "cube50"

    !------------------------------------------------------------------------------------
    !-- read the alya mesh fesh files in GMSH/ALYA FORMAT
    call read_alya_mesh_files(gmsh_file_path,gmsh_file_name)

    !------------------------------------------------------------------------------------
    !----- ok ! now it's my time to do cgns stuff :)
    !------------------------------------------------------------------------------------

    !let's do mesh partitioning
    call do_mesh_partitioning()

    !NOW I CAN PUT ALL THIS SHIT IN THE ALYA_MESH_MOD
#if 1
    call create_CGNSmesh_par(cgns_file_path,cgns_file_name)

    !----------------------------------------------------------------------------------------
    call init_comms()

    allocate(res_field(numNodesRankPar))

    !res_field(:) = (mpi_rank+1)*10.d0
    res_field(:) = 10.d0

    call add_write_doubleField_CGNSmesh_vertexSolution("res",res_field,index_res)

    call update_and_comm_dfield(res_field)

    call write_doubleField_CGNSmesh_vertexSolution(index_res,res_field)

    call close_CGNSmesh_par()
#endif
    call end_mpi()

end program test_alya_mesh