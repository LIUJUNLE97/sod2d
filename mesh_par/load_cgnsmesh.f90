program test_alya_mesh
    use mod_mpi
    use mod_cgns_mesh
    use mod_comms
    implicit none

    character(500) :: cgns_file_path, cgns_file_name

    real(8), dimension(:), allocatable :: res_field, aux_comm_field_s, aux_comm_field_r
    character(128) :: field_name
    integer :: index_res

!------------------------------------------------------------------------------------------------------
    call init_mpi()

    write(cgns_file_path,*) ""
    write(cgns_file_name,*) "cube50"

    call open_CGNSmesh_par(cgns_file_path,cgns_file_name)
    call init_comms()


    allocate(res_field(numNodesRankPar))
    !res_field(:) = (mpi_rank+1)*10.d0
    res_field(:) = 5.d0

    !index_sol_V=1
    !index_res=3

    field_name = "res"
    index_res = getFieldId_V(field_name)

    call write_doubleField_CGNSmesh_vertexSolution(index_res,res_field)

    call update_and_comm_dfield(res_field)

    call write_doubleField_CGNSmesh_vertexSolution(index_res,res_field)

    call close_CGNSmesh_par()

    call end_mpi()

end program test_alya_mesh