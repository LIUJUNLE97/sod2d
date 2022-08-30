program test_alya_mesh
    use mod_mpi
    use mod_mpi_mesh
    use mod_comms
    use mod_hdf5
    implicit none

    character(500) :: gmsh_file_path,gmsh_file_name
    character(500) :: mesh_h5_file_path,mesh_h5_file_name

    real(8), dimension(:), allocatable :: d_field
    real(4), dimension(:), allocatable :: f_field
    logical :: periodic,loadMesh

!------------------------------------------------------------------------------------------------------
    call init_mpi()

    if(mpi_rank.eq.0) write(*,*) ' #Implementing cgns format...'

    write(gmsh_file_path,*) "./mesh/"
    !write(gmsh_file_path,*) "./mesh_cube50/"
    !write(gmsh_file_path,*) "./mesh2/"
    write(gmsh_file_name,*) "cube" ! Nsys

    write(mesh_h5_file_path,*) ""
    write(mesh_h5_file_name,*) "cube" ! Nsys

    !periodic = .false.
    periodic = .true.

    loadMesh = .false.
    !loadMesh = .true.
    
    if(loadMesh) then
        call load_hdf5_meshfile(mesh_h5_file_path,mesh_h5_file_name)
    else
        call read_alyaMesh_part_and_create_hdf5Mesh(gmsh_file_path,gmsh_file_name,periodic,&
                                                    mesh_h5_file_path,mesh_h5_file_name)
    end if

    !------------------------------------------------------------------------------------
    !-- read the alya mesh fesh files in GMSH/ALYA FORMAT
    !call read_alya_mesh_files(gmsh_file_path,gmsh_file_name,periodic)

    !------------------------------------------------------------------------------------
    !----- Do mesh partitioning!
    !------------------------------------------------------------------------------------
    !call do_mesh_partitioning()

    !------------------------------------------------------------------------------------
    !----- Create HDF5 File
    !------------------------------------------------------------------------------------
    !call create_hdf5_meshfile(h5_file_path,h5_file_name)

#if 0
    !NOW I CAN PUT ALL THIS SHIT IN THE ALYA_MESH_MOD

    !------------------------------------------------------------------------------------
    !----- Do COMMS
    !------------------------------------------------------------------------------------

    ALL THIS PART HAS TO BE REWRITTEN IF WANT TO BE USED!

    !----------------------------------------------------------------------------------------
    call init_comms(useIntInComms,useFloatInComms,useDoubleInComms)

    allocate(d_field(numNodesRankPar))
    allocate(f_field(numNodesRankPar))

    !res_field(:) = (mpi_rank+1)*10.d0
    d_field(:) = 10.d0
    f_field(:) = 5.d0
!--------------------------------------------------------------------------------
    call update_and_comm_doubleField(d_field)

    call print_csv_file_dfield('dfield',d_field)

!--------------------------------------------------------------------------------
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
!--------------------------------------------------------------------------------
    call init_shared_mem_window()

    call update_and_comm_shared_mem_floatField(f_field)

    call close_shared_mem_windows()
    call print_csv_file_ffield('ffield',f_field)
!--------------------------------------------------------------------------------
#endif
    call end_mpi()

end program test_alya_mesh