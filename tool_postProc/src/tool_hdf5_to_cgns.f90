program tool_hdf5_to_cgns
    use mod_mpi
    use mod_mpi_mesh
    use mod_hdf5
    use mod_cgns_mesh
    implicit none

    character(512) :: mesh_h5_file_path,mesh_h5_file_name
    character(512) :: cgns_file_path,cgns_file_name
    character(512) :: results_h5_file_path,results_h5_file_name
    !logical :: periodic,loadMesh

    integer :: load_step
    real(rp) :: time
    real(rp), allocatable :: rho(:),pr(:),E(:),eta(:),csound(:),machno(:),divU(:),Qcrit(:)
    real(rp), allocatable :: envit(:),mut(:),mu_fluid(:)
    real(rp), allocatable :: u(:,:),gradRho(:,:),curlU(:,:)

!------------------------------------------------------------------------------------------------------

    write(mesh_h5_file_path,*) ""
    write(mesh_h5_file_name,*) "cube10"

    write(cgns_file_path,*) ""
    write(cgns_file_name,*) "cube10"

    write(results_h5_file_path,*) ""
    write(results_h5_file_name,*) "results"

    load_step = 501

!------------------------------------------------------------------------------------------------------
    call init_mpi()

    if(mpi_rank.eq.0) write(*,*) '## CONVERSION TOOL HDF5 -> CGNS ##'

    call init_hdf5_interface(mesh_h5_file_path,mesh_h5_file_name,results_h5_file_path,results_h5_file_name)

    if(mpi_rank.eq.0) write(*,*) '## LOADING HDF5 MESH FILE... ##'
    call load_hdf5_meshfile(mesh_h5_file_path,mesh_h5_file_name)

    allocate(rho(numNodesRankPar))
    allocate(pr(numNodesRankPar))
    allocate(E(numNodesRankPar))
    allocate(eta(numNodesRankPar))
    allocate(csound(numNodesRankPar))
    allocate(machno(numNodesRankPar))
    allocate(divU(numNodesRankPar))
    allocate(Qcrit(numNodesRankPar))
    allocate(envit(numNodesRankPar))
    allocate(mut(numNodesRankPar))
    allocate(mu_fluid(numNodesRankPar))
    allocate(u(numNodesRankPar,ndime))
    allocate(gradRho(numNodesRankPar,ndime))
    allocate(curlU(numNodesRankPar,ndime))

    if(mpi_rank.eq.0) write(*,*) '## LOADING HDF5 RESULTS FILE... ##'

    call load_hdf5_resultsFile_allArrays(load_step,time,rho,u,pr,E,eta,csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,envit,mut)

    if(mpi_rank.eq.0) write(*,*) '## CREATING CGNS FILE... ##'
    call create_CGNSmesh_par(cgns_file_path,cgns_file_name)
    !call create_hdf5_meshfile(h5_file_path,h5_file_name)

    call add_write_floatField_CGNSmesh_vertexSolution('rho',rho)
    call add_write_floatField_CGNSmesh_vertexSolution('pr',pr)

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
    ! End hdf5 interface
    call end_hdf5_interface()

    call end_mpi()

end program tool_hdf5_to_cgns