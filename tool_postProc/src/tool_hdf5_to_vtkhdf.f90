
program tool_hdf5_to_cgns
    use mod_mpi
    use mod_read_inputFile
    use mod_mpi_mesh
    use mod_hdf5
    !use mod_cgns_mesh
    implicit none

    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    character(512) :: results_h5_filePath,results_h5_fileName
    character(999) :: full_fileName,input_file,read_sec,read_val
    character(256) :: parameter2read
    integer :: lineCnt

    logical :: do_averages=.false.
    integer :: first_step,last_step,nstep,iStep,numAvgSteps,stat,aux_int, i, j
    real(rp) :: time

    !---------------------- ARRAYS -------------------------
    !-- For average 
    real(rp), allocatable :: avrho(:),avpre(:),avmueff(:)
    real(rp), allocatable :: avvel(:,:),avve2(:,:)
    real(rp), allocatable :: favrho(:),favpre(:),favmueff(:)
    real(rp), allocatable :: favvel(:,:),favve2(:,:)
    !-- For inst 
    real(rp), allocatable :: rho(:),pr(:),E(:),eta(:),csound(:),machno(:),divU(:),Qcrit(:)
    real(rp), allocatable :: envit(:),mut(:),mu_fluid(:)
    real(rp), allocatable :: u(:,:),gradRho(:,:),curlU(:,:)

!------------------------------------------------------------------------------------------------------

    call init_mpi()

    if(mpi_rank.eq.0) write(*,*) '## CONVERSION TOOL HDF5 -> VTKHDF ##'

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, input_file)
        if(mpi_rank.eq.0) write(*,*) 'Reading input file: ',trim(adjustl(input_file))
    else
        if(mpi_rank.eq.0) write(*,*) 'You must call this amazing tool with an input file!!!'
        call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
    endif
    !------------------------------------------------------------------------------
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

    !3. results_h5_filePath--------------------------------------------------------------
    parameter2read = 'results_h5_filePath'
    call read_inputFile_string(lineCnt,parameter2read,results_h5_filePath)

    !4. results_h5_fileName--------------------------------------------------------------
    parameter2read = 'results_h5_fileName'
    call read_inputFile_string(lineCnt,parameter2read,results_h5_fileName)

    !5. do_averages--------------------------------------------------------------------------
    parameter2read = 'do_averages'
    call read_inputFile_logical(lineCnt,parameter2read,do_averages)

    !6. first_step--------------------------------------------------------------------------
    parameter2read = 'first_step'
    call read_inputFile_integer(lineCnt,parameter2read,first_step)

    !7. last_step--------------------------------------------------------------------------
    parameter2read = 'last_step'
    call read_inputFile_integer(lineCnt,parameter2read,last_step)

    !8. nstep--------------------------------------------------------------------------
    parameter2read = 'nstep'
    call read_inputFile_integer(lineCnt,parameter2read,nstep)

    close(99)
    if(mpi_rank.eq.0) write(*,*) '## End of Reading input file: ',trim(adjustl(input_file))

!---------------------------------------------------------------------------------------------------------

    call init_hdf5_interface()
    call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName)
    call set_hdf5_baseResultsFile_name(results_h5_filePath,results_h5_fileName,mesh_h5_fileName)

    if(mpi_rank.eq.0) write(*,*) '# Loading HDF5 mesh file...'
    !call load_hdf5_meshfile(mesh_h5_filePath,mesh_h5_fileName)
    call load_hdf5_meshfile()

    if(do_averages) then
        numAvgSteps = 0
        allocate(avrho(numNodesRankPar))
        allocate(avpre(numNodesRankPar))
        allocate(avmueff(numNodesRankPar))
        allocate(avvel(numNodesRankPar,ndime))
        allocate(avve2(numNodesRankPar,ndime))
        allocate(favrho(numNodesRankPar))
        allocate(favpre(numNodesRankPar))
        allocate(favmueff(numNodesRankPar))
        allocate(favvel(numNodesRankPar,ndime))
        allocate(favve2(numNodesRankPar,ndime))
        do i = 1, numNodesRankPar 
           favrho(i)   = 0.0_rp
           favpre(i)   = 0.0_rp
           favmueff(i) = 0.0_rp
           do j = 1, ndime
              favvel(i,j) = 0.0_rp
              favve2(i,j) = 0.0_rp
           end do
        end do
    else
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
    end if

    do iStep=first_step,last_step,nstep
        if(mpi_rank.eq.0) write(*,*) '## Doing iStep',iStep   

        if(mpi_rank.eq.0) write(*,*) '# Loading HDF5 Results file...'

        if(do_averages) then
            numAvgSteps = numAvgSteps + 1
            call load_hdf5_avgResultsFile(iStep,avvel,avve2,avrho,avpre,avmueff)
        else
            call load_hdf5_resultsFile_allArrays(iStep,time,rho,u,pr,E,eta,csound,&
                                machno,gradRho,curlU,divU,Qcrit,mu_fluid,envit,mut)
        end if

        if(mpi_rank.eq.0) write(*,*) '# Creating VTKHDF file...'

        if(do_averages) then

            call save_vtkhdf_avgResultsFile(iStep,avrho,avpre,avmueff,avvel,avve2)

            favrho(:)   = favrho(:)   + avrho(:)
            favpre(:)   = favpre(:)   + avpre(:)
            favmueff(:) = favmueff(:) + avmueff(:)
            favvel(:,1:3) = favvel(:,1:3) + avvel(:,1:3)
            favve2(:,1:3) = favve2(:,1:3) + avve2(:,1:3)
        else
            call save_vtkhdf_instResultsFile(iStep,rho,pr,E,eta,csound,machno,divU,Qcrit,&
                                            envit,mut,mu_fluid,u,gradRho,curlU)
        endif
    end do

    if(do_averages) then
        if(mpi_rank.eq.0) write(*,*) '# Doing final Avg...'

        do i = 1, numNodesRankPar 
           favrho(i)   = favrho(i)   / real(numAvgSteps, rp)
           favpre(i)   = favpre(i)   / real(numAvgSteps, rp)
           favmueff(i) = favmueff(i) / real(numAvgSteps, rp)
           do j = 1, ndime
              favvel(i,j) = favvel(i,j) / real(numAvgSteps, rp)
              favve2(i,j) = favve2(i,j) / real(numAvgSteps, rp)
           end do
        end do

        call save_vtkhdf_finalAvgResultsFile(favrho,favpre,favmueff,favvel,favve2)
    end if

    call end_hdf5_interface()

    call end_mpi()

end program tool_hdf5_to_cgns
