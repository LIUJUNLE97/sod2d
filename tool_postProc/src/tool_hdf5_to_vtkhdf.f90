#define average 0

module str2int_mod
contains 

  elemental subroutine str2int(str,int,stat)
    implicit none
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat

    read(str,*,iostat=stat)  int
  end subroutine str2int

end module

program tool_hdf5_to_cgns
    use str2int_mod
    use mod_mpi
    use mod_mpi_mesh
    use mod_hdf5
    !use mod_cgns_mesh
    implicit none

    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    character(512) :: results_h5_filePath,results_h5_fileName
    character(999) :: full_fileName,input_file,read_sec,read_val
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
    ! Reading the parameters
    open(99,file=input_file,status="old")

    !1. mesh_h5_filePath--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'mesh_h5_filePath') then
        mesh_h5_filePath = read_val
        if(mpi_rank.eq.0) write(*,*) 'mesh_h5_filePath: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 1 must be mesh_h5_filePath value'
        call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
    endif

    !2. mesh_h5_fileName--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'mesh_h5_fileName') then
        mesh_h5_fileName = read_val
        if(mpi_rank.eq.0) write(*,*) 'mesh_h5_fileName: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 2 must be mesh_h5_fileName value'
        call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
    endif

    !3. results_h5_filePath--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'results_h5_filePath') then
        results_h5_filePath = read_val
        if(mpi_rank.eq.0) write(*,*) 'results_h5_filePath: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 3 must be results_h5_filePath value'
        call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
    endif

    !4. results_h5_fileName--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'results_h5_fileName') then
        results_h5_fileName = read_val
        if(mpi_rank.eq.0) write(*,*) 'results_h5_fileName: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 4 must be results_h5_fileName value'
        call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
    endif

    !5. do_averages--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'do_averages') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            if(aux_int.eq.0) then
                do_averages = .false.
            elseif(aux_int .eq. 1) then
                do_averages = .true.
            else
                if(mpi_rank.eq.0) then
                    write(*,*) 'ERROR! do_averages must be 0(false)/1(true)',&
                               ' | wrong value ',trim(adjustl(read_val))
                    call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
                end if
            end if
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! do_averages must be 0(false)/1(true)',&
                           ' | wrong value ',trim(adjustl(read_val))
                call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'do_averages: ',do_averages
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 5 must be do_averages 0/1'
    endif

    !6. first_step--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'first_step') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            first_step = aux_int
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! first_step must be an integer',&
                           ' | wrong value ',trim(adjustl(read_val))
                call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'first_step: ',first_step
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 6 must be first_step value'
    endif

    !7. last_step--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'last_step') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            last_step = aux_int
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! last_step must be an integer',&
                           ' | wrong value ',trim(adjustl(read_val))
                call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'last_step: ',last_step
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 7 must be last_step value'
    endif

    !8. nstep--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'nstep') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            nstep = aux_int
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! nstep must be an integer',&
                           ' | wrong value ',read_val
                call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'nstep: ',nstep
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 8 must be nstep value'
    endif

    close(99)
    if(mpi_rank.eq.0) write(*,*) '## End of Reading input file: ',trim(adjustl(input_file))

!---------------------------------------------------------------------------------------------------------

    call init_hdf5_interface(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName)

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
