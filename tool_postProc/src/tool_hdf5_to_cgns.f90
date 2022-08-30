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
    use mod_cgns_mesh
    implicit none

    character(512) :: mesh_h5_filePath,mesh_h5_fileName
    character(512) :: cgns_filePath,cgns_fileName
    character(512) :: results_h5_filePath,results_h5_fileName
    character(128) :: output_cgns_fileName
    character(999) :: full_fileName,input_file,read_sec,read_val
    logical :: do_averages=.false.

    integer :: first_step,last_step,nstep,iStep,numAvgSteps,stat,aux_int
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

    if(mpi_rank.eq.0) write(*,*) '## CONVERSION TOOL HDF5 -> CGNS ##'

    !------------------------------------------------------------------------------
    ! Reading input file
    if(command_argument_count() .eq. 1) then
        call get_command_argument(1, input_file)
        !read(cnsteps,'(i)') numNodesSrl
        if(mpi_rank.eq.0) write(*,*) 'Reading input file: ',trim(adjustl(input_file))
    else
        stop 0
    endif
    !------------------------------------------------------------------------------
    ! Reading the parameters
    open(99,file=input_file,status="old")

    !mesh_h5_filePath--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'mesh_h5_filePath') then
        mesh_h5_filePath = read_val
        if(mpi_rank.eq.0) write(*,*) 'mesh_h5_filePath: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 1 must be mesh_h5_filePath value'
        stop 0
    endif

    !mesh_h5_fileName--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'mesh_h5_fileName') then
        mesh_h5_fileName = read_val
        if(mpi_rank.eq.0) write(*,*) 'mesh_h5_fileName: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 2 must be mesh_h5_fileName value'
        stop 0
    endif

    !cgns_filePath--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'cgns_filePath') then
        cgns_filePath = read_val
        if(mpi_rank.eq.0) write(*,*) 'cgns_filePath: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 3 must be cgns_filePath value'
        stop 0
    endif

    !cgns_fileName--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'cgns_fileName') then
        cgns_fileName = read_val
        if(mpi_rank.eq.0) write(*,*) 'cgns_fileName: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 4 must be cgns_fileName value'
        stop 0
    endif

    !results_h5_filePath--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'results_h5_filePath') then
        results_h5_filePath = read_val
        if(mpi_rank.eq.0) write(*,*) 'results_h5_filePath: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 5 must be results_h5_filePath value'
        stop 0
    endif

    !results_h5_fileName--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'results_h5_fileName') then
        results_h5_fileName = read_val
        if(mpi_rank.eq.0) write(*,*) 'results_h5_fileName: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 6 must be results_h5_fileName value'
        stop 0
    endif

    !output_cgns_fileName--------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'output_cgns_fileName') then
        output_cgns_fileName = read_val
        if(mpi_rank.eq.0) write(*,*) 'output_cgns_fileName: ',trim(adjustl(read_val))
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 7 must be output_cgns_fileName value'
        stop 0
    endif

    !do_averages--------------------------------------------------------------------------
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
                    stop 0
                end if
            end if
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! do_averages must be 0(false)/1(true)',&
                           ' | wrong value ',trim(adjustl(read_val))
                stop 0
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'do_averages: ',do_averages
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 8 must be do_averages 0/1'
    endif

    !first_step--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'first_step') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            first_step = aux_int
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! first_step must be 0(false)/1(true)',&
                           ' | wrong value ',trim(adjustl(read_val))
                stop 0
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'first_step: ',first_step
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 9 must be first_step value'
    endif

    !last_step--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'last_step') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            last_step = aux_int
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! last_step must be 0(false)/1(true)',&
                           ' | wrong value ',trim(adjustl(read_val))
                stop 0
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'last_step: ',last_step
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 10 must be last_step value'
    endif

    !nstep--------------------------------------------------------------------------
    read(99,*) read_sec,read_val! Section header
    if(read_sec.eq.'nstep') then
        call str2int(read_val,aux_int,stat)
        if(stat.eq.0) then
            nstep = aux_int
        else
            if(mpi_rank.eq.0) then
                write(*,*) 'ERROR! nstep must be 0(false)/1(true)',&
                           ' | wrong value ',read_val
                stop 0
            end if
        end if
        if(mpi_rank.eq.0) write(*,*) 'nstep: ',nstep
    else
        if(mpi_rank.eq.0) write(*,*) 'Error! Line 11 must be nstep value'
    endif

    close(99)
    if(mpi_rank.eq.0) write(*,*) '## End of Reading input file: ',trim(adjustl(input_file))

!---------------------------------------------------------------------------------------------------------

    call init_hdf5_interface(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName)

    if(mpi_rank.eq.0) write(*,*) '# Loading HDF5 mesh file...'
    call load_hdf5_meshfile(mesh_h5_filePath,mesh_h5_fileName)

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

    call init_CGNSmesh_arrays()

    !here the loop!
    do iStep=first_step,last_step,nstep
        if(mpi_rank.eq.0) write(*,*) '## Doing iStep',iStep   

        call set_CGNS_full_fileName(full_fileName,cgns_filePath,cgns_fileName,output_cgns_fileName,iStep)

        if(mpi_rank.eq.0) write(*,*) '# Loading HDF5 Results file...'

        if(do_averages) then
            numAvgSteps = numAvgSteps + 1
            call load_hdf5_avgResultsFile(iStep,avvel,avve2,avrho,avpre,avmueff)
        else
            call load_hdf5_resultsFile_allArrays(iStep,time,rho,u,pr,E,eta,csound,&
                                machno,gradRho,curlU,divU,Qcrit,mu_fluid,envit,mut)
        end if

        if(mpi_rank.eq.0) write(*,*) '# Creating CGNS file...'
        call create_CGNSmesh_par(full_fileName)

        if(do_averages) then
            call add_write_floatField_CGNSmesh_vertexSolution('avrho',avrho)
            call add_write_floatField_CGNSmesh_vertexSolution('avpre',avpre)
            call add_write_floatField_CGNSmesh_vertexSolution('avmueff',avmueff)
            call add_write_floatField_CGNSmesh_vertexSolution('avvelX',avvel(:,1))
            call add_write_floatField_CGNSmesh_vertexSolution('avvelY',avvel(:,2))
            call add_write_floatField_CGNSmesh_vertexSolution('avvelZ',avvel(:,3))
            call add_write_floatField_CGNSmesh_vertexSolution('avve2X',avve2(:,1))
            call add_write_floatField_CGNSmesh_vertexSolution('avve2Y',avve2(:,2))
            call add_write_floatField_CGNSmesh_vertexSolution('avve2Z',avve2(:,3))

            favrho(:)   = favrho(:)   + avrho(:)
            favpre(:)   = favpre(:)   + avpre(:)
            favmueff(:) = favmueff(:) + avmueff(:)
            favvel(:,:) = favvel(:,:) + avvel(:,:)
            favve2(:,:) = favve2(:,:) + avve2(:,:)
        else
            call add_write_floatField_CGNSmesh_vertexSolution('rho',rho)
            call add_write_floatField_CGNSmesh_vertexSolution('VelocityX',u(:,1))
            call add_write_floatField_CGNSmesh_vertexSolution('VelocityY',u(:,2))
            call add_write_floatField_CGNSmesh_vertexSolution('VelocityZ',u(:,3))
            call add_write_floatField_CGNSmesh_vertexSolution('E',E)
            call add_write_floatField_CGNSmesh_vertexSolution('pr',pr)
            call add_write_floatField_CGNSmesh_vertexSolution('eta',eta)
            call add_write_floatField_CGNSmesh_vertexSolution('csound',csound)
            call add_write_floatField_CGNSmesh_vertexSolution('machno',machno)
            call add_write_floatField_CGNSmesh_vertexSolution('gradRhoX',gradRho(:,1))
            call add_write_floatField_CGNSmesh_vertexSolution('gradRhoY',gradRho(:,2))
            call add_write_floatField_CGNSmesh_vertexSolution('gradRhoZ',gradRho(:,3))
            call add_write_floatField_CGNSmesh_vertexSolution('curlUX',curlU(:,1))
            call add_write_floatField_CGNSmesh_vertexSolution('curlUY',curlU(:,2))
            call add_write_floatField_CGNSmesh_vertexSolution('curlUZ',curlU(:,3))
            call add_write_floatField_CGNSmesh_vertexSolution('divU',divU)
            call add_write_floatField_CGNSmesh_vertexSolution('Qcrit',Qcrit)
            call add_write_floatField_CGNSmesh_vertexSolution('mu_fluid',mu_fluid)
            call add_write_floatField_CGNSmesh_vertexSolution('envit',envit)
            call add_write_floatField_CGNSmesh_vertexSolution('mut',mut)
        endif

        call close_CGNSmesh_par()
    end do

    if(do_averages) then
        if(mpi_rank.eq.0) write(*,*) '# Doing final Avg...'

        favrho(:)   = favrho(:)   / real(numAvgSteps, rp)
        favpre(:)   = favpre(:)   / real(numAvgSteps, rp)
        favmueff(:) = favmueff(:) / real(numAvgSteps, rp)
        favvel(:,:) = favvel(:,:) / real(numAvgSteps, rp)
        favve2(:,:) = favve2(:,:) / real(numAvgSteps, rp)
    
        write(output_cgns_fileName,*) "results_finalAVG"
        iStep = 0
        call set_CGNS_full_fileName(full_fileName,cgns_filePath,cgns_fileName,output_cgns_fileName,iStep)

        call create_CGNSmesh_par(full_fileName)

        call add_write_floatField_CGNSmesh_vertexSolution('favrho',favrho)
        call add_write_floatField_CGNSmesh_vertexSolution('favpre',favpre)
        call add_write_floatField_CGNSmesh_vertexSolution('favmueff',favmueff)
        call add_write_floatField_CGNSmesh_vertexSolution('favvelX',favvel(:,1))
        call add_write_floatField_CGNSmesh_vertexSolution('favvelY',favvel(:,2))
        call add_write_floatField_CGNSmesh_vertexSolution('favvelZ',favvel(:,3))
        call add_write_floatField_CGNSmesh_vertexSolution('favve2X',favve2(:,1))
        call add_write_floatField_CGNSmesh_vertexSolution('favve2Y',favve2(:,2))
        call add_write_floatField_CGNSmesh_vertexSolution('favve2Z',favve2(:,3))

        call close_CGNSmesh_par()
    end if

    call end_CGNSmesh_arrays()

    call end_hdf5_interface()

    call end_mpi()

end program tool_hdf5_to_cgns
