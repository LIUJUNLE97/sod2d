module mod_comms_boundaries
    use mod_mpi_mesh
    implicit none

    !---- for the comms---
    integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
    !integer(KIND=MPI_ADDRESS_KIND) :: memPos_t
    integer :: bnd_worldGroup,bnd_commGroup

    integer(4), dimension(:), allocatable :: aux_bnd_intField_s, aux_bnd_intField_r
    real(4), dimension(:), allocatable :: aux_bnd_floatField_s, aux_bnd_floatField_r
    real(8), dimension(:), allocatable :: aux_bnd_doubleField_s, aux_bnd_doubleField_r

    integer :: window_id_bnd_int,window_id_bnd_float,window_id_bnd_double

    logical :: bnd_isInt,bnd_isFloat,bnd_isDouble

contains

    subroutine init_comms_bnd(useInt,useFloat,useDouble)
        implicit none
        logical, intent(in) :: useInt,useFloat,useDouble
        logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier

        bnd_isInt=.false.
        bnd_isFloat=.false.
        bnd_isDouble=.false.

        if(useInt) then
            bnd_isInt = .true.

            allocate(aux_bnd_intField_s(bnd_numNodesToComm))
            allocate(aux_bnd_intField_r(bnd_numNodesToComm))
            !$acc enter data create(aux_bnd_intField_s(:))
            !$acc enter data create(aux_bnd_intField_r(:))

            !call init_window_intField()
        end if

        if(useFloat) then
            bnd_isFloat = .true.

            allocate(aux_bnd_floatField_s(bnd_numNodesToComm))
            allocate(aux_bnd_floatField_r(bnd_numNodesToComm))
            !$acc enter data create(aux_bnd_floatField_s(:))
            !$acc enter data create(aux_bnd_floatField_r(:))

            !call init_window_floatField()
        end if

        if(useDouble) then
            bnd_isDouble = .true.

            allocate(aux_bnd_doubleField_s(numNodesToComm))
            allocate(aux_bnd_doubleField_r(numNodesToComm))
            !$acc enter data create(aux_bnd_doubleField_s(:))
            !$acc enter data create(aux_bnd_doubleField_r(:))       

            !call init_window_doubleField()
        end if

        call MPI_Comm_group(MPI_COMM_WORLD,bnd_worldGroup,mpi_err)
	     call MPI_Group_incl(bnd_worldGroup,bnd_numRanksWithComms,bnd_ranksToComm,bnd_commGroup,mpi_err);

        !useFenceFlags=.false. !by default
        !useAssertNoCheckFlags=.true. !by default
        !isPSCWBarrier=.true.!si faig molts loops amb aquesta opció a false la comm queda bloquejada
        !useLockBarrier=.true.!.false. with false it fails!
        !call setFenceFlags(useFenceFlags) 
        !call setPSCWAssertNoCheckFlags(useAssertNoCheckFlags)
        !call setLockBarrier(useLockBarrier)

    end subroutine init_comms_bnd

    subroutine end_comms_bnd()
        implicit none

        if(bnd_isInt) then
            !$acc exit data delete(aux_bnd_intField_s(:))
            !$acc exit data delete(aux_bnd_intField_r(:))
            deallocate(aux_bnd_intField_s)
            deallocate(aux_bnd_intField_r)

            call close_window_intField_bnd()            
        end if

        if(bnd_isFloat) then
           !$acc exit data delete(aux_bnd_floatField_s(:))
           !$acc exit data delete(aux_bnd_floatField_r(:))
            deallocate(aux_bnd_floatField_s)
            deallocate(aux_bnd_floatField_r)

            call close_window_floatField_bnd()
        end if

        if(bnd_isDouble) then
            !$acc exit data delete(aux_bnd_doubleField_s(:))
            !$acc exit data delete(aux_bnd_doubleField_r(:))
            deallocate(aux_bnd_doubleField_s)
            deallocate(aux_bnd_doubleField_r)

            call close_window_doubleField_bnd()
        end if

    end subroutine end_comms_bnd
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_intField_bnd()
        implicit none

        window_buffer_size = mpi_integer_size*bnd_numNodesToComm
        call MPI_Win_create(aux_bnd_intField_r,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_bnd_int,mpi_err)
    end subroutine init_window_intField_bnd

    subroutine close_window_intField_bnd()
        implicit none
        
        call MPI_Win_free(window_id_bnd_int,mpi_err)
    end subroutine close_window_intField_bnd
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_floatField_bnd()
        implicit none

        window_buffer_size = mpi_float_size*bnd_numNodesToComm
        call MPI_Win_create(aux_bnd_floatField_r,window_buffer_size,mpi_float_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_bnd_float,mpi_err)
    end subroutine init_window_floatField_bnd

    subroutine close_window_floatField_bnd()
        implicit none
        
        call MPI_Win_free(window_id_bnd_float,mpi_err)
    end subroutine close_window_floatField_bnd
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_doubleField_bnd()
        implicit none

        window_buffer_size = mpi_double_size*bnd_numNodesToComm
        call MPI_Win_create(aux_bnd_doubleField_r,window_buffer_size,mpi_double_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_bnd_double,mpi_err)
    end subroutine init_window_doubleField_bnd

    subroutine close_window_doubleField_bnd()
        implicit none
        
        call MPI_Win_free(window_id_bnd_double,mpi_err)
    end subroutine close_window_doubleField_bnd
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
 
!-----------------------------------------------------------------------------------------------------------------------
    subroutine fill_boundary_sendBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_matrixCommScheme(i,1)
            aux_bnd_intField_s(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_bnd_intField_r(:)=0.
        !$acc end kernels
    end subroutine fill_boundary_sendBuffer_int
!-------------------------------------------------------------------------
    subroutine fill_boundary_sendBuffer_float(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_matrixCommScheme(i,1)
            aux_bnd_floatField_s(i) = floatField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_bnd_floatField_r(:)=0.
        !$acc end kernels
    end subroutine fill_boundary_sendBuffer_float
!-------------------------------------------------------------------------
    subroutine fill_boundary_sendBuffer_double(doubleField)
        implicit none
        real(8), intent(inout) :: doubleField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_matrixCommScheme(i,1)
            aux_bnd_doubleField_s(i) = doubleField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_bnd_doubleField_r(:)=0.
        !$acc end kernels
    end subroutine fill_boundary_sendBuffer_double
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
    subroutine copy_from_min_boundary_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_matrixCommScheme(i,1)
            !$acc atomic update
            intField(iNodeL) = min(intField(iNodeL),aux_bnd_intField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_min_boundary_rcvBuffer_int
!-------------------------------------------------------------------------
    subroutine copy_from_max_boundary_rcvBuffer_float(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_matrixCommScheme(i,1)
            !$acc atomic update
            floatField(iNodeL) = max(floatField(iNodeL),aux_bnd_floatField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_max_boundary_rcvBuffer_float

!------------- min SEND/RECV -------------------------------------------
    ! INTEGER ---------------------------------------------------
    subroutine mpi_halo_min_boundary_update_int_iSendiRcv(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*bnd_numRanksWithComms)

        call fill_boundary_sendBuffer_int(intField)

        ireq=0
        do i=1,bnd_numRanksWithComms
            ngbRank  = bnd_ranksToComm(i)
            tagComm  = 0
            memPos_l = bnd_commsMemPosInLoc(i)
            memSize  = bnd_commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_bnd_intfield_r(mempos_l),memSize,MPI_INTEGER,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_bnd_intfield_s(mempos_l),memSize,MPI_INTEGER,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_min_boundary_rcvBuffer_int(intField)
    end subroutine mpi_halo_min_boundary_update_int_iSendiRcv
!------------- max SEND/RECV -------------------------------------------
    ! FLOAT ---------------------------------------------------
    subroutine mpi_halo_max_boundary_update_float_iSendiRcv(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*bnd_numRanksWithComms)

        call fill_boundary_sendBuffer_float(floatField)

        ireq=0
        do i=1,bnd_numRanksWithComms
            ngbRank  = bnd_ranksToComm(i)
            tagComm  = 0
            memPos_l = bnd_commsMemPosInLoc(i)
            memSize  = bnd_commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_bnd_floatfield_r(mempos_l),memSize,MPI_FLOAT,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_bnd_floatfield_s(mempos_l),memSize,MPI_FLOAT,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_max_boundary_rcvBuffer_float(floatField)
    end subroutine mpi_halo_max_boundary_update_float_iSendiRcv


end module mod_comms_boundaries
