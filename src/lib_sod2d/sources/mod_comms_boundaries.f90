module mod_comms_boundaries
    use mod_mpi_mesh

!-- Select type of communication for mpi_boundary_atomic_updates
#define _ISENDIRCV_ 1

    implicit none

    !---- for the comms---
    integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
    integer :: bnd_worldGroup,bnd_commGroup

    integer(4), dimension(:), allocatable :: aux_bnd_intField_s, aux_bnd_intField_r
    real(rp), dimension(:), allocatable :: aux_bnd_realField_s, aux_bnd_realField_r

    integer :: window_id_bnd_int,window_id_bnd_real

    logical :: bnd_isInt,bnd_isReal

contains

    subroutine init_comms_bnd(useInt,useReal)
        implicit none
        logical, intent(in) :: useInt,useReal
        logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier

#if _ISENDIRCV_
        write(*,*) "--| Boundary Comm. scheme: iSend-iRecv"
#endif

        bnd_isInt=.false.
        bnd_isReal=.false.

        if(useInt) then
            bnd_isInt = .true.

            allocate(aux_bnd_intField_s(bnd_numNodesToComm))
            allocate(aux_bnd_intField_r(bnd_numNodesToComm))
            !$acc enter data create(aux_bnd_intField_s(:))
            !$acc enter data create(aux_bnd_intField_r(:))

            !call init_window_intField()
        end if

        if(useReal) then
            bnd_isReal = .true.

            allocate(aux_bnd_realField_s(bnd_numNodesToComm))
            allocate(aux_bnd_realField_r(bnd_numNodesToComm))
            !$acc enter data create(aux_bnd_realField_s(:))
            !$acc enter data create(aux_bnd_realField_r(:))

            !call init_window_realField()
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

            !call close_window_intField_bnd()            
        end if

        if(bnd_isReal) then
           !$acc exit data delete(aux_bnd_realField_s(:))
           !$acc exit data delete(aux_bnd_realField_r(:))
            deallocate(aux_bnd_realField_s)
            deallocate(aux_bnd_realField_r)

            !call close_window_realField_bnd()
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
    subroutine init_window_realField_bnd()
        implicit none

        window_buffer_size = mpi_real_size*bnd_numNodesToComm
        call MPI_Win_create(aux_bnd_realField_r,window_buffer_size,mpi_real_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_bnd_real,mpi_err)
    end subroutine init_window_realField_bnd

    subroutine close_window_realField_bnd()
        implicit none
        
        call MPI_Win_free(window_id_bnd_real,mpi_err)
    end subroutine close_window_realField_bnd

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!   Standard fill boundary buffers
!-------------------------------------------------------------------------------------
    subroutine fill_boundary_sendBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            aux_bnd_intField_s(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_bnd_intField_r(:)=0
        !$acc end kernels
    end subroutine fill_boundary_sendBuffer_int
!-------------------------------------------------------------------------------------
    subroutine fill_boundary_sendBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            aux_bnd_realField_s(i) = realField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_bnd_realField_r(:)=0.
        !$acc end kernels
    end subroutine fill_boundary_sendBuffer_real

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!   Standard copy boundary buffers
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine copy_from_boundary_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            !$acc atomic update
            intField(iNodeL) = intField(iNodeL) + aux_bnd_intField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_boundary_rcvBuffer_int
!-------------------------------------------------------------------------
    subroutine copy_from_boundary_rcvBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = realField(iNodeL) + aux_bnd_realField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_boundary_rcvBuffer_real

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!   Specialized copy boundary buffers
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine copy_from_min_boundary_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            !$acc atomic update
            intField(iNodeL) = min(intField(iNodeL),aux_bnd_intField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_min_boundary_rcvBuffer_int
!-------------------------------------------------------------------------
    subroutine copy_from_max_boundary_rcvBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = max(realField(iNodeL),aux_bnd_realField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_max_boundary_rcvBuffer_real
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

    subroutine mpi_halo_boundary_atomic_update_int(intField)
        implicit none
        integer, intent(inout) :: intField(:)

#if _ISENDIRCV_
        call mpi_halo_boundary_atomic_update_int_iSendiRcv(intField)
#endif

    end subroutine mpi_halo_boundary_atomic_update_int

    subroutine mpi_halo_boundary_atomic_update_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)

#if _ISENDIRCV_
        call mpi_halo_boundary_atomic_update_real_iSendiRcv(realField)
#endif

    end subroutine mpi_halo_boundary_atomic_update_real

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

!------------- ISEND/IRECV -------------------------------------------
    !INTEGER ---------------------------------------------------------
    subroutine mpi_halo_boundary_atomic_update_int_iSendiRcv(intField)
        implicit none
        integer, intent(inout) :: intField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*bnd_numRanksWithComms)

        call fill_boundary_sendBuffer_int(intField)

        ireq=0
        !$acc host_data use_device (aux_bnd_intfield_r(:),aux_bnd_intfield_s(:))
        do i=1,bnd_numRanksWithComms
            ngbRank  = bnd_ranksToComm(i)
            tagComm  = 0
            memPos_l = bnd_commsMemPosInLoc(i)
            memSize  = bnd_commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_bnd_intfield_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_bnd_intfield_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_boundary_rcvBuffer_int(intField)
    end subroutine mpi_halo_boundary_atomic_update_int_iSendiRcv
    !REAL ---------------------------------------------------------
    subroutine mpi_halo_boundary_atomic_update_real_iSendiRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*numRanksWithComms)

        call fill_boundary_sendBuffer_real(realField)

        ireq=0
        !$acc host_data use_device (aux_bnd_realfield_r(:),aux_bnd_realfield_s(:))
        do i=1,bnd_numRanksWithComms
            ngbRank  = bnd_ranksToComm(i)
            tagComm  = 0
            memPos_l = bnd_commsMemPosInLoc(i)
            memSize  = bnd_commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_bnd_realfield_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_bnd_realfield_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_boundary_rcvBuffer_real(realField)
    end subroutine mpi_halo_boundary_atomic_update_real_iSendiRcv

!------------- MIN iSENDiRECV -------------------------------------------
    ! INTEGER ---------------------------------------------------
    subroutine mpi_halo_min_boundary_update_int_iSendiRcv(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*bnd_numRanksWithComms)

        call fill_boundary_sendBuffer_int(intField)

        ireq=0
        !$acc host_data use_device (aux_bnd_intfield_r(:),aux_bnd_intfield_s(:))
        do i=1,bnd_numRanksWithComms
            ngbRank  = bnd_ranksToComm(i)
            tagComm  = 0
            memPos_l = bnd_commsMemPosInLoc(i)
            memSize  = bnd_commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_bnd_intfield_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_bnd_intfield_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_min_boundary_rcvBuffer_int(intField)
    end subroutine mpi_halo_min_boundary_update_int_iSendiRcv
!------------- MAX iSENDiRECV -------------------------------------------
    ! REAL ---------------------------------------------------
    subroutine mpi_halo_max_boundary_update_real_iSendiRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*bnd_numRanksWithComms)

        call fill_boundary_sendBuffer_real(realField)

        ireq=0
        !$acc host_data use_device (aux_bnd_realfield_r(:),aux_bnd_realfield_s(:))
        do i=1,bnd_numRanksWithComms
            ngbRank  = bnd_ranksToComm(i)
            tagComm  = 0
            memPos_l = bnd_commsMemPosInLoc(i)
            memSize  = bnd_commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_bnd_realfield_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_bnd_realfield_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_max_boundary_rcvBuffer_real(realField)
    end subroutine mpi_halo_max_boundary_update_real_iSendiRcv


end module mod_comms_boundaries
