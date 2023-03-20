module mod_comms
    use mod_mpi_mesh

!-- Select type of communication for mpi_atomic_updates
#define _SENDRCV_ 0
#define _ISENDIRCV_ 1
#define _PUTFENCE_ 0

    implicit none

    !---- for the comms---
    integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
    integer(KIND=MPI_ADDRESS_KIND) :: memPos_t
    integer :: worldGroup,commGroup

    integer(4),dimension(:),allocatable :: aux_intField_s,  aux_intField_r
    real(rp),dimension(:),allocatable   :: aux_realField_s, aux_realField_r
    !real(4), dimension(:), allocatable :: aux_floatField_s, aux_floatField_r
    !real(8), dimension(:), allocatable :: aux_doubleField_s, aux_doubleField_r
    !real(4), dimension(:), allocatable :: aux_floatField_5s, aux_floatField_5r
    !using two buffers because if we use only one maybe we send the info after other proc has added something in me!
    !then i'll duplicate the info :S

    integer :: window_id_int,window_id_real!,window_id_float,window_id_double
    !integer :: window_id_float5
    integer :: window_id_sm
    integer :: beginFence=0,endFence=0
    integer :: startAssert=0,postAssert=0

    logical :: isInt,isReal!isFloat,isDouble,isFloat5
    logical :: isLockBarrier,isPSCWBarrier

    !integer :: ms_rank,ms_size,ms_newComm
    type(c_ptr) :: c_ms_ptr

contains

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
    subroutine init_comms(useInt,useReal)
        implicit none
        logical, intent(in) :: useInt,useReal!Float,useDouble
        logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier
        !logical :: useFloat5=.false. !i think that will dissapear... but...

#if _SENDRCV_
        write(111,*) "--| Comm. scheme: Send-Recv"
#endif
#if _ISENDIRCV_
        write(111,*) "--| Comm. scheme: iSend-iRecv"
#endif
#if _PUTFENCE_
        write(111,*) "--| Comm. scheme: Put(Fence)"
#endif

        isInt=.false.
        isReal=.false.
        !isFloat=.false.
        !isDouble=.false.
        !isFloat5=.false. !i think that will dissapear... but...

        if(useInt) then
            isInt = .true.

            allocate(aux_intField_s(numNodesToComm))
            allocate(aux_intField_r(numNodesToComm))
            !$acc enter data create(aux_intField_s(:))
            !$acc enter data create(aux_intField_r(:))

            call init_window_intField()
        end if

        if(useReal) then
            isReal = .true.

            allocate(aux_RealField_s(numNodesToComm))
            allocate(aux_RealField_r(numNodesToComm))
            !$acc enter data create(aux_RealField_s(:))
            !$acc enter data create(aux_RealField_r(:))

            call init_window_realField()
        end if
#if 0
        if(useFloat) then
            isFloat = .true.

            allocate(aux_floatField_s(numNodesToComm))
            allocate(aux_floatField_r(numNodesToComm))
            !$acc enter data create(aux_floatField_s(:))
            !$acc enter data create(aux_floatField_r(:))

            call init_window_floatField()
        end if

        if(useDouble) then
            isDouble = .true.

            allocate(aux_doubleField_s(numNodesToComm))
            allocate(aux_doubleField_r(numNodesToComm))
            !$acc enter data create(aux_doubleField_s(:))
            !$acc enter data create(aux_doubleField_r(:))       

            call init_window_doubleField()
        end if

        if(useFloat5) then
            isFloat5 = .true.

            allocate(aux_floatField_5s(5*numNodesToComm))
            allocate(aux_floatField_5r(5*numNodesToComm))
            !$acc enter data create(aux_floatField_5s(:))
            !$acc enter data create(aux_floatField_5r(:))

            call init_window_floatField5()
        end if
#endif
        call MPI_Comm_group(MPI_COMM_WORLD,worldGroup,mpi_err)
	    call MPI_Group_incl(worldGroup,numRanksWithComms,ranksToComm,commGroup,mpi_err);

        useFenceFlags=.false. !by default
        useAssertNoCheckFlags=.true. !by default
        isPSCWBarrier=.true.!si faig molts loops amb aquesta opció a false la comm queda bloquejada
        useLockBarrier=.true.!.false. with false it fails!
        call setFenceFlags(useFenceFlags) 
        call setPSCWAssertNoCheckFlags(useAssertNoCheckFlags)
        call setLockBarrier(useLockBarrier)

    end subroutine init_comms

    subroutine end_comms()
        implicit none

        if(isInt) then
            !$acc exit data delete(aux_intField_s(:))
            !$acc exit data delete(aux_intField_r(:))
            deallocate(aux_intField_s)
            deallocate(aux_intField_r)

            call close_window_intField()            
        end if

        if(isreal) then
           !$acc exit data delete(aux_realField_s(:))
           !$acc exit data delete(aux_realField_r(:))
            deallocate(aux_realField_s)
            deallocate(aux_realField_r)

            call close_window_realField()
        end if

    end subroutine end_comms
    
    subroutine setFenceFlags(useFenceFlags)
        implicit none
        logical,intent(in) :: useFenceFlags

        if(useFenceFlags) then
            beginFence = MPI_MODE_NOPRECEDE
	        endFence   = IOR(MPI_MODE_NOSTORE,IOR(MPI_MODE_NOPUT,MPI_MODE_NOSUCCEED))
	        !endFence   = IOR(MPI_MODE_NOSTORE,MPI_MODE_NOPUT)
        else
            beginFence = 0
	        endFence   = 0
        end if
    end subroutine

    subroutine setPSCWAssertNoCheckFlags(useAssertNoCheckFlags)
        implicit none
        logical,intent(in) :: useAssertNoCheckFlags

        if(useAssertNoCheckFlags) then
	        postAssert  = MPI_MODE_NOCHECK
	        startAssert = MPI_MODE_NOCHECK
        else
	        postAssert  = 0
	        startAssert = 0
        endif
    end subroutine

    subroutine setLockBarrier(useLockBarrier)
        implicit none
        logical,intent(in) :: useLockBarrier

        isLockBarrier = useLockBarrier
    end subroutine
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_intField()
        implicit none

        window_buffer_size = mpi_integer_size*numNodesToComm
        call MPI_Win_create(aux_intField_r,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_int,mpi_err)
    end subroutine init_window_intField

    subroutine close_window_intField()
        implicit none
        
        call MPI_Win_free(window_id_int,mpi_err)
    end subroutine close_window_intField
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_realField()
        implicit none

        window_buffer_size = mpi_real_size*numNodesToComm
        call MPI_Win_create(aux_realField_r,window_buffer_size,mpi_real_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_real,mpi_err)
    end subroutine init_window_realField

    subroutine close_window_realField()
        implicit none
        
        call MPI_Win_free(window_id_real,mpi_err)
    end subroutine close_window_realField
#if 0
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_floatField()
        implicit none

        window_buffer_size = mpi_float_size*numNodesToComm
        call MPI_Win_create(aux_floatField_r,window_buffer_size,mpi_float_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_float,mpi_err)
    end subroutine init_window_floatField

    subroutine close_window_floatField()
        implicit none
        
        call MPI_Win_free(window_id_float,mpi_err)
    end subroutine close_window_floatField
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_doubleField()
        implicit none

        window_buffer_size = mpi_double_size*numNodesToComm
        call MPI_Win_create(aux_doubleField_r,window_buffer_size,mpi_double_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_double,mpi_err)
    end subroutine init_window_doubleField

    subroutine close_window_doubleField()
        implicit none
        
        call MPI_Win_free(window_id_double,mpi_err)
    end subroutine close_window_doubleField
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_floatField5()
        implicit none

        window_buffer_size = 5*mpi_float_size*numNodesToComm
        call MPI_Win_create(aux_floatField_5r,window_buffer_size,mpi_float_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_float5,mpi_err)
    end subroutine init_window_floatField5

    subroutine close_window_floatField5()
        implicit none
        
        call MPI_Win_free(window_id_float5,mpi_err)
    end subroutine close_window_floatField5
!-------------------------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!    SUBROUTINES COPY/FROM SEND/RCV BUFFERS
!-----------------------------------------------------------------------------------------------------------------------
    subroutine fill_sendBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_intField_s(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_intField_r(:)=0
        !$acc end kernels
    end subroutine fill_sendBuffer_int
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_real(realField)
        implicit none
        real(rp),intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_realField_s(i) = realField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_realField_r(:)=0.
        !$acc end kernels
    end subroutine fill_sendBuffer_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_get_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_intField_r(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_intField_s(:)=0
        !$acc end kernels
    end subroutine fill_sendBuffer_get_int
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_get_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_realField_r(i) = realField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_realField_s(:)=0.
        !$acc end kernels
    end subroutine fill_sendBuffer_get_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            intField(iNodeL) = intField(iNodeL) + aux_intField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_int
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = realField(iNodeL) + aux_realField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_get_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            intField(iNodeL) = intField(iNodeL) + aux_intField_s(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_get_int
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_get_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = realField(iNodeL) + aux_realField_s(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_get_real
!-----------------------------------------------------------------------------------------------------------------------
    subroutine copy_from_conditional_ave_rcvBuffer_real(cond,realField)
        implicit none
        real(rp),intent(in) :: cond
        real(rp),intent(inout) :: realField(:)
        integer :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
           iNodeL = nodesToComm(i)
           if(abs(aux_realField_r(i)).gt. cond) then 
              if(abs(realField(iNodeL)).gt. cond) then 
                 !$acc atomic update
                 realField(iNodeL) = realField(iNodeL)+aux_realField_r(i)
                 !$acc end atomic
              else
                 !$acc atomic write
                 realField(iNodeL) = aux_realField_r(i)
                 !$acc end atomic
              end if
           end if
        end do
        !$acc end parallel loop
    end subroutine copy_from_conditional_ave_rcvBuffer_real
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

    subroutine mpi_halo_atomic_update_int(intField)
        implicit none
        integer, intent(inout) :: intField(:)
#if _SENDRCV_
        call mpi_halo_atomic_update_int_sendRcv(intField)
#endif
#if _ISENDIRCV_
        call mpi_halo_atomic_update_int_iSendiRcv(intField)
#endif
#if _PUTFENCE_
        call mpi_halo_atomic_update_int_put_fence(intField)
#endif
    end subroutine mpi_halo_atomic_update_int

    subroutine mpi_halo_atomic_update_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)

#if _SENDRCV_
        call mpi_halo_atomic_update_real_sendRcv(realField)
#endif
#if _ISENDIRCV_
        call mpi_halo_atomic_update_real_iSendiRcv(realField)
#endif
#if _PUTFENCE_
        call mpi_halo_atomic_update_real_put_fence(realField)
#endif

    end subroutine mpi_halo_atomic_update_real

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!------------- SEND/RECV -------------------------------------------
    ! INTEGER ---------------------------------------------------
    subroutine mpi_halo_atomic_update_int_sendRcv(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ngbRank,tagComm
        integer :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_intfield_s(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              aux_intfield_r(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        end do

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_sendRcv
    ! REAL ---------------------------------------------------
    subroutine mpi_halo_atomic_update_real_sendRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank,tagComm
        integer :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_realfield_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_realfield_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        end do

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_sendRcv

!------------- ISEND/IRECV -------------------------------------------
    !INTEGER ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_iSendiRcv(intField)
        implicit none
        integer, intent(inout) :: intField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*numRanksWithComms)

        call fill_sendBuffer_int(intField)

        ireq=0
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_intfield_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_intfield_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_iSendiRcv
    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_iSendiRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*numRanksWithComms)

        call fill_sendBuffer_real(realField)

        ireq=0
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_realfield_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realfield_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_iSendiRcv

!------------- conditional average ISEND/IRECV -------------------------------------------
    ! REAL ---------------------------------------------------
    subroutine mpi_halo_conditional_ave_update_real_iSendiRcv(cond,realField)
        implicit none
        real(rp),intent(in) :: cond
        real(rp),intent(inout) :: realField(:)
        integer :: i,ireq,ngbRank,tagComm
        integer :: memPos_l,memSize
        integer :: requests(2*numRanksWithComms)

        call fill_sendBuffer_real(realField)

        ireq=0
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_realfield_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realfield_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,MPI_COMM_WORLD,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_conditional_ave_rcvBuffer_real(cond,realField)
    end subroutine mpi_halo_conditional_ave_update_real_iSendiRcv
!------------- PUT FENCE -------------------------------------------
    !INT
    subroutine mpi_halo_atomic_update_int_put_fence(intField)
        implicit none
        integer, intent(inout) :: intField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        call MPI_Win_fence(beginFence,window_id_int,mpi_err)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do

        call MPI_Win_fence(endFence,window_id_int,mpi_err)

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_put_fence
    !REAL
    subroutine mpi_halo_atomic_update_real_put_fence(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        call MPI_Win_fence(beginFence,window_id_real,mpi_err)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do

        call MPI_Win_fence(endFence,window_id_real,mpi_err)

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_fence
!------------- PUT PSCW -------------------------------------------
    !INTEGER--------------------------------------------------
    subroutine mpi_halo_atomic_update_int_put_pscw(intField)
        implicit none
        integer, intent(inout) :: intField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        if(isPSCWBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err);

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do

	    call MPI_Win_complete(window_id_int,mpi_err);
	    call MPI_Win_wait(window_id_int,mpi_err);

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_put_pscw
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_put_pscw(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        if(isPSCWBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err);

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do

	    call MPI_Win_complete(window_id_real,mpi_err);
	    call MPI_Win_wait(window_id_real,mpi_err);

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_pscw

!------------- PUT LOCK -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_put_lock(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)
		
            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
        end do

        if(isLockBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_put_lock
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_put_lock(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)
		
            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
        end do

        if(isLockBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_lock

!------------- GET FENCE -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_get_fence(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_get_int(intField)

        call MPI_Win_fence(beginFence,window_id_int,mpi_err)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do

        call MPI_Win_fence(endFence,window_id_int,mpi_err)

        call copy_from_rcvBuffer_get_int(intField)
    end subroutine mpi_halo_atomic_update_int_get_fence
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_get_fence(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        call MPI_Win_fence(beginFence,window_id_real,mpi_err)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do

        call MPI_Win_fence(endFence,window_id_real,mpi_err)

        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_fence

!------------- GET PSCW -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_get_pscw(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_get_int(intField)

        if(isPSCWBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err);

        do i=1,numRanksWithComms
            ngbRank   = ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize   = commsMemSize(i)

            call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do

	    call MPI_Win_complete(window_id_int,mpi_err);
	    call MPI_Win_wait(window_id_int,mpi_err);

        call copy_from_rcvBuffer_get_int(intField)
    end subroutine mpi_halo_atomic_update_int_get_pscw
    !FLOAT-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_get_pscw(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        if(isPSCWBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err);

        do i=1,numRanksWithComms
            ngbRank   = ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize   = commsMemSize(i)

            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do

	    call MPI_Win_complete(window_id_real,mpi_err);
	    call MPI_Win_wait(window_id_real,mpi_err);

        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_pscw

!------------- GET LOCK -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_get_lock(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_get_int(intField)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)
		
            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
        end do

        if(isLockBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

        call copy_from_rcvBuffer_get_int(intField)
    end subroutine mpi_halo_atomic_update_int_get_lock
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_get_lock(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer :: i,ngbRank
        integer :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)
		
            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
        end do

        if(isLockBarrier) call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_lock

!------------- ONLY BUFFERS -------------------------------------------
!for testing and devel stuff
    !INTEGER ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_onlybuffers(intfield) 
        implicit none
        integer(4), intent(inout) :: intfield(:)
        integer :: i,ngbrank,tagcomm
        integer :: mempos_l,memsize

        call fill_sendbuffer_int(intfield)
        call copy_from_rcvbuffer_int(intfield)

    end subroutine mpi_halo_atomic_update_int_onlybuffers
    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_onlybuffers(realfield) 
        implicit none
        real(rp), intent(inout) :: realfield(:)
        integer :: i,ngbrank,tagcomm
        integer :: mempos_l,memsize

        call fill_sendbuffer_real(realfield)
        call copy_from_rcvbuffer_real(realfield)

    end subroutine mpi_halo_atomic_update_real_onlybuffers

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!------------- OLD FUNCTIONS & DEPRECATED - TO BE DELETED --------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
#if 0
    subroutine update_and_comm_doubleField(doubleField)
        implicit none
        real(8), intent(inout) :: doubleField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_doubleField_s(i) = doubleField(iNodeL)
        end do
        aux_doubleField_r(:)=0.

        !window_buffer_size = mpi_double_size*numNodesToComm
        !call MPI_Win_create(aux_doubleField_r,window_buffer_size,mpi_double_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_double,mpi_err)
        call MPI_Win_fence(0,window_id_double,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            !write(*,*) '[',mpi_rank,']iRank->',iRank
            !write(*,*) '[',mpi_rank,']memPos_l->',memPos_l
            !write(*,*) '[',mpi_rank,']memPos_t->',memPos_t
            !write(*,*) '[',mpi_rank,']memSize->',memSize

            call MPI_Accumulate(aux_doubleField_s(memPos_l),memSize,mpi_datatype_real8,&
                                iRank,memPos_t,memSize,mpi_datatype_real8,MPI_SUM,window_id_double,mpi_err)
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_double,mpi_err)
        !call MPI_Win_free(window_id_double,mpi_err)

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !$acc atomic update
            doubleField(iNodeL) = doubleField(iNodeL) + aux_doubleField_r(i) 
            !$acc end atomic
        end do

    end subroutine update_and_comm_doubleField

    subroutine update_and_comm_floatField(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize
 
        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = floatField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_floatField_r(:)=0.
        !$acc end kernels

        call MPI_Win_fence(0,window_id_float,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            !write(*,*) '[',mpi_rank,']iRank->',iRank
            !write(*,*) '[',mpi_rank,']memPos_l->',memPos_l
            !write(*,*) '[',mpi_rank,']memPos_t->',memPos_t
            !write(*,*) '[',mpi_rank,']memSize->',memSize

            call MPI_Put(aux_floatField_s(memPos_l),memSize,mpi_datatype_real,iRank,memPos_t,memSize,mpi_datatype_real,window_id_float,mpi_err)
            !call MPI_Accumulate(aux_floatField_s(memPos_l),memSize,mpi_datatype_real,&
            !                    iRank,memPos_t,memSize,mpi_datatype_real,MPI_SUM,window_id_float,mpi_err)
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_float,mpi_err)

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !$acc atomic update
            floatField(iNodeL) = floatField(iNodeL) + aux_floatField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop

    end subroutine update_and_comm_floatField

    subroutine update_and_comm_floatField_all(numNodesInRank,fmass,fener,fmomx,fmomy,fmomz)
    !TESTEJAR FER LA COMUNICACIO DE TOTS ELS FIELDS AQUI A LA VEGADA, AMB UN SOL PUT
    !I NO FER 6 COMUNICACIONS DIFERENTS
    !POT SER ALGO HARDCODEJAT PERO POT VALER LA PENA...
    !MIRAR BE PERQUE AIXO TE MIGA SI ES VOL FER ELS 5 ARRAYS DE COP!
    !CREC QUE NO ES TAN STRAIGHTFORWARD....
        implicit none
        integer :: numNodesInRank
        real(4), intent(inout) :: fmass(numNodesInRank),fener(numNodesInRank)
        real(4), intent(inout) :: fmomx(numNodesInRank),fmomy(numNodesInRank),fmomz(numNodesInRank)
        integer :: i,iRank,iNodeL,ngbRank,tagComm
        integer :: memPos_l,memSize,numFields,origMemPos,newMemPos
        
        numFields=5

        !$acc kernels
        do iRank=1,numRanksWithComms
            ngbRank=ranksToComm(iRank)
            
            origMemPos = commsMemPosInLoc(iRank)-1
            newMemPos  = (commsMemPosInLoc(iRank)-1)*numFields
            memSize    = commsMemSize(iRank)

            do i=1,memSize
                iNodeL = matrixCommScheme(origMemPos+i,1)

                aux_floatField_5s(newMemPos+i)           = fmass(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*1) = fener(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*2) = fmomx(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*3) = fmomy(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*4) = fmomz(iNodeL)
            end do 
        end do

        aux_floatField_5r(:)=0.
        !$acc end kernels

        call MPI_Win_fence(0,window_id_float5,mpi_err)

        do i=1,numRanksWithComms
            ngbRank=ranksToComm(i)

            memPos_l  = (commsMemPosInLoc(i)-1)*numFields + 1
            memPos_t  = (commsMemPosInNgb(i)-1)*numFields
            memSize   = commsMemSize(i)*numFields

            call MPI_Put(aux_floatField_5s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_float5,mpi_err)
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_float5,mpi_err)

        !$acc kernels
        do iRank=1,numRanksWithComms
            origMemPos = commsMemPosInLoc(iRank)-1
            newMemPos  = (commsMemPosInLoc(iRank)-1)*numFields
            memSize    = commsMemSize(iRank)

            do i=1,memSize
                iNodeL = matrixCommScheme(origMemPos+i,1)
                !$acc atomic update
                fmass(iNodeL) = fmass(iNodeL) + aux_floatField_5r(newMemPos+i)
                !$acc end atomic
                !$acc atomic update
                fener(iNodeL) = fener(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*1)
                !$acc end atomic
                !$acc atomic update
                fmomx(iNodeL) = fmomx(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*2)
                !$acc end atomic
                !$acc atomic update
                fmomy(iNodeL) = fmomy(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*3)
                !$acc end atomic
                !$acc atomic update
                fmomz(iNodeL) = fmomz(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*4)
                !$acc end atomic
            end do 
        end do
        !$acc end kernels

    end subroutine update_and_comm_floatField_all

    subroutine sendRcv_floatField(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL,ngbRank,tagComm
        integer :: memPos_l,memSize
 
        !$acc data present(aux_floatField_s(:),aux_floatField_r(:))
        !$acc parallel loop 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = floatField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_floatField_r(:)=0.
        !$acc end kernels

        do i=1,numRanksWithComms
            ngbRank=ranksToComm(i)
            tagComm=0
            memPos_l  = commsMemPosInLoc(i)
            memSize = commsMemSize(i)
            !memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            !write(*,*) '[',mpi_rank,']iRank->',iRank
            !write(*,*) '[',mpi_rank,']memPos_l->',memPos_l
            !write(*,*) '[',mpi_rank,']memSize->',memSize
          !$acc host_data use_device (aux_floatField_s,aux_floatField_r)
            call MPI_Sendrecv(aux_floatfield_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_floatfield_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
!            call MPI_Sendrecv(aux_floatField_s, memSize, mpi_datatype_real, ngbRank, tagComm, &
!                              aux_floatField_r, memSize, mpi_datatype_real, ngbRank, tagComm, &
!                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
          !$acc end host_data
        end do

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !$acc atomic update
            floatField(iNodeL) = floatField(iNodeL) + aux_floatField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
        !$acc end data

    end subroutine sendRcv_floatField

    subroutine sendRcv_floatField_noGPU(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL,ngbRank,tagComm
        integer :: memPos_l,memSize
 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = floatField(iNodeL)
        end do
        aux_floatField_r(:)=0.

        do i=1,numRanksWithComms
            ngbRank=ranksToComm(i)
            tagComm=0
            memPos_l  = commsMemPosInLoc(i)
            memSize = commsMemSize(i)

            call MPI_Sendrecv(aux_floatfield_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_floatfield_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        end do

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            floatField(iNodeL) = floatField(iNodeL) + aux_floatField_r(i)
        end do

    end subroutine sendRcv_floatField_noGPU

    subroutine sendRcv_floatField_devel(floatField)
        implicit none
        real(4), intent(inout) :: floatField(numNodesRankPar)
        real(4) :: aux_ff_s(numNodesToComm), aux_ff_r(numNodesToComm)
        integer :: i,iNodeL,ngbRank,tagComm
        integer :: memPos_l,memSize
 
        !$acc parallel loop copyin(floatField(:),matrixCommScheme(:,:)) copyout(aux_ff_s(:))
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_ff_s(i) = floatField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels copyout(aux_ff_r(:))
        aux_ff_r(:)=0.0_rp
        !$acc end kernels
        call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
        do i=1,numRanksWithComms
            ngbRank=ranksToComm(i)
            tagComm=0
            memPos_l  = commsMemPosInLoc(i)
            memSize = commsMemSize(i)

            call MPI_Sendrecv(aux_ff_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_ff_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        end do
        call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
        !$acc parallel loop copyin(aux_ff_r(:),matrixCommScheme(:,:)) copy(floatField(:)) 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !$acc atomic update
            floatField(iNodeL) = floatField(iNodeL) + aux_ff_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine sendRcv_floatField_devel

    subroutine sendRcv_floatField_all(numNodesInRank,fmass,fener,fmomx,fmomy,fmomz)
        implicit none
        integer :: numNodesInRank
        real(4), intent(inout) :: fmass(numNodesInRank),fener(numNodesInRank)
        real(4), intent(inout) :: fmomx(numNodesInRank),fmomy(numNodesInRank),fmomz(numNodesInRank)
        integer :: i,iRank,iNodeL,ngbRank,tagComm
        integer :: memPos_l,memSize,numFields,origMemPos,newMemPos
 
        numFields=5

        !$acc kernels
        do iRank=1,numRanksWithComms
            ngbRank=ranksToComm(iRank)
            
            origMemPos = commsMemPosInLoc(iRank)-1
            newMemPos  = (commsMemPosInLoc(iRank)-1)*numFields
            memSize    = commsMemSize(iRank)

            do i=1,memSize
                iNodeL = matrixCommScheme(origMemPos+i,1)

                aux_floatField_5s(newMemPos+i)           = fmass(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*1) = fener(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*2) = fmomx(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*3) = fmomy(iNodeL)
                aux_floatField_5s(newMemPos+i+memSize*4) = fmomz(iNodeL)
            end do 
        end do

        aux_floatField_5r(:)=0.
        !$acc end kernels

        do i=1,numRanksWithComms
            ngbRank=ranksToComm(i)
            tagComm=0
            memPos_l  = (commsMemPosInLoc(i)-1)*numFields + 1
            memSize   = commsMemSize(i)*numFields

            call MPI_Sendrecv(aux_floatField_5s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_floatField_5r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        end do

        !$acc kernels
        do iRank=1,numRanksWithComms
            origMemPos = commsMemPosInLoc(iRank)-1
            newMemPos  = (commsMemPosInLoc(iRank)-1)*numFields
            memSize    = commsMemSize(iRank)

            do i=1,memSize
                iNodeL = matrixCommScheme(origMemPos+i,1)
                !$acc atomic update
                fmass(iNodeL) = fmass(iNodeL) + aux_floatField_5r(newMemPos+i)
                !$acc end atomic
                !$acc atomic update
                fener(iNodeL) = fener(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*1)
                !$acc end atomic
                !$acc atomic update
                fmomx(iNodeL) = fmomx(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*2)
                !$acc end atomic
                !$acc atomic update
                fmomy(iNodeL) = fmomy(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*3)
                !$acc end atomic
                !$acc atomic update
                fmomz(iNodeL) = fmomz(iNodeL) + aux_floatField_5r(newMemPos+i+memSize*4)
                !$acc end atomic
            end do 
        end do
        !$acc end kernels

    end subroutine sendRcv_floatField_all

    subroutine sendRcv_intField(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL,ngbRank,tagComm
        integer :: memPos_l,memSize
 
        !$acc data present(aux_intField_s(:),aux_intField_r(:))
        !$acc parallel loop 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_intField_s(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels
        aux_intField_r(:)=0.
        !$acc end kernels

        do i=1,numRanksWithComms
            ngbRank=ranksToComm(i)
            tagComm=0
            memPos_l  = commsMemPosInLoc(i)
            memSize = commsMemSize(i)
          !$acc host_data use_device (aux_floatField_s,aux_floatField_r)
            call MPI_Sendrecv(aux_intField_s(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              aux_intField_r(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
          !$acc end host_data
        end do

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !$acc atomic update
            intField(iNodeL) = intField(iNodeL) + aux_intField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
        !$acc end data

    end subroutine sendRcv_intField

    subroutine update_and_comm_intField(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize
 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_intField_s(i) = intField(iNodeL)
        end do
        aux_intField_r(:)=0

        call MPI_Win_fence(0,window_id_int,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,iRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_int,mpi_err)

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !$acc atomic update
            intField(iNodeL) = intField(iNodeL) + aux_intField_r(i)
            !$acc end atomic
        end do

    end subroutine update_and_comm_intField

    subroutine init_shared_mem_window()
        implicit none

        !1.creating new world with shared mem
        !if all procs have shared mem we do not need it but... keep it for now
        !call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,ms_newComm,mpi_err)
        !call mpi_comm_rank(ms_newComm, ms_rank, mpi_err)
        !call mpi_comm_size(ms_newComm, ms_size, mpi_err)

        !if (ms_newComm .ne. MPI_COMM_NULL) then
        !    write(*,*) 'New comm created! old rank', mpi_rank, ' ms_rank ', ms_rank,' ms_size ',ms_size
        !end if
        
        window_buffer_size = mpi_float_size*numNodesToComm
        call MPI_Win_allocate_shared(window_buffer_size,mpi_float_size,MPI_INFO_NULL,&
                                     smNode_comm,c_ms_ptr,window_id_sm,mpi_err)

    end subroutine init_shared_mem_window

    subroutine close_shared_mem_windows()
        implicit none
        call MPI_Win_free(window_id_sm, mpi_err)
    end subroutine close_shared_mem_windows

    subroutine update_and_comm_shared_mem_floatField(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,j,pos_l,pos_t,iNodeL,iRank
        integer :: memPos_l,memSize,memDisp,size_p
        integer(KIND=MPI_ADDRESS_KIND) :: ssize
        real(4), pointer :: f_ms_ptr(:)
        real(4) :: addVal

        size_p = numNodesToComm*mpi_size

        call MPI_Win_shared_query(window_id_sm,mpi_rank,ssize,memDisp,c_ms_ptr,mpi_err)
        call c_f_pointer(c_ms_ptr, f_ms_ptr, SHAPE = [size_p])

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = floatField(iNodeL)
            f_ms_ptr(i) = 0.
        end do

        call MPI_Win_fence(0,window_id_sm,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) !- 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            call MPI_Win_shared_query(window_id_sm,iRank,ssize,memDisp,c_ms_ptr,mpi_err)
            !call MPI_Win_shared_query(window_id_sm,0,ssize,memDisp,c_ms_ptr,mpi_err)
            !write(*,*) 'ms[',mpi_rank,'] ssize',ssize,' memDisp ',memDisp
            call c_f_pointer(c_ms_ptr, f_ms_ptr, SHAPE = [size_p])

            do j=1,memSize
                pos_l=j+memPos_l-1
                pos_t=j+memPos_t-1
                f_ms_ptr(pos_t) = f_ms_ptr(pos_t)+aux_floatField_s(pos_l)
            end do
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_sm,mpi_err)

        call MPI_Win_shared_query(window_id_sm,mpi_rank,ssize,memDisp,c_ms_ptr,mpi_err)
        call c_f_pointer(c_ms_ptr, f_ms_ptr, SHAPE = [size_p])
        !write(*,*) 'ms[',mpi_rank,'] ssize',ssize,' memDisp ',memDisp
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            floatField(iNodeL) = floatField(iNodeL) + f_ms_ptr(i)
        end do

    end subroutine update_and_comm_shared_mem_floatField
#endif

!-------------------------------------------------------------------------------
! ####      OTHER FUNCS      ---------------------------------------------------
#if 0
    subroutine test_multiupdate_and_comm_floatField(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize
 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = (mpi_rank+1)*0.1!floatField(iNodeL)
        end do
        aux_floatField_r(:)=0.

        call MPI_Win_fence(0,window_id_float,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = 1!commsMemPosInLoc(i)
            memPos_t  = 0!commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = numNodesToComm !commsMemSize(i)

            call MPI_Accumulate(aux_floatField_s(memPos_l),memSize,mpi_datatype_real,&
                                iRank,memPos_t,memSize,mpi_datatype_real,MPI_SUM,window_id_float,mpi_err)
        end do


        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_float,mpi_err)

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            floatField(iNodeL) = aux_floatField_r(i) + aux_floatField_s(i)
        end do


    end subroutine test_multiupdate_and_comm_floatField

    subroutine test_new_feature()
        integer :: newRank, newSize, mpi_newcomm
        integer :: window_id,iRank,disp
        type(c_ptr) :: c_window_ptr
        integer, pointer :: f_window_ptr(:)
        integer(KIND=MPI_ADDRESS_KIND) :: ssize
        !call MPI_COMM_TYPE_SHARED(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,mpi_key,MPI_INFO_NULL,mpi_newcomm,mpi_err)
        !call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,mpi_newcomm,mpi_err)
        
        call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,mpi_newcomm,mpi_err)
        
        call mpi_comm_rank(mpi_newcomm, newRank, mpi_err)
        call mpi_comm_size(mpi_newcomm, newSize, mpi_err)

        if (mpi_newcomm .ne. MPI_COMM_NULL) then
            write(*,*) 'NEW COMM NOT NULL! mpi_newcomm ',mpi_newcomm,' old_comm ',MPI_COMM_WORLD
            write(*,*) 'old rank', mpi_rank, ' new_rank ', newRank,' new_size ',newSize
        end if

        window_buffer_size = 1*mpi_integer_size
        call MPI_Win_allocate_shared(window_buffer_size,mpi_integer_size,MPI_INFO_NULL,&
                                     mpi_newcomm,c_window_ptr,window_id,mpi_err)

        call c_f_pointer(c_window_ptr, f_window_ptr, SHAPE = [mpi_size])

        f_window_ptr(0) = mpi_rank+1
        write(*,*) 'pre val[',mpi_rank,'] ',f_window_ptr(0)
        if(mpi_rank.eq.0) f_window_ptr(1) = 3
        if(mpi_rank.eq.2) f_window_ptr(-2) = 5

        write(*,*) 'post val[',mpi_rank,'] ',f_window_ptr(0)


        iRank=0
        if(mpi_rank.ne.iRank) then
            call MPI_Win_shared_query(window_id,iRank,ssize,disp,c_window_ptr,mpi_err)
            !write(*,*) 't[',mpi_rank,'] ssize',ssize,' disp ',disp
        end if
        call c_f_pointer(c_window_ptr, f_window_ptr, SHAPE = [mpi_size])

        write(*,*) 'post t',mpi_rank,'] ',f_window_ptr(0)


        call MPI_Win_free(window_id, mpi_err)

    end subroutine test_new_feature
#endif

end module mod_comms
