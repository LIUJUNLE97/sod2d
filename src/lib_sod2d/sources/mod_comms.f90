module mod_comms
    use mod_mpi_mesh
    use mod_custom_types
    use mod_nvtx
#ifdef NCCL_COMMS
    use cudafor
    use nccl
#endif

!-- Select type of communication for mpi_atomic_updates
#define _SENDRCV_ 0
#define _ISENDIRCV_ 1
#define _PUTFENCE_ 0
#define _NOCUDAAWARE_ 0

#ifdef AUTOCOMP
 #define _NOCUDAAWARE_ 1
 #define _SENDRCV_ 0
 #define _ISENDIRCV_ 0
 #define _PUTFENCE_ 0
#endif


    implicit none

    !---- for the comms---
    integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
    integer(KIND=MPI_ADDRESS_KIND) :: memPos_t
    integer(4) :: worldGroup,commGroup

    integer(4),dimension(:),allocatable :: aux_intField_s,  aux_intField_r
    real(rp),dimension(:),allocatable   :: aux_realField_s, aux_realField_r

    integer(4) :: maxIntBufferArraySize=1,maxRealBufferArraySize=1
    integer(4) :: window_id_int,window_id_real
    integer(4) :: window_id_sm
    integer(4) :: beginFence=0,endFence=0
    integer(4) :: startAssert=0,postAssert=0

    logical :: isInt=.false.,isReal=.false.,isIntWindow=.false.,isRealWindow=.false.
    logical :: isLockBarrier,isPSCWBarrier

#ifdef NCCL_COMMS
    integer                        :: cuda_stat
    type(ncclResult)               :: nccl_stat
    type(ncclUniqueId)             :: nccl_uid
    type(ncclComm)                 :: nccl_comm
    integer(kind=cuda_stream_kind) :: nccl_stream
#endif

contains

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
    subroutine init_comms(useInt,useReal,intBufferMult,realBufferMult,initWindowsArg)
        implicit none
        logical,intent(in) :: useInt,useReal
        integer(4),intent(in),optional :: intBufferMult,realBufferMult
        logical,intent(in),optional :: initWindowsArg
        logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier,initWindows=.false.

        if(present(initWindowsArg)) then
           initWindows = initWindowsArg
        end if

        maxIntBufferArraySize  = 1
        maxRealBufferArraySize = 1

        if(present(intBufferMult)) then
            maxIntBufferArraySize = intBufferMult
        end if

        if(present(realBufferMult)) then
            maxRealBufferArraySize = realBufferMult
        end if

#if _SENDRCV_
        if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: Send-Recv"
#endif
#if _ISENDIRCV_
        if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: iSend-iRecv"
#endif
#if _PUTFENCE_
        if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: Put(Fence)"
        initWindows = .true. !forcing initWindows to true
#endif
#if _NOCUDAAWARE_
        if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: NO CUDA-Aware iSend-iRecv"
#endif
#ifdef NCCL_COMMS
        if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: NCCL"
#endif

        isInt        =.false.
        isReal       =.false.
        isIntWindow  =.false.
        isRealWindow =.false.

        if(useInt) then
            isInt = .true.

            allocate(aux_intField_s(maxIntBufferArraySize*numNodesToComm))
            allocate(aux_intField_r(maxIntBufferArraySize*numNodesToComm))
            !$acc enter data create(aux_intField_s(:))
            !$acc enter data create(aux_intField_r(:))

            if(initWindows) call init_window_intField()
        end if

        if(useReal) then
            isReal = .true.

            allocate(aux_realField_s(maxRealBufferArraySize*numNodesToComm))
            allocate(aux_realField_r(maxRealBufferArraySize*numNodesToComm))
            !$acc enter data create(aux_realField_s(:))
            !$acc enter data create(aux_realField_r(:))

            if(initWindows) call init_window_realField()
        end if

        call MPI_Comm_group(app_comm,worldGroup,mpi_err)
	    call MPI_Group_incl(worldGroup,numRanksWithComms,ranksToComm,commGroup,mpi_err);

        useFenceFlags=.false. !by default
        useAssertNoCheckFlags=.true. !by default
        isPSCWBarrier=.true.!si faig molts loops amb aquesta opció a false la comm queda bloquejada
        useLockBarrier=.true.!.false. with false it fails!
        call setFenceFlags(useFenceFlags)
        call setPSCWAssertNoCheckFlags(useAssertNoCheckFlags)
        call setLockBarrier(useLockBarrier)

#ifdef NCCL_COMMS
        if (mpi_rank == 0) then
            nccl_stat = ncclGetUniqueId(nccl_uid)
        end if
        call MPI_Bcast(nccl_uid, int( sizeof(ncclUniqueId), kind = 4 ), MPI_BYTE, 0, app_comm, mpi_err)
        nccl_stat = ncclCommInitRank(nccl_comm, mpi_size, nccl_uid, mpi_rank);
        cuda_stat = cudaStreamCreate(nccl_stream);
#endif

    end subroutine init_comms

    subroutine end_comms()
        implicit none

        if(isInt) then
            !$acc exit data delete(aux_intField_s(:))
            !$acc exit data delete(aux_intField_r(:))
            deallocate(aux_intField_s)
            deallocate(aux_intField_r)

            if(isIntWindow) call close_window_intField()
        end if

        if(isReal) then
           !$acc exit data delete(aux_realField_s(:))
           !$acc exit data delete(aux_realField_r(:))
            deallocate(aux_realField_s)
            deallocate(aux_realField_r)

            if(isRealWindow) call close_window_realField()
        end if

#ifdef NCCL_COMMS
        nccl_stat = ncclCommDestroy(nccl_comm)
        cuda_stat = cudaStreamDestroy(nccl_stream)
#endif
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

        isIntWindow=.true.

        window_buffer_size = mpi_integer_size*numNodesToComm
        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        call MPI_Win_create(aux_intField_r,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_int,mpi_err)
        !$acc end host_data
    end subroutine init_window_intField

    subroutine close_window_intField()
        implicit none

        call MPI_Win_free(window_id_int,mpi_err)
    end subroutine close_window_intField
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_realField()
        implicit none

        isRealWindow=.true.

        window_buffer_size = mpi_real_size*numNodesToComm
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        call MPI_Win_create(aux_realField_r,window_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_real,mpi_err)
        !$acc end host_data
    end subroutine init_window_realField

    subroutine close_window_realField()
        implicit none

        call MPI_Win_free(window_id_real,mpi_err)
    end subroutine close_window_realField

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!    SUBROUTINES COPY/FROM SEND/RCV BUFFERS
!-----------------------------------------------------------------------------------------------------------------------
    subroutine fill_sendBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop async(1)
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_intField_s(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_intField_r(:)=0
        !$acc end kernels

        !$acc wait
    end subroutine fill_sendBuffer_int
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_real(realField)
        implicit none
        real(rp),intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop async(1)
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_realField_s(i) = realField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_realField_r(:)=0.0_rp
        !$acc end kernels

        !$acc wait
    end subroutine fill_sendBuffer_real

    subroutine fill_sendBuffer_arrays_real(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp),intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,iArray,iNodeL

        call nvtxStartRange("fillBuffer_arrays")
        !$acc parallel loop async(1)
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            do iArray = 1,numArrays
                aux_realField_s((i-1)*numArrays+iArray) = arrays2comm(iNodeL,iArray)
            end do
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_realField_r(:)=0.0_rp
        !$acc end kernels

        !$acc wait
        call nvtxEndRange
    end subroutine fill_sendBuffer_arrays_real

    subroutine fill_sendBuffer_mass_ener_momentum_real(mass,ener,momentum)
        implicit none
        real(rp),intent(in) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,idime,iNodeL
        integer(4),parameter :: nArrays=5

        call nvtxStartRange("fillBuffer_async")
        !$acc parallel loop async(1)
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            do idime = 1,ndime
                aux_realField_s((i-1)*nArrays+idime) = momentum(iNodeL,idime)
            end do
            aux_realField_s((i-1)*nArrays+4) = mass(iNodeL)
            aux_realField_s((i-1)*nArrays+5) = ener(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_realField_r(:)=0.0_rp
        !$acc end kernels

        !$acc wait
        call nvtxEndRange
    end subroutine fill_sendBuffer_mass_ener_momentum_real

    subroutine fill_sendBuffer_arraysPtr_real(numArrays,arraysPtr2comm)
        implicit none
        integer(4), intent(in) :: numArrays
        type(ptr_array1d_rp),intent(inout) :: arraysPtr2comm(numArrays)
        integer(4) :: iNode,iArray,iNodeL

        call nvtxStartRange("fillBuffer_arrays")
        !$acc parallel loop present(arraysPtr2comm(:)) async(1)
        do iNode=1,numNodesToComm
            iNodeL = nodesToComm(iNode)
            do iArray = 1,numArrays
                aux_realField_s((iNode-1)*numArrays+iArray) = arraysPtr2comm(iArray)%ptr(iNodeL) 
            end do
        end do
        !$acc end parallel loop

        !$acc kernels async(2)
        aux_realField_r(:)=0.0_rp
        !$acc end kernels

        !$acc wait
        call nvtxEndRange
    end subroutine fill_sendBuffer_arraysPtr_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_get_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop async(1)
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_intField_r(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_intField_s(:)=0
        !$acc end kernels

        !$acc wait
    end subroutine fill_sendBuffer_get_int
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_get_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop async(1)
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            aux_realField_r(i) = realField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_realField_s(:)=0.0_rp
        !$acc end kernels

        !$acc wait
    end subroutine fill_sendBuffer_get_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

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
    subroutine copy_from_min_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            intField(iNodeL) = min(intField(iNodeL), aux_intField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_min_rcvBuffer_int
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = realField(iNodeL) + aux_realField_r(i)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_real

    subroutine copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp),intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,iArray,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            do iArray = 1,numArrays
                !$acc atomic update
                arrays2comm(iNodeL,iArray) = arrays2comm(iNodeL,iArray) + aux_realField_r((i-1)*numArrays+iArray)
                !$acc end atomic
            end do
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_arrays_real

    subroutine copy_from_rcvBuffer_mass_ener_momentum_real(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,idime,iNodeL
        integer(4),parameter :: nArrays=5

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            do idime = 1,ndime
                !$acc atomic update
                momentum(iNodeL,idime) = momentum(iNodeL,idime) + aux_realField_r((i-1)*nArrays+idime)
                !$acc end atomic
            end do
            !$acc atomic update
            mass(iNodeL) = mass(iNodeL) + aux_realField_r((i-1)*nArrays+4)
            !$acc end atomic
            !$acc atomic update
            ener(iNodeL) = ener(iNodeL) + aux_realField_r((i-1)*nArrays+5)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_mass_ener_momentum_real

    subroutine copy_from_rcvBuffer_arraysPtr_real(numArrays,arraysPtr2comm)
        implicit none
        integer(4), intent(in) :: numArrays
        type(ptr_array1d_rp),intent(inout) :: arraysPtr2comm(numArrays)
        integer(4) :: iNode,iArray,iNodeL

        !$acc parallel loop present(arraysPtr2comm(:))
        do iNode=1,numNodesToComm
            iNodeL = nodesToComm(iNode)
            do iArray = 1,numArrays
                !$acc atomic update
                arraysPtr2comm(iArray)%ptr(iNodeL) = arraysPtr2comm(iArray)%ptr(iNodeL) + aux_realField_r((iNode-1)*numArrays+iArray)
                !$acc end atomic
            end do
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_arraysPtr_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine copy_from_min_rcvBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = min(realField(iNodeL) , aux_realField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_min_rcvBuffer_real
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_get_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

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
        integer(4) :: i,iNodeL

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
        integer(4) :: i,iNodeL

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
!-------------------------------------------------------------------------------------
    subroutine mpi_halo_atomic_update_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
#if _SENDRCV_
        call mpi_halo_atomic_update_int_sendRcv(intField)
#endif
#if _ISENDIRCV_
        call mpi_halo_atomic_update_int_iSendiRcv(intField)
#endif
#if _PUTFENCE_
        call mpi_halo_atomic_update_int_put_fence(intField)
#endif
#if _NOCUDAAWARE_
        call mpi_halo_atomic_update_int_iSendiRcv_noCudaAware(intField)
#endif
    end subroutine mpi_halo_atomic_update_int
!-------------------------------------------------------------------------------------
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
#if _NOCUDAAWARE_
        call mpi_halo_atomic_update_real_iSendiRcv_noCudaAware(realField)
#endif
    end subroutine mpi_halo_atomic_update_real
!-------------------------------------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_mass_ener_momentum(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)

#if _SENDRCV_
        call mpi_halo_atomic_update_real_mass_ener_momentum_SendRcv(mass,ener,momentum)
#endif
#if _ISENDIRCV_
        call mpi_halo_atomic_update_real_mass_ener_momentum_iSendiRcv(mass,ener,momentum)
#endif
#if _PUTFENCE_
        write(*,*) 'Crashing. Not Implemented for PUTFENCE option'
        call MPI_Abort(app_comm,-1,mpi_err)
#endif
#if _NOCUDAAWARE_
        call mpi_halo_atomic_update_real_mass_ener_mom_iSendiRcv_noCudaAware(mass,ener,momentum)
#endif

    end subroutine mpi_halo_atomic_update_real_mass_ener_momentum
!-------------------------------------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_arrays(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp), intent(inout) :: arrays2comm(:,:)

#if _SENDRCV_
        call mpi_halo_atomic_update_real_arrays_SendRcv(numArrays,arrays2comm)
#endif
#if _ISENDIRCV_
        call mpi_halo_atomic_update_real_arrays_iSendiRcv(numArrays,arrays2comm)
#endif
#if _PUTFENCE_
        write(*,*) 'Crashing. Not Implemented for PUTFENCE option'
        call MPI_Abort(app_comm,-1,mpi_err)
#endif
#if _NOCUDAAWARE_
        call mpi_halo_atomic_update_real_arrays_iSendiRcv_noCudaAware(numArrays,arrays2comm)
#endif

    end subroutine mpi_halo_atomic_update_real_arrays

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!------------- SEND/RECV -------------------------------------------
    ! INTEGER ---------------------------------------------------
    subroutine mpi_halo_atomic_update_int_sendRcv(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_intField_s(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              aux_intField_r(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do
        !$acc end host_data

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_sendRcv
    ! REAL ---------------------------------------------------
    subroutine mpi_halo_atomic_update_real_sendRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do
        !$acc end host_data

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_sendRcv

<<<<<<< HEAD

 !INTEGER ---------------------------------------------------
    subroutine mpi_halo_atomic_min_update_int_sendRcv(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_intField_s(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              aux_intField_r(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do
        !$acc end host_data

        call copy_from_min_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_min_update_int_sendRcv


     ! REAL ---------------------------------------------------
=======
    ! REAL MIN ---------------------------------------------------
>>>>>>> origin/master
    subroutine mpi_halo_atomic_min_update_real_sendRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do
        !$acc end host_data

        call copy_from_min_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_min_update_real_sendRcv

    ! MASS,ENER,MOMENTUM ---------------------------------------------------
    subroutine mpi_halo_atomic_update_real_mass_ener_momentum_sendRcv(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4),parameter :: nArrays=5

        call fill_sendBuffer_mass_ener_momentum_real(mass,ener,momentum)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*nArrays+1
            memSize  = commsMemSize(i)*nArrays

            call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do
        !$acc end host_data

        call copy_from_rcvBuffer_mass_ener_momentum_real(mass,ener,momentum)
    end subroutine mpi_halo_atomic_update_real_mass_ener_momentum_sendRcv

    ! REAL N-ARRAYS ---------------------------------------------------
    subroutine mpi_halo_atomic_update_real_arrays_sendRcv(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp), intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
            memSize  = commsMemSize(i)*numArrays

            call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do
        !$acc end host_data

        call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)
    end subroutine mpi_halo_atomic_update_real_arrays_sendRcv

!------------- ISEND/IRECV -------------------------------------------
    !INTEGER ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_iSendiRcv(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_int(intField)

        ireq=0
        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_iSendiRcv
    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_iSendiRcv(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_real(realField)

#if NCCL_COMMS
        nccl_stat = ncclGroupStart()
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            nccl_stat = ncclRecv(aux_realField_r(mempos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
            nccl_stat = ncclSend(aux_realField_s(mempos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
        end do
        !$acc end host_data
        nccl_stat = ncclGroupEnd()
        cuda_stat = cudaStreamSynchronize(nccl_stream)
#else
        ireq=0
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
#endif

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_iSendiRcv

    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_mass_ener_momentum_iSendiRcv(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)
        integer(4),parameter :: nArrays=5

        call fill_sendBuffer_mass_ener_momentum_real(mass,ener,momentum)

#if NCCL_COMMS
        nccl_stat = ncclGroupStart()
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = (commsMemPosInLoc(i)-1)*nArrays+1
            memSize  = commsMemSize(i)*nArrays

            nccl_stat = ncclRecv(aux_realField_r(memPos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
            nccl_stat = ncclSend(aux_realField_s(memPos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
        end do
        !$acc end host_data
        nccl_stat = ncclGroupEnd()
        cuda_stat = cudaStreamSynchronize(nccl_stream)
#else
        ireq=0
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*nArrays+1
            memSize  = commsMemSize(i)*nArrays

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
#endif

        call copy_from_rcvBuffer_mass_ener_momentum_real(mass,ener,momentum)
    end subroutine mpi_halo_atomic_update_real_mass_ener_momentum_iSendiRcv


    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp), intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

#if NCCL_COMMS
        nccl_stat = ncclGroupStart()
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
            memSize  = commsMemSize(i)*numArrays

            nccl_stat = ncclRecv(aux_realField_r(memPos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
            nccl_stat = ncclSend(aux_realField_s(memPos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
        end do
        !$acc end host_data
        nccl_stat = ncclGroupEnd()
        cuda_stat = cudaStreamSynchronize(nccl_stream)
#else
        ireq=0
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
            memSize  = commsMemSize(i)*numArrays

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
#endif

        call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)
    end subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv

    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_arraysPtr_iSendiRcv(numArrays,arraysPtr2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        type(ptr_array1d_rp),intent(inout) :: arraysPtr2comm(numArrays)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_arraysPtr_real(numArrays,arraysPtr2comm)

#if NCCL_COMMS
        nccl_stat = ncclGroupStart()
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
            memSize  = commsMemSize(i)*numArrays

            nccl_stat = ncclRecv(aux_realField_r(memPos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
            nccl_stat = ncclSend(aux_realField_s(memPos_l), memSize, ncclFloat, ngbRank, nccl_comm, nccl_stream)
        end do
        !$acc end host_data
        nccl_stat = ncclGroupEnd()
        cuda_stat = cudaStreamSynchronize(nccl_stream)
#else
        ireq=0
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
            memSize  = commsMemSize(i)*numArrays

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
#endif

        call copy_from_rcvBuffer_arraysPtr_real(numArrays,arraysPtr2comm)

    end subroutine mpi_halo_atomic_update_real_arraysPtr_iSendiRcv


!------------- conditional average ISEND/IRECV -------------------------------------------
    ! REAL ---------------------------------------------------
    subroutine mpi_halo_conditional_ave_update_real_iSendiRcv(cond,realField)
        implicit none
        real(rp),intent(in) :: cond
        real(rp),intent(inout) :: realField(:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_real(realField)

        ireq=0
        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_realfield_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realfield_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do
        !$acc end host_data

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        call copy_from_conditional_ave_rcvBuffer_real(cond,realField)
    end subroutine mpi_halo_conditional_ave_update_real_iSendiRcv
!------------- PUT FENCE -------------------------------------------
    !INT
    subroutine mpi_halo_atomic_update_int_put_fence(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        call MPI_Win_fence(beginFence,window_id_int,mpi_err)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do
        !$acc end host_data

        call MPI_Win_fence(endFence,window_id_int,mpi_err)

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_put_fence
    !REAL
    subroutine mpi_halo_atomic_update_real_put_fence(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        call MPI_Win_fence(beginFence,window_id_real,mpi_err)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do
        !$acc end host_data

        call MPI_Win_fence(endFence,window_id_real,mpi_err)

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_fence
!------------- PUT PSCW -------------------------------------------
    !INTEGER--------------------------------------------------
    subroutine mpi_halo_atomic_update_int_put_pscw(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err);

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do
        !$acc end host_data

	    call MPI_Win_complete(window_id_int,mpi_err);
	    call MPI_Win_wait(window_id_int,mpi_err);

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_put_pscw
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_put_pscw(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err);

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do
        !$acc end host_data

	    call MPI_Win_complete(window_id_real,mpi_err);
	    call MPI_Win_wait(window_id_real,mpi_err);

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_pscw

!------------- PUT LOCK -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_put_lock(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_int(intField)

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
        end do
        !$acc end host_data

        if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_put_lock
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_put_lock(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
        end do
        !$acc end host_data

        if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_lock

!------------- GET FENCE -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_get_fence(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_int(intField)

        call MPI_Win_fence(beginFence,window_id_int,mpi_err)

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do
        !$acc end host_data

        call MPI_Win_fence(endFence,window_id_int,mpi_err)

        call copy_from_rcvBuffer_get_int(intField)
    end subroutine mpi_halo_atomic_update_int_get_fence
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_get_fence(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        call MPI_Win_fence(beginFence,window_id_real,mpi_err)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do
        !$acc end host_data

        call MPI_Win_fence(endFence,window_id_real,mpi_err)

        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_fence

!------------- GET PSCW -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_get_pscw(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_int(intField)

        if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err);

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank   = ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize   = commsMemSize(i)

            call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
        end do
        !$acc end host_data

	    call MPI_Win_complete(window_id_int,mpi_err);
	    call MPI_Win_wait(window_id_int,mpi_err);

        call copy_from_rcvBuffer_get_int(intField)
    end subroutine mpi_halo_atomic_update_int_get_pscw
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_get_pscw(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
	    call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err);
	    call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err);

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank   = ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize   = commsMemSize(i)

            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do
        !$acc end host_data

	    call MPI_Win_complete(window_id_real,mpi_err);
	    call MPI_Win_wait(window_id_real,mpi_err);

        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_pscw

!------------- GET LOCK -------------------------------------------
    !INTEGER-------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_get_lock(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_int(intField)

        !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Get(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
        end do
        !$acc end host_data

        if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

        call copy_from_rcvBuffer_get_int(intField)
    end subroutine mpi_halo_atomic_update_int_get_lock
    !REAL-------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_get_lock(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
            call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
        end do
        !$acc end host_data

        if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_lock

!---------------------------------------------------------------------------------
!    FOR TESTING AND DEVEL STUFF!
!---------------------------------------------------------------------------------

!------------- ONLY BUFFERS -------------------------------------------
    !INTEGER ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_int_onlybuffers(intfield)
        implicit none
        integer(4), intent(inout) :: intfield(:)
        integer(4) :: i,ngbrank,tagcomm
        integer(4) :: mempos_l,memsize

        call fill_sendbuffer_int(intfield)
        call copy_from_rcvbuffer_int(intfield)

    end subroutine mpi_halo_atomic_update_int_onlybuffers
    !REAL ---------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_onlybuffers(realfield)
        implicit none
        real(rp), intent(inout) :: realfield(:)
        integer(4) :: i,ngbrank,tagcomm
        integer(4) :: mempos_l,memsize

        call fill_sendbuffer_real(realfield)
        call copy_from_rcvbuffer_real(realfield)

    end subroutine mpi_halo_atomic_update_real_onlybuffers
!---------------------------------------------------------------------------------
    subroutine mpi_halo_atomic_update_real_iSendiRcv_noCudaAware(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_real(realField)

        ireq=0

        !$acc update host(aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        !$acc update device(aux_realField_r(:))
        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_iSendiRcv_noCudaAware

    subroutine mpi_halo_atomic_update_real_sendRcv_noCudaAware(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank,tagComm
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        !$acc update host(aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                              app_comm, MPI_STATUS_IGNORE, mpi_err)
        end do

        !$acc update device(aux_realField_r(:))
        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_sendRcv_noCudaAware

    !REAL
    subroutine mpi_halo_atomic_update_real_put_fence_noCudaAware(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_real(realField)

        call MPI_Win_fence(beginFence,window_id_real,mpi_err)

        !$acc update host(aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do

        call MPI_Win_fence(endFence,window_id_real,mpi_err)

        !$acc update device(aux_realField_r(:))
        call copy_from_rcvBuffer_real(realField)
    end subroutine mpi_halo_atomic_update_real_put_fence_noCudaAware
 
    subroutine mpi_halo_atomic_update_real_get_fence_noCudaAware(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,ngbRank
        integer(4) :: memPos_l,memSize

        call fill_sendBuffer_get_real(realField)

        call MPI_Win_fence(beginFence,window_id_real,mpi_err)

        !$acc update host(aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            memPos_l = commsMemPosInLoc(i)
            memPos_t = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement
            memSize  = commsMemSize(i)

            call MPI_Get(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
        end do

        call MPI_Win_fence(endFence,window_id_real,mpi_err)

        !$acc update device(aux_realField_s(:))
        call copy_from_rcvBuffer_get_real(realField)
    end subroutine mpi_halo_atomic_update_real_get_fence_noCudaAware

    subroutine mpi_halo_atomic_update_real_mass_ener_mom_iSendiRcv_noCudaAware(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)
        integer(4),parameter :: nArrays=5

        call fill_sendBuffer_mass_ener_momentum_real(mass,ener,momentum)

        ireq=0
        !$acc update host(aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*nArrays+1
            memSize  = commsMemSize(i)*nArrays

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        !$acc update device(aux_realField_r(:))
        call copy_from_rcvBuffer_mass_ener_momentum_real(mass,ener,momentum)
    end subroutine mpi_halo_atomic_update_real_mass_ener_mom_iSendiRcv_noCudaAware

    subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv_noCudaAware(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp), intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

        ireq=0
        !$acc update host(aux_realField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
            memSize  = commsMemSize(i)*numArrays

            ireq = ireq+1
            call MPI_Irecv(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

        !$acc update device(aux_realField_r(:))
        call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)
    end subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv_noCudaAware

    subroutine mpi_halo_atomic_update_int_iSendiRcv_noCudaAware(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,ireq,ngbRank,tagComm
        integer(4) :: memPos_l,memSize
        integer(4) :: requests(2*numRanksWithComms)

        call fill_sendBuffer_int(intField)

        ireq=0
        !$acc update host(aux_intField_s(:))
        do i=1,numRanksWithComms
            ngbRank  = ranksToComm(i)
            tagComm  = 0
            memPos_l = commsMemPosInLoc(i)
            memSize  = commsMemSize(i)

            ireq = ireq+1
            call MPI_Irecv(aux_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
            ireq = ireq+1
            call MPI_ISend(aux_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
        end do

        call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
        !$acc update device(aux_intField_r(:))
        call copy_from_rcvBuffer_int(intField)
    end subroutine mpi_halo_atomic_update_int_iSendiRcv_noCudaAware

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

end module mod_comms
