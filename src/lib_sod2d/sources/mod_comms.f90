module mod_comms
   use mod_mpi_mesh
   use mod_custom_types
   use mod_nvtx
#ifdef NCCL_SUPPORT
   use cudafor
   use nccl
#endif

!-- Select type of communication for mpi_atomic_updates
#define _SENDRCV_ 0
#define _ISENDIRCV_ 1
#define _PUTFENCE_ 0
#define _GETFENCE_ 0
#define _NOCUDAAWARE_ 0
#define _NCCL_ 0

#ifdef AUTOCOMP
 #define _NOCUDAAWARE_ 1
 #define _SENDRCV_ 0
 #define _ISENDIRCV_ 0
 #define _PUTFENCE_ 0
 #define _GETFENCE_ 0
 #define _NCCL_ 0
#endif

!-----------------------------------

   implicit none

   !---- for the comms---
   integer(KIND=MPI_ADDRESS_KIND) :: memPos_t
   integer(4) :: worldGroup,commGroup

   integer(4),dimension(:),allocatable :: aux_intField_s,  aux_intField_r
   real(rp),dimension(:),allocatable   :: aux_realField_s, aux_realField_r

   integer(4) :: window_id_int,window_id_real
   integer(4) :: window_id_sm
   integer(4) :: beginFence=0,endFence=0
   integer(4) :: startAssert=0,postAssert=0

   logical :: isInt=.false.,isReal=.false.,isIntWindow=.false.,isRealWindow=.false.
   logical :: isLockBarrier,isPSCWBarrier

#ifdef NCCL_SUPPORT
   integer                        :: cuda_stat
   type(ncclResult)               :: nccl_stat
   type(ncclUniqueId)             :: nccl_uid
   type(ncclComm)                 :: nccl_comm
   integer(kind=cuda_stream_kind) :: nccl_stream
#endif

contains

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
   subroutine init_comms(useInt,useReal,intBufferMultIn,realBufferMultIn,initWinModeIn)
      implicit none
      logical,intent(in) :: useInt,useReal
      integer(4),intent(in),optional :: intBufferMultIn,realBufferMultIn,initWinModeIn
      integer(4) :: maxIntBufferArraySize=1,maxRealBufferArraySize=1,winMode
      logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier,initWindows=.false.

      initWindows = .false.
      winMode = 0

      if(present(initWinModeIn)) then
         initWindows = .true.
         winMode = initWinModeIn

         if((winMode.ne.1).and.(winMode.ne.2)) then
            write(*,*) 'Crashing in init_comms! Bad winMode:',winMode,' (must be 1 (PUT) or 2 (GET))'
            call MPI_Abort(app_comm,-1,mpi_err)
         end if
      end if

      maxIntBufferArraySize  = 1
      maxRealBufferArraySize = 1

      if(present(intBufferMultIn)) then
         maxIntBufferArraySize = intBufferMultIn
      end if

      if(present(realBufferMultIn)) then
         maxRealBufferArraySize = realBufferMultIn
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
      winMode = 1 !forcing winMode 1 (PUT)
#endif
#if _GETFENCE_
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: Get(Fence)"
      initWindows = .true. !forcing initWindows to true
      winMode = 2 !forcing winMode 2 (GET)
#endif
#if _NCCL_
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: NCCL"
#endif
#if _NOCUDAAWARE_
      if(mpi_rank.eq.0) write(111,*) "--| Comm. scheme: NO CUDA-Aware iSend-iRecv"
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

         if(initWindows) call init_window_intField(winMode,maxIntBufferArraySize)
      end if

      if(useReal) then
         isReal = .true.

         allocate(aux_realField_s(maxRealBufferArraySize*numNodesToComm))
         allocate(aux_realField_r(maxRealBufferArraySize*numNodesToComm))
         !$acc enter data create(aux_realField_s(:))
         !$acc enter data create(aux_realField_r(:))

         if(initWindows) call init_window_realField(winMode,maxRealBufferArraySize)
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

#ifdef NCCL_SUPPORT
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

#ifdef NCCL_SUPPORT
      nccl_stat = ncclCommDestroy(nccl_comm)
      cuda_stat = cudaStreamDestroy(nccl_stream)
#endif

      !setting false all general flags
      isInt=.false.
      isReal=.false.
      isIntWindow=.false.
      isRealWindow=.false.
      isLockBarrier=.false.
      isPSCWBarrier=.false.

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
    subroutine init_window_intField(windowMode,maxIntBufferArraySize)
        implicit none
        integer(4),intent(in) :: windowMode,maxIntBufferArraySize
        integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
        !------------------------------------------

        window_buffer_size = mpi_integer_size*(maxIntBufferArraySize*numNodesToComm)

        if(windowMode.eq.1) then !For PUT Comms, the buffer exposed is aux_intField_r
            !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
            call MPI_Win_create(aux_intField_r,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_int,mpi_err)
            !$acc end host_data
        else if(windowMode.eq.2) then !For GET Comms, the buffer exposed is aux_intField_s
            !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
            call MPI_Win_create(aux_intField_s,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_int,mpi_err)
            !$acc end host_data
        else
            write(*,*) 'Crashing. Bad windowsMode in init_window_intField! Must be 1 (PUT) or 2 (GET)'
            call MPI_Abort(app_comm,-1,mpi_err)
        end if

        isIntWindow=.true.

    end subroutine init_window_intField

    subroutine close_window_intField()
        implicit none

        call MPI_Win_free(window_id_int,mpi_err)
    end subroutine close_window_intField
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine init_window_realField(windowMode,maxRealBufferArraySize)
        implicit none
        integer(4),intent(in) :: windowMode,maxRealBufferArraySize
        integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
        !------------------------------------------

        window_buffer_size = mpi_real_size*(maxRealBufferArraySize*numNodesToComm)

        if(windowMode.eq.1) then !For PUT Comms, the buffer exposed is aux_intField_r
            !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
            call MPI_Win_create(aux_realField_r,window_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_real,mpi_err)
            !$acc end host_data
        else if(windowMode .eq. 2) then !For GET Comms, the buffer exposed is aux_intField_s
            !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
            call MPI_Win_create(aux_realField_s,window_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_real,mpi_err)
            !$acc end host_data
        else
            write(*,*) 'Crashing. Bad windowsMode in init_window_realField! Must be 1 (PUT) or 2 (GET)'
            call MPI_Abort(app_comm,-1,mpi_err)
        end if

        isRealWindow=.true.

    end subroutine init_window_realField

    subroutine close_window_realField()
        implicit none

        call MPI_Win_free(window_id_real,mpi_err)
    end subroutine close_window_realField

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!-------- Fill Buffer Subroutines ----------------
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
!-------------------------------------------------------------------------
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
!-------------------------------------------------------------------------
    subroutine fill_sendBuffer_massEnerMom_real(mass,ener,momentum)
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
    end subroutine fill_sendBuffer_massEnerMom_real
!-------------------------------------------------------------------------
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
!-------- Copy Buffer Subroutines ----------------
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
!-------------------------------------------------------------------------
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
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)
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
    end subroutine copy_from_rcvBuffer_massEnerMom_real
!-------------------------------------------------------------------------
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
!-------------------------------------------------------------------------
    subroutine copy_from_min_rcvBuffer_get_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = min(realField(iNodeL), aux_realField_s(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_min_rcvBuffer_get_real
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_arrays_get_real(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp),intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,iArray,iNodeL

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            do iArray = 1,numArrays
                !$acc atomic update
                arrays2comm(iNodeL,iArray) = arrays2comm(iNodeL,iArray) + aux_realField_s((i-1)*numArrays+iArray)
                !$acc end atomic
            end do
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_arrays_get_real
!-------------------------------------------------------------------------
    subroutine copy_from_rcvBuffer_massEnerMom_get_real(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,idime,iNodeL
        integer(4),parameter :: nArrays=5

        !$acc parallel loop
        do i=1,numNodesToComm
            iNodeL = nodesToComm(i)
            do idime = 1,ndime
                !$acc atomic update
                momentum(iNodeL,idime) = momentum(iNodeL,idime) + aux_realField_s((i-1)*nArrays+idime)
                !$acc end atomic
            end do
            !$acc atomic update
            mass(iNodeL) = mass(iNodeL) + aux_realField_s((i-1)*nArrays+4)
            !$acc end atomic
            !$acc atomic update
            ener(iNodeL) = ener(iNodeL) + aux_realField_s((i-1)*nArrays+5)
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_rcvBuffer_massEnerMom_get_real

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
#if _NCCL_
      call mpi_halo_atomic_update_int_ncclSendRcv(intField)
#endif
#if _PUTFENCE_
      call mpi_halo_atomic_update_int_put_fence(intField)
#endif
#if _GETFENCE_
      call mpi_halo_atomic_update_int_get_fence(intField)
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
#if _NCCL_
      call mpi_halo_atomic_update_real_ncclSendRcv(realField)
#endif
#if _PUTFENCE_
      call mpi_halo_atomic_update_real_put_fence(realField)
#endif
#if _GETFENCE_
      call mpi_halo_atomic_update_real_get_fence(intField)
#endif
#if _NOCUDAAWARE_
      call mpi_halo_atomic_update_real_iSendiRcv_noCudaAware(realField)
#endif
   end subroutine mpi_halo_atomic_update_real
!-------------------------------------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)

#if _SENDRCV_
      call mpi_halo_atomic_update_real_massEnerMom_SendRcv(mass,ener,momentum)
#endif
#if _ISENDIRCV_
      call mpi_halo_atomic_update_real_massEnerMom_iSendiRcv(mass,ener,momentum)
#endif
#if _NCCL_
      call mpi_halo_atomic_update_real_massEnerMom_ncclSendRcv(mass,ener,momentum)
#endif
#if _PUTFENCE_
      call mpi_halo_atomic_update_real_massEnerMom_put_fence(mass,ener,momentum)
#endif
#if _GITFENCE_
      call mpi_halo_atomic_update_real_massEnerMom_get_fence(mass,ener,momentum)
#endif
#if _NOCUDAAWARE_
      call mpi_halo_atomic_update_real_massEnerMom_iSendiRcv_noCudaAware(mass,ener,momentum)
#endif
    end subroutine mpi_halo_atomic_update_real_massEnerMom
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
#if _NCCL_
      call mpi_halo_atomic_update_real_arrays_ncclSendRcv(numArrays,arrays2comm)
#endif
#if _PUTFENCE_
      call mpi_halo_atomic_update_real_arrays_put_fence(numArrays,arrays2comm)
#endif
#if _GETFENCE_
      call mpi_halo_atomic_update_real_arrays_get_fence(numArrays,arrays2comm)
#endif
#if _NOCUDAAWARE_
      call mpi_halo_atomic_update_real_arrays_iSendiRcv_noCudaAware(numArrays,arrays2comm)
#endif
   end subroutine mpi_halo_atomic_update_real_arrays
!-----------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------
!----------------------- MPI BASE COMMS -----------------------
!--------------------------------------------------------------

!------------- Send / Recv -------------

   subroutine mpi_base_comms_int_sendRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,tagComm,memPos_l,memSize
      !-------------------------------------------------

      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = commsMemSize(i)*numArrays

         call MPI_Sendrecv(aux_intField_s(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                           aux_intField_r(mempos_l), memSize, mpi_datatype_int, ngbRank, tagComm, &
                           app_comm, MPI_STATUS_IGNORE, mpi_err)
      end do
      !$acc end host_data
   end subroutine mpi_base_comms_int_sendRcv

   subroutine mpi_base_comms_real_sendRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,tagComm,memPos_l,memSize
      !-------------------------------------------------

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
   end subroutine mpi_base_comms_real_sendRcv

   subroutine mpi_base_comms_real_sendRcv_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,tagComm,memPos_l,memSize
      !-------------------------------------------------

      !$acc update host(aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = commsMemSize(i)*numArrays

         call MPI_Sendrecv(aux_realField_s(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                           aux_realField_r(mempos_l), memSize, mpi_datatype_real, ngbRank, tagComm, &
                           app_comm, MPI_STATUS_IGNORE, mpi_err)
      end do
   end subroutine mpi_base_comms_real_sendRcv_noCudaAware

!------------- iSend / iRecv -------------

   subroutine mpi_base_comms_int_iSendiRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)
      !-------------------------------------------------

      ireq=0
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = commsMemSize(i)*numArrays

         ireq = ireq+1
         call MPI_Irecv(aux_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      enddo
      !$acc end host_data

      call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

   end subroutine mpi_base_comms_int_iSendiRcv

   subroutine mpi_base_comms_int_iSendiRcv_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)
      !-------------------------------------------------

      ireq=0
      !$acc update host(aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = commsMemSize(i)*numArrays

         ireq = ireq+1
         call MPI_Irecv(aux_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do

      call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)
   end subroutine mpi_base_comms_int_iSendiRcv_noCudaAware

   subroutine mpi_base_comms_real_iSendiRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)
      !-------------------------------------------------

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
   end subroutine mpi_base_comms_real_iSendiRcv

   subroutine mpi_base_comms_real_iSendiRcv_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)
      !-------------------------------------------------

      ireq=0
      !$acc update host(aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         tagComm  = 0
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = commsMemSize(i)*numArrays

         ireq = ireq+1
         call MPI_Irecv(aux_realField_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_realField_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      call MPI_Waitall((2*numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

   end subroutine mpi_base_comms_real_iSendiRcv_noCudaAware

!------------- NCCL Send / Recv -------------

   subroutine mpi_base_comms_int_ncclSendRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !-------------------------------------------------

#ifdef NCCL_SUPPORT
      nccl_stat = ncclGroupStart()
      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = commsMemSize(i)*numArrays

         nccl_stat = ncclRecv(aux_intField_r(memPos_l), memSize, ncclInt, ngbRank, nccl_comm, nccl_stream)
         nccl_stat = ncclSend(aux_intField_s(memPos_l), memSize, ncclInt, ngbRank, nccl_comm, nccl_stream)
      end do
      !$acc end host_data
      nccl_stat = ncclGroupEnd()
      cuda_stat = cudaStreamSynchronize(nccl_stream)
#else
      write(*,*) 'Crashing in mpi_base_comms_int_ncclSendRcv! SOD2D compiled without NCCL support'
      call MPI_Abort(app_comm,-1,mpi_err)
#endif
   end subroutine mpi_base_comms_int_ncclSendRcv

   subroutine mpi_base_comms_real_ncclSendRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !-------------------------------------------------

#ifdef NCCL_SUPPORT
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
      write(*,*) 'Crashing in mpi_base_comms_real_ncclSendRcv! SOD2D compiled without NCCL support'
      call MPI_Abort(app_comm,-1,mpi_err)
#endif
   end subroutine mpi_base_comms_real_ncclSendRcv

!------------- Put -------------
! --- Put :: Fence ---

   subroutine mpi_base_comms_int_put_fence(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      call MPI_Win_fence(beginFence,window_id_int,mpi_err)
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data
      call MPI_Win_fence(endFence,window_id_int,mpi_err)

   end subroutine mpi_base_comms_int_put_fence

   subroutine mpi_base_comms_real_put_fence(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      call MPI_Win_fence(beginFence,window_id_real,mpi_err)
      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data
      call MPI_Win_fence(endFence,window_id_real,mpi_err)
   end subroutine mpi_base_comms_real_put_fence

   subroutine mpi_base_comms_real_put_fence_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

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

   end subroutine mpi_base_comms_real_put_fence_noCudaAware

! --- Put :: Post/Start/Complete/Wait ---

   subroutine mpi_base_comms_int_put_pscw(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err)
      call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err)

      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data

      call MPI_Win_complete(window_id_int,mpi_err)
      call MPI_Win_wait(window_id_int,mpi_err)

   end subroutine mpi_base_comms_int_put_pscw

   subroutine mpi_base_comms_real_put_pscw(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err)
      call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err)

      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
          ngbRank  = ranksToComm(i)
          memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
          memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
          memSize  = commsMemSize(i)*numArrays

          call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data

      call MPI_Win_complete(window_id_real,mpi_err)
      call MPI_Win_wait(window_id_real,mpi_err)

   end subroutine mpi_base_comms_real_put_pscw

! --- Put :: Lock ---

   subroutine mpi_base_comms_int_put_lock(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Put(aux_intField_s(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
      end do
      !$acc end host_data

      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

   end subroutine mpi_base_comms_int_put_lock

   subroutine mpi_base_comms_real_put_lock(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Put(aux_realField_s(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
      end do
      !$acc end host_data

      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

   end subroutine mpi_base_comms_real_put_lock

!------------- Get -------------
! --- Get :: Fence ---

   subroutine mpi_base_comms_int_get_fence(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      call MPI_Win_fence(beginFence,window_id_int,mpi_err)

      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Get(aux_intField_r(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data

      call MPI_Win_fence(endFence,window_id_int,mpi_err)

   end subroutine mpi_base_comms_int_get_fence

   subroutine mpi_base_comms_real_get_fence(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      call MPI_Win_fence(beginFence,window_id_real,mpi_err)

      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Get(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data

      call MPI_Win_fence(endFence,window_id_real,mpi_err)

   end subroutine mpi_base_comms_real_get_fence

   subroutine mpi_base_comms_real_get_fence_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      call MPI_Win_fence(beginFence,window_id_real,mpi_err)

      !$acc update host(aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Get(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do

      call MPI_Win_fence(endFence,window_id_real,mpi_err)

   end subroutine mpi_base_comms_real_get_fence_noCudaAware

! --- Get :: Post/Start/Complete/Wait ---

   subroutine mpi_base_comms_int_get_pscw(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call MPI_Win_post(commGroup,postAssert,window_id_int,mpi_err)
      call MPI_Win_start(commGroup,startAssert,window_id_int,mpi_err)

      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank   = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Get(aux_intField_r(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
      end do
      !$acc end host_data

      call MPI_Win_complete(window_id_int,mpi_err)
      call MPI_Win_wait(window_id_int,mpi_err)

   end subroutine mpi_base_comms_int_get_pscw

   subroutine mpi_base_comms_real_get_pscw(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      if(isPSCWBarrier) call MPI_Barrier(app_comm,mpi_err)
      call MPI_Win_post(commGroup,postAssert,window_id_real,mpi_err)
      call MPI_Win_start(commGroup,startAssert,window_id_real,mpi_err)

      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank   = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Get(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
      end do
      !$acc end host_data

      call MPI_Win_complete(window_id_real,mpi_err)
      call MPI_Win_wait(window_id_real,mpi_err)

   end subroutine mpi_base_comms_real_get_pscw

! --- Get :: Lock ---

   subroutine mpi_base_comms_int_get_lock(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      !$acc host_data use_device(aux_intField_r(:),aux_intField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_int,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Get(aux_intField_r(memPos_l),memSize,mpi_datatype_int,ngbRank,memPos_t,memSize,mpi_datatype_int,window_id_int,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_int,mpi_err)
      end do
      !$acc end host_data

      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

   end subroutine mpi_base_comms_int_get_lock

   subroutine mpi_base_comms_real_get_lock(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ngbRank,memPos_l,memSize
      !---------------------------------------

      !$acc host_data use_device(aux_realField_r(:),aux_realField_s(:))
      do i=1,numRanksWithComms
         ngbRank  = ranksToComm(i)
         memPos_l = (commsMemPosInLoc(i)-1)*numArrays+1
         memPos_t = (commsMemPosInNgb(i)-1)*numArrays !not adding +1 is because this value is the target displacement
         memSize  = commsMemSize(i)*numArrays

         call MPI_Win_lock(MPI_LOCK_SHARED,ngbRank,0,window_id_real,mpi_err); !May we can try MPI_LOCK_EXCLSUIVE
         call MPI_Get(aux_realField_r(memPos_l),memSize,mpi_datatype_real,ngbRank,memPos_t,memSize,mpi_datatype_real,window_id_real,mpi_err)
         call MPI_Win_unlock(ngbRank,window_id_real,mpi_err)
      end do
      !$acc end host_data

      if(isLockBarrier) call MPI_Barrier(app_comm,mpi_err)

   end subroutine mpi_base_comms_real_get_lock

!---------------------------- INTEGER -------------------------------------------

   ! INTEGER :: Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_int_sendRcv(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_sendRcv(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_sendRcv

   !INTEGER-MIN :: Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_int_sendRcv(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_sendRcv(numArrays)

      call copy_from_min_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_min_update_int_sendRcv

   !INTEGER :: iSend/iRecv ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_iSendiRcv(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_iSendiRcv(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_iSendiRcv

   !INTEGER-MIN :: iSend/iRecv ---------------------------------------------------------
   subroutine mpi_halo_atomic_min_update_int_iSendiRcv(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_iSendiRcv(numArrays)

      call copy_from_min_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_min_update_int_iSendiRcv

   !INTEGER :: NCCL Send/Recv ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_ncclSendRcv(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_ncclSendRcv(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_ncclSendRcv

   !INTEGER :: Put-Fence ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_put_fence(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_put_fence(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_put_fence

   !INTEGER :: Put-PSCW ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_put_pscw(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_put_pscw(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_put_pscw

   !INTEGER :: Put-Lock ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_put_lock(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_put_lock(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_put_lock

   !INTEGER :: Get-Fence ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_get_fence(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_get_fence(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_get_fence

   !INTEGER :: Get-PSCW ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_get_pscw(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_get_pscw(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_get_pscw

   !INTEGER :: Get-Lock ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_get_lock(intField)
      implicit none
      integer(4),intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_get_lock(numArrays)

      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_get_lock

!---------------------------- REAL -------------------------------------------

   ! REAL :: Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_sendRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_sendRcv(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_sendRcv

   ! REAL MASS-ENER-MOM :: Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_sendRcv(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_sendRcv(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_sendRcv

   ! REAL N-ARRAYS :: Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_sendRcv(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_sendRcv(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_sendRcv

   ! REAL MIN :: Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_sendRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_sendRcv(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_sendRcv

   ! REAL :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_iSendiRcv(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_iSendiRcv

   ! REAL MASS-ENER-MOM :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_iSendiRcv(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_iSendiRcv

   ! REAL N-ARRAYS :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv

   ! REAL MIN :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_iSendiRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_iSendiRcv

   ! REAL CONDITIONAL :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_conditional_ave_update_real_iSendiRcv(cond,realField)
      implicit none
      real(rp),intent(in) :: cond
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_conditional_ave_rcvBuffer_real(cond,realField)

   end subroutine mpi_halo_conditional_ave_update_real_iSendiRcv

   ! REAL ARRAYS-PTR :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arraysPtr_iSendiRcv(numArrays,arraysPtr2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      type(ptr_array1d_rp),intent(inout) :: arraysPtr2comm(numArrays)
      !----------------------------------------------------------------

      call fill_sendBuffer_arraysPtr_real(numArrays,arraysPtr2comm)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_rcvBuffer_arraysPtr_real(numArrays,arraysPtr2comm)

   end subroutine mpi_halo_atomic_update_real_arraysPtr_iSendiRcv

   ! REAL :: NCCL Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_ncclSendRcv(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_ncclSendRcv(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_ncclSendRcv

   ! REAL MASS-ENER-MOM :: NCCL Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_ncclSendRcv(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_ncclSendRcv(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_ncclSendRcv

   ! REAL N-ARRAYS :: NCCL Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_ncclSendRcv(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_ncclSendRcv(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_ncclSendRcv

   ! REAL MIN :: NCCL Send/Recv ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_ncclSendRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_iSendiRcv(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_ncclSendRcv

!------------- PUT FENCE -------------------------------------------

   ! REAL :: Put-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_put_fence(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_fence(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_put_fence

   ! REAL MASS-ENER-MOM :: Put-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_put_fence(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_put_fence(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_put_fence

   ! REAL N-ARRAYS :: Put-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_put_fence(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_put_fence(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_put_fence

   ! REAL MIN :: Put-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_put_fence(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_fence(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_put_fence

!------------ PUT PSCW -------------------------------------------

   !REAL :: Put-pscw -------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_put_pscw(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_pscw(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_put_pscw

   ! REAL MASS-ENER-MOM :: Put-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_put_pscw(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_put_pscw(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_put_pscw

   ! REAL N-ARRAYS :: Put-pscw ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_put_pscw(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_put_pscw(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_put_pscw

   ! REAL MIN :: Put-pscw ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_put_pscw(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_pscw(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_put_pscw

!------------ PUT LOCK -------------------------------------------

   !REAL :: Put-lock -------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_put_lock(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_lock(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_put_lock

   ! REAL MASS-ENER-MOM :: Put-lock ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_put_lock(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_put_lock(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_put_lock

   ! REAL N-ARRAYS :: Put-lock ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_put_lock(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_put_lock(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_put_lock

   ! REAL MIN :: Put-lock ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_put_lock(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_lock(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_put_lock

!------------ GET FENCE -------------------------------------------

   !REAL :: Get-fence -------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_get_fence(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_get_fence(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_get_fence

   ! REAL MASS-ENER-MOM :: get-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_massEnerMom_get_fence(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !--------------------------------------------------------

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_get_fence(numArrays)

      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_get_fence

   ! REAL N-ARRAYS :: get-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_update_real_arrays_get_fence(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !--------------------------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_get_fence(numArrays)

      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_atomic_update_real_arrays_get_fence

   ! REAL MIN :: get-fence ---------------------------------------------------
   subroutine mpi_halo_atomic_min_update_real_get_fence(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !---------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_get_fence(numArrays)

      call copy_from_min_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_min_update_real_get_fence

!------------ GET PSCW -------------------------------------------

   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_get_pscw(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_get_pscw(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_get_pscw

!------------ GET LOCK -------------------------------------------

   !REAL-------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_get_lock(realField)
      implicit none
      real(rp),intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_get_lock(numArrays)

      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_get_lock

!---------------------------------------------------------------------------------
!    FOR TESTING AND DEVEL STUFF!
!---------------------------------------------------------------------------------

!------------- ONLY BUFFERS -------------------------------------------
   !INTEGER ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_int_onlybuffers(intfield)
      implicit none
      integer(4), intent(inout) :: intfield(:)

      call fill_sendbuffer_int(intfield)
      call copy_from_rcvbuffer_int(intfield)

   end subroutine mpi_halo_atomic_update_int_onlybuffers
   !REAL ---------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_onlybuffers(realfield)
      implicit none
      real(rp), intent(inout) :: realfield(:)

      call fill_sendbuffer_real(realfield)
      call copy_from_rcvbuffer_real(realfield)

   end subroutine mpi_halo_atomic_update_real_onlybuffers
!---------------------------------------------------------------------------------
   subroutine mpi_halo_atomic_update_real_iSendiRcv_noCudaAware(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_iSendiRcv_noCudaAware(numArrays)

      !$acc update device(aux_realField_r(:))
      call copy_from_rcvBuffer_real(realField)

    end subroutine mpi_halo_atomic_update_real_iSendiRcv_noCudaAware

   subroutine mpi_halo_atomic_update_real_sendRcv_noCudaAware(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_sendRcv_noCudaAware(numArrays)

      !$acc update device(aux_realField_r(:))
      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_sendRcv_noCudaAware

   subroutine mpi_halo_atomic_update_real_put_fence_noCudaAware(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_put_fence_noCudaAware(numArrays)

      !$acc update device(aux_realField_r(:))
      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_put_fence_noCudaAware

   subroutine mpi_halo_atomic_update_real_get_fence_noCudaAware(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_real(realField)

      call mpi_base_comms_real_get_fence_noCudaAware(numArrays)

      !$acc update device(aux_realField_s(:))
      call copy_from_rcvBuffer_real(realField)

   end subroutine mpi_halo_atomic_update_real_get_fence_noCudaAware

   subroutine mpi_halo_atomic_update_real_massEnerMom_iSendiRcv_noCudaAware(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5

      call fill_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_comms_real_iSendiRcv_noCudaAware(numArrays)

      !$acc update device(aux_realField_r(:))
      call copy_from_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_atomic_update_real_massEnerMom_iSendiRcv_noCudaAware

   subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv_noCudaAware(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !-----------------------------------------

      call fill_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_comms_real_iSendiRcv_noCudaAware(numArrays)

      !$acc update device(aux_realField_r(:))
      call copy_from_rcvBuffer_arrays_real(numArrays,arrays2comm)

    end subroutine mpi_halo_atomic_update_real_arrays_iSendiRcv_noCudaAware

   subroutine mpi_halo_atomic_update_int_iSendiRcv_noCudaAware(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_sendBuffer_int(intField)

      call mpi_base_comms_int_iSendiRcv_noCudaAware(numArrays)

      !$acc update device(aux_intField_r(:))
      call copy_from_rcvBuffer_int(intField)

   end subroutine mpi_halo_atomic_update_int_iSendiRcv_noCudaAware

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

end module mod_comms
