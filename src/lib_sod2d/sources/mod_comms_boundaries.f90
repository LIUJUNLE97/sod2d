module mod_comms_boundaries
   use mod_mpi_mesh
#ifdef NCCL_SUPPORT
   use cudafor
   use nccl
#endif

!-- Select type of communication for mpi_bnd_atomic
#define _ISENDIRCV_ 1
#define _NOCUDAAWARE_ 0

#ifdef AUTOCOMP
 #define _NOCUDAAWARE_ 1
 #define _ISENDIRCV_ 0
#endif

   implicit none

   !---- for the comms---
   integer(4) :: bnd_worldGroup,bnd_commGroup

   integer(4), dimension(:), allocatable :: aux_bnd_intField_s, aux_bnd_intField_r
   real(rp), dimension(:), allocatable :: aux_bnd_realField_s, aux_bnd_realField_r

   integer(4) :: window_id_bnd_int,window_id_bnd_real

   logical :: bnd_isInt,bnd_isReal

#ifdef NCCL_SUPPORT
   integer                        :: cuda_bnd_stat
   type(ncclResult)               :: nccl_bnd_stat
   type(ncclComm)                 :: nccl_bnd_comm
   integer(kind=cuda_stream_kind) :: nccl_bnd_stream
#endif

contains

   subroutine init_comms_bnd(useInt,useReal,intBndBufferMultIn,realBndBufferMultIn)
      implicit none
      logical, intent(in) :: useInt,useReal
      integer(4),intent(in),optional :: intBndBufferMultIn,realBndBufferMultIn
      integer(4) :: intBndBufferArraySize,realBndBufferArraySize,intBndBufferMult=1,realBndBufferMult=1
      integer(kind=mpi_address_kind) :: win_buffer_size
      logical :: useFenceFlags,useAssertNoCheckFlags,useLockBarrier

      intBndBufferMult  = 1
      realBndBufferMult = 1

      if(present(intBndBufferMultIn)) then
         intBndBufferMult = intBndBufferMultIn
      end if

      if(present(realBndBufferMultIn)) then
         realBndBufferMult = realBndBufferMultIn
      end if

#if _ISENDIRCV_
      if(mpi_rank.eq.0) write(111,*) "--| Boundary Comm. scheme: iSend-iRecv"
#endif
#if _NOCUDAAWARE_
      if(mpi_rank.eq.0) write(111,*) "--| Boundary Comm. scheme: NO CUDA-Aware iSend-iRecv"
#endif

      bnd_isInt=.false.
      bnd_isReal=.false.

      if(useInt) then
         bnd_isInt = .true.

         intBndBufferArraySize = intBndBufferMult*bnd_numNodesToComm 

         allocate(aux_bnd_intField_s(intBndBufferArraySize))
         allocate(aux_bnd_intField_r(intBndBufferArraySize))
         !$acc enter data create(aux_bnd_intField_s(:))
         !$acc enter data create(aux_bnd_intField_r(:))

         !call init_window_intField()
      end if

      if(useReal) then
         bnd_isReal = .true.

         realBndBufferArraySize = realBndBufferMult*bnd_numNodesToComm 

         allocate(aux_bnd_realField_s(realBndBufferArraySize))
         allocate(aux_bnd_realField_r(realBndBufferArraySize))
         !$acc enter data create(aux_bnd_realField_s(:))
         !$acc enter data create(aux_bnd_realField_r(:))

         !call init_window_realField()
      end if

      call MPI_Comm_group(app_comm,bnd_worldGroup,mpi_err)
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
      end if

      if(bnd_isReal) then
         !$acc exit data delete(aux_bnd_realField_s(:))
         !$acc exit data delete(aux_bnd_realField_r(:))
          deallocate(aux_bnd_realField_s)
          deallocate(aux_bnd_realField_r)
      end if

#ifdef NCCL_SUPPORT
      nccl_bnd_stat = ncclCommDestroy(nccl_bnd_comm)
      cuda_bnd_stat = cudaStreamDestroy(nccl_bnd_stream)
#endif

   end subroutine end_comms_bnd
!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------

   subroutine init_bnd_put_window_intField(intBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: intBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      win_buffer_size = mpi_integer_size*intBndBufferArraySize

      !$acc host_data use_device(aux_bnd_intField_r(:))
      call MPI_Win_create(aux_bnd_intField_r,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_bnd_int,mpi_err)
      !$acc end host_data

   end subroutine init_bnd_put_window_intField

   subroutine init_bnd_put_window_intField_noCudaAware(intBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: intBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      win_buffer_size = mpi_integer_size*intBndBufferArraySize

      call MPI_Win_create(aux_bnd_intField_r,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_bnd_int,mpi_err)

   end subroutine init_bnd_put_window_intField_noCudaAware

   subroutine init_bnd_get_window_intField(intBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: intBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      win_buffer_size = mpi_integer_size*intBndBufferArraySize

      !$acc host_data use_device(aux_bnd_intField_s(:))
      call MPI_Win_create(aux_bnd_intField_s,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_bnd_int,mpi_err)
      !$acc end host_data

   end subroutine init_bnd_get_window_intField

   subroutine init_bnd_get_window_intField_noCudaAware(intBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: intBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      win_buffer_size = mpi_integer_size*intBndBufferArraySize

      call MPI_Win_create(aux_bnd_intField_s,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id_bnd_int,mpi_err)

   end subroutine init_bnd_get_window_intField_noCudaAware

   subroutine close_bnd_window_intField()
      implicit none

      call MPI_Win_free(window_id_bnd_int,mpi_err)
   end subroutine close_bnd_window_intField

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------

   subroutine init_bnd_put_window_realField(realBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: realBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------
      
      win_buffer_size = mpi_real_size*realBndBufferArraySize

      !$acc host_data use_device(aux_bnd_realField_r(:))
      call MPI_Win_create(aux_bnd_realField_r,win_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_bnd_real,mpi_err)
      !$acc end host_data

   end subroutine init_bnd_put_window_realField

   subroutine init_bnd_put_window_realField_noCudaAware(realBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: realBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      call MPI_Win_create(aux_bnd_realField_r,win_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_bnd_real,mpi_err)

   end subroutine init_bnd_put_window_realField_noCudaAware

   subroutine init_bnd_get_window_realField(realBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: realBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      !$acc host_data use_device(aux_bnd_realField_s(:))
      call MPI_Win_create(aux_bnd_realField_s,win_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_bnd_real,mpi_err)
      !$acc end host_data

   end subroutine init_bnd_get_window_realField

   subroutine init_bnd_get_window_realField_noCudaAware(realBndBufferArraySize)
      implicit none
      integer(4),intent(in) :: realBndBufferArraySize
      integer(kind=mpi_address_kind) :: win_buffer_size
      !------------------------------------------------------------

      call MPI_Win_create(aux_bnd_realField_s,win_buffer_size,mpi_real_size,MPI_INFO_NULL,app_comm,window_id_bnd_real,mpi_err)

   end subroutine init_bnd_get_window_realField_noCudaAware

   subroutine close_bnd_window_realField()
      implicit none

       call MPI_Win_free(window_id_bnd_real,mpi_err)
   end subroutine close_bnd_window_realField

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!   Standard fill boundary buffers
!-------------------------------------------------------------------------------------
    subroutine fill_bnd_sendBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop async(1)
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            aux_bnd_intField_s(i) = intField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_bnd_intField_r(:)=0
        !$acc end kernels

        !$acc wait
    end subroutine fill_bnd_sendBuffer_int
!-------------------------------------------------------------------------------------
    subroutine fill_bnd_sendBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop async(1)
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            aux_bnd_realField_s(i) = realField(iNodeL)
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_bnd_realField_r(:)=0.
        !$acc end kernels

        !$acc wait
    end subroutine fill_bnd_sendBuffer_real
!-------------------------------------------------------------------------------------
    subroutine fill_bnd_sendBuffer_arrays_real(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp), intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,iNodeL,iArray

        !$acc parallel loop async(1)
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            do iArray = 1,numArrays
                aux_bnd_realField_s((i-1)*numArrays+iArray) = arrays2comm(iNodeL,iArray)
            end do
        end do
        !$acc end parallel loop
        !$acc kernels async(2)
        aux_bnd_realField_r(:)=0.
        !$acc end kernels

        !$acc wait
    end subroutine fill_bnd_sendBuffer_arrays_real
!-------------------------------------------------------------------------------------
    subroutine fill_bnd_sendBuffer_massEnerMom_real(mass,ener,momentum)
        implicit none
        real(rp),intent(in) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,idime,iNodeL
        integer(4),parameter :: nArrays=5

        !$acc parallel loop async(1)
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            do idime = 1,ndime
                aux_bnd_realField_s((i-1)*nArrays+idime) = momentum(iNodeL,idime)
            end do
            aux_bnd_realField_s((i-1)*nArrays+4) = mass(iNodeL)
            aux_bnd_realField_s((i-1)*nArrays+5) = ener(iNodeL)
        end do
        !$acc end parallel loop

        !$acc kernels async(2)
        aux_bnd_realField_r(:)=0.
        !$acc end kernels

        !$acc wait
    end subroutine fill_bnd_sendBuffer_massEnerMom_real

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!   Standard copy boundary buffers
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
   subroutine copy_from_bnd_rcvBuffer_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4) :: i,iNodeL

      !$acc parallel loop
      do i=1,bnd_numNodesToComm
         iNodeL = bnd_nodesToComm(i)
         !$acc atomic update
         intField(iNodeL) = intField(iNodeL) + aux_bnd_intField_r(i)
         !$acc end atomic
      end do
      !$acc end parallel loop
   end subroutine copy_from_bnd_rcvBuffer_int
!-------------------------------------------------------------------------
   subroutine copy_from_bnd_rcvBuffer_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4) :: i,iNodeL

      !$acc parallel loop
      do i=1,bnd_numNodesToComm
         iNodeL = bnd_nodesToComm(i)
         !$acc atomic update
         realField(iNodeL) = realField(iNodeL) + aux_bnd_realField_r(i)
         !$acc end atomic
      end do
      !$acc end parallel loop
   end subroutine copy_from_bnd_rcvBuffer_real
!-------------------------------------------------------------------------
   subroutine copy_from_bnd_rcvBuffer_arrays_real(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp),intent(inout) :: arrays2comm(:,:)
      integer(4) :: i,iArray,iNodeL

      !$acc parallel loop
      do i=1,bnd_numNodesToComm
         iNodeL = bnd_nodesToComm(i)
         do iArray = 1,numArrays
            !$acc atomic update
            arrays2comm(iNodeL,iArray) = arrays2comm(iNodeL,iArray) + aux_bnd_realField_r((i-1)*numArrays+iArray)
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop
   end subroutine copy_from_bnd_rcvBuffer_arrays_real

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!   Specialized copy boundary buffers
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
    subroutine copy_from_min_bnd_rcvBuffer_int(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            !$acc atomic update
            intField(iNodeL) = min(intField(iNodeL),aux_bnd_intField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_min_bnd_rcvBuffer_int
!-------------------------------------------------------------------------
    subroutine copy_from_max_bnd_rcvBuffer_real(realField)
        implicit none
        real(rp), intent(inout) :: realField(:)
        integer(4) :: i,iNodeL

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            !$acc atomic update
            realField(iNodeL) = max(realField(iNodeL),aux_bnd_realField_r(i))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_max_bnd_rcvBuffer_real

    subroutine copy_from_max_bnd_rcvBuffer_arrays_real(numArrays,arrays2comm)
        implicit none
        integer(4),intent(in) :: numArrays
        real(rp), intent(inout) :: arrays2comm(:,:)
        integer(4) :: i,iNodeL,iArray

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            do iArray = 1,numArrays
                !$acc atomic update
                arrays2comm(iNodeL,iArray) = max(arrays2comm(iNodeL,iArray), aux_bnd_realField_r((i-1)*numArrays+iArray))
                !$acc end atomic
            end do
        end do
        !$acc end parallel loop
    end subroutine copy_from_max_bnd_rcvBuffer_arrays_real

    subroutine copy_from_max_bnd_rcvBuffer_massEnerMom_real(mass,ener,momentum)
        implicit none
        real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
        integer(4) :: i,idime,iNodeL
        integer(4),parameter :: numArrays=5

        !$acc parallel loop
        do i=1,bnd_numNodesToComm
            iNodeL = bnd_nodesToComm(i)
            do idime = 1,ndime
                !$acc atomic update
                momentum(iNodeL,idime) = max(momentum(iNodeL,idime), aux_bnd_realField_r((i-1)*numArrays+idime))
                !$acc end atomic
            end do
            !$acc atomic update
            mass(iNodeL) = max(mass(iNodeL), aux_bnd_realField_r((i-1)*numArrays+4))
            !$acc end atomic
            !$acc atomic update
            ener(iNodeL) = max(ener(iNodeL), aux_bnd_realField_r((i-1)*numArrays+5))
            !$acc end atomic
        end do
        !$acc end parallel loop
    end subroutine copy_from_max_bnd_rcvBuffer_massEnerMom_real

!-----------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------

   subroutine mpi_halo_bnd_atomic_int(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)

#if _ISENDIRCV_
      call mpi_halo_bnd_atomic_int_iSendiRcv(intField)
#endif
#if _NOCUDAAWARE_
      call mpi_halo_bnd_atomic_int_iSendiRcv_noCudaAware(intField)
#endif

   end subroutine mpi_halo_bnd_atomic_int

   subroutine mpi_halo_bnd_atomic_real(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)

#if _ISENDIRCV_
      call mpi_halo_bnd_atomic_real_iSendiRcv(realField)
#endif
#if _NOCUDAAWARE_
      call mpi_halo_bnd_atomic_real_iSendiRcv_noCudaAware(realField)
#endif

   end subroutine mpi_halo_bnd_atomic_real

!-----------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------
!---------------- MPI BASE BOUNDARY COMMS ---------------------
!--------------------------------------------------------------

   subroutine mpi_base_bnd_comms_int_iSendiRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*bnd_numRanksWithComms)
      !-------------------------------------------------

      ireq=0
      !$acc host_data use_device (aux_bnd_intField_r(:),aux_bnd_intField_s(:))
      do i=1,bnd_numRanksWithComms
         ngbRank  = bnd_ranksToComm(i)
         tagComm  = 0
         memPos_l = (bnd_commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = bnd_commsMemSize(i)*numArrays

         ireq = ireq+1
         call MPI_Irecv(aux_bnd_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_bnd_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      !$acc end host_data

      call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

   end subroutine mpi_base_bnd_comms_int_iSendiRcv

   subroutine mpi_base_bnd_comms_int_iSendiRcv_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*bnd_numRanksWithComms)
      !-------------------------------------------------

      !$acc update host(aux_bnd_intField_s(:))

      ireq=0
      do i=1,bnd_numRanksWithComms
         ngbRank  = bnd_ranksToComm(i)
         tagComm  = 0
         memPos_l = (bnd_commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = bnd_commsMemSize(i)*numArrays

         ireq = ireq+1
         call MPI_Irecv(aux_bnd_intField_r(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_bnd_intField_s(mempos_l),memSize,mpi_datatype_int,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

      !$acc update device(aux_bnd_intField_r(:))

   end subroutine mpi_base_bnd_comms_int_iSendiRcv_noCudaAware

   subroutine mpi_base_bnd_comms_real_iSendiRcv(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)
      !-------------------------------------------------

      ireq=0
      !$acc host_data use_device (aux_bnd_realField_r(:),aux_bnd_realField_s(:))
      do i=1,bnd_numRanksWithComms
         ngbRank  = bnd_ranksToComm(i)
         tagComm  = 0
         memPos_l = (bnd_commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = bnd_commsMemSize(i)*numArrays


         ireq = ireq+1
         call MPI_Irecv(aux_bnd_realField_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_bnd_realField_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      !$acc end host_data

      call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

   end subroutine mpi_base_bnd_comms_real_iSendiRcv


   subroutine mpi_base_bnd_comms_real_iSendiRcv_noCudaAware(numArrays)
      implicit none
      integer(4),intent(in) :: numArrays
      integer(4) :: i,ireq,ngbRank,tagComm,memPos_l,memSize
      integer(4) :: requests(2*numRanksWithComms)
      !-------------------------------------------------

      !$acc update host(aux_bnd_realField_s(:))

      ireq=0
      do i=1,bnd_numRanksWithComms
         ngbRank  = bnd_ranksToComm(i)
         tagComm  = 0
         memPos_l = (bnd_commsMemPosInLoc(i)-1)*numArrays+1
         memSize  = bnd_commsMemSize(i)*numArrays

         ireq = ireq+1
         call MPI_Irecv(aux_bnd_realField_r(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
         ireq = ireq+1
         call MPI_ISend(aux_bnd_realField_s(mempos_l),memSize,mpi_datatype_real,ngbRank,tagComm,app_comm,requests(ireq),mpi_err)
      end do
      call MPI_Waitall((2*bnd_numRanksWithComms),requests,MPI_STATUSES_IGNORE,mpi_err)

      !$acc update device(aux_bnd_realField_r(:))

   end subroutine mpi_base_bnd_comms_real_iSendiRcv_noCudaAware

!--------------------------------------------------------------------------------
!---------------------------- INTEGER -------------------------------------------
   
   !INTEGER :: iSend/iRcv ---------------------------------------------------------
   subroutine mpi_halo_bnd_atomic_int_iSendiRcv(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_int(intField)

      call mpi_base_bnd_comms_int_iSendiRcv(numArrays) 

      call copy_from_bnd_rcvBuffer_int(intField)

   end subroutine mpi_halo_bnd_atomic_int_iSendiRcv

   !INTEGER-MIN :: iSend/iRcv ---------------------------------------------------------
   subroutine mpi_halo_bnd_atomic_min_int_iSendiRcv(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_int(intField)

      call mpi_base_bnd_comms_int_iSendiRcv(numArrays)

      call copy_from_min_bnd_rcvBuffer_int(intField)

   end subroutine mpi_halo_bnd_atomic_min_int_iSendiRcv

!---------------------------- REAL -------------------------------------------

   ! REAL :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_real_iSendiRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_real(realField)

      call mpi_base_bnd_comms_real_iSendiRcv(numArrays)

      call copy_from_bnd_rcvBuffer_real(realField)

   end subroutine mpi_halo_bnd_atomic_real_iSendiRcv

   ! REAL N-ARRAYS :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_real_arrays_iSendiRcv(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !-----------------------------------------

      call fill_bnd_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_bnd_comms_real_iSendiRcv(numArrays)

      call copy_from_bnd_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_bnd_atomic_real_arrays_iSendiRcv

   ! REAL-MAX :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_max_real_iSendiRcv(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_real(realField)

      call mpi_base_bnd_comms_real_iSendiRcv(numArrays)

      call copy_from_max_bnd_rcvBuffer_real(realField)

   end subroutine mpi_halo_bnd_atomic_max_real_iSendiRcv

   ! REAL-MAX N-ARRAYS :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_max_real_arrays_iSendiRcv(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !-----------------------------------------

      call fill_bnd_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_bnd_comms_real_iSendiRcv(numArrays)

      call copy_from_max_bnd_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_bnd_atomic_max_real_arrays_iSendiRcv

   ! REAL-MAX MASS-ENER-MOM :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_max_real_massEnerMom_iSendiRcv(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !-----------------------------------------

      call fill_bnd_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_bnd_comms_real_iSendiRcv(numArrays)

      call copy_from_max_bnd_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_bnd_atomic_max_real_massEnerMom_iSendiRcv

!----------------------------------------------------------------------------------
!        NO-CUDA-Aware versions
!----------------------------------------------------------------------------------

   !INTEGER :: iSend/iRcv ---------------------------------------------------------
   subroutine mpi_halo_bnd_atomic_int_iSendiRcv_noCudaAware(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_int(intField)

      call mpi_base_bnd_comms_int_iSendiRcv_noCudaAware(numArrays) 

      call copy_from_bnd_rcvBuffer_int(intField)
      
   end subroutine mpi_halo_bnd_atomic_int_iSendiRcv_noCudaAware

   !INTEGER-MIN :: iSend/iRcv ---------------------------------------------------------
   subroutine mpi_halo_bnd_atomic_min_int_iSendiRcv_noCudaAware(intField)
      implicit none
      integer(4), intent(inout) :: intField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_int(intField)

      call mpi_base_bnd_comms_int_iSendiRcv_noCudaAware(numArrays)

      call copy_from_min_bnd_rcvBuffer_int(intField)

   end subroutine mpi_halo_bnd_atomic_min_int_iSendiRcv_noCudaAware

!---------------------------- REAL -------------------------------------------

   ! REAL :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_real_iSendiRcv_noCudaAware(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_real(realField)

      call mpi_base_bnd_comms_real_iSendiRcv_noCudaAware(numArrays)

      call copy_from_bnd_rcvBuffer_real(realField)

   end subroutine mpi_halo_bnd_atomic_real_iSendiRcv_noCudaAware

   ! REAL-MAX :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_max_real_iSendiRcv_noCudaAware(realField)
      implicit none
      real(rp), intent(inout) :: realField(:)
      integer(4),parameter :: numArrays=1
      !-----------------------------------------

      call fill_bnd_sendBuffer_real(realField)

      call mpi_base_bnd_comms_real_iSendiRcv_noCudaAware(numArrays)

      call copy_from_max_bnd_rcvBuffer_real(realField)

   end subroutine mpi_halo_bnd_atomic_max_real_iSendiRcv_noCudaAware

   ! REAL-MAX N-ARRAYS :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_max_real_arrays_iSendiRcv_noCudaAware(numArrays,arrays2comm)
      implicit none
      integer(4),intent(in) :: numArrays
      real(rp), intent(inout) :: arrays2comm(:,:)
      !-----------------------------------------

      call fill_bnd_sendBuffer_arrays_real(numArrays,arrays2comm)

      call mpi_base_bnd_comms_real_iSendiRcv_noCudaAware(numArrays)

      call copy_from_max_bnd_rcvBuffer_arrays_real(numArrays,arrays2comm)

   end subroutine mpi_halo_bnd_atomic_max_real_arrays_iSendiRcv_noCudaAware

   ! REAL-MAX MASS-ENER-MOM :: iSend/iRecv ---------------------------------------------------
   subroutine mpi_halo_bnd_atomic_max_real_massEnerMom_iSendiRcv_noCudaAware(mass,ener,momentum)
      implicit none
      real(rp),intent(inout) :: mass(:),ener(:),momentum(:,:)
      integer(4),parameter :: numArrays=5
      !-----------------------------------------

      call fill_bnd_sendBuffer_massEnerMom_real(mass,ener,momentum)

      call mpi_base_bnd_comms_real_iSendiRcv_noCudaAware(numArrays)

      call copy_from_max_bnd_rcvBuffer_massEnerMom_real(mass,ener,momentum)

   end subroutine mpi_halo_bnd_atomic_max_real_massEnerMom_iSendiRcv_noCudaAware

end module mod_comms_boundaries
