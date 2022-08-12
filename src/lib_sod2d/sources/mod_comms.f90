module mod_comms
    use mod_mpi_mesh
    use mod_parTimer
    implicit none

    !---- for the comms---
    integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
    integer(KIND=MPI_ADDRESS_KIND) :: memPos_t

    integer(4), dimension(:), allocatable :: aux_intField_s, aux_intField_r
    real(4), dimension(:), allocatable :: aux_floatField_s, aux_floatField_r
    real(8), dimension(:), allocatable :: aux_doubleField_s, aux_doubleField_r

    integer :: window_id_int,window_id_float,window_id_double
    integer :: window_id_sm

    logical :: isInt,isFloat,isDouble

    type(parTimer) :: timer_d1,timer_d2,timer_d3
    type(parTimer) :: timer_f1,timer_f2,timer_f3
    type(parTimer) :: timer_sm1,timer_sm2,timer_sm3

    !integer :: ms_rank,ms_size,ms_newComm
    type(c_ptr) :: c_ms_ptr



contains

    subroutine init_comms(useInt,useFloat,useDouble)
        implicit none
        logical, intent(in) :: useInt, useFloat, useDouble
        !using two buffers because if we use only one maybe we send the info after other proc has added something in me!
        !then i'll duplicate the info :S

        isInt=.false.
        isFloat=.false.
        isDouble=.false.

        if(useInt) then
            isInt = .true.

            allocate(aux_intField_s(numNodesToComm))
            allocate(aux_intField_r(numNodesToComm))
            !$acc enter data create(aux_intField_s(:))
            !$acc enter data create(aux_intField_r(:))

            call init_window_intField()
        end if

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

        if(isFloat) then
           !$acc exit data delete(aux_floatField_s(:))
           !$acc exit data delete(aux_floatField_r(:))
            deallocate(aux_floatField_s)
            deallocate(aux_floatField_r)

            call close_window_floatField()
        end if

        if(isDouble) then
            !$acc exit data delete(aux_doubleField_s(:))
            !$acc exit data delete(aux_doubleField_r(:))
            deallocate(aux_doubleField_s)
            deallocate(aux_doubleField_r)

            call close_window_doubleField()
        end if

    end subroutine end_comms

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


    subroutine update_and_comm_doubleField(doubleField)
        implicit none
        real(8), intent(inout) :: doubleField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize

        !call timer_d1%start_timer()

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_doubleField_s(i) = doubleField(iNodeL)
        end do
        aux_doubleField_r(:)=0.

        !call timer_d2%start_timer()

        !window_buffer_size = mpi_double_size*numNodesToComm
        !call MPI_Win_create(aux_doubleField_r,window_buffer_size,mpi_double_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id_double,mpi_err)
        call MPI_Win_fence(0,window_id_double,mpi_err)

        call timer_d3%start_timer()

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            !write(*,*) '[',mpi_rank,']iRank->',iRank
            !write(*,*) '[',mpi_rank,']memPos_l->',memPos_l
            !write(*,*) '[',mpi_rank,']memPos_t->',memPos_t
            !write(*,*) '[',mpi_rank,']memSize->',memSize

            call MPI_Accumulate(aux_doubleField_s(memPos_l),memSize,MPI_DOUBLE,&
                                iRank,memPos_t,memSize,MPI_DOUBLE,MPI_SUM,window_id_double,mpi_err)
        end do

        !call timer_d3%stop_timer()

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_double,mpi_err)
        !call MPI_Win_free(window_id_double,mpi_err)

        !call timer_d2%stop_timer()

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            !doubleField(iNodeL) = aux_doubleField_r(i) + aux_doubleField_s(i)
            doubleField(iNodeL) = doubleField(iNodeL) + aux_doubleField_r(i) 
        end do

        !call timer_d1%stop_timer()
    end subroutine update_and_comm_doubleField

    subroutine update_and_comm_floatField(floatField)
        implicit none
        real(4), intent(inout) :: floatField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize
 
        !call timer_f1%start_timer()
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = floatField(iNodeL)
        end do
        aux_floatField_r(:)=0.
        !call timer_f2%start_timer()

        call MPI_Win_fence(0,window_id_float,mpi_err)

        !call timer_f3%start_timer()
        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            !write(*,*) '[',mpi_rank,']iRank->',iRank
            !write(*,*) '[',mpi_rank,']memPos_l->',memPos_l
            !write(*,*) '[',mpi_rank,']memPos_t->',memPos_t
            !write(*,*) '[',mpi_rank,']memSize->',memSize

            call MPI_Put(aux_floatField_s(memPos_l),memSize,MPI_FLOAT,iRank,memPos_t,memSize,MPI_FLOAT,window_id_float,mpi_err)
            !call MPI_Accumulate(aux_floatField_s(memPos_l),memSize,MPI_FLOAT,&
            !                    iRank,memPos_t,memSize,MPI_FLOAT,MPI_SUM,window_id_float,mpi_err)
        end do

        !call timer_f3%stop_timer()

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_float,mpi_err)

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            floatField(iNodeL) = floatField(iNodeL) + aux_floatField_r(i)
        end do
        !call timer_f1%stop_timer()

    end subroutine update_and_comm_floatField

    subroutine update_and_comm_intField(intField)
        implicit none
        integer(4), intent(inout) :: intField(:)
        integer :: i,iNodeL,iRank
        integer :: memPos_l,memSize
 
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_intField_s(i) = intField(iNodeL)
        end do
        aux_intField_r(:)=0.

        call MPI_Win_fence(0,window_id_int,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            call MPI_Put(aux_intField_s(memPos_l),memSize,MPI_INTEGER,iRank,memPos_t,memSize,MPI_INTEGER,window_id_int,mpi_err)
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_int,mpi_err)

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            intField(iNodeL) = intField(iNodeL) + aux_intField_r(i)
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

        !call timer_sm1%start_timer()
 
        call MPI_Win_shared_query(window_id_sm,mpi_rank,ssize,memDisp,c_ms_ptr,mpi_err)
        call c_f_pointer(c_ms_ptr, f_ms_ptr, SHAPE = [size_p])

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_floatField_s(i) = floatField(iNodeL)
            f_ms_ptr(i) = 0.
        end do

        !call timer_sm2%start_timer()

        call MPI_Win_fence(0,window_id_sm,mpi_err)

        !call timer_sm3%start_timer()

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

        !call timer_sm3%stop_timer()

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id_sm,mpi_err)

        !call timer_sm2%stop_timer()

        call MPI_Win_shared_query(window_id_sm,mpi_rank,ssize,memDisp,c_ms_ptr,mpi_err)
        call c_f_pointer(c_ms_ptr, f_ms_ptr, SHAPE = [size_p])
        !write(*,*) 'ms[',mpi_rank,'] ssize',ssize,' memDisp ',memDisp
        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            floatField(iNodeL) = floatField(iNodeL) + f_ms_ptr(i)
        end do

        !call timer_sm1%stop_timer()

    end subroutine update_and_comm_shared_mem_floatField

!-------------------------------------------------------------------------------
! ####       TIMERS      -------------------------------------------------------

    subroutine print_dtimers()
        real(8) :: time_t1,time_t2,time_t3

        time_t1 = timer_d1%get_totalTime()
        time_t2 = timer_d2%get_totalTime()
        time_t3 = timer_d3%get_totalTime()

        write(*,*) 'd: t1 ',time_t1,' t2 ',time_t2,' t3 ',time_t3
    end subroutine print_dtimers

    subroutine print_ftimers()
        real(8) :: time_t1,time_t2,time_t3

        time_t1 = timer_f1%get_totalTime()
        time_t2 = timer_f2%get_totalTime()
        time_t3 = timer_f3%get_totalTime()

        write(*,*) 'f: t1 ',time_t1,' t2 ',time_t2,' t3 ',time_t3
    end subroutine print_ftimers

    subroutine print_smtimers()
        real(8) :: time_t1,time_t2,time_t3

        time_t1 = timer_sm1%get_totalTime()
        time_t2 = timer_sm2%get_totalTime()
        time_t3 = timer_sm3%get_totalTime()

        write(*,*) 'sm: t1 ',time_t1,' t2 ',time_t2,' t3 ',time_t3
    end subroutine print_smtimers

    subroutine init_dtimers()
        implicit none

        call timer_d1%init_timer()
        call timer_d2%init_timer()
        call timer_d3%init_timer()
    end subroutine init_dtimers

    subroutine init_ftimers()
        implicit none

        call timer_f1%init_timer()
        call timer_f2%init_timer()
        call timer_f3%init_timer()
    end subroutine init_ftimers

    subroutine init_smtimers()
        implicit none

        call timer_sm1%init_timer()
        call timer_sm2%init_timer()
        call timer_sm3%init_timer()
    end subroutine init_smtimers



!-------------------------------------------------------------------------------
! ####      OTHER FUNCS      ---------------------------------------------------

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

            call MPI_Accumulate(aux_floatField_s(memPos_l),memSize,MPI_FLOAT,&
                                iRank,memPos_t,memSize,MPI_FLOAT,MPI_SUM,window_id_float,mpi_err)
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
    

    !TO DO: TRY TO USE the MPI_Win_allocate_shared
    !read: https://stackoverflow.com/quest ons/37221134/mpi-how-to-use-mpi-win-allocate-shared-properly
! or the older form: INCLUDE ’mpif.h’
!MPI_COMM_SPLIT_TYPE(COMM, SPLIT_TYPE, KEY, INFO, NEWCOMM, IERROR)
!    INTEGER    COMM, SPLIT_TYPE, KEY, INFO, NEWCOMM, IERROR


end module mod_comms