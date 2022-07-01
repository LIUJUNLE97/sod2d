module mod_comms
    use mod_mpi_mesh
    implicit none

    !---- for the comms---
    integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
    integer(KIND=MPI_ADDRESS_KIND) :: memPos_t

    real(8), dimension(:), allocatable :: aux_dfield_s, aux_dfield_r

contains

    subroutine init_comms()
        implicit none
        !using two buffers because if we use only one maybe we send the info after other proc has added something in me!
        !then i'll duplicate the info :S
        allocate(aux_dfield_s(numNodesToComm))
        allocate(aux_dfield_r(numNodesToComm))
    end subroutine init_comms

    subroutine update_and_comm_dfield(dfield)
        implicit none
        real(8), intent(inout) :: dfield(:)
        integer :: i,iNodeL,iRank
        integer :: window_id,memPos_l,memSize

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            aux_dfield_s(i) = dfield(iNodeL)
        end do

        window_buffer_size = mpi_double_size*numNodesToComm
        call MPI_Win_create(aux_dfield_r,window_buffer_size,mpi_double_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
        call MPI_Win_fence(0,window_id,mpi_err)

        do i=1,numRanksWithComms
            iRank=ranksToComm(i)
            memPos_l  = commsMemPosInLoc(i)
            memPos_t  = commsMemPosInNgb(i) - 1 !the -1 is because this value is the target displacement

            memSize = commsMemSize(i)

            !write(*,*) '[',mpi_rank,']memPos_l->',memPos_l
            !write(*,*) '[',mpi_rank,']memPos_t->',memPos_t
            !write(*,*) '[',mpi_rank,']memSize->',memSize

            call MPI_Accumulate(aux_dfield_s(memPos_l),memSize,MPI_DOUBLE,&
                                iRank,memPos_t,memSize,MPI_DOUBLE,MPI_SUM,window_id,mpi_err)
        end do

        !! Wait for the MPI_Get issued to complete before going any further
        call MPI_Win_fence(0,window_id,mpi_err)
        call MPI_Win_free(window_id,mpi_err)

        do i=1,numNodesToComm
            iNodeL = matrixCommScheme(i,1)
            dfield(iNodeL) = aux_dfield_r(i) + aux_dfield_s(i)
        end do




    end subroutine update_and_comm_dfield

    !TO DO: TRY TO USE the MPI_Win_allocate_shared
    !read: https://stackoverflow.com/questions/37221134/mpi-how-to-use-mpi-win-allocate-shared-properly
! or the older form: INCLUDE ’mpif.h’
!MPI_COMM_SPLIT_TYPE(COMM, SPLIT_TYPE, KEY, INFO, NEWCOMM, IERROR)
!    INTEGER    COMM, SPLIT_TYPE, KEY, INFO, NEWCOMM, IERROR


end module mod_comms