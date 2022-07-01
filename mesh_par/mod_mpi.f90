module mod_mpi
    use mpi
    implicit none

    integer :: mpi_rank, mpi_size, mpi_err   
    integer :: mpi_integer_size,mpi_double_size,mpi_float_size

    contains

    subroutine init_mpi()
        call mpi_init(mpi_err)
        call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
        call mpi_comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)

        call MPI_Type_size(MPI_INTEGER, mpi_integer_size, mpi_err)
        call MPI_Type_size(MPI_FLOAT, mpi_float_size, mpi_err)
        call MPI_Type_size(MPI_DOUBLE, mpi_double_size, mpi_err)

    end subroutine init_mpi

    subroutine end_mpi()
        call mpi_finalize(mpi_err)
    end subroutine end_mpi

end module mod_mpi