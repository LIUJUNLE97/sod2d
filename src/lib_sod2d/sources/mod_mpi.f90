module mod_mpi
   use mpi
#ifndef NOACC
   use openacc
#endif
   implicit none

   integer :: mpi_rank, mpi_size, mpi_err   
   integer :: mpi_integer_size,mpi_double_size,mpi_float_size

   integer :: smNode_comm,smNode_rank,smNode_size
   integer :: num_devices, id_device

   contains

   subroutine init_mpi()
      implicit none

      call mpi_init(mpi_err)
      call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, mpi_err)
      call mpi_comm_size(MPI_COMM_WORLD, mpi_size, mpi_err)

      call MPI_Type_size(MPI_INTEGER, mpi_integer_size, mpi_err)
      call MPI_Type_size(MPI_FLOAT, mpi_float_size, mpi_err)
      call MPI_Type_size(MPI_DOUBLE, mpi_double_size, mpi_err)

      call init_sharedMemoryNode_comm()

   end subroutine init_mpi

   subroutine init_sharedMemoryNode_comm()
      implicit none
      call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,smNode_comm,mpi_err)

      call mpi_comm_rank(smNode_comm, smNode_rank, mpi_err)
      call mpi_comm_size(smNode_comm, smNode_size, mpi_err)
#ifndef NOACC
      num_devices = acc_get_num_devices(acc_device_nvidia);
      if(num_devices.ge.1) then
         id_device   = mod(smNode_rank,num_devices)
         call acc_set_device_num(id_device, acc_device_nvidia);
         write(*,*) 'rank ',smNode_rank,' num_devices ',num_devices,' id_device ',id_device
      else
         id_device = 0
         write(*,*) 'NO GPU FOUND IN THIS NODE!'
      end if
#else
      num_devices = 0
      id_device = 0
#endif

    end subroutine init_sharedMemoryNode_comm

    subroutine end_mpi()
        call mpi_finalize(mpi_err)
    end subroutine end_mpi

end module mod_mpi
