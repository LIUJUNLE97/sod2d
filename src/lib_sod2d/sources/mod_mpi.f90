module mod_mpi
   use mpi
   use mod_constants
#ifndef NOACC
   use openacc
#endif
   implicit none

   integer :: mpi_rank, mpi_size, mpi_err
   integer :: mpi_world_rank, mpi_world_size
   integer :: mpi_integer_size,mpi_int4_size,mpi_int8_size
   integer :: mpi_real_size,mpi_real4_size,mpi_real8_size
   integer :: mpi_datatype_int,mpi_datatype_int4,mpi_datatype_int8
   integer :: mpi_datatype_real,mpi_datatype_real4,mpi_datatype_real8

   integer :: smNode_comm,smNode_rank,smNode_size
   integer :: num_devices, id_device
   integer :: app_comm

   contains

   subroutine init_mpi()
      integer :: num_args, equal_pos, iarg, ierr, color
      character(len=16) :: arg, color_str="0"

      ! Get color for communicator if --tag is passed as command line argument
      num_args = command_argument_count()
      do iarg = 1, num_args
         call get_command_argument(iarg, arg)
         equal_pos = scan(adjustl(trim(arg)), "=")
         if (adjustl(trim(arg(:equal_pos-1))) .eq. "--tag") then
            color_str = trim(adjustl(arg(equal_pos+1:)))
         end if
      end do
      read(color_str, *) color

      ! Init MPI_COMM_WORLD communicator (shared accros all applications launched with mpirun ... : ...)
      call mpi_init(mpi_err)
      call mpi_comm_rank(MPI_COMM_WORLD, mpi_world_rank, mpi_err)
      call mpi_comm_size(MPI_COMM_WORLD, mpi_world_size, mpi_err)
      ! Split MPI_COMM_WORLD and create a communicator for this app only (color must be unique for each application)
      call MPI_Comm_split (MPI_COMM_WORLD, color, mpi_world_rank, app_comm, mpi_err)
      call MPI_Comm_rank (app_comm, mpi_rank, mpi_err)
      call MPI_Comm_size (app_comm, mpi_size , mpi_err)

      mpi_datatype_int = MPI_INTEGER
      mpi_datatype_int4 = MPI_INTEGER4
      mpi_datatype_int8 = MPI_INTEGER8
      mpi_datatype_real4 = MPI_REAL4
      mpi_datatype_real8 = MPI_REAL8

      if(rp.eq.4) then
         mpi_datatype_real = MPI_REAL4
      else if(rp.eq.8) then
         mpi_datatype_real = MPI_REAL8
      else
         write(*,*) 'Fatal error in init_mpi()! rp is not 4 or 8 >> CRASH!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      call MPI_Type_size(mpi_datatype_int,mpi_integer_size, mpi_err)
      call MPI_Type_size(mpi_datatype_int4,mpi_int4_size, mpi_err)
      call MPI_Type_size(mpi_datatype_int8,mpi_int8_size, mpi_err)

      call MPI_Type_size(mpi_datatype_real,mpi_real_size,mpi_err)
      call MPI_Type_size(mpi_datatype_real4,mpi_real4_size,mpi_err)
      call MPI_Type_size(mpi_datatype_real8,mpi_real8_size,mpi_err)

      call init_sharedMemoryNode_comm()

   end subroutine init_mpi

   subroutine init_sharedMemoryNode_comm()
      implicit none
      call MPI_Comm_split_type(app_comm,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,smNode_comm,mpi_err)

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
