module mod_smartredis
#if SMARTREDIS

   include "enum_fortran.inc"
   use iso_c_binding, only: c_bool
   use mod_constants
   use mod_numerical_params, only: do_smartredis
   use mod_mpi
   use smartredis_client, only: client_type
   implicit none

   type(client_type) :: client ! Client instance of SmartRedis to communicate with Redis database
   integer, dimension(:), allocatable, public :: state_sizes, state_displs, actions_sizes, actions_displs

   contains

   ! Initialise SmartRedis client
   ! Since each MPI rank communicates with the database, all of them need to initialise the client
   subroutine init_smartredis(client, state_local_size, actions_local_size)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: state_local_size, actions_local_size
      integer :: state_counter, actions_counter, error, i

      allocate(state_sizes(mpi_size))
      allocate(state_displs(mpi_size))
      allocate(actions_sizes(mpi_size))
      allocate(actions_displs(mpi_size))
      !$acc enter data create(state_sizes(:))
      !$acc enter data create(state_displs(:))
      !$acc enter data create(actions_sizes(:))
      !$acc enter data create(actions_displs(:))

      ! TODO: compute size of witness points per each rank and their displacements
      ! https://gist.github.com/jnvance/7b8cabebb06f91e2c1e788334f5de6c7
      ! 1. Gather the individual sizes to get total size and offsets in root process (0)
      call mpi_gather( &
         state_local_size, 1, mpi_integer,   &  ! everyone sends 1 int from state_local_size
         state_sizes, 1, mpi_integer,        &  ! root receives 1 int from each proc into state_sizes
         0, mpi_comm_world, error            &  ! rank 0 is root
      )
      call mpi_gather( &
         actions_local_size, 1, mpi_integer, &  ! everyone sends 1 int from actions_local_size
         actions_sizes, 1, mpi_integer,      &  ! root receives 1 int from each proc into s
         0, mpi_comm_world, error            &  ! rank 0 is root
      )
      ! 2. Compute displacements
      if (mpi_rank .eq. 0) then
         state_counter = 0
         actions_counter = 0
         do i = 1, mpi_size
            state_displs(i) = state_counter
            actions_displs(i) = actions_counter
            state_counter = state_counter + state_sizes(i)
            actions_counter = actions_counter + actions_sizes(i)
         end do
      end if

      error = client%initialize(do_smartredis)
      if (error /= SRNoError) stop 'Error in SmartRedis client initialization'
   end subroutine init_smartredis

   ! Destroy SmartRedis client
   subroutine end_smartredis(client)
      type(client_type), intent(inout) :: client
      integer :: error

      error = client%destructor()
      if (error /= SRNoError) stop 'Error in SmartRedis client destruction'
   end subroutine end_smartredis

   ! Write witness points state into DB
   subroutine write_state(client, state_local_size, state_global_size, state_local, key)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: state_local_size               ! local number of witness points
      integer, intent(in) :: state_global_size              ! total number of witness points
      real(rp), intent(in) :: state_local(state_local_size) ! local witness points state values
      character(len=*), intent(in) :: key                   ! state name to write to database

      integer :: error
      real(rp) :: state_global(state_global_size)

      ! gather the local states into a global state
      call mpi_gatherv( &
         state_local, state_local_size, mpi_datatype_real, &           ! everyone sends state_local data
         state_global, state_sizes, state_displs, mpi_datatype_real, & ! root receives it into state_global
         0, mpi_comm_world, error &                                    ! rank 0 is root
      )

      ! write global state into DB
      if (mpi_rank .eq. 0) then
         error = client%put_tensor(key, state_global, shape(state_global))
         if (error /= SRNoError) stop 'Error during SmartRedis state writting.'
      end if
   end subroutine write_state

   ! Read actions from DB
   subroutine read_actions(client, actions_local_size, actions_global_size, actions_local, key)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: actions_local_size                 ! local number of actions
      integer, intent(in) :: actions_global_size                ! total number of actions
      real(rp), intent(in) :: actions_local(actions_local_size) ! local actions values
      character(len=*), intent(in) :: key                       ! actions name to read from database

      integer, parameter :: interval = 10              ! polling interval in milliseconds
      integer, parameter :: tries = huge(1)            ! infinite number of polling tries
      logical(kind=c_bool) :: exists                   ! receives whether the tensor exists
      integer :: found, error
      real(rp) :: actions_global(actions_global_size)

      ! wait (poll) until the actions array is found in the DB, then read
      if (mpi_rank .eq. 0) then
         found = client%poll_tensor(key, interval, tries, exists)
         if (found /= 0) stop 'Error in SmartRedis actions reading. Actions array not found.'
         error = client%unpack_tensor(key, actions_global, shape(actions_global))
         if (error /= SRNoError) stop 'Error in SmartRedis actions reading. Tensor could not be deleted.'
      end if
      call mpi_barrier(mpi_comm_world, error) ! all processes wait for root to read actions

      ! scatter the global actions into local actions
      call mpi_scatterv( &
         actions_global, actions_sizes, actions_displs, mpi_datatype_real, &  ! actions_global is scattered according to actions_sizes and displs
         actions_local, actions_local_size, mpi_datatype_real, &              ! actions_local to receive the data
         0, mpi_comm_world, error &                                           ! rank 0 is the one sending data
      )
   end subroutine read_actions

#endif
end module mod_smartredis