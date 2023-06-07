module mod_smartredis
#ifdef SMARTREDIS

   use mod_constants
   use mod_mpi
   use smartredis_client, only: client_type
   implicit none

   type(client_type) :: client ! Client instance of SmartRedis to communicate with Redis database
   integer, dimension(:), allocatable, public :: state_sizes, state_displs, action_sizes, action_displs
   integer :: state_local_size, state_global_size, action_local_size, action_global_size
   contains

   ! Initialise SmartRedis client
   ! Since each MPI rank communicates with the database, all of them need to initialise the client
   subroutine init_smartredis(client, state_local_size2, action_local_size2, db_clustered)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: state_local_size2, action_local_size2
      logical, intent(in) :: db_clustered
      integer :: state_counter, action_counter, error, i
      logical :: is_error

      allocate(state_sizes(mpi_size))
      allocate(state_displs(mpi_size))
      allocate(action_sizes(mpi_size))
      allocate(action_displs(mpi_size))
      !$acc enter data create(state_sizes(:))
      !$acc enter data create(state_displs(:))
      !$acc enter data create(action_sizes(:))
      !$acc enter data create(action_displs(:))

      ! TODO: compute size of witness points per each rank and their displacements
      ! https://gist.github.com/jnvance/7b8cabebb06f91e2c1e788334f5de6c7
      ! 1. Gather the individual sizes to get total size and offsets in root process (0)
      call mpi_gather( &
         state_local_size2, 1, mpi_integer,   &  ! everyone sends 1 int from state_local_size
         state_sizes, 1, mpi_integer,        &  ! root receives 1 int from each proc into state_sizes
         0, mpi_comm_world, error            &  ! rank 0 is root
      )
      call mpi_gather( &
         action_local_size2, 1, mpi_integer, &  ! everyone sends 1 int from actions_local_size
         action_sizes, 1, mpi_integer,      &  ! root receives 1 int from each proc into s
         0, mpi_comm_world, error            &  ! rank 0 is root
      )
      ! 2. Compute displacements
      if (mpi_rank .eq. 0) then
         state_counter = 0
         action_counter = 0
         do i = 1, mpi_size
            state_displs(i) = state_counter
            action_displs(i) = action_counter
            state_counter = state_counter + state_sizes(i)
            action_counter = action_counter + action_sizes(i)
         end do
      end if

      ! Save vars
      state_local_size = state_local_size2
      action_local_size = action_local_size2
      call mpi_allreduce(state_local_size, state_global_size, 1, mpi_integer, mpi_sum, mpi_comm_world, mpi_err)
      call mpi_allreduce(action_local_size, action_global_size, 1, mpi_integer, mpi_sum, mpi_comm_world, mpi_err)

      ! Init client
      error = client%initialize(db_clustered)
      is_error = client%SR_error_parser(error)
      if (error /= 0) stop 'Error in SmartRedis client initialization'
   end subroutine init_smartredis

   ! Destroy SmartRedis client
   subroutine end_smartredis(client)
      type(client_type), intent(inout) :: client
      integer :: error
      logical :: is_error

      error = client%destructor()
      is_error = client%SR_error_parser(error)
      if (error /= 0) stop 'Error in SmartRedis client destruction'
   end subroutine end_smartredis

   ! Write witness points state into DB
   subroutine write_state(client, state_local, key)
      type(client_type), intent(inout) :: client
      real(rp), intent(in) :: state_local(state_local_size) ! local witness points state values
      character(len=*), intent(in) :: key                        ! state name to write to database
      integer :: error
      logical :: is_error
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
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error during SmartRedis state writting.'
      end if
   end subroutine write_state

   ! Read actions from DB
   ! subroutine read_actions(client, actions_local, key)
   !    type(client_type), intent(inout) :: client
   !    real(rp), intent(in) :: actions_local(actions_local_size) ! local actions values
   !    character(len=*), intent(in) :: key                       ! actions name to read from database

   !    integer, parameter :: interval = 10              ! polling interval in milliseconds
   !    integer, parameter :: tries = huge(1)            ! infinite number of polling tries
   !    logical(1) :: exists                             ! receives whether the tensor exists
   !    logical :: is_error
   !    integer :: found, error
   !    real(rp) :: actions_global(actions_global_size)

   !    ! wait (poll) until the actions array is found in the DB, then read, then delete
   !    if (mpi_rank .eq. 0) then
   !       found = client%poll_tensor(key, interval, tries, exists)
   !       is_error = client%SR_error_parser(found)
   !       if (found /= 0) stop 'Error in SmartRedis actions reading. Actions array not found.'
   !       error = client%unpack_tensor(key, actions_global, shape(actions_global))
   !       is_error = client%SR_error_parser(error)
   !       if (error /= 0) stop 'Error in SmartRedis actions reading. Tensor could not be unpacked.'
   !       error = client%delete_tensor(key)
   !       is_error = client%SR_error_parser(error)
   !       if (error /= 0) stop 'Error in SmartRedis actions reading. Tensor could not be deleted.'
   !    end if
   !    call mpi_barrier(mpi_comm_world, error) ! all processes wait for root to read actions

   !    ! scatter the global actions into local actions
   !    call mpi_scatterv( &
   !       actions_global, action_sizes, action_displs, mpi_datatype_real, &  ! actions_global is scattered according to action_sizes and displs
   !       actions_local, actions_local_size, mpi_datatype_real, &              ! actions_local to receive the data
   !       0, mpi_comm_world, error &                                           ! rank 0 is the one sending data
   !    )
   ! end subroutine read_actions

   ! Indicate environment time step status -> 0: init time step. 1: mid time step. 2: end time step
   subroutine write_step_type(client, step_type, key)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: step_type(1)
      character(len=*), intent(in) :: key

      integer :: error
      logical :: is_error

      if (mpi_rank .eq. 0) then
         error = client%put_tensor(key, step_type, shape(step_type))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis step_type writing.'
      end if
   end subroutine write_step_type

#endif
end module mod_smartredis