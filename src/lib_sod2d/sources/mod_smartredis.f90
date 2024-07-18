module mod_smartredis
#ifdef SMARTREDIS

   use mod_constants
   use mod_mpi
   use smartredis_client, only: client_type
   implicit none

   type(client_type) :: client ! Client instance of SmartRedis to communicate with Redis database
   integer, dimension(:), allocatable :: state_sizes, state_displs
   real(rp), dimension(:), allocatable :: action_global, action_global_previous
   integer :: state_local_size, state_global_size, action_global_size, step_type_mod, n_pseudo_envs
   contains

   ! Initialise SmartRedis client
   ! State is stored in arrays of different sizes on each MPI rank. Actions is a global array living in all processes.
   subroutine init_smartredis(client, state_local_size2, action_global_size2, n_pseudo_envs2, tag, db_clustered)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: state_local_size2, action_global_size2, n_pseudo_envs2
      character(len=*), intent(in) :: tag
      logical, intent(in) :: db_clustered

      integer :: state_counter, error, i, state_global_size_tensor(1), action_global_size_tensor(1)
      logical :: is_error

      allocate(state_sizes(mpi_size))
      allocate(state_displs(mpi_size))
      allocate(action_global(action_global_size2))
      allocate(action_global_previous(action_global_size2))
      !$acc enter data create(state_sizes(:))
      !$acc enter data create(state_displs(:))
      !$acc enter data create(action_global(:))
      !$acc enter data create(action_global_previous(:))

      action_global(:) = 0.0_rp
      action_global_previous(:) = 0.0_rp
      !$acc update device(action_global(:))
      !$acc update device(action_global_previous(:))

      ! https://gist.github.com/jnvance/7b8cabebb06f91e2c1e788334f5de6c7
      ! 1. Gather the individual sizes to get total size and offsets in root process (0)
      call mpi_gather( &
         state_local_size2, 1, mpi_integer,  &  ! everyone sends 1 int from state_local_size
         state_sizes, 1, mpi_integer,        &  ! root receives 1 int from each proc into state_sizes
         0, app_comm, error                  &  ! rank 0 is root
      )

      ! 2. Compute displacements
      if (mpi_rank .eq. 0) then
         state_counter = 0
         do i = 1, mpi_size
            state_displs(i) = state_counter
            state_counter = state_counter + state_sizes(i)
         end do
      end if

      ! Store in module variables
      n_pseudo_envs = n_pseudo_envs2
      action_global_size = action_global_size2
      state_local_size = state_local_size2
      call mpi_allreduce(state_local_size, state_global_size, 1, mpi_integer, mpi_sum, app_comm, mpi_err)

      ! Init client (only root process!) and write global state and action sizes into DB.
      if (mpi_rank .eq. 0) then
         error = client%initialize(db_clustered)
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis client initialization'

         ! Write global size of state and action into DB
         state_global_size_tensor(1) = state_global_size
         action_global_size_tensor(1) = action_global_size

         ! if (tag == "0") then
            error = client%put_tensor("state_size", state_global_size_tensor, shape(state_global_size_tensor))
            is_error = client%SR_error_parser(error)
            if (error /= 0) stop 'Error during SmartRedis state size writting.'

            error = client%put_tensor("action_size", action_global_size_tensor, shape(action_global_size_tensor))
            is_error = client%SR_error_parser(error)
            if (error /= 0) stop 'Error during SmartRedis state size writting.'
         ! end if
      end if
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
      character(len=*), intent(in) :: key ! state name to write to database
      integer :: error
      logical :: is_error
      real(rp) :: state_global(state_global_size)

      ! gather the local states into a global state
      call mpi_gatherv( &
         state_local, state_local_size, mpi_datatype_real, &           ! everyone sends state_local data
         state_global, state_sizes, state_displs, mpi_datatype_real, & ! root receives it into state_global
         0, app_comm, error &                                          ! rank 0 is root
      )

      ! write global state into DB
      if (mpi_rank .eq. 0) then
         error = client%put_tensor(key, state_global, shape(state_global))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error during SmartRedis write_state.'
      end if
   end subroutine write_state

   ! Read actions from DB
   subroutine read_action(client, key)
      type(client_type), intent(inout) :: client
      character(len=*), intent(in) :: key ! actions name to read from database

      integer, parameter :: interval = 100 ! polling interval in milliseconds
      integer, parameter :: tries = 1000 ! huge(1) ! infinite number of polling tries
      ! logical(1) :: exists ! receives whether the tensor exists
      logical :: exists ! receives whether the tensor exists
      logical :: is_error
      integer :: found, error

      ! wait (poll) until the actions array is found in the DB, then read, then delete
      if (mpi_rank .eq. 0) then
         found = client%poll_tensor(trim(adjustl(key)), interval, tries, exists) ! wait indefinitely for new actions to appear
         is_error = client%SR_error_parser(found)
         if (found /= 0) stop 'Error in SmartRedis read_action. Actions array not found.'
         error = client%unpack_tensor(trim(adjustl(key)), action_global, shape(action_global))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis read_action. Tensor could not be unpacked.'
         error = client%delete_tensor(trim(adjustl(key)))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis read_action. Tensor could not be deleted.'
      end if

      ! broadcast rank 0 global action array to all processes
      call mpi_bcast(action_global, action_global_size, mpi_datatype_real, 0, app_comm, error)
      !$acc update device(action_global(:))
   end subroutine read_action

   ! Writes the reward values: wall shear stress integral for both positive and negative values
   subroutine write_reward(client, Ftau_neg, key)
      type(client_type), intent(inout) :: client
      real(rp), intent(in) :: Ftau_neg(n_pseudo_envs)
      character(len=*), intent(in) :: key

      integer :: error
      logical :: is_error

      if (mpi_rank .eq. 0) then
         error = client%put_tensor(key, Ftau_neg, shape(Ftau_neg))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis write_reward.'
      end if
   end subroutine write_reward

   ! Indicate environment time step status -> 1: init time step. 2: mid time step. 0: end time step
   subroutine write_step_type(client, step_type, key)
      type(client_type), intent(inout) :: client
      integer, intent(in) :: step_type
      character(len=*), intent(in) :: key

      integer :: error
      logical :: is_error

      if (mpi_rank .eq. 0) then
         error = client%put_tensor(key, [step_type], shape([step_type]))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis write_step_type.'
      end if
      step_type_mod = step_type
   end subroutine write_step_type

   ! Write flow time
   subroutine write_time(client, time, key)
      type(client_type), intent(inout) :: client
      real(rp), intent(in) :: time
      character(len=*), intent(in) :: key

      integer :: error
      logical :: is_error
      real(rp) :: time_tensor(1)

      if (mpi_rank .eq. 0) then
         time_tensor(1) = time
         error = client%put_tensor(key, time_tensor, shape(time_tensor))
         is_error = client%SR_error_parser(error)
         if (error /= 0) stop 'Error in SmartRedis write_time.'
      end if
   end subroutine write_time

#endif
end module mod_smartredis