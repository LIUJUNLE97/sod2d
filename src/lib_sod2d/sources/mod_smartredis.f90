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

   contains

   ! Initialise SmartRedis client
   ! Since each MPI rank communicates with the database, all of them need to initialise the client
   subroutine init_smartredis(client)
      type(client_type), intent(inout) :: client
      integer :: return_code

      return_code = client%initialize(do_smartredis)
      if (return_code /= SRNoError) stop 'Error in SmartRedis client initialization'
   end subroutine init_smartredis

   ! Destroy SmartRedis client
   subroutine end_smartredis(client)
      type(client_type), intent(inout) :: client
      integer :: return_code

      return_code = client%destructor()
      if (return_code /= SRNoError) stop 'Error in SmartRedis client destruction'
   end subroutine end_smartredis

   ! Write witness points state into DB
   subroutine write_witness_state(client, state_size, state, key_name)
      type(client_type), intent(inout) :: client
      integer(4), intent(in) :: state_size         ! witness points size
      real(rp), intent(in) :: state(state_size)    ! witness points state values
      character(len=*), intent(in) :: key_name     ! array name to write to database
      character(len=8) :: key_mpi_prefix           ! mpi rank prefix
      character(len=32) :: key                     ! name + mpi rank prefix
      integer :: return_code

      write(key_mpi_prefix,'(I6.6,A)') mpi_rank, '_'
      key = key_mpi_prefix // trim(key_name)

      return_code = client%put_tensor(key, state, shape(state))
      if (return_code /= SRNoError) stop 'Error during SmartRedis state writting.'
   end subroutine write_witness_state

   ! Read actions from DB
   subroutine read_actions(client, actions_size, actions, key_name)
      type(client_type), intent(inout) :: client
      integer(4), intent(inout) :: actions_size        ! actions array buffer
      real(rp), intent(inout) :: actions(actions_size) ! actions array buffer
      character(len=*), intent(in) :: key_name         ! actions name to read from database
      integer, parameter :: interval = 10              ! polling interval in milliseconds
      integer, parameter :: tries = huge(1)            ! infinite number of polling tries
      character(len=8) :: key_mpi_prefix               ! mpi rank prefix
      character(len=32) :: key                         ! name + mpi rank prefix
      logical(kind=c_bool) :: exists                   ! receives whether the tensor exists
      logical :: found
      integer :: return_code

      write(key_mpi_prefix,'(I6.6,A)') mpi_rank, '_'
      key = trim(key_mpi_prefix) // trim(key_name)

      found = client%poll_tensor(key, interval, tries, exists)
      if (.not. found) stop 'Error in SmartRedis actions reading. Actions array not found.'
      return_code = client%unpack_tensor(key, actions, shape(actions))
      if (return_code /= SRNoError) stop 'Error during SmartRedis actions reading.'
      return_code = client%delete_tensor(key)
      if (return_code /= SRNoError) stop 'Error in SmartRedis actions reading. Tensor could not be deleted.'
   end subroutine read_actions

#endif
end module mod_smartredis