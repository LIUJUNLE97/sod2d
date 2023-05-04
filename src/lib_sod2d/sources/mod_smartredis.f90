module mod_smartredis
#ifdef SMARTREDIS
   use smartredis_client, only: client_type
   use mod_numerical_params, only: do_smartredis, do_smartredis
   implicit none

   type(client_type) :: client ! Client instance of SmartRedis to communicate with Redis database

contains

   ! Initialise SmartRedis client
   ! Since each MPI rank communicates with the database, all of them need to initialise the client
   subroutine init_smartredis(client)
      type(client_type), intent(inout) :: client
      integer :: return_code

      if(do_smartredis) then
         return_code = client%initialize(do_smartredis)
      endif
   end subroutine init_smartredis

   ! Destroy SmartRedis client
   subroutine end_smartredis(client)
      type(client_type), intent(inout) :: client
      integer :: return_code

      if(do_smartredis) then
         return_code = client%destructor()
      endif
   end subroutine end_smartredis
#endif

end module mod_smartredis