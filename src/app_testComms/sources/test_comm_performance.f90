program test_comm_perfomance
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms
   use mod_comms_peformance
   implicit none

   character(10) :: cnsteps
   integer :: numNodesSrl,numNodesB_1r,numIters
   integer :: iter
   
   call init_mpi()

   !---------------------------------------------------------------------
   !Num nodes in Serial (total number of nodes)
   numNodesSrl = 100
   if (command_argument_count() > 0) then
     call get_command_argument(1, cnsteps)
     read(cnsteps,'(i)') numNodesSrl
   endif
   if (numNodesSrl <= 0) numNodesSrl = 100
   !----------------------------------------------
   !Num nodes in boundary in one side, (total num in each rank is x2)
   numNodesB_1r = 1 
   if (command_argument_count() > 1) then
     call get_command_argument(2, cnsteps)
     read(cnsteps,'(i)') numNodesB_1r
   endif
   if (numNodesB_1r <= 0) numNodesB_1r = 1
   !----------------------------------------------
   !Number of iterations
   numIters = 1
   if (command_argument_count() > 2) then
     call get_command_argument(3, cnsteps)
     read(cnsteps,'(i)') numIters
   endif
   if (numIters <= 0) numIters = 1
   !----------------------------------------------

   if(mpi_rank.eq.0) write(*,*) 'numNodesSrl ',numNodesSrl,' numNodesB_1r ',numNodesB_1r,' numIters ',numIters
   !---------------------------------------------------------------------

   call init_comms_peformance(numNodesSrl,numNodesB_1r)
        
   !call nvtxStartRange("saxpy loop")
   !call do_saxpy_loop(numIters)
   !call nvtxEndRange

   !call test_mpi_cudaware(numNodesSrl,numIters)

   call do_crazy_mpi_test(numNodesSrl,numIters)

   call end_comms()
   call end_mpi()

end program test_comm_perfomance
