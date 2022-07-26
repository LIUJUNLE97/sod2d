! main.f90
program main
   use mod_constants
   !use ThermalChannelFlowSolver_mod
   use BluffBodySolver_mod

   implicit none
   type(BluffBodySolver)  :: bluff
   !type(ThermalChannelFlowSolver)  :: channel

   call bluff%run()
   !call channel%run()

end program main

