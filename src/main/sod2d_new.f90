


! main.f90
program main
   use mod_constants
   use ThermalChannelFlowSolver_mod

   implicit none
   type(ThermalChannelFlowSolver)  :: channel

   call channel%run()

end program main

