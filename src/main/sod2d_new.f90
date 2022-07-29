! main.f90
program main
   use mod_constants
   !use TGVSolver_mod
   use ThermalChannelFlowSolver_mod
   !use BluffBodySolver_mod

   implicit none
   !type(BluffBodySolver)  :: bluff
   type(ThermalChannelFlowSolver)  :: channel
   !type(TGVSolver)  :: tgv

   !call bluff%run()
   call channel%run()
   !call tgv%run()

end program main

