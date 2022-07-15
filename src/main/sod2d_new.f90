
module ChannelFlowSolver_mod
   use mod_constants
   use CFDSolverBase_mod
   implicit none
   private

   type, public, extends(CFDSolverBase) :: ChannelFlowSolver

   contains
      procedure, public :: printDt => ChannelFlowSolver_printDt
   end type ChannelFlowSolver
contains

   subroutine ChannelFlowSolver_printDt(this)
      class(ChannelFlowSolver), intent(inout) :: this

      write(*,*) " Channel Dt ",this%dt
   end subroutine ChannelFlowSolver_printDt
end module ChannelFlowSolver_mod


! main.f90
program main
   use mod_constants
   use CFDSolverBase_mod
   use ChannelFlowSolver_mod

   implicit none
   type(CFDSolverBase)      :: cas
   type(ChannelFlowSolver)  :: channel

   cas%dt = 1e-4
   channel%dt = 1e-3
   call cas%printDt()
   call channel%printDt()


   call cas%printAll()
   call channel%printAll()

end program main

