! main.f90

program main

   use mod_numerical_params
   use mod_arrays
   use json_module
   use CFDSolverBase_mod
   use TGVSolver_mod
   use TGVCompSolver_mod
   use TGVSolverIncomp_mod
   use ChannelFlowSolver_mod
   use ChannelFlowSolverIncomp_mod
   use BluffBodySolver_mod
   use BluffBodySolverIncomp_mod
   use BluffBody3DSolver_mod
   use BluffBody3DSolverIncomp_mod
   use BLFlowSolver_mod
   use ABlFlowSolverIncomp_mod
   implicit none

   logical :: found
   character(len=:) , allocatable :: value
   character(len=100) :: json_filename
   class(CFDSolverBase), pointer :: solver

   ! Get the name of the JSON file
   call get_command_argument(1, json_filename)

   ! Append the extension
   json_filename = trim(json_filename) // ".json"
   write(*,*), "Reading the JSON file : ", json_filename

   call json%initialize()
   call json%load_file(json_filename)
   if (json%failed()) then 
      write(*,*) " There is a syntax error of the JSON file "
      stop 1
   end if

   call json%get("type", value, found)

   if(value .eq. "TGVSolver") then
      allocate(TGVSolver::solver) 
   else if(value .eq. "TGVCompSolver") then
      allocate(TGVCompSolver::solver) 
   else if(value .eq. "TGVSolverIncomp") then
      allocate(TGVSolverIncomp::solver) 
   else if(value .eq. "ChannelFlowSolver") then
      allocate(ChannelFlowSolver::solver) 
   else if(value .eq. "ChannelFlowSolverIncomp") then
      allocate(ChannelFlowSolverIncomp::solver) 
   else if(value .eq. "BluffBodySolver") then
      allocate(BluffBodySolver::solver) 
   else if(value .eq. "BluffBodySolverIncomp") then
      allocate(BluffBodySolverIncomp::solver) 
   else if(value .eq. "BluffBody3DSolver") then
      allocate(BluffBody3DSolver::solver) 
   else if(value .eq. "BluffBody3DSolverIncomp") then
      allocate(BluffBody3DSolverIncomp::solver) 
   else if(value .eq. "BLFlowSolver") then
      allocate(BLFlowSolver::solver) 
   else if(value .eq. "ABlFlowSolverIncomp") then
      allocate(ABlFlowSolverIncomp::solver) 
   else
      write(*,*) " Solver not implemented in SOD2D : ",value
      stop 1
   end if

   call solver%run()

end program main

