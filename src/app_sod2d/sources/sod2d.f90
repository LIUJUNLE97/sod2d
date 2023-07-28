! main.f90

#define _tgv_ 1
#define _tgv_multi 0
#define _channel_ 0
#define _bluff_ 0
#define _bluff3d_ 0
#define _bl_ 0

program main
   use mod_numerical_params
#if _tgv_
   use TGVSolver_mod
#endif
#if _tgv_multi
   use TGVMultiSolver_mod
#endif
#if _channel_
   use ChannelFlowSolver_mod
   !use ThermalChannelFlowSolver_mod
#endif
#if _bluff_
   use BluffBodySolver_mod
#endif
#if _bluff3d_
   use BluffBody3DSolver_mod
#endif
#if _bl_
   use BLFlowSolver_mod
#endif
   implicit none

#if _tgv_
   type(TGVSolver)  :: tgv
   call tgv%run()
#endif
#if _tgv_multi
   type(TGVMultiSolver)  :: tgv_multi
   call tgv_multi%run()
#endif
#if _channel_
   !type(ThermalChannelFlowSolver) :: channel
   type(ChannelFlowSolver)  :: channel
   call channel%run()
#endif
#if _bluff_
   type(BluffBodySolver)  :: bluff
   call bluff%run()
#endif
#if _bluff3d_
   type(BluffBody3DSolver)  :: bluff3d
   call bluff3d%run()
#endif

#if _bl_
   type(BLFlowSolver)  :: blflow
   call blflow%run()
#endif




end program main

