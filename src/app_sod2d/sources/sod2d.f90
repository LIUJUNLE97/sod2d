! main.f90

#define _tgv_ 0
#define _channel_ 0
#define _bluff_ 0
#define _bluff3d_ 0
#define _bl_ 0
#define _tsb_ 1

program main
   use mod_numerical_params
#if _tgv_
   use TGVSolver_mod
#endif
#if _channel_
   use ChannelFlowSolver_mod
   !use ThermalChannelFlowSolver_mod
#endif
#if _bluff_
   use BluffBodySolver_mod
   use DLRSolver_mod
#endif
#if _bluff3d_
   use BluffBody3DSolver_mod
#endif
#if _bl_
   !use BLFlowSolver_mod
   use BLAPGFlowSolver_mod
#endif
#if _tsb_
   use BLTSBFlowSolver_mod
   use BLTSBDRLFlowSolver_mod
#endif
   implicit none

#if _tgv_
   type(TGVSolver)  :: tgv
   call tgv%run()
#endif
#if _channel_
   !type(ThermalChannelFlowSolver) :: channel
   type(ChannelFlowSolver)  :: channel
   call channel%run()
#endif
#if _bluff_
   !type(BluffBodySolver)  :: bluff
   type(DLRSolver)  :: bluff
   call bluff%run()
#endif
#if _bluff3d_
   type(BluffBody3DSolver)  :: bluff3d
   call bluff3d%run()
#endif

#if _bl_
   !type(BLFlowSolver)    :: blflow
   type(BLAPGFlowSolver) :: blflow
   call blflow%run()
#endif

#if _tsb_
   ! type(BLTSBFlowSolver)  :: bltsbflow
   type(BLTSBDRLFlowSolver)  :: bltsbdrlflow
   ! call bltsbflow%run()
   call bltsbdrlflow%run()
#endif




end program main

