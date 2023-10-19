! main.f90

#define _tgv_ 0
#define _tgv_multi 0
#define _tgv_comp 0
#define _tgv_incomp 0
#define _channel_ 0
#define _channel_incomp 0
#define _bluff_ 0
#define _bluff_incomp 0
#define _bluff3d_ 0
#define _bluff3d_incomp 0
#define _bl_ 0
#define _bl_incomp 1
#define _tsb_ 0

program main
   use mod_numerical_params
#if _tgv_
   use TGVSolver_mod
#endif
#if _tgv_multi
   use TGVMultiSolver_mod
#endif
#if _tgv_comp
   use TGVCompSolver_mod
#endif
#if _tgv_incomp
   use TGVSolverIncomp_mod
#endif
#if _channel_
   use ChannelFlowSolver_mod
   !use ThermalChannelFlowSolver_mod
#endif
#if _channel_incomp
   use ChannelFlowSolverIncomp_mod
#endif
#if _bluff_
   use BluffBodySolver_mod
   use DLRSolver_mod
#endif
#if _bluff_incomp
   use BluffBodySolverIncomp_mod
#endif
#if _bluff3d_
   use BluffBody3DSolver_mod
#endif
#if _bluff3d_incomp
   use BluffBody3DSolverIncomp_mod
#endif
#if _bl_
   ! use BLFlowSolver_mod
   use BLFlowSolverTest_mod
   ! use BLAPGFlowSolver_mod
#endif
#if _bl_incomp
   ! use BLFlowSolverIncomp_mod
   ! use BLFlowSolverIncompTest_mod
   use BLFlowSolverIncompTest2_mod
   ! use BLMARLFlowSolverIncomp_mod
#endif
#if _tsb_
   ! use BLTSBFlowSolver_mod
   use BLTSBDRLFlowSolver_mod
#endif
   implicit none

#if _tgv_
   type(TGVSolver)  :: tgv
   call tgv%run()
#endif
#if _tgv_comp
   type(TGVCompSolver)  :: tgv
   call tgv%run()
#endif
#if _tgv_incomp
   type(TGVSolverIncomp)  :: tgv
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
#if _channel_incomp
   type(ChannelFlowSolverIncomp)  :: channel
   call channel%run()
#endif
#if _bluff_
   !type(BluffBodySolver)  :: bluff
   type(DLRSolver)  :: bluff
   call bluff%run()
#endif
#if _bluff_incomp
   type(BluffBodySolverIncomp)  :: bluff
   call bluff%run()
#endif
#if _bluff3d_
   type(BluffBody3DSolver)  :: bluff3d
   call bluff3d%run()
#endif
#if _bluff3d_incomp
   type(BluffBody3DSolverIncomp)  :: bluff3d
   call bluff3d%run()
#endif
#if _bl_
   ! type(BLFlowSolver)    :: blflow
   type(BLFlowSolverTest)    :: blflow
   ! type(BLAPGFlowSolver) :: blflow
   call blflow%run()
#endif
#if _bl_incomp
   !  type(BLFlowSolverIncomp)    :: blflow
   !  type(BLFlowSolverIncompTest)    :: blflow
    type(BLFlowSolverIncompTest2)    :: blflow
   !  type(BLMARLFlowSolverIncomp)    :: blflow
   call blflow%run()
#endif
#if _tsb_
   !  type(BLTSBFlowSolver)  :: bltsbflow
   type(BLTSBDRLFlowSolver)  :: bltsbdrlflow
   !  call bltsbflow%run()
   call bltsbdrlflow%run()
#endif




end program main

