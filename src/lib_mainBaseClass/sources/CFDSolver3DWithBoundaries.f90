module CFDSolver3DWithBoundaries_mod
   use mod_arrays
   use mod_nvtx
#ifndef NOACC
   use cudafor
#endif   
   use mod_veclen

   use elem_qua
   use elem_hex
   use jacobian_oper
   use quadrature_rules
   use mesh_reader
   use inicond_reader
   use mass_matrix
   use mod_geom
   use mod_output
   use mod_period
   use time_integ
   use mod_analysis
   use mod_constants
   use mod_time_ops
   use mod_fluid_viscosity
   use mod_postpro
   use mod_aver
   use CFDSolverBase_mod
   implicit none
   private

   type, public, extends(CFDSolverBase) :: CFDSolver3DWithBoundaries

   contains
      procedure, public :: fillBCTypes             =>CFDSolver3DWithBoundaries_fill_BC_Types
      procedure, public :: callTimeIntegration     =>CFDSolver3DWithBoundaries_callTimeIntegration
      procedure, public :: saveAverages            =>CFDSolver3DWithBoundaries_saveAverages
      procedure, public :: savePosprocessingFields =>CFDSolver3DWithBoundaries_savePosprocessingFields
   end type CFDSolver3DWithBoundaries
contains

   subroutine CFDSolver3DWithBoundaries_fill_BC_Types(this)
      class(CFDSolver3DWithBoundaries), intent(inout) :: this
   
      if(mpi_rank.eq.0) write(111,*) "--| Boundary types must be defined "
      STOP(1)
   end subroutine CFDSolver3DWithBoundaries_fill_BC_Types

   subroutine CFDSolver3DWithBoundaries_callTimeIntegration(this)
      class(CFDSolver3DWithBoundaries), intent(inout) :: this

      call rk_4_main(0,0,numElemsInRank,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,point2elem,lnbn,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
         1,connecParWork,Ngp,dNgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
         rho,u,u_wall,rho_wall,mu_wall,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor, &
         ndofRankPar,numBoundaryNodesRankPar,ldofPar,lbnodesPar,boundPar,bouCodesPar,bouCodesNodesPar,numBoundCodes,wgp_b,boundNormalPar,bouCodes2WallModel,coordPar,normalsAtNodes) ! Optional args

   end subroutine CFDSolver3DWithBoundaries_callTimeIntegration

   subroutine CFDSolver3DWithBoundaries_savePosprocessingFields(this,istep)
      class(CFDSolver3DWithBoundaries), intent(inout) :: this
      integer(4)        , intent(in)   :: istep

      if(save_hdf5) then
         call save_hdf5_resultsFile(istep,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),eta(:,2),csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs)
      end if

   end subroutine CFDSolver3DWithBoundaries_savePosprocessingFields

   subroutine CFDSolver3DWithBoundaries_saveAverages(this,istep)
      class(CFDSolver3DWithBoundaries), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      call eval_average_window(isMeshPeriodic,numNodesRankPar,numElemsInRank,acuvel,acuve2,acurho,acupre,acumueff,this%acutim,&
											avvel,avve2,avrho,avpre,avmueff,nPerRankPar,masSlaRankPar)

      if(save_hdf5) then
         call save_hdf5_avgResultsFile(istep,avvel,avve2,avrho,avpre,avmueff)
      end if

   end subroutine CFDSolver3DWithBoundaries_saveAverages
end module CFDSolver3DWithBoundaries_mod
