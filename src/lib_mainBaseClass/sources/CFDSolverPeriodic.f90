module CFDSolverPeriodic_mod
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
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use CFDSolverBase_mod
   implicit none
   private

   type, public, extends(CFDSolverBase) :: CFDSolverPeriodic

   contains
      procedure, public :: callTimeIntegration     =>CFDSolverPeriodic_callTimeIntegration
      procedure, public :: saveAverages            =>CFDSolverPeriodic_saveAverages
      procedure, public :: savePosprocessingFields =>CFDSolverPeriodic_savePosprocessingFields
   end type CFDSolverPeriodic
contains

   subroutine CFDSolverPeriodic_callTimeIntegration(this)
      class(CFDSolverPeriodic), intent(inout) :: this

      this%noBoundaries = .true.
      call rk_4_main(this%noBoundaries,this%isWallModelOn,0,0,numElemsInRank,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbn,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
         1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor)

   end subroutine CFDSolverPeriodic_callTimeIntegration

   subroutine CFDSolverPeriodic_saveAverages(this,istep)
      class(CFDSolverPeriodic), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      call eval_average_window(isMeshPeriodic,numNodesRankPar,numElemsInRank,acuvel,acuve2,acuvex,acurho,acupre,acumueff,this%acutim,&
											avvel,avve2,avvex,avrho,avpre,avmueff,nPerRankPar,masSlaRankPar)

      if(save_vtk) then
         call write_vtkAVG_binary(istep,numNodesRankPar,numElemsInRank,coordPar,connecVTK,avvel,avve2,avrho,avpre,avmueff)
      end if

      if(save_hdf5) then
         call save_hdf5_avgResultsFile(istep,avvel,avve2,avvex,avrho,avpre,avmueff)
      end if

   end subroutine CFDSolverPeriodic_saveAverages

   subroutine CFDSolverPeriodic_savePosprocessingFields(this,istep)
      class(CFDSolverPeriodic), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      if(save_vtk) then
         call write_vtk_binary(isMeshPeriodic,istep,numNodesRankPar,numElemsInRank,coordPar,connecVTK, &
            rho(:,2),u(:,:,2),pr(:,2),E(:,2),csound,machno, &
            gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs,nPerRankPar,masSlaRankPar)
      end if

      if(save_hdf5) then
         call save_hdf5_resultsFile(istep,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),eta(:,2),csound,machno,gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs)
      end if

   end subroutine CFDSolverPeriodic_savePosprocessingFields

end module CFDSolverPeriodic_mod
