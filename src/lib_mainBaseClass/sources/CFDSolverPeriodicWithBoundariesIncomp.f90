module CFDSolverPeriodicWithBoundariesIncomp_mod
   use mod_arrays
   use mod_nvtx
#ifndef NOACC
   use cudafor
#endif
   

   use elem_qua
   use elem_hex
   use jacobian_oper
   use quadrature_rules
   use mod_inicond_reader
   use mass_matrix
   use mod_geom
   use time_integ_incomp
   use mod_analysis
   use mod_numerical_params
   use mod_time_ops
   use mod_fluid_viscosity
   use mod_postpro
   use mod_aver
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use CFDSolverPeriodicIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicIncomp) :: CFDSolverPeriodicWithBoundariesIncomp

   contains
      procedure, public :: callTimeIntegration              => CFDSolverPeriodicWithBoundariesIncomp_callTimeIntegration
   end type CFDSolverPeriodicWithBoundariesIncomp
contains
  
   subroutine CFDSolverPeriodicWithBoundariesIncomp_callTimeIntegration(this,istep)
      class(CFDSolverPeriodicWithBoundariesIncomp), intent(inout) :: this
      integer(4) , intent(in)    :: istep

      this%noBoundaries = .false.
      call ab_main_incomp(istep,this%save_logFile_next,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbnNodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
         1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor,mue_l, &
         ndofRankPar,numBoundaryNodesRankPar,ldofPar,lbnodesPar,boundPar,bouCodesPar,bouCodesNodesPar, & ! Optional args
         listBoundsWallModel,wgp_b,boundNormalPar,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)       ! Optional args

   end subroutine CFDSolverPeriodicWithBoundariesIncomp_callTimeIntegration

end module CFDSolverPeriodicWithBoundariesIncomp_mod
