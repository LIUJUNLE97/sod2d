module CFDSolver3DWithBoundariesIncomp_mod
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
   use CFDSolverPeriodicIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicIncomp) :: CFDSolver3DWithBoundariesIncomp

   contains
      procedure, public :: fillBCTypes             =>CFDSolver3DWithBoundariesIncomp_fill_BC_Types
      procedure, public :: callTimeIntegration     =>CFDSolver3DWithBoundariesIncomp_callTimeIntegration
   end type CFDSolver3DWithBoundariesIncomp
contains

   subroutine CFDSolver3DWithBoundariesIncomp_fill_BC_Types(this)
      class(CFDSolver3DWithBoundariesIncomp), intent(inout) :: this
   
      if(mpi_rank.eq.0) write(111,*) "--| Boundary types must be defined "
      stop 1
   end subroutine CFDSolver3DWithBoundariesIncomp_fill_BC_Types

   subroutine CFDSolver3DWithBoundariesIncomp_callTimeIntegration(this,istep)
      class(CFDSolver3DWithBoundariesIncomp), intent(inout) :: this
      integer(4)        , intent(in)   :: istep


      this%noBoundaries = .false.
      call ab_main_incomp(istep,this%save_logFile_next,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbn,lnbnNodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor, &
            mue_l,convertIJK,al_weights,am_weights,an_weights, &
            ndofRankPar,numBoundaryNodesRankPar,ldofPar,lbnodesPar,boundPar,bouCodesPar,bouCodesNodesPar, & ! Optional args
            listBoundsWallModel,wgp_b,boundNormalPar,normalsAtNodes,u_buffer,tauw,walave_u)                   ! Optional args
   end subroutine CFDSolver3DWithBoundariesIncomp_callTimeIntegration

end module CFDSolver3DWithBoundariesIncomp_mod
