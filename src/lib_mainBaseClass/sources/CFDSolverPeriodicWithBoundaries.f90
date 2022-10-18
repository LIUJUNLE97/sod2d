module CFDSolverPeriodicWithBoundaries_mod
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
   use CFDSolverPeriodic_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodic) :: CFDSolverPeriodicWithBoundaries

   contains
      procedure, public :: fillBCTypes             =>CFDSolverPeriodicWithBoundaries_fill_BC_Types
      procedure, public :: callTimeIntegration     =>CFDSolverPeriodicWithBoundaries_callTimeIntegration
   end type CFDSolverPeriodicWithBoundaries
contains

   subroutine CFDSolverPeriodicWithBoundaries_fill_BC_Types(this)
      class(CFDSolverPeriodicWithBoundaries), intent(inout) :: this
   
      if(mpi_rank.eq.0) write(111,*) "--| Boundary types must be defined "
      stop 1
   end subroutine CFDSolverPeriodicWithBoundaries_fill_BC_Types

   subroutine CFDSolverPeriodicWithBoundaries_callTimeIntegration(this)
      class(CFDSolverPeriodicWithBoundaries), intent(inout) :: this

      this%noBoundaries = .false.
      call rk_4_main(this%noBoundaries,0,0,numElemsInRank,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,point2elem,lnbn,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
         1,connecParWork,Ngp,dNgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
         rho,u,u_wall,rho_wall,mu_wall,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor, &
         ndofRankPar,numBoundaryNodesRankPar,ldofPar,lbnodesPar,boundPar,bouCodesPar,bouCodesNodesPar,numBoundCodes,wgp_b,boundNormalPar,bouCodes2WallModel,coordPar,normalsAtNodes,source_term) ! Optional args

   end subroutine CFDSolverPeriodicWithBoundaries_callTimeIntegration

end module CFDSolverPeriodicWithBoundaries_mod
