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
   use mod_period
   use time_integ
   use mod_analysis
   use mod_numerical_params
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
   end type CFDSolverPeriodic
contains

   subroutine CFDSolverPeriodic_callTimeIntegration(this,istep)
      class(CFDSolverPeriodic), intent(inout) :: this
      integer(4)              , intent(in)   :: istep

      this%noBoundaries = .true.
      if(flag_implicit == 1) then
         if(implicit_solver == implicit_solver_bdf2_rk10) then
            call rk_implicit_bdf2_rk10_main(istep,this%save_logFile_next,this%currentNonLinealIter,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbn,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
               1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
               rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor)
         endif
      else
         call rk_4_main(this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbn,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor)
      end if

   end subroutine CFDSolverPeriodic_callTimeIntegration

end module CFDSolverPeriodic_mod
