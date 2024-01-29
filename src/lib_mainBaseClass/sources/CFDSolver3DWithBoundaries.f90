module CFDSolver3DWithBoundaries_mod
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
   use time_integ
   use time_integ_imex
   use mod_analysis
   use mod_numerical_params
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
   end type CFDSolver3DWithBoundaries
contains

   subroutine CFDSolver3DWithBoundaries_fill_BC_Types(this)
      class(CFDSolver3DWithBoundaries), intent(inout) :: this
   
      if(mpi_rank.eq.0) write(111,*) "--| Boundary types must be defined "
      stop 1
   end subroutine CFDSolver3DWithBoundaries_fill_BC_Types

   subroutine CFDSolver3DWithBoundaries_callTimeIntegration(this,istep)
      class(CFDSolver3DWithBoundaries), intent(inout) :: this
      integer(4)        , intent(in)   :: istep


      this%noBoundaries = .false.
      if(flag_implicit == 1) then        
         if (implicit_solver == implicit_solver_imex) then 
            call imex_main(istep,this%save_logFile_next,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbnNodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                        1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
                        rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor,mue_l, &
                        ndofRankPar,numBoundaryNodesRankPar,ldofPar,lbnodesPar,boundPar,bouCodesPar,bouCodesNodesPar, & ! Optional args
                        listBoundsWallModel,wgp_b,boundNormalPar,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)         
         endif
      else
         call rk_4_main(this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,mu_fluid,mu_factor,mue_l, &
            ndofRankPar,numBoundaryNodesRankPar,ldofPar,lbnodesPar,boundPar,bouCodesPar,bouCodesNodesPar, & ! Optional args
            listBoundsWallModel,wgp_b,boundNormalPar,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)                   ! Optional args
      end if

   end subroutine CFDSolver3DWithBoundaries_callTimeIntegration

end module CFDSolver3DWithBoundaries_mod
