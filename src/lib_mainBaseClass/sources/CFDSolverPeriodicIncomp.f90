module CFDSolverPeriodicIncomp_mod
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
   use time_integ_species_imex
   use time_integ_species
   use mod_analysis
   use mod_numerical_params
   use mod_time_ops
   use mod_fluid_viscosity
   use mod_postpro
   use mod_aver
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_saveFields
   use CFDSolverBase_mod
   implicit none
   private

   type, public, extends(CFDSolverBase) :: CFDSolverPeriodicIncomp

   contains
      procedure, public :: callTimeIntegration              => CFDSolverPeriodicIncomp_callTimeIntegration
      procedure, public :: initializeDefaultSaveFields      => CFDSolverPeriodicIncomp_initializeDefaultSaveFields
      procedure, public :: initializeDefaultCustomSettings  => CFDSolverPeriodicIncomp_initializeDefaultCustomSettings
      procedure, public :: evalInitialDt                    => CFDSolverPeriodicIncomp_evalInitialDt
      procedure, public :: evalDt                           => CFDSolverPeriodicIncomp_evalDt
      procedure, public :: initNSSolver                     => CFDSolverPeriodicIncomp_initNSSolver
      procedure, public :: endNSSolver                      => CFDSolverPeriodicIncomp_endNSSolver
      procedure, public :: eval_vars_after_load_hdf5_resultsFile => CFDSolverPeriodicIncomp_eval_vars_after_load_hdf5_resultsFile 
   end type CFDSolverPeriodicIncomp
contains
   subroutine CFDSolverPeriodicIncomp_eval_vars_after_load_hdf5_resultsFile(this)
      implicit none
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      integer :: iNodeL,idime

      !values loaded -> rho,u,pr,E,mu_e,mu_sgs

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         eta(iNodeL,2) =0.5_rp*dot_product(u(iNodeL,1:ndime,2),u(iNodeL,1:ndime,2))

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         u(iNodeL,1:ndime,3) = u(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         eta(iNodeL,3) =  eta(iNodeL,2)
         pr(iNodeL,3) =  pr(iNodeL,2)  

         q(iNodeL,1:ndime,4) = q(iNodeL,1:ndime,2)
         u(iNodeL,1:ndime,4) = u(iNodeL,1:ndime,2)
         rho(iNodeL,4) = rho(iNodeL,2)
         eta(iNodeL,4) =  eta(iNodeL,2)
         pr(iNodeL,4) =  pr(iNodeL,2)  

         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop


      !$acc kernels
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels

   end subroutine CFDSolverPeriodicIncomp_eval_vars_after_load_hdf5_resultsFile

   subroutine CFDSolverPeriodicIncomp_initNSSolver(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      
      call init_rk4_solver_incomp(numNodesRankPar) 
      if(flag_use_species .eqv. .true.) then
         call init_imex_species_solver(numNodesRankPar,numElemsRankPar)
         !call init_rk4_ls_species_solver(numNodesRankPar,numElemsRankPar)
      end if

   end subroutine CFDSolverPeriodicIncomp_initNSSolver

   subroutine CFDSolverPeriodicIncomp_endNSSolver(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      
       call end_rk4_solver_incomp()

   end subroutine CFDSolverPeriodicIncomp_endNSSolver

   subroutine CFDSolverPeriodicIncomp_evalInitialDt(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this

      !$acc kernels
      csound(:) = 0.0_rp
      !$acc end kernels

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating initial dt..."
      call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
      if(mpi_rank.eq.0) write(111,*) "--| Initial time-step dt := ",this%dt,"s"

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverPeriodicIncomp_evalInitialDt

   subroutine CFDSolverPeriodicIncomp_evalDt(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      
      !$acc kernels
      csound(:) = 0.0_rp
      !$acc end kernels

      call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)

   end subroutine CFDSolverPeriodicIncomp_evalDt

   subroutine CFDSolverPeriodicIncomp_initializeDefaultSaveFields(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this

      call initializeDefaultSaveFieldsIncompressible()

   end subroutine CFDSolverPeriodicIncomp_initializeDefaultSaveFields

   subroutine CFDSolverPeriodicIncomp_initializeDefaultCustomSettings(this)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this

      flag_real_diff=1
      flag_diff_suth=0
      flag_walave = .true.

   end subroutine CFDSolverPeriodicIncomp_initializeDefaultCustomSettings

   subroutine CFDSolverPeriodicIncomp_callTimeIntegration(this,istep)
      class(CFDSolverPeriodicIncomp), intent(inout) :: this
      integer(4)                    , intent(in)    :: istep
      integer(4) :: ispc

      this%noBoundaries = .true.

      if(flag_use_species .eqv. .true.) then
         do ispc = 1,nspecies
            call imex_species_main(ispc,istep,this%local_step,this%save_logFile_next,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,invMl,gpvol,this%dt,helem,helem_l,this%Cp,this%Prt, &
            rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,workingNodesPar,mu_fluid,mue_l)   
            !call rk_4_ls_species_main(ispc,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numBoundsWMRankPar,point2elem,lnbnNodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            !1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,invMl,gpvol,this%dt,helem,helem_l,this%Cp,this%Prt, &
            !rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,workingNodesPar,mu_fluid,mue_l)   
         end do
      end if

      call bdfext3_main_incomp(istep,this%local_step,this%save_logFile_next,this%noBoundaries,this%isWallModelOn,numElemsRankPar,numBoundsRankPar,numNodesRankPar,numWorkingNodesRankPar,numOwnedNodesRankPar,numBoundsWMRankPar,&
            point2elem,lnbnNodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,this%leviCivi,&
            1,connecParWork,Ngp,dNgp,coordPar,wgp,He,Ml,invMl,gpvol,this%dt,helem,helem_l,this%Rgas,this%gamma_gas,this%Cp,this%Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,workingNodesPar,ownedNodesPar,mu_fluid,mu_factor,mue_l)
      

   end subroutine CFDSolverPeriodicIncomp_callTimeIntegration

end module CFDSolverPeriodicIncomp_mod
