module ChannelFlowSolver_mod
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
   use CFDSolverPeriodicWithBoundaries_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: ChannelFlowSolver

      real(rp) , public  :: vo, M, delta, U0, rho0, Retau, Re, utau, to, po, mu

   contains
      procedure, public :: fillBCTypes           => ChannelFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters  => ChannelFlowSolver_initializeParameters
      procedure, public :: initializeSourceTerms => ChannelFlowSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => ChannelFlowSolver_evalInitialConditions
   end type ChannelFlowSolver
contains

   subroutine ChannelFlowSolver_fill_BC_Types(this)
      class(ChannelFlowSolver), intent(inout) :: this

      !bouCodes2BCType(1) = bc_type_non_slip_adiabatic
      bouCodes2BCType(1) = bc_type_slip_wall_model
      bouCodes2WallModel(1) = 1

   end subroutine ChannelFlowSolver_fill_BC_Types

   subroutine ChannelFlowSolver_initializeSourceTerms(this)
      class(ChannelFlowSolver), intent(inout) :: this

        allocate(source_term(ndime))
        source_term(1) = (this%utau*this%utau*this%rho0/this%delta)
        source_term(2) = 0.00_rp
        source_term(3) = 0.00_rp

   end subroutine ChannelFlowSolver_initializeSourceTerms

   subroutine ChannelFlowSolver_initializeParameters(this)
      class(ChannelFlowSolver), intent(inout) :: this
      real(rp) :: mur

      write(this%gmsh_file_path,*) "./mesh/"
      write(this%gmsh_file_name,*) "channel"

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "channel"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      this%isPeriodic = .true.
      this%loadMesh = .true.

!      this%loadResults = .true.
!      this%continue_oldLogs = .false.
!      this%load_step = 400001

      this%nstep = 9000000 
      this%cfl_conv = 1.5_rp
      this%cfl_diff = 1.5_rp
      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1
      this%nleap = 50000 ! Saving interval, TODO: input
      this%nleap2 = 50  ! Saving interval, TODO: input
      this%nleapAVG = 50000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 1.0_rp
      this%M  = 0.2_rp
      this%delta  = 1.0_rp
      this%U0     = 1.0_rp
      this%rho0   = 1.0_rp
      this%Retau  = 950.0_rp
      this%gamma_gas = 1.40_rp

      this%Re     = exp((1.0_rp/0.88_rp)*log(this%Retau/0.09_rp))
      this%mu    = (this%rho0*2.0_rp*this%delta*this%vo)/this%Re
      this%utau   = (this%Retau*this%mu)/(this%delta*this%rho0)
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = this%mu/mur
      write(1,*) " Gp ", this%utau*this%utau*this%rho0/this%delta
      nscbc_rho_inf = this%rho0
      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas
      nscbc_T_C = this%to

   end subroutine ChannelFlowSolver_initializeParameters

   subroutine ChannelFlowSolver_evalInitialConditions(this)
      class(ChannelFlowSolver), intent(inout) :: this
      integer(4) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, ti(3), yp
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      logical :: readFiles

      readFiles = .true.

      this%interpInitialResults = .true.

      if(readFiles) then
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         call read_veloc_from_file_Par(numNodesRankPar,totalNumNodesSrl,this%gmsh_file_path,u(:,:,2),matGidSrlOrdered)

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = this%po
            rho(iNodeL,2) = this%po/this%Rgas/this%to
            e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
            Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
            E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
            q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
            csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
            eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         end do
         !$acc end parallel loop
      else
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         auxCnt = 1
         !!$acc parallel loop
         do iLine = 1,totalNumNodesSrl
            call random_number(ti)
            if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
               iNodeL = matGidSrlOrdered(auxCnt,1)
               auxCnt=auxCnt+1
               if(coordPar(iNodeL,2)<this%delta) then
                  yp = coordPar(iNodeL,2)*this%utau*this%rho0/this%mu
               else
                  yp = abs(coordPar(iNodeL,2)-2.0_rp*this%delta)*this%utau*this%rho0/this%mu
               end if

               velo = this%utau*((1.0_rp/0.41_rp)*log(1.0_rp+0.41_rp*yp)+7.8_rp*(1.0_rp-exp(-yp/11.0_rp)-(yp/11.0_rp)*exp(-yp/3.0_rp))) 

               u(iNodeL,1,2) = velo*(1.0_rp + 0.1_rp*(ti(1) -0.5_rp))
               u(iNodeL,2,2) = velo*(0.1_rp*(ti(2) -0.5_rp))
               u(iNodeL,3,2) = velo*(0.1_rp*(ti(3) -0.5_rp))
            end if
         end do
         !!$acc end parallel loop
         !!$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = this%po
            rho(iNodeL,2) = this%po/this%Rgas/this%to
            e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
            Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
            E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
            q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
            csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
            eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         end do
         !!$acc end parallel loop
      end if

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)
      end do
      !$acc end parallel loop

      !$acc kernels
      mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels
      call nvtxEndRange

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop
   end subroutine ChannelFlowSolver_evalInitialConditions

end module ChannelFlowSolver_mod
