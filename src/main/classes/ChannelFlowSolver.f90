module ChannelFlowSolver_mod
   use mod_arrays
   use mod_nvtx
   use cudafor
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
   use CFDSolverPeriodicWithBoundaries_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: ChannelFlowSolver

      real(rp) , public  :: vo, M, delta, U0, rho0, Retau, Re, utau, to, po

   contains
      procedure, public :: initializeParameters  => ChannelFlowSolver_initializeParameters
      procedure, public :: initializeSourceTerms => ChannelFlowSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => ChannelFlowSolver_evalInitialConditions
   end type ChannelFlowSolver
contains

   subroutine ChannelFlowSolver_initializeSourceTerms(this)
      class(ChannelFlowSolver), intent(inout) :: this

        allocate(source_term(ndime))
        source_term(1) = (this%utau*this%utau*this%rho0/this%delta)
        source_term(2) = 0.00_rp
        source_term(3) = 0.00_rp

   end subroutine ChannelFlowSolver_initializeSourceTerms

   subroutine ChannelFlowSolver_initializeParameters(this)
      class(ChannelFlowSolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%file_name,*) "channel" 

      this%nstep = 90000000 
      this%cfl_conv = 2.2_rp
      this%cfl_diff = 2.2_rp
      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1
      this%nleap = 20000 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 10  ! Saving interval, TODO: input
      this%nleapAVG = 20000
      this%isPeriodic = 1 ! TODO: make it a read parameter (0 if not periodic, 1 if periodic)
      this%nper = 32131 ! TODO: if periodic, request number of periodic nodes

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
      mul    = (this%rho0*2.0_rp*this%delta*this%vo)/this%Re
      this%utau   = (this%Retau*mul)/(this%delta*this%rho0)
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = mul/mur
      write(1,*) " Gp ", this%utau*this%utau*this%rho0/this%delta
      nscbc_p_inf = this%po

   end subroutine ChannelFlowSolver_initializeParameters

   subroutine ChannelFlowSolver_evalInitialConditions(this)
      class(ChannelFlowSolver), intent(inout) :: this
      integer(rp) :: ipoin

      call read_veloc(this%npoin,this%file_path,u(:,:,2))
      !$acc parallel loop
      do ipoin = 1,this%npoin
         pr(ipoin,2) = this%po
         rho(ipoin,2) = this%po/this%Rgas/this%to
         e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(this%gamma_gas-1.0_rp))
         Tem(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*this%Rgas)
         E(ipoin,2) = rho(ipoin,2)*(0.5_rp*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
         q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
         csound(ipoin) = sqrt(this%gamma_gas*pr(ipoin,2)/rho(ipoin,2))
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do ipoin = 1,this%npoin
         machno(ipoin) = dot_product(u(ipoin,:,2),u(ipoin,:,2))/csound(ipoin)
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
      do ipoin = 1,this%npoin
         mu_factor(ipoin) = flag_mu_factor
      end do
      !$acc end parallel loop
   end subroutine ChannelFlowSolver_evalInitialConditions

end module ChannelFlowSolver_mod
