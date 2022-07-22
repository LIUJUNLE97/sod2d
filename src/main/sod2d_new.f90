

module ThermalChannelFlowSolver_mod
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
   use CFDSolverBase_mod
   implicit none
   private

   type, public, extends(CFDSolverBase) :: ThermalChannelFlowSolver

      real(rp) , public  ::  delta, rho0, Retau, utau, to, po, mul

   contains
      procedure, public :: initializeParameters  => ThermalChannelFlowSolver_initializeParameters
      procedure, public :: initializeSourceTerms => ThermalChannelFlowSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => ThermalChannelFlowSolver_evalInitialConditions
   end type ThermalChannelFlowSolver
contains

   subroutine ThermalChannelFlowSolver_initializeSourceTerms(this)
      class(ThermalChannelFlowSolver), intent(inout) :: this

        allocate(source_term(ndime))
        source_term(1) = (this%utau*this%utau*this%rho0/this%delta)
        source_term(2) = 0.00_rp
        source_term(3) = 0.00_rp

   end subroutine ThermalChannelFlowSolver_initializeSourceTerms

   subroutine ThermalChannelFlowSolver_initializeParameters(this)
      class(ThermalChannelFlowSolver), intent(inout) :: this

      write(this%file_name,*) "channel_sem" 

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
      this%nper = 58561 ! TODO: if periodic, request number of periodic nodes

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%to = 0.5_rp*(293.0_rp+586.0_rp)
      this%delta  = 0.0015_rp*2.0_rp
      this%rho0   = 0.80396_rp
      this%Retau  = 400.0_rp
      this%gamma_gas = 1.40_rp

      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%po   = this%rho0*this%Rgas*this%to
      this%mul = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      this%utau = (this%Retau*this%mul)/(this%delta*this%rho0)
      flag_mu_factor = 1.0_rp
      write(1,*) " Gp ", this%utau*this%utau*this%rho0/this%delta, " utau ",this%utau
      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas

   end subroutine ThermalChannelFlowSolver_initializeParameters

   subroutine ThermalChannelFlowSolver_evalInitialConditions(this)
      class(ThermalChannelFlowSolver), intent(inout) :: this
      integer(4) :: ipoin, readFiles
      real(rp) :: velo, ti(3), yp

      readFiles = 1

      if(readFiles .eq. 1) then
         call read_veloc(this%npoin,this%file_path,u(:,:,2))
         call read_densi(this%npoin,this%file_path,rho(:,2))
         call read_temper(this%npoin,this%file_path,Tem(:,2))
         !!$acc parallel loop
         do ipoin = 1,this%npoin
            pr(ipoin,2) = Tem(ipoin,2)*this%Rgas*rho(ipoin,2)
            e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(this%gamma_gas-1.0_rp))
            E(ipoin,2) = rho(ipoin,2)*(0.5_rp*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
            q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
            csound(ipoin) = sqrt(this%gamma_gas*pr(ipoin,2)/rho(ipoin,2))
         end do
         !!$acc end parallel loop
      else
         !!$acc parallel loop
         do ipoin = 1,this%npoin
            if(coord(ipoin,2)<this%delta) then
               yp = coord(ipoin,2)*this%utau*this%rho0/this%mul
            else
               yp = abs(coord(ipoin,2)-2.0_rp*this%delta)*this%utau*this%rho0/this%mul
            end if

            velo = this%utau*((1.0_rp/0.41_rp)*log(1.0_rp+0.41_rp*yp)+7.8_rp*(1.0_rp-exp(-yp/11.0_rp)-(yp/11.0_rp)*exp(-yp/3.0_rp))) 

            call random_number(ti)


            u(ipoin,1,2) = velo*(1.0_rp + 0.1_rp*(ti(1) -0.5_rp))
            u(ipoin,2,2) = velo*(0.1_rp*(ti(2) -0.5_rp))
            u(ipoin,3,2) = velo*(0.1_rp*(ti(3) -0.5_rp))
            pr(ipoin,2) = this%po
            rho(ipoin,2) = this%rho0
            e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(this%gamma_gas-1.0_rp))
            Tem(ipoin,2) = this%to
            E(ipoin,2) = rho(ipoin,2)*(0.5_rp*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
            q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
            csound(ipoin) = sqrt(this%gamma_gas*pr(ipoin,2)/rho(ipoin,2))
         end do
         !!$acc end parallel loop
      end if


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
   end subroutine ThermalChannelFlowSolver_evalInitialConditions

end module ThermalChannelFlowSolver_mod

! main.f90
program main
   use mod_constants
   use ThermalChannelFlowSolver_mod

   implicit none
   type(ThermalChannelFlowSolver)  :: channel

   call channel%run()

end program main

