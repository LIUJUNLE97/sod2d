
module ThermalChannelFlowSolver_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: ThermalChannelFlowSolver

      real(rp) , public  ::  delta, rho, Retau, utau, tC,tH, po,to, mu,tauw

   contains
      procedure, public :: initializeParameters  => ThermalChannelFlowSolver_initializeParameters
      procedure, public :: initializeSourceTerms => ThermalChannelFlowSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => ThermalChannelFlowSolver_evalInitialConditions
      procedure, public :: afterDt => ThermalChannelFlowSolver_afterDt
   end type ThermalChannelFlowSolver
contains

   subroutine ThermalChannelFlowSolver_initializeSourceTerms(this)
      class(ThermalChannelFlowSolver), intent(inout) :: this

        allocate(source_term(ndime))
        source_term(1) = this%tauw/this%delta
        source_term(2) = 0.00_rp
        source_term(3) = 0.00_rp

   end subroutine ThermalChannelFlowSolver_initializeSourceTerms

   subroutine ThermalChannelFlowSolver_initializeParameters(this)
      class(ThermalChannelFlowSolver), intent(inout) :: this

      write(this%gmsh_file_path,*) "./mesh/"
      write(this%gmsh_file_name,*) "channel_sem" 

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "channel_sem" ! Nsys

      this%isPeriodic = .true.
      this%loadMesh = .true.
      this%loadResults = .true.
      this%continue_oldLogs = .true.
      this%load_step = 3880001
      this%initial_istep=3880001

      this%nstep = 90000000
      this%cfl_conv = 1.5_rp
      this%cfl_diff = 1.5_rp
      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1
      this%nleap = 40000 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 50  ! Saving interval, TODO: input
      this%nleapAVG = 40000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%tC = 293.0_rp
      this%tH = 586.0_rp
      this%delta  = 0.0015_rp*2.0_rp
      this%gamma_gas = 1.40_rp
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = 0.5_rp*(this%tC+this%tH)
      this%po  = 101325.0_rp 
      this%rho = this%po/(this%Rgas*this%to)

      this%mu = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)

      this%Retau = 400.0_rp

      this%utau = (this%Retau*this%mu)/(this%delta*this%rho)

      this%tauw = this%rho*this%utau*this%utau

      flag_mu_factor = 1.0_rp
      write(111,*) " Gp ", this%tauw/this%delta
      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas

   end subroutine ThermalChannelFlowSolver_initializeParameters

   subroutine ThermalChannelFlowSolver_evalInitialConditions(this)
      class(ThermalChannelFlowSolver), intent(inout) :: this
      real(rp) :: iniU(totalNumNodesSrl,ndime), iniRho(totalNumNodesSrl), iniT(totalNumNodesSrl)
      integer(4) :: iNodeL,iNodeGSrl
      logical :: readFiles
      real(rp) :: velo, ti(3), yp

      readFiles = .true.

      if(readFiles) then
         !call read_veloc(this%npoin,this%gmsh_file_path,u(:,:,2))
         !call read_densi(this%npoin,this%gmsh_file_path,rho(:,2))
         !call read_temper(this%npoin,this%gmsh_file_path,Tem(:,2))
         call read_veloc_from_file_Srl(totalNumNodesSrl,this%gmsh_file_path,iniU)
         call read_densi_from_file_Srl(totalNumNodesSrl,this%gmsh_file_path,iniRho)
         call read_temper_from_file_Srl(totalNumNodesSrl,this%gmsh_file_path,iniT)

         do iNodeL=1,numNodesRankPar
            iNodeGSrl=globalIdSrl(iNodeL)
            u(iNodeL,1,2) = iniU(iNodeGSrl,1)
            u(iNodeL,2,2) = iniU(iNodeGSrl,2)
            u(iNodeL,3,2) = iniU(iNodeGSrl,3)
            rho(iNodeL,2) = iniRho(iNodeGSrl)
            Tem(iNodeL,2)  = iniT(iNodeGSrl)
         end do

         !!$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = Tem(iNodeL,2)*this%Rgas*rho(iNodeL,2)
            e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
            E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
            q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
            csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         end do
         !!$acc end parallel loop
      else
         !!$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if(coordPar(iNodeL,2)<this%delta) then
               yp = coordPar(iNodeL,2)*this%utau*this%rho/this%mu
            else
               yp = abs(coordPar(iNodeL,2)-2.0_rp*this%delta)*this%utau*this%rho/this%mu
            end if

            velo = this%utau*((1.0_rp/0.41_rp)*log(1.0_rp+0.41_rp*yp)+7.8_rp*(1.0_rp-exp(-yp/11.0_rp)-(yp/11.0_rp)*exp(-yp/3.0_rp))) 

            call random_number(ti)


            u(iNodeL,1,2) = velo*(1.0_rp + 0.1_rp*(ti(1) -0.5_rp))
            u(iNodeL,2,2) = velo*(0.1_rp*(ti(2) -0.5_rp))
            u(iNodeL,3,2) = velo*(0.1_rp*(ti(3) -0.5_rp))
            pr(iNodeL,2) = this%po
            rho(iNodeL,2) = this%rho
            e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
            Tem(iNodeL,2) = this%to
            E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
            q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
            csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
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
   end subroutine ThermalChannelFlowSolver_evalInitialConditions

   subroutine ThermalChannelFlowSolver_afterDt(this,istep)
      class(ThermalChannelFlowSolver), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      integer(4) :: codeH, codeC
      real(rp) :: area,tw,RetauC,RetauH,Retau,rhoC,rhoH,muH,muC,twH,twC,utauH,utauC

      if(istep == this%nsave2) then
         codeH = 6
         codeC = 7
         area = this%delta*2.0*v_pi*this%delta*v_pi
         twH = Ftau(codeH,1)/area
         twC = Ftau(codeC,1)/area

         rhoH=this%po/(this%Rgas*this%tH)
         rhoC=this%po/(this%Rgas*this%tC)

         utauH = sqrt(twH/rhoH)
         utauC = sqrt(twC/rhoC)

         muH = 0.000001458_rp*(this%tH**1.50_rp)/(this%tH+110.40_rp)
         muC = 0.000001458_rp*(this%tC**1.50_rp)/(this%tC+110.40_rp)

         RetauH = this%delta*utauH*rhoH/muH
         RetauC = this%delta*utauC*rhoC/muC

         Retau = 0.5_rp*(RetauH+RetauC)

         if(Retau .le. this%Retau) then 
            source_term(1) = source_term(1)*1.05_rp
         else  
            source_term(1) = source_term(1)*0.95_rp
         end if
        source_term(2) = 0.00_rp
        source_term(3) = 0.00_rp
      end if

   end subroutine ThermalChannelFlowSolver_afterDt

end module ThermalChannelFlowSolver_mod
