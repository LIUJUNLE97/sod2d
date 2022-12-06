#define ABL 0

module BluffBody3DSolver_mod
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
   use CFDSolver3DWithBoundaries_mod
   implicit none
   private

   type, public, extends(CFDSolver3DWithBoundaries) :: BluffBody3DSolver

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po

   contains
      procedure, public :: fillBCTypes           =>BluffBody3DSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBody3DSolver_initializeParameters
      procedure, public :: evalInitialConditions => BluffBody3DSolver_evalInitialConditions
      procedure, public :: evalViscosityFactor=>BluffBody3DSolver_evalViscosityFactor
   end type BluffBody3DSolver
contains

   subroutine BluffBody3DSolver_fill_BC_Types(this)
      class(BluffBody3DSolver), intent(inout) :: this

#if ABL
      bouCodes2BCType(1) = bc_type_inlet
      bouCodes2BCType(2) = bc_type_outlet
      bouCodes2BCType(3) = bc_type_inlet
      bouCodes2BCType(4) = bc_type_slip_wall_model
      bouCodes2BCType(5) = bc_type_slip_wall_model
#else
      bouCodes2BCType(1) = bc_type_inlet
      bouCodes2BCType(2) = bc_type_slip_wall_model
      bouCodes2BCType(3) = bc_type_non_slip_adiabatic
      bouCodes2BCType(4) = bc_type_slip_wall_model
      bouCodes2BCType(5) = bc_type_slip_adiabatic
      bouCodes2BCType(6) = bc_type_outlet
#endif

   end subroutine BluffBody3DSolver_fill_BC_Types

   subroutine BluffBody3DSolver_initializeParameters(this)
      class(BluffBody3DSolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%gmsh_file_path,*) "./mesh/"
      write(this%gmsh_file_name,*) "auto" 

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "auto"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      this%isPeriodic = .false.

      this%loadMesh = .true.
      this%loadResults = .true.

      this%continue_oldLogs = .false.
      this%load_step = 140001

      this%nstep = 800000001 !250001
#if ABL
      this%cfl_conv = 1.5_rp !0.1_rp
      this%cfl_diff = 1.5_rp !0.1_rp
#else
      this%cfl_conv = 0.95_rp !0.1_rp
      this%cfl_diff = 0.95_rp !0.1_rp
#endif

      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1
      this%nleap = 20000!25 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 25  ! Saving interval, TODO: input
      this%nleapAVG = 20000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 1.0_rp
      this%M  = 0.2_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp
#if ABL
      this%Re     =  37250.0_rp
#else
      this%Re     =  2900000.0_rp
#endif      

      mul    = (this%rho0*this%delta*this%vo)/this%Re
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = mul/mur

      nscbc_u_inf = this%vo
      nscbc_p_inf = this%po
      nscbc_rho_inf = this%rho0
      nscbc_gamma_inf = this%gamma_gas
      nscbc_c_inf = sqrt(this%gamma_gas*this%po/this%rho0)

   end subroutine BluffBody3DSolver_initializeParameters

   subroutine BluffBody3DSolver_evalInitialConditions(this)
      class(BluffBody3DSolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      readFiles = .false.

      if(readFiles) then
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         call read_veloc_from_file_Par(numElemsInRank,numNodesRankPar,totalNumNodesSrl,this%gmsh_file_path,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)
      else
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
              u(iNodeL,1,2) = this%vo
              u(iNodeL,2,2) = 0.0_rp
              u(iNodeL,3,2) = 0.0_rp
         end do
         !$acc end parallel loop
      end if

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

   end subroutine BluffBody3DSolver_evalInitialConditions

   subroutine BluffBody3DSolver_evalViscosityFactor(this)
      class(BluffBody3DSolver), intent(inout) :: this
      integer(4) :: iNodeL

      ! set out of the buffer zone
      ! remember that the mu_factor field has to we filled at least with the
      ! flag_mu_factor
#if ABL
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
         if(coordPar(iNodeL,1)>15.0_rp) then
            mu_factor(iNodeL) = flag_mu_factor*10000.0_rp
         end if
      end do
      !$acc end parallel loop
#else
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
         if(coordPar(iNodeL,1)<-4.5_rp) then
            mu_factor(iNodeL) = flag_mu_factor*10000.0_rp
         end if
         !if(coordPar(iNodeL,2)>0.9_rp) then
         !   mu_factor(iNodeL) = flag_mu_factor*10000.0_rp
         !end if
         !if(coordPar(iNodeL,3)<-8.0_rp) then
         !   mu_factor(iNodeL) = flag_mu_factor*10.0_rp
         !end if
         !if(coordPar(iNodeL,3)>1.2_rp) then
         !   mu_factor(iNodeL) = flag_mu_factor*10000.0_rp
         !end if
         !if(coordPar(iNodeL,2)<-0.9_rp) then
         !   mu_factor(iNodeL) = flag_mu_factor*10000.0_rp
         !end if
         if(coordPar(iNodeL,1)>5.5_rp) then
            mu_factor(iNodeL) = flag_mu_factor*10000.0_rp
         end if
      end do
      !$acc end parallel loop
#endif      

   end subroutine BluffBody3DSolver_evalViscosityFactor
end module BluffBody3DSolver_mod
