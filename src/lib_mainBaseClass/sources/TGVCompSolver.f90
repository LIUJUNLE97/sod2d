module TGVCompSolver_mod
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
   use CFDSolverPeriodic_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodic) :: TGVCompSolver

      real(rp) , public  ::  rho0,Re,to, po

   contains
      procedure, public :: initializeParameters  => TGVCompSolver_initializeParameters
      procedure, public :: evalInitialConditions => TGVCompSolver_evalInitialConditions
   end type TGVCompSolver
contains

   subroutine TGVCompSolver_initializeParameters(this)
      class(TGVCompSolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%gmsh_file_path,*) "./mesh_cube/"
      write(this%gmsh_file_name,*) "cube" 

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "cube" ! Nsys

      this%isPeriodic = .true.
      this%doGlobalAnalysis = .true.
      this%loadMesh = .false.

      this%nstep = 10000 
      this%cfl_conv = 0.95_rp
      this%cfl_diff = 0.95_rp
      this%nsave  = 1  ! First step to save, TODO: input
      this%nsave2 = 1   ! First step to save, TODO: input
      this%nsaveAVG = 1
      this%nleap = 250 ! Saving interval, TODO: input
      this%tleap = 0.5_rp ! Saving interval, TODO: input
      this%nleap2 = 10  ! Saving interval, TODO: input
      this%nleapAVG = 2000000000

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%to  = 1.0_rp
      this%Re = 1600.0_rp
      this%gamma_gas = 1.40_rp

      this%Rgas = 1.0_rp*1.0_rp/(1.4_rp*1.0_rp*1.25_rp*1.25_rp)
      this%rho0 = (1.0_rp/(1.4_rp*1.25*1.25))/(this%Rgas*this%to)

      mul    = (this%rho0*1.0_rp*1.0_rp)/this%Re
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = mul/mur

      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas

   end subroutine TGVCompSolver_initializeParameters

   subroutine TGVCompSolver_evalInitialConditions(this)
      class(TGVCompSolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL

      call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
      call read_veloc_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,this%gmsh_file_path,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)
      call read_densi_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,this%gmsh_file_path,rho(:,2),connecParOrig,Ngp_l,matGidSrlOrdered)
      call read_press_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,this%gmsh_file_path,pr(:,2),connecParOrig,Ngp_l,matGidSrlOrdered)

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
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

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop
   end subroutine TGVCompSolver_evalInitialConditions

end module TGVCompSolver_mod
