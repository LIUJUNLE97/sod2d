#define CRM 1

module BluffBody3DSolver_mod
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
   use mod_analysis
   use mod_numerical_params
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

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po, aoa

   contains
      procedure, public :: fillBCTypes           =>BluffBody3DSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBody3DSolver_initializeParameters
      procedure, public :: evalInitialConditions => BluffBody3DSolver_evalInitialConditions
      procedure, public :: initialBuffer => BluffBody3DSolver_initialBuffer
   end type BluffBody3DSolver
contains

   subroutine BluffBody3DSolver_fill_BC_Types(this)
      class(BluffBody3DSolver), intent(inout) :: this

#if CRM
      
       bouCodes2BCType(1) = bc_type_slip_wall_model
       bouCodes2BCType(2) = bc_type_far_field 
       bouCodes2BCType(3) = bc_type_far_field
      bouCodes2BCType(4) = bc_type_far_field
      bouCodes2BCType(5) = bc_type_far_field


#else
     ! bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! floor
     ! bouCodes2BCType(2) = bc_type_far_field ! top wall
     ! bouCodes2BCType(3) = bc_type_far_field ! inlet
     ! bouCodes2BCType(4) = bc_type_far_field ! outlet
     ! bouCodes2BCType(5) = bc_type_non_slip_adiabatic !back
     ! bouCodes2BCType(6) = bc_type_non_slip_adiabatic !pins
     ! bouCodes2BCType(7) = bc_type_non_slip_adiabatic !car

      bouCodes2BCType(1) = bc_type_slip_wall_model
      bouCodes2BCType(2) = bc_type_far_field
      bouCodes2BCType(3) = bc_type_far_field
      bouCodes2BCType(4) = bc_type_far_field
      bouCodes2BCType(5) = bc_type_non_slip_adiabatic
      bouCodes2BCType(6) = bc_type_slip_wall_model
#endif
      !$acc update device(bouCodes2BCType(:))

   end subroutine BluffBody3DSolver_fill_BC_Types

      subroutine BluffBody3DSolver_initialBuffer(this)
      class(BluffBody3DSolver), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
#if CRM
            u_buffer(iNodeL,1) = this%vo*cos(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = this%vo*sin(this%aoa*v_pi/180.0_rp)    
#else
            u_buffer(iNodeL,1) = this%vo
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
#endif
      end do
      !$acc end parallel loop

   end subroutine BluffBody3DSolver_initialBuffer

   subroutine BluffBody3DSolver_initializeParameters(this)
      class(BluffBody3DSolver), intent(inout) :: this
      real(rp) :: mul, mur

#if CRM
      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "crm"
#else
      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "windsor"
      !write(this%mesh_h5_file_name,*) "auto"
#endif

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 800000001

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 5000
      this%save_resultsFile_step = 5000

      this%save_restartFile_first = 5000
      this%save_restartFile_step = 5000
      this%loadRestartFile = .true.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%initial_avgTime = 0.0_rp
      this%saveAvgFile = .false.
      this%loadAvgFile = .false.
      this%saveSurfaceResults = .true.
      this%saveInitialField = .false.
      
      this%saveSurfaceResults = .true.
      !----------------------------------------------

      ! numerical params
      flag_les = 1
      flag_implicit = 1
      flag_rk_order=4
      !flag_total_enthalpy = .true.

      maxIter = 20
      tol = 1e-3

      period_walave   = 1.0_rp
      flag_walave     = .true.

#if CRM
      !this%dt = 1e-3
      !flag_use_constant_dt = 1
      this%cfl_conv = 0.95_rp 
      this%cfl_diff = 100.0_rp 
#else  
      !this%dt = 5e-3
      !flag_use_constant_dt = 1 
      this%cfl_conv = 0.5_rp
      this%cfl_diff = 0.5_rp
#endif
      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 1.0_rp
      this%M  = 0.2_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp
#if CRM
      this%Re  =  5600000.0_rp
      this%aoa = 11.00_rp
#else
      this%Re     =  2900000.0_rp
      this%aoa = 0.0_rp
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
      nscbc_Rgas_inf = this%Rgas


      flag_buffer_on = .true.
#if CRM
! Case1 Vangelis
      flag_buffer_on_east = .true.
      flag_buffer_e_min = 200.0_rp
      flag_buffer_e_size = 36.0_rp 

      flag_buffer_on_west = .true.
      flag_buffer_w_min = -200.0_rp
      flag_buffer_w_size = 36.0_rp 

      flag_buffer_on_north = .true.
      flag_buffer_n_min = 200.0_rp
      flag_buffer_n_size = 36.0_rp 
      
      flag_buffer_on_top = .true.
      flag_buffer_t_min = 200.0_rp
      flag_buffer_t_size = 36.0_rp

      flag_buffer_on_bottom = .true.
      flag_buffer_b_min = -200.0_rp
      flag_buffer_b_size = 36.0_rp

! Xevi
!       flag_buffer_e_min = 80.0_rp
!       flag_buffer_e_size = 30.0_rp 
! 
!       flag_buffer_on_west = .true.
!       flag_buffer_w_min = -80.0_rp
!       flag_buffer_w_size = 30.0_rp 
! 
!       flag_buffer_on_north = .true.
!       flag_buffer_n_min = 80.0_rp
!       flag_buffer_n_size = 30.0_rp 
!       
!       flag_buffer_on_top = .true.
!       flag_buffer_t_min = 80.0_rp
!       flag_buffer_t_size = 30.0_rp
! 
!       flag_buffer_on_bottom = .true.
!       flag_buffer_b_min = -80.0_rp
!       flag_buffer_b_size = 30.0_rp
#else 
      !!windsor
      flag_buffer_on_east = .true.
      flag_buffer_e_min = 5.5_rp
      flag_buffer_e_size = 1.0_rp

      flag_buffer_on_west = .true.
      flag_buffer_w_min = -3.5_rp
      flag_buffer_w_size = 1.0_rp

      flag_buffer_on_north = .true.
      flag_buffer_n_min = 1.5_rp
      flag_buffer_n_size = 0.5_rp

      flag_buffer_on_south = .true.
      flag_buffer_s_min = -1.5_rp
      flag_buffer_s_size = 0.5_rp

      flag_buffer_on_top = .true.
      flag_buffer_t_min = 2.0_rp
      flag_buffer_t_size = 0.5_rp
#endif
      

   end subroutine BluffBody3DSolver_initializeParameters

   subroutine BluffBody3DSolver_evalInitialConditions(this)
      class(BluffBody3DSolver), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles
      character(512) :: initialField_filePath

      readFiles = .false.

      if(readFiles) then
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         write(initialField_filePath,*) ""
         call read_veloc_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)
      else
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
#if CRM
            u(iNodeL,1,2) = this%vo*cos(this%aoa*v_pi/180.0_rp)
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)    
#else
            u(iNodeL,1,2) = this%vo
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = 0.0_rp  
#endif
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
         
         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) = eta(iNodeL,2) 
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

end module BluffBody3DSolver_mod
