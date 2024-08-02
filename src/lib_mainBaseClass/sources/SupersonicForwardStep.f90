module SupersonicForwardStep_mod
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
   use CFDSolverPeriodicWithBoundaries_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: SupersonicForwardStep

      real(rp) , public  :: rho0,to, po

   contains
      procedure, public :: fillBCTypes           =>SupersonicForwardStep_fill_BC_Types
      procedure, public :: initializeParameters  => SupersonicForwardStep_initializeParameters
      procedure, public :: evalInitialConditions => SupersonicForwardStep_evalInitialConditions
      procedure, public :: initialBuffer => SupersonicForwardStep_initialBuffer
   end type SupersonicForwardStep
contains

   subroutine SupersonicForwardStep_fill_BC_Types(this)
      class(SupersonicForwardStep), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine SupersonicForwardStep_fill_BC_Types

   subroutine SupersonicForwardStep_initializeParameters(this)
      use json_module
      implicit none
      class(SupersonicForwardStep), intent(inout) :: this
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value
      real(rp) :: mul, mur

      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)

      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"cylin"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value

      call json%get("results_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%results_h5_file_path,*) value
      call json%get("results_h5_file_name",value, found,"results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value

      !  --------------  I/O params -------------
      call json%get("final_istep",this%final_istep, found,1000001); call this%checkFound(found,found_aux)

      call json%get("save_logFile_first",this%save_logFile_first, found, 1); call this%checkFound(found,found_aux)
      call json%get("save_logFile_step",this%save_logFile_step, found, 10); call this%checkFound(found,found_aux)

      call json%get("save_resultsFile_first",this%save_resultsFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_resultsFile_step" ,this%save_resultsFile_step, found,10000); call this%checkFound(found,found_aux)

      call json%get("save_restartFile_first",this%save_restartFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_restartFile_step" ,this%save_restartFile_step, found,10000); call this%checkFound(found,found_aux)


      call json%get("loadRestartFile" ,this%loadRestartFile, found, .false.); call this%checkFound(found,found_aux)
      call json%get("restartFile_to_load" ,this%restartFile_to_load, found,1); call this%checkFound(found,found_aux)

      call json%get("continue_oldLogs" ,this%continue_oldLogs, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveAvgFile" ,this%saveAvgFile, found, .true.); call this%checkFound(found,found_aux)
      call json%get("loadAvgFile" ,this%loadAvgFile, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)

      call json%get("doTimerAnalysis",this%doTimerAnalysis, found,.false.)

      ! numerical params
      call json%get("flag_implicit",flag_implicit, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("flag_rk_ls",flag_rk_ls, found,.true.); 
      call json%get("flag_rk_ls_stages",flag_rk_ls_stages, found,5); 
      call json%get("flag_high_mach",flag_high_mach, found,.true.); call this%checkFound(found,found_aux)

      this%saveInitialField = .false.

      !ce = 100.0_rp
      this%maxPhysTime = 5
      flag_real_diff = 0
      flag_force_2D = .true.
      flag_normalise_entropy = 1
      !flag_use_constant_dt = 1
      !this%dt = 0.0001_rp
      !flag_rk_ls = .false.

      ! fixed by the type of base class parameters
      !this%Rgas = 0.714285_rp
      !this%gamma_gas = this%Rgas/1.78571_rp + 1.0_rp
      !this%Cp = this%Rgas + 1.78571_rp 
      !this%to = 1.0_rp
      !this%po = 1.0_rp
      
      this%to = 1.0_rp
      this%po = 1.0_rp
      this%rho0 = 1.4_rp
      this%Rgas = this%po/(this%rho0*this%to)
      this%gamma_gas = 1.40_rp
      this%Cp = this%Rgas*this%gamma_gas/(this%gamma_gas-1.0_rp)
      this%rho0 = this%po/this%Rgas/this%to

      mul    = real(18.0e-6,rp)
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = 0.0_rp*mul/mur
      !this%Prt = this%Cp*real(18.0e-6,rp)/real(32.3e-6,rp)
      this%Prt = 0.71_rp

      nscbc_c_inf = sqrt(this%gamma_gas*this%Rgas*this%to)
      nscbc_u_inf = 3.0_rp*nscbc_c_inf
      nscbc_p_inf = this%po
      nscbc_rho_inf = this%po/this%Rgas/this%to
      nscbc_gamma_inf = this%gamma_gas
      nscbc_Rgas_inf = this%Rgas

      call this%readJSONBuffer()

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'
   end subroutine SupersonicForwardStep_initializeParameters

   subroutine SupersonicForwardStep_evalInitialConditions(this)
      class(SupersonicForwardStep), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = 3.0_rp*nscbc_c_inf
         u(iNodeL,2,2) = 0.0_rp
         u(iNodeL,3,2) = 0.0_rp
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%rho0
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = this%to
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         !eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(rho(iNodeL,2)*e_int(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = sqrt(dot_product(u(iNodeL,:,2),u(iNodeL,:,2)))/csound(iNodeL)

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) = eta(iNodeL,2) 
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

   end subroutine SupersonicForwardStep_evalInitialConditions


   subroutine SupersonicForwardStep_initialBuffer(this)
      class(SupersonicForwardStep), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = 3.0_rp*nscbc_c_inf
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine SupersonicForwardStep_initialBuffer

end module SupersonicForwardStep_mod
