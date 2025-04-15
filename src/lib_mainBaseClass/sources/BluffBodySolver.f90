module BluffBodySolver_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BluffBodySolver

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po, aoa

   contains
      procedure, public :: fillBCTypes           =>BluffBodySolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolver_initializeParameters
      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditions
      procedure, public :: initialBuffer => BluffBodySolver_initialBuffer
   end type BluffBodySolver
contains

   subroutine BluffBodySolver_fill_BC_Types(this)
      class(BluffBodySolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BluffBodySolver_fill_BC_Types

   subroutine BluffBodySolver_initializeParameters(this)
      use json_module
      implicit none
      class(BluffBodySolver), intent(inout) :: this
      real(rp) :: mul, mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

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
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("flag_implicit",flag_implicit, found,0); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("flag_imex_stages",flag_imex_stages, found,4); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("flag_high_mach",flag_high_mach, found,.true.); call this%checkFound(found,found_aux)

      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,1.0_rp); call this%checkFound(found,found_aux)

      call json%get("period_walave",period_walave, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("wmles_walex",wmles_walex, found,0.1_rp); !optional depending of the model

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("Cp",this%Cp, found,1004.0_rp); call this%checkFound(found,found_aux)
      call json%get("Prt",this%Prt, found,0.71_rp); call this%checkFound(found,found_aux)
      call json%get("v0",this%vo, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("M",this%M, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("delta",this%delta, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Re",this%Re, found,10000.0_rp); call this%checkFound(found,found_aux)
      call json%get("gamma_gas",this%gamma_gas, found,1.4_rp); call this%checkFound(found,found_aux)
      call json%get("aoa",this%aoa, found,0.0_rp); call this%checkFound(found,found_aux)

      call json%get("flag_rk_ls",flag_rk_ls, found,.true.); 
      call json%get("flag_rk_ls_stages",flag_rk_ls_stages, found,5); 
      call json%get("flag_rk_ls_n",flag_rk_ls_n, found, 1)
      call json%get("c_sgs",c_sgs, found,0.025_rp); 

      call json%get("flag_lps_stab",flag_lps_stab, found,.true.); call this%checkFound(found,found_aux)

      call json%get("flag_trip_element",flag_trip_element, found,.false.); call this%checkFound(found,found_aux)
      call json%get("x_trip_o",x_trip_o, found,0.0_rp); call this%checkFound(found,found_aux)
      call json%get("y_trip_o",y_trip_o, found,0.0_rp); call this%checkFound(found,found_aux)
      call json%get("l_trip_x",l_trip_x, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("l_trip_y",l_trip_y, found,1.0_rp); call this%checkFound(found,found_aux)


      ! fixed by the type of base class parameters
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

      !Witness points parameters
      call json%get("have_witness",this%have_witness, found,.false.)
      if(this%have_witness .eqv. .true.) then
         call json%get("witness_inp_file_name",value, found,"witness.txt"); call this%checkFound(found,found_aux)
         write(this%witness_inp_file_name,*) value
         call json%get("witness_h5_file_name",value, found,"resultwit.h5"); call this%checkFound(found,found_aux)
         write(this%witness_h5_file_name,*) value

         call json%get("leapwit",this%leapwit, found,1); call this%checkFound(found,found_aux)
         call json%get("nwit",this%nwit, found,17986); call this%checkFound(found,found_aux)
         call json%get("wit_save_u_i",this%wit_save_u_i, found,.true.); call this%checkFound(found,found_aux)
         call json%get("wit_save_pr",this%wit_save_pr, found,.true.); call this%checkFound(found,found_aux)
         call json%get("wit_save_rho",this%wit_save_rho, found,.true.); call this%checkFound(found,found_aux)
         call json%get("continue_witness",this%continue_witness, found,.false.); call this%checkFound(found,found_aux)
      end if  

      call this%readJSONBuffer()
      call this%readJSONWMTypes()


      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'
   end subroutine BluffBodySolver_initializeParameters

   subroutine BluffBodySolver_evalInitialConditions(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      call nvtxStartRange("BluffBodySolver Init")
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo*cos(this%aoa*v_pi/180.0_rp)
         u(iNodeL,2,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)
         u(iNodeL,3,2) = 0.0_rp
      end do
      !$acc end parallel loop

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
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)

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

   end subroutine BluffBodySolver_evalInitialConditions


   subroutine BluffBodySolver_initialBuffer(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo*cos(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolver_initialBuffer

end module BluffBodySolver_mod
