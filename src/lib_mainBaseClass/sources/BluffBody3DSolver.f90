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

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po, aoa_alpha, aoa_beta


   contains
      procedure, public :: fillBCTypes           =>BluffBody3DSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBody3DSolver_initializeParameters
      procedure, public :: evalInitialConditions => BluffBody3DSolver_evalInitialConditions
      procedure, public :: initialBuffer => BluffBody3DSolver_initialBuffer
   end type BluffBody3DSolver
contains

   subroutine BluffBody3DSolver_fill_BC_Types(this)
      class(BluffBody3DSolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BluffBody3DSolver_fill_BC_Types

   subroutine BluffBody3DSolver_initialBuffer(this)
      class(BluffBody3DSolver), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u_buffer(iNodeL,1) = this%vo*cos(this%aoa_alpha*v_pi/180.0_rp)*cos(this%aoa_beta*v_pi/180.0_rp)
         u_buffer(iNodeL,2) = this%vo*sin(this%aoa_beta*v_pi/180.0_rp) 
         u_buffer(iNodeL,3) = this%vo*sin(this%aoa_alpha*v_pi/180.0_rp)    
      end do
      !$acc end parallel loop

   end subroutine BluffBody3DSolver_initialBuffer

   subroutine BluffBody3DSolver_initializeParameters(this)
      use json_module
      implicit none
      class(BluffBody3DSolver), intent(inout) :: this
      real(rp) :: mul, mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)

      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"cetaceo"); call this%checkFound(found,found_aux)
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

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.true.); call this%checkFound(found,found_aux)
      !----------------------------------------------

      ! numerical params

      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("flag_implicit",flag_implicit, found,0); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)

      call json%get("period_walave",period_walave, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("flag_walave",flag_walave, found,.true.); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("Cp",this%Cp, found,1004.0_rp); call this%checkFound(found,found_aux)
      call json%get("Prt",this%Prt, found,0.71_rp); call this%checkFound(found,found_aux)
      call json%get("v0",this%vo, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("M",this%M, found,0.2_rp); call this%checkFound(found,found_aux)
      call json%get("delta",this%delta, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Re",this%Re, found,5600000.0_rp); call this%checkFound(found,found_aux)
      call json%get("gamma_gas",this%gamma_gas, found,1.4_rp); call this%checkFound(found,found_aux)
      call json%get("aoa_alpha",this%aoa_alpha, found,11.0_rp); call this%checkFound(found,found_aux)
      call json%get("aoa_beta",this%aoa_beta, found,0.0_rp); call this%checkFound(found,found_aux)

      !Witness points parameters
      call json%get("have_witness",this%have_witness, found,.false.); call this%checkFound(found,found_aux)
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

      call this%readJSONBuffer()

      call json%destroy()
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
            u(iNodeL,1,2) = this%vo*cos(this%aoa_alpha*v_pi/180.0_rp)*cos(this%aoa_beta*v_pi/180.0_rp)
            u(iNodeL,2,2) = this%vo*sin(this%aoa_beta*v_pi/180.0_rp) 
            u(iNodeL,3,2) = this%vo*sin(this%aoa_alpha*v_pi/180.0_rp)    
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
