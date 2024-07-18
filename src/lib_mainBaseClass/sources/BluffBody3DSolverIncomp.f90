module BluffBody3DSolverIncomp_mod
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
   use CFDSolver3DWithBoundariesIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolver3DWithBoundariesIncomp) :: BluffBody3DSolverIncomp

      real(rp) , public  :: vo,delta, rho0, Re, aoa_alpha, aoa_beta

   contains
      procedure, public :: fillBCTypes           =>BluffBody3DSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBody3DSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => BluffBody3DSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer => BluffBody3DSolverIncomp_initialBuffer
   end type BluffBody3DSolverIncomp
contains

   subroutine BluffBody3DSolverIncomp_fill_BC_Types(this)
      class(BluffBody3DSolverIncomp), intent(inout) :: this
      
      call this%readJSONBCTypes()

   end subroutine BluffBody3DSolverIncomp_fill_BC_Types

   subroutine BluffBody3DSolverIncomp_initialBuffer(this)
      class(BluffBody3DSolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar

         u_buffer(iNodeL,1) = this%vo*cos(this%aoa_alpha*v_pi/180.0_rp)*cos(this%aoa_beta*v_pi/180.0_rp)
         u_buffer(iNodeL,2) = this%vo*sin(this%aoa_beta*v_pi/180.0_rp) 
         u_buffer(iNodeL,3) = this%vo*sin(this%aoa_alpha*v_pi/180.0_rp)    
         !just a momentary trick
         pr(iNodeL,2) = 0.0_rp 
         rho(iNodeL,2) = this%rho0            

         rho(iNodeL,3) = rho(iNodeL,2)
         pr(iNodeL,3) =  pr(iNodeL,2)
      end do
      !$acc end parallel loop

   end subroutine BluffBody3DSolverIncomp_initialBuffer

   subroutine BluffBody3DSolverIncomp_initializeParameters(this)
      use json_module
      implicit none
      class(BluffBody3DSolverIncomp), intent(inout) :: this
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

      !----------------------------------------------
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
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("period_walave",period_walave, found,1.0_rp); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("v0",this%vo, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("delta",this%delta, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Re",this%Re, found,5600000.0_rp); call this%checkFound(found,found_aux)
      call json%get("aoa_alpha",this%aoa_alpha, found,11.0_rp); call this%checkFound(found,found_aux)
      call json%get("aoa_beta",this%aoa_beta, found,0.0_rp); call this%checkFound(found,found_aux)

      call json%get("c_sgs",c_sgs, found,0.025_rp); 

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

      ! fixed by the type of base class parameters

      mul    = (this%rho0*this%delta*this%vo)/this%Re
      incomp_viscosity = mul
      flag_mu_factor = 1.0_rp

      nscbc_u_inf = this%vo
      nscbc_p_inf = 0.0_rp
      nscbc_rho_inf = this%rho0

      call this%readJSONBuffer()

      call json%destroy()

   end subroutine BluffBody3DSolverIncomp_initializeParameters

   subroutine BluffBody3DSolverIncomp_evalInitialConditions(this)
      class(BluffBody3DSolverIncomp), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo*cos(this%aoa_alpha*v_pi/180.0_rp)*cos(this%aoa_beta*v_pi/180.0_rp)
         u(iNodeL,2,2) = this%vo*sin(this%aoa_beta*v_pi/180.0_rp) 
         u(iNodeL,3,2) = this%vo*sin(this%aoa_alpha*v_pi/180.0_rp)    
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         rho(iNodeL,2) = this%rho0

         eta(iNodeL,2) =0.5_rp*dot_product(u(iNodeL,1:ndime,2),u(iNodeL,1:ndime,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         u(iNodeL,1:ndime,3) = u(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         eta(iNodeL,3) =  eta(iNodeL,2)
         pr(iNodeL,3) =  pr(iNodeL,2)  

         q(iNodeL,1:ndime,4) = q(iNodeL,1:ndime,2)
         u(iNodeL,1:ndime,4) = u(iNodeL,1:ndime,2)
         rho(iNodeL,4) = rho(iNodeL,2)
         eta(iNodeL,4) =  eta(iNodeL,2)
         pr(iNodeL,4) =  pr(iNodeL,2)  
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

   end subroutine BluffBody3DSolverIncomp_evalInitialConditions

end module BluffBody3DSolverIncomp_mod
