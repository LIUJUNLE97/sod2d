module TGVSolverIncomp_mod
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
   use mod_saveFields
   use CFDSolverPeriodicIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicIncomp) :: TGVSolverIncomp

      real(rp) , public  ::  Re

   contains
      procedure, public :: initializeParameters  => TGVSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => TGVSolverIncomp_evalInitialConditions
   end type TGVSolverIncomp
contains

   subroutine TGVSolverIncomp_initializeParameters(this)
      use json_module
      implicit none
      class(TGVSolverIncomp), intent(inout) :: this
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)
      
      ! get(label,target,is found?, default value)

      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"cube"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value
      
      call json%get("results_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%results_h5_file_path,*) value
      call json%get("results_h5_file_name",value, found,"results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value

      call json%get("save_logFile_first",this%save_logFile_first, found, 1); call this%checkFound(found,found_aux)
      call json%get("save_logFile_step",this%save_logFile_step, found, 10); call this%checkFound(found,found_aux)

      call json%get("save_resultsFile_first",this%save_resultsFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_resultsFile_step" ,this%save_resultsFile_step, found,200); call this%checkFound(found,found_aux)

      call json%get("save_restartFile_first",this%save_restartFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_restartFile_step" ,this%save_restartFile_step, found,200); call this%checkFound(found,found_aux)

      call json%get("loadRestartFile" ,this%loadRestartFile, found, .false.); call this%checkFound(found,found_aux)
      call json%get("restartFile_to_load" ,this%restartFile_to_load, found,1); call this%checkFound(found,found_aux)

      call json%get("continue_oldLogs" ,this%continue_oldLogs, found, .false.); call this%checkFound(found,found_aux)
      call json%get("saveAvgFile" ,this%saveAvgFile, found, .false.); call this%checkFound(found,found_aux)

      call json%get("loadAvgFile" ,this%loadAvgFile, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)

      call json%get("saveInitialField",this%saveInitialField, found,.true.); call this%checkFound(found,found_aux)

      call json%get("doGlobalAnalysis",this%doGlobalAnalysis, found, .false.); call this%checkFound(found,found_aux)
      call json%get("doTimerAnalysis",this%doTimerAnalysis, found, .true.); call this%checkFound(found,found_aux)

      call json%get("final_istep",this%final_istep, found,500001); call this%checkFound(found,found_aux)
      call json%get("maxPhysTime",this%maxPhysTime, found,20.0_rp); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.5_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.5_rp); call this%checkFound(found,found_aux)

      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.01d0); call this%checkFound(found,found_aux)

      call json%get("Re",this%Re, found,1600.0_rp); call this%checkFound(found,found_aux)

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
      incomp_viscosity = 1.0_rp/this%Re
      flag_fs_fix_pressure = .false.
      flag_mu_factor = 1.0_rp
      nscbc_p_inf = 0.0_rp

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'
   end subroutine TGVSolverIncomp_initializeParameters

   subroutine TGVSolverIncomp_evalInitialConditions(this)
      class(TGVSolverIncomp), intent(inout) :: this
      real(4) :: V0,L
      real(4) :: x,y,z
      integer(4) :: iNodeL,iNodeGSrl

      if(mpi_rank.eq.0) write(*,*) "--| TGV - Setting Initial Conditions..."

      call nvtxStartRange("TGVSolver_incomp Init")
      V0 = 1.0_rp
      L  = 1.0_rp

      !$acc parallel loop
      do iNodeL=1,numNodesRankPar
         x = coordPar(iNodeL,1)
         y = coordPar(iNodeL,2)
         z = coordPar(iNodeL,3)

         u(iNodeL,1,2) =  V0*sin(x/(L))*cos(y/(L))*cos(z/(L))
         u(iNodeL,2,2) = -V0*cos(x/(L))*sin(y/(L))*cos(z/(L))
         u(iNodeL,3,2) = 0.0

         pr(iNodeL,2)  = nscbc_p_inf+((1.0_rp*V0*V0)/(16.0_rp))*(cos(2.0_rp*x/L)+cos(2.0_rp*y/L))*(cos(2.0_rp*z/L)+2.0_rp)
         rho(iNodeL,2) = 1.0_rp
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
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

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop
   end subroutine TGVSolverIncomp_evalInitialConditions

end module TGVSolverIncomp_mod
