module WindFarmSolverIncomp_mod
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
   use CFDSolverPeriodicWithBoundariesIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: WindFarmSolverIncomp

      real(rp) , public  :: rough,vinf,Lhub,rho0,Lz,ustar

   contains
      procedure, public :: fillBCTypes           => WindFarmSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => WindFarmSolverIncomp_initializeParameters
      procedure, public :: initializeSourceTerms => WindFarmSolverIncomp_initializeSourceTerms
      procedure, public :: evalInitialConditions => WindFarmSolverIncomp_evalInitialConditions
      procedure, public :: eval_vars_after_load_hdf5_resultsFile => WindFarmSolverIncomp_eval_vars_after_load_hdf5_resultsFile 
   end type WindFarmSolverIncomp
contains

   subroutine WindFarmSolverIncomp_eval_vars_after_load_hdf5_resultsFile(this)
      implicit none
      class(WindFarmSolverIncomp), intent(inout) :: this
      integer :: iNodeL,idime

      !values loaded -> rho,u,pr,E,mu_e,mu_sgs

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         eta(iNodeL,2) =0.5_rp*dot_product(u(iNodeL,1:ndime,2),u(iNodeL,1:ndime,2))

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

         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop


      !$acc kernels
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      zo(:) = this%rough
      !$acc end kernels

   end subroutine WindFarmSolverIncomp_eval_vars_after_load_hdf5_resultsFile

   subroutine WindFarmSolverIncomp_fill_BC_Types(this)
      class(WindFarmSolverIncomp), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine WindFarmSolverIncomp_fill_BC_Types

   subroutine WindFarmSolverIncomp_initializeSourceTerms(this)
      class(WindFarmSolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL

      allocate(source_term(numNodesRankPar,ndime))
      !$acc enter data create(source_term(:,:))

      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) =  this%rho0*this%ustar**2/this%Lz
         source_term(iNodeL,2) = 0.00_rp
         source_term(iNodeL,3) = 0.00_rp
      end do
      !$acc end parallel loop

   end subroutine WindFarmSolverIncomp_initializeSourceTerms

   subroutine WindFarmSolverIncomp_initializeParameters(this)
      use json_module
      implicit none
      class(WindFarmSolverIncomp), intent(inout) :: this
      real(rp) :: mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)
      
      ! get(label,target,is found?, default value)

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "wf"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      !----------------------------------------------
      !  --------------  I/O params -------------
      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"channel"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value

      call json%get("results_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%results_h5_file_path,*) value
      call json%get("results_h5_file_name",value, found,"results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value

      !  --------------  I/O params -------------

      call json%get("final_istep",this%final_istep, found,5000001); call this%checkFound(found,found_aux)

      call json%get("save_logFile_first",this%save_logFile_first, found, 1); call this%checkFound(found,found_aux)
      call json%get("save_logFile_step",this%save_logFile_step, found, 10); call this%checkFound(found,found_aux)

      call json%get("save_resultsFile_first",this%save_resultsFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_resultsFile_step" ,this%save_resultsFile_step, found,10000); call this%checkFound(found,found_aux)

      call json%get("save_restartFile_first",this%save_restartFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_restartFile_step" ,this%save_restartFile_step, found,10000); call this%checkFound(found,found_aux)


      call json%get("loadRestartFile" ,this%loadRestartFile, found, .true.); call this%checkFound(found,found_aux)
      call json%get("restartFile_to_load" ,this%restartFile_to_load, found,1); call this%checkFound(found,found_aux)

      call json%get("continue_oldLogs" ,this%continue_oldLogs, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveAvgFile" ,this%saveAvgFile, found, .true.); call this%checkFound(found,found_aux)
      call json%get("loadAvgFile" ,this%loadAvgFile, found, .false.); call this%checkFound(found,found_aux)

      call json%get("initial_avgTime",this%initial_avgTime, found,3600.0_rp); call this%checkFound(found,found_aux)

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)

      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("flag_walave",flag_walave, found,.true.); call this%checkFound(found,found_aux)
      call json%get("period_walave",period_walave, found,3600.0_rp); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Lz",this%Lz, found,1500.0_rp); call this%checkFound(found,found_aux)
      call json%get("Lhub",this%Lhub, found,90.0_rp); call this%checkFound(found,found_aux)
      call json%get("vinf",this%vinf, found,8.0_rp); call this%checkFound(found,found_aux)
      call json%get("rough",this%rough, found,0.1682_rp); call this%checkFound(found,found_aux)

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
      flag_fs_fix_pressure = .false.
      flag_type_wmles = wmles_type_abl
      nscbc_p_inf = 0.0_rp

      this%ustar = this%vinf*0.41_rp/log(1.0_rp+this%Lhub/this%rough)
      incomp_viscosity = 1.81e-5
      flag_mu_factor = 1.0_rp

   end subroutine WindFarmSolverIncomp_initializeParameters

   subroutine WindFarmSolverIncomp_evalInitialConditions(this)
      class(WindFarmSolverIncomp), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, rti(3), zp,velo_aux1
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      logical :: readFiles
      character(512) :: initialField_filePath

      readFiles = .false.

      if(readFiles) then
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         write(initialField_filePath,*) ""
         call read_veloc_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = 0.0_rp 
            rho(iNodeL,2) = this%rho0            

            rho(iNodeL,3) = rho(iNodeL,2)
            pr(iNodeL,3) =  pr(iNodeL,2)
         end do
         !$acc end parallel loop
      else
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         auxCnt = 1
         !!!!$acc parallel loop 
         serialLoop : do iLine = 1,totalNumNodesSrl
            call random_number(rti)
            if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
               iNodeL = matGidSrlOrdered(auxCnt,1)
               auxCnt=auxCnt+1

               zp = coordPar(iNodeL,3)
               velo =  this%ustar*log(1.0_rp+zp/this%rough)/0.41_rp

               u(iNodeL,1,2) = velo*(1.0_rp + 0.1_rp*(rti(1) -0.5_rp))
               u(iNodeL,2,2) = velo*(0.1_rp*(rti(2) -0.5_rp))
               u(iNodeL,3,2) = velo*(0.1_rp*(rti(3) -0.5_rp))
            end if
            if(auxCnt.gt.numNodesRankPar) then
               exit serialLoop
            end if
         end do serialLoop
         !!!!$acc end parallel loop

         !$acc update device(u(:,:,:))

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            pr(iNodeL,2) = 0.0_rp 
            rho(iNodeL,2) = this%rho0            
         end do
         !$acc end parallel loop
      end if

      !$acc kernels
      mu_e(:,:) = 0.0_rp ! Element syabilization viscosity
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      zo(:) = this%rough
      !$acc end kernels
      call nvtxEndRange

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

      !
      ! Initialize exponential averaging for wall law 
      !
      call nvtxStartRange("Wall Average init")
      if(flag_walave .eqv. .true.) then
         !$acc kernels
         walave_u(:,:) = u(:,:,2)
         !$acc end kernels
      end if
      call nvtxEndRange

   end subroutine WindFarmSolverIncomp_evalInitialConditions

end module WindFarmSolverIncomp_mod
