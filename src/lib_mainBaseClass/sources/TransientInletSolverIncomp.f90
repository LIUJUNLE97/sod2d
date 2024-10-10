module TransientInletSolverIncomp_mod
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

   integer(8), allocatable, dimension(:,:)   :: matGidSrlOrdered
   !real(rp), allocatable, dimension(:,:)   :: inletDB
   real(rp), allocatable, dimension(:,:) :: uInletDB, vInletDB, wInletDB
   integer(4), allocatable, dimension(:)   :: ipoinInletDB
   integer(4), allocatable, dimension(:)   :: myPointsDB
   integer(4), public :: npoinDB, nStepsDB    

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: TransientInletSolverIncomp

      real(rp) , public  :: vo, delta, rho0, Re, Tp, deltaT          

   contains
      procedure, public :: fillBCTypes           =>TransientInletSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => TransientInletSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => TransientInletSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer => TransientInletSolverIncomp_initialBuffer
      procedure, public :: beforeTimeIteration =>TransientInletSolverIncomp_beforeTimeIteration
   end type TransientInletSolverIncomp
contains

   subroutine TransientInletSolverIncomp_fill_BC_Types(this)
      class(TransientInletSolverIncomp), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine TransientInletSolverIncomp_fill_BC_Types

   subroutine TransientInletSolverIncomp_initializeParameters(this)
      use json_module
      implicit none
      class(TransientInletSolverIncomp), intent(inout) :: this
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
      !----------------------------------------------
      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("v0",this%vo, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("delta",this%delta, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Re",this%Re, found,10000.0_rp); call this%checkFound(found,found_aux)

      !----------------------------------------------
      ! Inlet Database file
      call json%get("inlet_hdf_file_name",value, found,"inflowTurb"); call this%checkFound(found,found_aux)
      write(this%inlet_hdf_file_name,*) value
      ! Read the period T of the inflow velocity
      call json%get("inlet_hdf_period",this%T, found,1.0_rp); call this%checkFound(found,found_aux)
      ! Read the time step of the inflow velocity
      call json%get("inlet_hdf_time_step",this%deltaT, found,1.0_rp); call this%checkFound(found,found_aux)

      ! fixed by the type of base class parameters
      incomp_viscosity = (this%rho0*this%delta*this%vo)/this%Re
      flag_mu_factor = 1.0_rp

      nscbc_u_inf = this%vo
      nscbc_p_inf = 0.0_rp
      nscbc_rho_inf = this%rho0

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

      call json%destroy()

   end subroutine TransientInletSolverIncomp_initializeParameters

   subroutine TransientInletSolverIncomp_evalInitialConditions(this)
      class(TransientInletSolverIncomp), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      call nvtxStartRange("BluffBody_incomp Init")
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo
         u(iNodeL,2,2) = 0.0_rp
         u(iNodeL,3,2) = 0.0_rp
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
      mu_e(:,:) = 0.0_rp 
      mu_sgs(:,:) = 0.0_rp
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      au(:,:) = 0.0_rp
      !$acc end kernels
      call nvtxEndRange

   end subroutine TransientInletSolverIncomp_evalInitialConditions


   subroutine TransientInletSolverIncomp_initialBuffer(this)
      class(TransientInletSolverIncomp), intent(inout) :: this
      integer(4) :: inode, ipoin,iLine,auxCnt,bcode

      allocate(matGidSrlOrdered(numNodesRankPar,2))
      !$acc enter data create(matGidSrlOrdered(:,:))
      call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
      !$acc update device(matGidSrlOrdered(:,:))

      !here we read the npoinDB, nStepsDB
      allocate(inletDB(npoinDB*nStepsDB,4),ipoinInletDB(npoinDB),myPointsDB(npoinDB))
      !$acc enter data create(inletDB(:,:),ipoinInletDB(:),myPointsDB(:))

      !$acc kernels
      myPointsDB(:) = 0
      !$acc end kernels

      ! you need to read your files and update to device

      auxCnt = 1
      do iLine = 1,totalNumNodesSrl
         if(iLine.eq.matGidSrlOrdered(ipoin,2)) then
            inode = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1
            bcode = bouCodesNodesPar(inode) 
            if(bcode .lt. max_num_bou_codes) then
               if (bcode == bc_type_unsteady_inlet) then
                  do ipoin=1,npoinDB
                     if(auxCnt .eq. ipoinInletDB(ipoin)) then
                        myPointsDB(ipoin) = inode
                     end if
                  end do
               end if
            end if
         end if
      end do
      !$acc update device(myPointsDB(:))


      ! you need to do your fancy interpolation
      !$acc parallel loop
      do inode = 1,npoinDB
         ipoin = myPointsDB(ipoin)
         if(ipoin .gt. 0) then
            u_buffer(ipoin,1) = this%vo
            u_buffer(ipoin,2) = 0.0_rp
            u_buffer(ipoin,3) = 0.0_rp  
         end if
      end do
      !$acc end parallel loop

   end subroutine TransientInletSolverIncomp_initialBuffer

   subroutine TransientInletSolverIncomp_beforeTimeIteration(this)
      class(TransientInletSolverIncomp), intent(inout) :: this
      
      integer(4) :: inode, ipoin

      ! you need to do your fancy interpolation
      !$acc parallel loop
      do inode = 1,npoinDB
         ipoin = myPointsDB(ipoin)
         if(ipoin .gt. 0) then
            u_buffer(ipoin,1) = this%vo
            u_buffer(ipoin,2) = 0.0_rp
            u_buffer(ipoin,3) = 0.0_rp  
         end if
      end do
      !$acc end parallel loop

   end subroutine TransientInletSolverIncomp_beforeTimeIteration

end module TransientInletSolverIncomp_mod
