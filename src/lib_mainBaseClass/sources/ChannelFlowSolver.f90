module ChannelFlowSolver_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: ChannelFlowSolver

      real(rp) , public  :: vo, M, delta, rho0, Retau, Re, utau, to, po, mu

   contains
      procedure, public :: fillBCTypes           => ChannelFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters  => ChannelFlowSolver_initializeParameters
      procedure, public :: initializeSourceTerms => ChannelFlowSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => ChannelFlowSolver_evalInitialConditions
   end type ChannelFlowSolver
contains

   subroutine ChannelFlowSolver_fill_BC_Types(this)
      class(ChannelFlowSolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine ChannelFlowSolver_fill_BC_Types

   subroutine ChannelFlowSolver_initializeSourceTerms(this)
      class(ChannelFlowSolver), intent(inout) :: this
      integer(4) :: iNodeL
      real(rp) :: source_x

      allocate(source_term(numNodesRankPar,ndime+2))
      !$acc enter data create(source_term(:,:))

      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_x = (this%utau*this%utau*this%rho0/this%delta)
         
         source_term(iNodeL,1) = 0.0_rp ! Mass
         source_term(iNodeL,2) = 0.0_rp ! Energy
         ! Momentum
         source_term(iNodeL,3) = source_x
         source_term(iNodeL,4) = 0.0_rp
         source_term(iNodeL,5) = 0.0_rp
      end do
      !$acc end parallel loop

   end subroutine ChannelFlowSolver_initializeSourceTerms

   subroutine ChannelFlowSolver_initializeParameters(this)
      use json_module
      implicit none
      class(ChannelFlowSolver), intent(inout) :: this
      real(rp) :: mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)

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

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)

      call json%get("doTimerAnalysis",this%doTimerAnalysis, found,.false.)
      !----------------------------------------------

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
       
      call json%get("period_walave",period_walave, found,20.0_rp); call this%checkFound(found,found_aux)
      call json%get("wmles_walex",wmles_walex, found,0.1_rp); !optional depending of the model

      call json%get("cfl_conv",this%cfl_conv, found,1.5_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,1.5_rp); call this%checkFound(found,found_aux)

      call json%get("Cp",this%Cp, found,1004.0_rp); call this%checkFound(found,found_aux)
      call json%get("Prt",this%Prt, found,0.71_rp); call this%checkFound(found,found_aux)
      call json%get("v0",this%vo, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("M",this%M, found,0.2_rp); call this%checkFound(found,found_aux)
      call json%get("delta",this%delta, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Retau",this%Retau, found,950.0_rp); call this%checkFound(found,found_aux)
      call json%get("gamma_gas",this%gamma_gas, found,1.4_rp); call this%checkFound(found,found_aux)

      call json%get("flag_rk_ls",flag_rk_ls, found,.true.); 
      call json%get("flag_rk_ls_stages",flag_rk_ls_stages, found,5); 
      call json%get("flag_rk_ls_n",flag_rk_ls_n, found, 1)
      call json%get("c_sgs",c_sgs, found,0.025_rp); 

      call json%get("flag_lps_stab",flag_lps_stab, found,.true.); call this%checkFound(found,found_aux)
      call json%get("flag_use_ducros",flag_use_ducros, found,.false.); call this%checkFound(found,found_aux)
      call json%get("ducros_min_val",ducros_min_val, found,0.1_rp); call this%checkFound(found,found_aux)

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
      this%Re     = exp((1.0_rp/0.88_rp)*log(this%Retau/0.09_rp))
      this%mu    = (this%rho0*2.0_rp*this%delta*this%vo)/this%Re
      this%utau   = (this%Retau*this%mu)/(this%delta*this%rho0)
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = this%mu/mur
      nscbc_rho_inf = this%rho0
      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas
      nscbc_T_C = this%to
      
      call this%readJSONWMTypes()

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'

   end subroutine ChannelFlowSolver_initializeParameters

   subroutine ChannelFlowSolver_evalInitialConditions(this)
      class(ChannelFlowSolver), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, rti(3), yp,velo_aux1
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      logical :: readFiles
      character(512) :: initialField_filePath

      readFiles = .false.

      call nvtxStartRange("ChannelFlow Init")
      if(readFiles) then
         call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
         write(initialField_filePath,*) ""
         call read_veloc_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)

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

            rho(iNodeL,3) = rho(iNodeL,2)
            E(iNodeL,3) =  E(iNodeL,2)
            eta(iNodeL,3) = eta(iNodeL,2) 
            q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
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
               if(coordPar(iNodeL,2)<this%delta) then
                  yp = coordPar(iNodeL,2)*this%utau*this%rho0/this%mu
               else
                  yp = abs(coordPar(iNodeL,2)-2.0_rp*this%delta)*this%utau*this%rho0/this%mu
               end if

               velo = this%utau*((1.0_rp/0.41_rp)*log(1.0_rp+0.41_rp*yp)+7.8_rp*(1.0_rp-exp(-yp/11.0_rp)-(yp/11.0_rp)*exp(-yp/3.0_rp))) 

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
      end if

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

      !
      ! Initialize exponential averaging for wall law 
      !
      call nvtxStartRange("Wall Average init")
      if(flag_walave) then
         !$acc kernels
         walave_u(:,:) = u(:,:,2)
         !$acc end kernels
      end if
      call nvtxEndRange

   end subroutine ChannelFlowSolver_evalInitialConditions

end module ChannelFlowSolver_mod
