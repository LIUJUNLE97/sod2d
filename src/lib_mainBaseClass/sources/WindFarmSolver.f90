
module WindFarmSolver_mod
   use mod_arrays
   use mod_arrays_wf
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: WindFarmSolver

      real(rp) , public  :: rough,vinf,Lhub,rho0,Lz,ustar,wind_alpha, lambda_vertical,gamma_free,delta_capping,gamma_inversion,inversion_height,s_ra,nu_ra,latitude,earth_omega,fc,gradP,Ug_x,Ug_y,T_wall,Lz_ra,Ug_alpha
      integer(4), public :: N_ad
      logical, public :: capping

   contains
      procedure, public :: fillBCTypes           => WindFarmSolver_fill_BC_Types
      procedure, public :: initializeParameters  => WindFarmSolver_initializeParameters
      procedure, public :: initializeSourceTerms => WindFarmSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => WindFarmSolver_evalInitialConditions
      procedure, public :: afterDt => WindFarmSolver_afterDt
   end type WindFarmSolver
contains

   subroutine WindFarmSolver_afterDt(this,istep)
      class(WindFarmSolver), intent(inout) :: this
      integer(4), intent(in) :: istep
      integer :: iNodeL,idime,iAD
      real(rp) :: x_ad,y_ad,z_ad,radius
      real(8) :: vol_T(this%N_ad),vol_T2(this%N_ad)
      real(rp) :: vz,N,temp


      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,3) = -rho(iNodeL,2)*this%Ug_y*this%fc - this%fc*q(iNodeL,2,2)
         source_term(iNodeL,4) = rho(iNodeL,2)*this%Ug_x*this%fc + this%fc*q(iNodeL,1,2)
         source_term(iNodeL,5) = 0.0_rp 
      end do
      !$acc end parallel loop


   end subroutine WindFarmSolver_afterDt


   subroutine WindFarmSolver_fill_BC_Types(this)
      class(WindFarmSolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine WindFarmSolver_fill_BC_Types

   subroutine WindFarmSolver_initializeSourceTerms(this)
      class(WindFarmSolver), intent(inout) :: this
      integer(4) :: iNodeL,idime
      real(rp) :: vz,N


      allocate(source_term(numNodesRankPar,ndime+2))
      !$acc enter data create(source_term(:,:))


      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar         
         source_term(iNodeL,3) = -rho(iNodeL,2)*this%Ug_y*this%fc - this%fc*q(iNodeL,2,2)
         source_term(iNodeL,4) =  rho(iNodeL,2)*this%Ug_x*this%fc + this%fc*q(iNodeL,1,2)
         source_term(iNodeL,5) =  0.0_rp 
      end do
      !$acc end parallel loop

   end subroutine WindFarmSolver_initializeSourceTerms

   subroutine WindFarmSolver_initializeParameters(this)
      use json_module
      implicit none
      class(WindFarmSolver), intent(inout) :: this
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

      call json%get("saveInitialField",this%saveInitialField, found,.true.); call this%checkFound(found,found_aux)

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

      call json%get("doTimerAnalysis",this%doTimerAnalysis, found,.false.)

      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("flag_implicit",flag_implicit, found,0); call this%checkFound(found,found_aux)
      call json%get("flag_imex_stages",flag_imex_stages, found,4); call this%checkFound(found,found_aux)
      call json%get("flag_walave",flag_walave, found,.true.); call this%checkFound(found,found_aux)
      call json%get("period_walave",period_walave, found,3600.0_rp); call this%checkFound(found,found_aux)
      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,300.0_rp); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.95_rp); call this%checkFound(found,found_aux)

      
      call json%get("Lz",this%Lz, found,440.0_rp); call this%checkFound(found,found_aux)
      call json%get("Lhub",this%Lhub, found,90.0_rp); call this%checkFound(found,found_aux)
      call json%get("vinf",this%vinf, found,8.0_rp); call this%checkFound(found,found_aux)
      call json%get("rough",this%rough, found,0.1_rp); call this%checkFound(found,found_aux)
      call json%get("wind_alpha",this%wind_alpha, found,270.0_rp); call this%checkFound(found,found_aux)
      flag_high_mach = .true.
      flag_bouyancy_effect = .true.

      this%wind_alpha = 270.0_rp-this%wind_alpha !Comming North is 0 , East is 90, South is 180 and West is 270 in a x-y axis
      this%Ug_alpha = this%wind_alpha !+ 20.0_rp


      ! fixed by the type of base class parameters
      flag_type_wmles = wmles_type_abl
      nscbc_u_inf = this%vinf
      flag_mu_factor = 1.0_rp
      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%gamma_gas = 1.40_rp
      nscbc_gamma_inf = this%gamma_gas
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      nscbc_Rgas_inf = this%Rgas
      nscbc_p_inf = 101325.0_rp 
      this%T_wall = 265.0_rp
      nscbc_T_ref = 265.0_rp 
      nscbc_rho_inf = nscbc_p_inf/(this%Rgas*this%T_wall)
      nscbc_T_C = nscbc_T_ref
      nscbc_g_x = 0.0_rp
      nscbc_g_y = 0.0_rp
      nscbc_g_z = 9.81_rp
      

      this%gamma_free = 0.005_rp
      this%delta_capping = 100.0_rp
      this%gamma_inversion = 0.01_rp
      this%inversion_height = 200.0_rp
      this%lambda_vertical = 7000.0_rp
      this%s_ra = 1.5_rp
      this%nu_ra = 3.0_rp
      this%latitude = 40.0_rp*v_pi/180.0_rp
      this%earth_omega = 0.00007272_rp
               
      !extra calc
      this%fc =  0.000139_rp !2.0_rp*this%earth_omega*sin(this%latitude)
      this%gradP = 2.0_rp*this%vinf*this%fc ! 10 is an input: gradP equal to mag(Ug)/fc                  
      this%Ug_x = this%vinf
      this%Ug_y = 0.0_rp
      this%Lz_ra = this%lambda_vertical*1.5_rp

      if(mpi_rank.eq.0) write(*,*) "--| gradP :", this%gradP
      if(mpi_rank.eq.0) write(*,*) "--| Ugx :", this%Ug_x
      if(mpi_rank.eq.0) write(*,*) "--| Ugy :", this%Ug_y
      if(mpi_rank.eq.0) write(*,*) "--| fc :", this%fc
      if(mpi_rank.eq.0) write(*,*) "--| Tref :", nscbc_T_ref
      if(mpi_rank.eq.0) write(*,*) "--| gz :", nscbc_g_z
      if(mpi_rank.eq.0) write(*,*) "--| rho :", nscbc_rho_inf
      if(mpi_rank.eq.0) write(*,*) "--| R :", this%Rgas

      call this%readJSONBuffer()

   end subroutine WindFarmSolver_initializeParameters

   subroutine WindFarmSolver_evalInitialConditions(this)
      class(WindFarmSolver), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, rti(3), zp,velo_aux1, veloMatias(2), ugMatias, theta, p100, gcp
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      character(512) :: initialField_filePath

      call nvtxStartRange("WindFarm Init")

      gcp = nscbc_g_z/this%Cp

      call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
      auxCnt = 1
      serialLoop : do iLine = 1,totalNumNodesSrl
         call random_number(rti)
         if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
            iNodeL = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1

            zp = coordPar(iNodeL,3)
            
            ! GABLS
            if (zp.lt.(100.0_rp)) then 
               u(iNodeL,1,2) = this%vinf*(1.0_rp + 0.01_rp*(rti(1) -0.5_rp))
            else
               u(iNodeL,1,2) = this%vinf
            end if
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = 0.0_rp
            
            p100 =  nscbc_p_inf*((1.0_rp-gcp*100.0_rp/this%T_wall)**(this%gamma_gas/(this%gamma_gas-1.0_rp)))
            
            ! GABLS1
            if (zp.lt.(100.0_rp)) then ! capping inversion region of 100m
               theta =  this%T_wall
               pr(iNodeL,2) = nscbc_p_inf*((1.0_rp-gcp*zp/this%T_wall)**(this%gamma_gas/(this%gamma_gas-1.0_rp)))    
               Tem(iNodeL,2) =  this%T_wall - gcp*zp           
            else
               theta =  this%T_wall   -  (zp-100.0_rp)*0.01_rp      
               pr(iNodeL,2) = p100*((1.0_rp-((0.01_rp+gcp)/(this%T_wall-gcp*100.0_rp))*(zp-100.0_rp))**(this%gamma_gas/(this%gamma_gas-1.0_rp)))
               Tem(iNodeL,2) =  this%T_wall - (0.01_rp+gcp)*(zp-100.0_rp)  - gcp*100.0_rp
            end if
         end if
         if(auxCnt.gt.numNodesRankPar) then
            exit serialLoop
         end if
      end do serialLoop

      !$acc update device(u(:,:,:))
      !$acc update device(Tem(:,:))
      !$acc update device(pr(:,:))

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         
         rho(iNodeL,2) = pr(iNodeL,2)/this%Rgas/Tem(inodeL,2)
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) =  eta(iNodeL,2)
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


   end subroutine WindFarmSolver_evalInitialConditions

end module WindFarmSolver_mod
