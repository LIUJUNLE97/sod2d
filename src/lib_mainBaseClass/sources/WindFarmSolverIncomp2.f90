
module WindFarmSolverIncomp2_mod
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
   use CFDSolverPeriodicWithBoundariesIncomp_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: WindFarmSolverIncomp2

      real(rp) , public  :: rough,vinf,Lhub,rho0,Lz,ustar,wind_alpha, lambda_vertical,gamma_free,delta_capping,gamma_inversion,inversion_height,s_ra,nu_ra,latitude,earth_omega,fc,gradP,Ug_x,Ug_y,T_wall,Lz_ra,Ug_alpha
      integer(4), public :: N_ad
      logical, public :: capping

   contains
      procedure, public :: fillBCTypes           => WindFarmSolverIncomp2_fill_BC_Types
      procedure, public :: initializeParameters  => WindFarmSolverIncomp2_initializeParameters
      procedure, public :: initializeSourceTerms => WindFarmSolverIncomp2_initializeSourceTerms
      procedure, public :: evalInitialConditions => WindFarmSolverIncomp2_evalInitialConditions
      procedure, public :: afterDt => WindFarmSolverIncomp2_afterDt
      procedure, public :: initialBuffer => WindFarmSolverIncomp2_initialBuffer
   end type WindFarmSolverIncomp2
contains

   subroutine WindFarmSolverIncomp2_initialBuffer(this)
      class(WindFarmSolverIncomp2), intent(inout) :: this
      integer :: iNodeL

      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         if (coordPar(iNodeL,3).lt.(100.0_rp)) then ! capping inversion region of 100m
            Yk_buffer(iNodeL,1) =  this%T_wall
         else
            Yk_buffer(iNodeL,1) =   this%T_wall  +  (coordPar(iNodeL,3)-100.0_rp)*0.01_rp 
         end if
         u_buffer(iNodeL,1)  = this%vinf
         u_buffer(iNodeL,2)  = 0.0_rp
         u_buffer(iNodeL,3)  = 0.0_rp
      end do
      !$acc end parallel loop
   end subroutine WindFarmSolverIncomp2_initialBuffer


   subroutine WindFarmSolverIncomp2_afterDt(this,istep)
      class(WindFarmSolverIncomp2), intent(inout) :: this
      integer(4), intent(in) :: istep
      integer :: iNodeL


      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) = -rho(iNodeL,2)*this%Ug_y*this%fc + this%fc*rho(iNodeL,2)*u(iNodeL,2,2)
         source_term(iNodeL,2) = rho(iNodeL,2)*this%Ug_x*this%fc - this%fc*rho(iNodeL,2)*u(iNodeL,1,2)
         source_term(iNodeL,3) = -rho(iNodeL,2)*((Yk(iNodeL,1,2)-nscbc_T_ref)/265.0_rp)*nscbc_g_z 

         if(coordPar(iNodeL,3) .lt. 100_rp) then
            Yk_buffer(iNodeL,1) = this%T_wall - 0.25_rp*(this%time/3600.0_rp)
         end if
      end do
      !$acc end parallel loop
   end subroutine WindFarmSolverIncomp2_afterDt


   subroutine WindFarmSolverIncomp2_fill_BC_Types(this)
      class(WindFarmSolverIncomp2), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine WindFarmSolverIncomp2_fill_BC_Types

   subroutine WindFarmSolverIncomp2_initializeSourceTerms(this)
      class(WindFarmSolverIncomp2), intent(inout) :: this
      integer(4) :: iNodeL,idime
      real(rp) :: vz,N


      allocate(source_term(numNodesRankPar,ndime))
      !$acc enter data create(source_term(:,:))


      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar         
         source_term(iNodeL,1) = -rho(iNodeL,2)*this%Ug_y*this%fc + this%fc*rho(iNodeL,2)*u(iNodeL,2,2)
         source_term(iNodeL,2) =  rho(iNodeL,2)*this%Ug_x*this%fc - this%fc*rho(iNodeL,2)*u(iNodeL,1,2)
         source_term(iNodeL,3) =  -nscbc_rho_inf*((Yk(iNodeL,1,2)-nscbc_T_ref)/265.0_rp)*nscbc_g_z 
      end do
      !$acc end parallel loop

   end subroutine WindFarmSolverIncomp2_initializeSourceTerms

   subroutine WindFarmSolverIncomp2_initializeParameters(this)
      use json_module
      implicit none
      class(WindFarmSolverIncomp2), intent(inout) :: this
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
      call json%get("flag_walave",flag_walave, found,.true.); call this%checkFound(found,found_aux)
      call json%get("period_walave",period_walave, found,3600.0_rp); call this%checkFound(found,found_aux)
      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,300.0_rp); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)

      
      call json%get("Lz",this%Lz, found,440.0_rp); call this%checkFound(found,found_aux)
      call json%get("Lhub",this%Lhub, found,90.0_rp); call this%checkFound(found,found_aux)
      call json%get("vinf",this%vinf, found,8.0_rp); call this%checkFound(found,found_aux)
      call json%get("rough",this%rough, found,0.1_rp); call this%checkFound(found,found_aux)
      call json%get("wind_alpha",this%wind_alpha, found,270.0_rp); call this%checkFound(found,found_aux)      

      this%maxPhysTime = 9.0_rp*3600.0_rp

      this%wind_alpha = 270.0_rp-this%wind_alpha !Comming North is 0 , East is 90, South is 180 and West is 270 in a x-y axis
      this%Ug_alpha = this%wind_alpha !+ 20.0_rp


      ! fixed by the type of base class parameters
      flag_use_species = .true.
      nspecies = 1

      flag_type_wmles = wmles_type_abl
      nscbc_u_inf = this%vinf
      incomp_viscosity = 1.81e-5
      flag_mu_factor = 1.0_rp
      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%T_wall = 1.0_rp !265.0_rp
      nscbc_T_ref = 1.0_rp !265.0_rp 
      nscbc_rho_inf = 1.0_rp
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

   end subroutine WindFarmSolverIncomp2_initializeParameters

   subroutine WindFarmSolverIncomp2_evalInitialConditions(this)
      class(WindFarmSolverIncomp2), intent(inout) :: this
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
               u(iNodeL,1,2) = this%vinf
            else
               u(iNodeL,1,2) = this%vinf
            end if
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = 0.0_rp            
            
            ! GABLS1
            if (zp.lt.(100.0_rp)) then ! capping inversion region of 100m
               theta =  this%T_wall*(1.0_rp + 0.001_rp*(rti(1) -0.5_rp))
               Yk(iNodeL,1,2) =  theta
            else
               theta =  this%T_wall  +  (zp-100.0_rp)*0.01_rp      
               Yk(iNodeL,1,2) =  theta
            end if
         end if
         if(auxCnt.gt.numNodesRankPar) then
            exit serialLoop
         end if
      end do serialLoop

      !$acc update device(u(:,:,:))
      !$acc update device(Yk(:,:,:))

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = 0.0_rp
         rho(iNodeL,2) = nscbc_rho_inf
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


   end subroutine WindFarmSolverIncomp2_evalInitialConditions

end module WindFarmSolverIncomp2_mod
