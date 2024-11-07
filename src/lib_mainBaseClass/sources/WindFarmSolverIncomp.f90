module mod_arrays_wf
   use mod_constants

   implicit none

   real(rp), allocatable :: CT(:), D(:),pos_x(:),pos_y(:),pos_z(:),vol_correct(:),ad_alpha(:),static_sources(:,:)

end module mod_arrays_wf

module WindFarmSolverIncomp_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: WindFarmSolverIncomp

      real(rp) , public  :: rough,vinf,Lhub,rho0,Lz,ustar,wind_alpha, lambda_vertical,gamma_free,delta_capping,gamma_inversion,inversion_height,s_ra,nu_ra,latitude,earth_omega,fc,gradP,Ug_x,Ug_y,T_wall,Lz_ra,Ug_alpha
      integer(4), public :: N_ad
      logical, public :: capping

   contains
      procedure, public :: fillBCTypes           => WindFarmSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => WindFarmSolverIncomp_initializeParameters
      procedure, public :: initializeSourceTerms => WindFarmSolverIncomp_initializeSourceTerms
      procedure, public :: evalInitialConditions => WindFarmSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer         => WindFarmSolverIncomp_initialBuffer
      procedure, public :: eval_vars_after_load_hdf5_resultsFile => WindFarmSolverIncomp_eval_vars_after_load_hdf5_resultsFile 
      procedure, public :: afterDt => WindFarmSolverIncomp_afterDt
      procedure, public :: readJSONAD => WindFarmSolverIncomp_readJSONAD
   end type WindFarmSolverIncomp
contains

subroutine WindFarmSolverIncomp_readJSONAD(this)
   use json_module
   implicit none
   class(WindFarmSolverIncomp), intent(inout) :: this
   logical :: found, found_aux = .false.
   type(json_file) :: json
   integer :: json_nad,iAD,id
   TYPE(json_core) :: jCore
   TYPE(json_value), pointer :: buffPointer, testPointer, p
   character(len=:) , allocatable :: value

   call json%initialize()
   call json%load_file(json_filename)

   call json%info("AD",n_children=json_nad)
   call json%get_core(jCore)
   call json%get('AD', buffPointer, found_aux)

   if(json_nad .gt. 0) then

      this%N_ad =  json_nad 

      allocate(CT(this%N_ad))
      allocate(D(this%N_ad))
      allocate(pos_x(this%N_ad))
      allocate(pos_y(this%N_ad))
      allocate(pos_z(this%N_ad))
      allocate(vol_correct(this%N_ad))
      allocate(ad_alpha(this%N_ad))
      !$acc enter data create(CT(:),D(:),pos_x(:),pos_y(:),pos_z(:),vol_correct(:),ad_alpha(:))

      
      do iAD=1, this%N_ad
         call jCore%get_child(buffPointer, iAD, testPointer, found)
         
         call jCore%get_child(testPointer, 'id', p, found)
         if(found) then
            call jCore%get(p,id)
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD id'
               stop 1      
            end if
         end if

         call jCore%get_child(testPointer, 'CT', p, found)
         if(found) then
            call jCore%get(p,CT(id))
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD CT'
               stop 1      
            end if
         end if

         call jCore%get_child(testPointer, 'D', p, found)
         if(found) then
            call jCore%get(p,D(id))
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD D'
               stop 1      
            end if
         end if

         call jCore%get_child(testPointer, 'x', p, found)
         if(found) then
            call jCore%get(p,pos_x(id))
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD x'
               stop 1      
            end if
         end if

         call jCore%get_child(testPointer, 'y', p, found)
         if(found) then
            call jCore%get(p,pos_y(id))
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD y'
               stop 1      
            end if
         end if

         call jCore%get_child(testPointer, 'z', p, found)
         if(found) then
            call jCore%get(p,pos_z(id))
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD z'
               stop 1      
            end if
         end if

         call jCore%get_child(testPointer, 'ad_alpha', p, found)
         if(found) then
            call jCore%get(p,ad_alpha(id))
            ad_alpha(id) = ad_alpha(id)-90.0_rp 
         else
            if(mpi_rank .eq. 0) then
               write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define the AD ad_alpha'
               stop 1      
            end if
         end if
      end do
   end if

   !$acc update device(CT(:),D(:),pos_x(:),pos_y(:),pos_z(:),ad_alpha(:))

   call json%destroy()

end subroutine WindFarmSolverIncomp_readJSONAD

   subroutine WindFarmSolverIncomp_afterDt(this,istep)
      class(WindFarmSolverIncomp), intent(inout) :: this
      integer(4), intent(in) :: istep
      integer :: iNodeL,idime,iAD
      real(rp) :: x_ad,y_ad,z_ad,radius
      real(8) :: vol_T(this%N_ad),vol_T2(this%N_ad)
      real(rp) :: vz,N

      if(istep == 1) then !esto va a dar problemas en continues
         allocate(static_sources(numNodesRankPar,3))
         !$acc enter data create(static_sources(:,:))
         if(this%N_ad .gt. 0) then
            vol_T(:) = 0.0d0
            !! Wind farm pre-processing
            !$acc parallel loop reduction(+:vol_T)
            do iNodeL = 1,numNodesRankPar
               x_ad =  coordPar(iNodeL,1)*cos(ad_alpha(iAD)*v_pi/180.0_rp)+coordPar(iNodeL,2)*sin(ad_alpha(iAD)*v_pi/180.0_rp)
               y_ad = -coordPar(iNodeL,1)*sin(ad_alpha(iAD)*v_pi/180.0_rp)+coordPar(iNodeL,2)*cos(ad_alpha(iAD)*v_pi/180.0_rp)
               z_ad =  coordPar(iNodeL,3)
               
               ad(iNodeL) = 0

               !$acc loop seq
               do iAD=1,this%N_ad
                  if((x_ad .gt. (pos_x(iAD)-(0.05_rp*D(iAD)))) .and. (x_ad .le. (pos_x(iAD)+(0.05_rp*D(iAD))))) then
                     radius = sqrt( (y_ad-pos_y(iAD))**2 + (z_ad-pos_z(iAD))**2 )
                     if(radius .le. (D(iAD)*0.5_rp)) then 
                        vol_T(iAD) = vol_T(iAD) + real(Ml(iNodeL),8)
                        ad(iNodeL) = iAD
                     end if
                  end if
               end do
            end do
            !$acc end parallel loop

            do iAD=1,this%N_ad            
               call MPI_Allreduce(vol_T(iAD),vol_T2(iAD),1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
            end do

            !$acc loop seq
            do iAD=1,this%N_ad
               vol_correct(iAD) = real(vol_T2(iAD),rp)/(v_pi*((D(iAD)*0.5_rp)**2)*(D(iAD)*0.1_rp))
               if(mpi_rank.eq.0)  write(111,*) " volT",vol_T2(iAD)," corr ",vol_correct(iAD)
               if(mpi_rank.eq.0)  write(111,*) " s prec",(this%rho0*this%ustar**2/this%Lz)," s ad ",0.5_rp*(CT(iAD)/(0.1_rp*D(iAD)))*(this%vinf**2)
            end do

            !$acc update device(vol_correct(:))

            !$acc parallel loop  
            do iNodeL = 1,numNodesRankPar
               source_term(iNodeL,1) = -this%Ug_y*this%fc
               source_term(iNodeL,2) = +this%Ug_x*this%fc
               source_term(iNodeL,3) = 0.00_rp
               if(ad(iNodeL) .ne. 0) then
                  iAD = ad(iNodeL)
                  source_term(iNodeL,1) = -0.5_rp*(CT(iAD)/(0.1_rp*D(iAD)))*(this%vinf**2)*vol_correct(iAD)*cos(ad_alpha(iAD)*v_pi/180.0_rp)
                  source_term(iNodeL,2) = -0.5_rp*(CT(iAD)/(0.1_rp*D(iAD)))*(this%vinf**2)*vol_correct(iAD)*sin(ad_alpha(iAD)*v_pi/180.0_rp)
                  source_term(iNodeL,3) = 0.0_rp
               end if
               static_sources(iNodeL,1) = source_term(iNodeL,1) 
               static_sources(iNodeL,2) = source_term(iNodeL,2) 
               static_sources(iNodeL,3) = source_term(iNodeL,3) 
            end do
            !$acc end parallel loop
         end if
      end if
      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) = static_sources(iNodeL,1) 
         source_term(iNodeL,2) = static_sources(iNodeL,2) 
         source_term(iNodeL,3) = static_sources(iNodeL,3) 
         if(this%capping .eqv. .true.) then
            source_term(iNodeL,1) = source_term(iNodeL,1) - this%fc*u(iNodeL,2,2)
            source_term(iNodeL,2) = source_term(iNodeL,2) + this%fc*u(iNodeL,1,2)
            source_term(iNodeL,3) = source_term(iNodeL,3) + nscbc_rho_inf*((Yk(iNodeL,1,2)-nscbc_T_ref)/nscbc_T_ref)*nscbc_g_z

            !if(coordPar(iNodeL,3) .gt. (this%Lz-this%Lz_ra)) then
            !   N = sqrt(nscbc_g_z*this%gamma_free/nscbc_T_ref)
            !   vz = this%nu_ra*N*(1.0_rp - cos((v_pi/this%s_ra)*((coordPar(iNodeL,3)-(this%Lz-this%Lz_ra))/this%Lz_ra)) )
               
            !   source_term(iNodeL,1) = source_term(iNodeL,1) + vz*(u(iNodeL,1,2)-this%Ug_x)  
            !   source_term(iNodeL,2) = source_term(iNodeL,2) + vz*(u(iNodeL,2,2)-this%Ug_y)  
            !   source_term(iNodeL,3) = source_term(iNodeL,2) + vz*(u(iNodeL,3,2))
            !end if
         end if
      end do
      !$acc end parallel loop

   end subroutine WindFarmSolverIncomp_afterDt

   subroutine  WindFarmSolverIncomp_initialBuffer(this)
      class(WindFarmSolverIncomp), intent(inout) :: this
      integer(4) :: iNode
      real(rp) :: velo,zp, veloMatias(2), ugMatias
      

      !$acc parallel loop
      do iNode = 1,numNodesRankPar         
         zp = coordPar(iNode,3)
         if(this%capping .eqv. .true.) then
            ugMatias = sqrt(this%Ug_x**2+this%Ug_y**2)
            !call veloc_abl_new(veloMatias,zp,this%inversion_height,ugMatias, this%Ug_alpha, this%rough)

            !u_buffer(iNode,1) = veloMatias(1)
            !u_buffer(iNode,2) = veloMatias(2)
            !u_buffer(iNode,3) = 0.0_rp

            u_buffer(iNode,1) = this%vinf
            u_buffer(iNode,2) = 0.0_rp
            u_buffer(iNode,3) = 0.0_rp
         else
            velo =  this%ustar*log(1.0_rp+zp/this%rough)/0.41_rp

            u_buffer(iNode,1) = velo*cos(this%wind_alpha*v_pi/180.0_rp)
            u_buffer(iNode,2) = velo*sin(this%wind_alpha*v_pi/180.0_rp)
            u_buffer(iNode,3) = 0.0_rp   
         end if
         
         if(this%capping .eqv. .true.) then
            !if(zp .lt. this%inversion_height) then
            !   Yk_buffer(iNode,1) = this%T_wall
            !else
              !if (zp.lt.(this%inversion_height+this%delta_capping)) then ! capping inversion region of 100m
               !   Yk_buffer(iNode,1) =  this%T_wall +  this%gamma_inversion*(zp-this%inversion_height)
               !else
               !   Yk_buffer(iNode,1) =  this%T_wall   + &
               !                         this%gamma_inversion*this%delta_capping + &
               !                         this%gamma_free*(zp-(this%inversion_height+this%delta_capping))
               !end if               
            !end if
            if (zp.lt.(100.0_rp)) then ! capping inversion region of 100m
               Yk_buffer(iNode,1) =  this%T_wall
            else
               Yk_buffer(iNode,1) =  this%T_wall   + &
                                       (zp-100.0_rp)*0.01_rp
            end if
         end if
      end do
      !$acc end parallel loop

   end subroutine  WindFarmSolverIncomp_initialBuffer

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
      integer(4) :: iNodeL,idime
      real(rp) :: vz,N


      allocate(source_term(numNodesRankPar,ndime))
      !$acc enter data create(source_term(:,:))


      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar         
         source_term(iNodeL,1) = -this%Ug_y*this%fc
         source_term(iNodeL,2) = +this%Ug_x*this%fc
         source_term(iNodeL,3) = 0.00_rp
         
         if(this%capping .eqv. .true.) then
               source_term(iNodeL,1) = source_term(iNodeL,1) - this%fc*u(iNodeL,2,2)
               source_term(iNodeL,2) = source_term(iNodeL,2) + this%fc*u(iNodeL,1,2)
               source_term(iNodeL,3) = source_term(iNodeL,3) + nscbc_rho_inf*((Yk(iNodeL,1,2)-nscbc_T_ref)/nscbc_T_ref)*nscbc_g_z
               
               !if(coordPar(iNodeL,3) .gt. (this%Lz-this%Lz_ra)) then
               !   N = sqrt(nscbc_g_z*this%gamma_free/nscbc_T_ref)
               !   vz = this%nu_ra*N*(1.0_rp - cos((v_pi/this%s_ra)*((coordPar(iNodeL,3)-(this%Lz-this%Lz_ra))/this%Lz_ra)) )
            
               !   source_term(iNodeL,1) = source_term(iNodeL,1) + vz*(u(iNodeL,1,2) - this%Ug_x ) 
               !   source_term(iNodeL,2) = source_term(iNodeL,2) + vz*(u(iNodeL,2,2) - this%Ug_y ) 
               !   source_term(iNodeL,3) = source_term(iNodeL,3) + vz*(u(iNodeL,3,2)) 
               !end if
         end if
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
      call json%get("flag_fs_fix_pressure",flag_fs_fix_pressure, found,.false.); call this%checkFound(found,found_aux)
      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,300.0_rp); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Lz",this%Lz, found,1500.0_rp); call this%checkFound(found,found_aux)
      call json%get("Lhub",this%Lhub, found,90.0_rp); call this%checkFound(found,found_aux)
      call json%get("vinf",this%vinf, found,8.0_rp); call this%checkFound(found,found_aux)
      call json%get("rough",this%rough, found,0.1682_rp); call this%checkFound(found,found_aux)
      call json%get("wind_alpha",this%wind_alpha, found,270.0_rp); call this%checkFound(found,found_aux)

      this%wind_alpha = 270.0_rp-this%wind_alpha !Comming North is 0 , East is 90, South is 180 and West is 270 in a x-y axis
      this%Ug_alpha = this%wind_alpha !+ 20.0_rp

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
      flag_type_wmles = wmles_type_abl
      nscbc_p_inf = 0.0_rp
      nscbc_u_inf = this%vinf
      nscbc_rho_inf = this%rho0

      this%ustar = this%vinf*0.41_rp/log(1.0_rp+this%Lhub/this%rough)
      incomp_viscosity = 1.81e-5
      flag_mu_factor = 1.0_rp

      call json%get("capping_inversion",this%capping, found,.false.)
      if(this%capping .eqv. .true.) then
         flag_use_species = .true.
         nspecies = 1
         this%Prt=0.71_rp
         this%Cp=1004.0_rp
         this%gamma_free = 0.005_rp
         this%delta_capping = 100.0_rp
         this%gamma_inversion = 0.01_rp
         this%inversion_height = 200.0_rp
         this%lambda_vertical = 7000.0_rp
         this%s_ra = 1.5_rp
         this%nu_ra = 3.0_rp
         this%latitude = 40.0_rp*v_pi/180.0_rp
         this%earth_omega = 0.00007272_rp
         this%T_wall = 265.0_rp
         nscbc_T_ref = 0.5_rp*(2.0_rp*265.0_rp+340.0_rp*0.01_rp) !this%T_wall!+0.5_rp*this%gamma_inversion*this%delta_capping

         !extra calc
         this%fc =  0.000139_rp !2.0_rp*this%earth_omega*sin(this%latitude)
         this%gradP = 2.0_rp*this%vinf*this%fc ! 10 is an input: gradP equal to mag(Ug)/fc
         
         
         !this%Ug_x = this%gradP*sin(this%wind_alpha*v_pi/180.0_rp)/this%fc
         !this%Ug_y = -this%gradP*cos(this%wind_alpha*v_pi/180.0_rp)/this%fc
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
      end if

      call this%readJSONBuffer()
      call this%readJSONAD()

   end subroutine WindFarmSolverIncomp_initializeParameters

   subroutine WindFarmSolverIncomp_evalInitialConditions(this)
      class(WindFarmSolverIncomp), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime
      real(rp) :: velo, rti(3), zp,velo_aux1, veloMatias(2), ugMatias
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      logical :: readFiles
      character(512) :: initialField_filePath

      readFiles = .false.

      call nvtxStartRange("WindFarm_incomp Init")
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
               if(this%capping .eqv. .true.) then
                  ugMatias = sqrt(this%Ug_x**2+this%Ug_y**2)
                  !call veloc_abl_new(veloMatias,zp,this%inversion_height,ugMatias, this%Ug_alpha, this%rough)
                  

                  !if(zp .lt. this%inversion_height) then
                  !   u(iNodeL,1,2) = veloMatias(1)*(1.0_rp + 0.05_rp*(rti(1) -0.5_rp))
                  !   u(iNodeL,2,2) = veloMatias(2)*(1.0_rp + 0.05_rp*(rti(2) -0.5_rp))
                  !   u(iNodeL,3,2) = sqrt(veloMatias(1)**2+veloMatias(2)**2)*(0.05_rp*(rti(3) -0.5_rp))
                  !else
                     !u(iNodeL,1,2) = veloMatias(1)
                     !u(iNodeL,2,2) = veloMatias(2)
                     !u(iNodeL,3,2) = 0.0_rp
                  !end if

                     ! GABLS
                     if (zp.lt.(100.0_rp)) then 
                        u(iNodeL,1,2) = this%vinf!*(1.0_rp + 0.01_rp*(rti(1) -0.5_rp))
                     else
                        u(iNodeL,1,2) = this%vinf
                     end if
                     u(iNodeL,2,2) = 0.0_rp
                     u(iNodeL,3,2) = 0.0_rp
               else
                  velo =  this%ustar*log(1.0_rp+zp/this%rough)/0.41_rp

                  u(iNodeL,1,2) = velo*(1.0_rp + 0.05_rp*(rti(1) -0.5_rp))*cos(this%wind_alpha*v_pi/180.0_rp)
                  u(iNodeL,2,2) = velo*( 0.05_rp*(rti(2) -0.5_rp))*sin(this%wind_alpha*v_pi/180.0_rp)
                  u(iNodeL,3,2) = velo*(0.05_rp*(rti(3) -0.5_rp))
               end if

               
               
               if(this%capping .eqv. .true.) then
                  !if(zp .lt. this%inversion_height) then
                  !   Yk(iNodeL,1,2) = this%T_wall
                  !else
                     !if (zp.lt.(this%inversion_height+this%delta_capping)) then ! capping inversion region of 100m
                     !   Yk(iNodeL,1,2) =  this%T_wall +  this%gamma_inversion*(zp-this%inversion_height)
                     !else
                     !   Yk(iNodeL,1,2) =  this%T_wall   + &
                     !                     this%gamma_inversion*this%delta_capping + &
                     !                     this%gamma_free*(zp-(this%inversion_height+this%delta_capping))
                     !end if

                     ! GABLS1
                     if (zp.lt.(100.0_rp)) then ! capping inversion region of 100m
                        Yk(iNodeL,1,2) =  this%T_wall
                     else
                        Yk(iNodeL,1,2) =  this%T_wall   + &
                                          (zp-100.0_rp)*0.01_rp
                     end if
                  !end if
                  Yk(iNodeL,1,1) =  Yk(iNodeL,1,2)
                  Yk(iNodeL,1,3) =  Yk(iNodeL,1,2)
                  Yk(iNodeL,1,4) =  Yk(iNodeL,1,2)
               end if
            end if
            if(auxCnt.gt.numNodesRankPar) then
               exit serialLoop
            end if
         end do serialLoop
         !!!!$acc end parallel loop

         !$acc update device(u(:,:,:))
         !$acc update device(Yk(:,:,:))
         !$acc update device(eta_Yk(:,:,:))

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            eta(iNodeL,2) =0.5_rp*dot_product(u(iNodeL,1:ndime,2),u(iNodeL,1:ndime,2))
            pr(iNodeL,2) = 0.0_rp 
            rho(iNodeL,2) = this%rho0     
            
            u(iNodeL,1:ndime,3) = u(iNodeL,1:ndime,2)
            rho(iNodeL,3) = rho(iNodeL,2)
            eta(iNodeL,3) =  eta(iNodeL,2)
            pr(iNodeL,3) =  pr(iNodeL,2)  

            u(iNodeL,1:ndime,4) = u(iNodeL,1:ndime,2)
            rho(iNodeL,4) = rho(iNodeL,2)
            eta(iNodeL,4) =  eta(iNodeL,2)
            pr(iNodeL,4) =  pr(iNodeL,2) 
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

   subroutine veloc_abl(velabl,y,hpbl,vegeo, rough)
      ! input : 
      !    y : height over terrain
      !    hpbl : height of PBL (capping inversion is there)
      !    vegeo(2):   Geostrophic wind speed
      !    rough: Roughness
      ! output velabl
      !returns velocity profile for ABL, such geostrophic velocity is in the initial condition
      implicit none
      real(rp), intent(in) ::  y, hpbl,vegeo(2), rough
      real(rp), intent(out)::  velabl(2)
      real(rp)  ux, uy, xi, delta, theta
      real(rp)  fu, fv, uxd,uyd,uxdr,uydr
      real(rp)  blend1, blend2,ustar

      ! rotation angle
      theta = 15.0_rp*v_pi/180.0_rp
      ! ustar
      ustar = sqrt(vegeo(1)*vegeo(1) + vegeo(2)*vegeo(2))*1.1_rp*0.41_rp/log(1.0 + 0.5_rp*hpbl/rough)
      ! 
      delta =0.4_rp*hpbl ! mergin region width
      xi = y/hpbl
      fu = 1.57_rp*xi -2.68_rp*xi*xi
      fv = 13.2_rp*xi -8.7_rp*xi*xi

      uxd = ustar/0.41_rp*( log(1.0_rp+y/rough)+fu)
      uyd = - ustar/0.41_rp*fv


      !rotation
      uxdr = uxd*cos(theta) - uyd*sin(theta)
      uydr = uxd*sin(theta) + uyd*cos(theta)

      blend1 = 0.5_rp*(1.0_rp-tanh((xi-0.7_rp)*2.0_rp*hpbl/delta)) ! between 1 and 0 
      blend2 = 0.5_rp*(1.0_rp+tanh((xi-0.7_rp)*2.0_rp*hpbl/delta)) ! between 0 and 1

      ux = uxdr*blend1 + vegeo(1)*blend2
      uy = uydr*blend1 + vegeo(2)*blend2

      velabl(1) = ux
      velabl(2) = uy

   end subroutine

   subroutine veloc_abl_new(velabl,y,hpbl,Mod_Ugeo, Dir_Geos,  rough)
      ! input : 
      !    y : height over terrain
      !    hpbl : height of PBL (capping inversion is there)
      !    vegeo(2):   Geostrophic wind speed
      !    rough: Roughness
      ! output velabl
      !returns velocity profile for ABL, such geostrophic velocity is in the initial condition
      implicit none
      real(rp), intent(in) ::  y, hpbl,mod_ugeo, dir_geos, rough
      real(rp), intent(out)::  velabl(2)
      real(rp)  ux, uy, xi, delta, theta
      real(rp)  fu, fv, uxd,uyd,uxdr,uydr
      real(rp)  blend1, blend2,pi, ustar, blend3
      pi =4.0*atan(1.0)

      ustar =  mod_ugeo*1.1*0.41/log(1.0 + 0.5*hpbl/rough)
    
      ! 
      delta =0.4*hpbl ! mergin region width
      xi = y/hpbl
      fu = 1.57*xi -2.68*xi*xi
      fv = 13.2*xi -11.7*xi*xi
      !fv = fv*0.4
      
      uxd = ustar/0.41*( log(1.0+y/rough)+fu)
      uyd = - ustar/0.41*fv
      
      ! rotation angle of velocity close to ground, uxdown function
      theta = 20.0*pi/180.0
      !rotation
      uxdr = uxd*cos(theta) - uyd*sin(theta)
      uydr = uxd*sin(theta) + uyd*cos(theta)
      
      blend1 = 0.5*(1.0-tanh((xi-0.5)*6.0))!  2.0*hpbl/delta)) ! between 1 and 0 
      blend2 = 0.5*(1.0+tanh((xi-0.5)*6.0))!  2.0*hpbl/delta)) ! between 0 and 1
      ! Blending for Uy function
      blend3 = 0.5*(1.0-tanh((xi-0.75)*15.0))

      uxd = uxdr*blend1 + mod_ugeo*blend2
      uyd = uydr*blend3 ! + vegeo(2)*blend2
      
     

      !rotation to match geos direction
      theta = (Dir_Geos) * pi/180.0
      ux = uxd*cos(theta) - uyd*sin(theta)
      uy = uxd*sin(theta) + uyd*cos(theta)
      
      velabl(1) = ux
      velabl(2) = uy
    
  end subroutine veloc_abl_new


end module WindFarmSolverIncomp_mod
