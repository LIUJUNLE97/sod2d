module mod_arrays_heli
   use mod_constants

   implicit none

   real(rp), allocatable :: CT(:), D(:),pos_x(:),pos_y(:),pos_z(:),vol_correct(:),heli_ad(:)

end module mod_arrays_heli

module HelicopterSolverIncomp_mod
   use mod_arrays
   use mod_arrays_heli
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

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: HelicopterSolverIncomp

   real(rp) , public  :: vo,delta, rho0, Re, aoa_alpha, aoa_beta
   integer(4), public :: N_ad
   contains
      procedure, public :: fillBCTypes           => HelicopterSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => HelicopterSolverIncomp_initializeParameters
      procedure, public :: initializeSourceTerms => HelicopterSolverIncomp_initializeSourceTerms
      procedure, public :: evalInitialConditions => HelicopterSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer         => HelicopterSolverIncomp_initialBuffer
      procedure, public :: afterDt => HelicopterSolverIncomp_afterDt
      procedure, public :: readJSONAD => HelicopterSolverIncomp_readJSONAD
   end type HelicopterSolverIncomp
contains

subroutine HelicopterSolverIncomp_readJSONAD(this)
   use json_module
   implicit none
   class(HelicopterSolverIncomp), intent(inout) :: this
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
      !$acc enter data create(CT(:),D(:),pos_x(:),pos_y(:),pos_z(:),vol_correct(:))

      
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
      end do
   end if

   !$acc update device(CT(:),D(:),pos_x(:),pos_y(:),pos_z(:))

   call json%destroy()

end subroutine HelicopterSolverIncomp_readJSONAD

   subroutine HelicopterSolverIncomp_afterDt(this,istep)
      class(HelicopterSolverIncomp), intent(inout) :: this
      integer(4), intent(in) :: istep
      integer :: iNodeL,idime,iAD
      real(rp) :: x_ad,y_ad,z_ad,radius
      real(8) :: vol_T(this%N_ad),vol_T2(this%N_ad)

      if(istep == 1) then !esto va a dar problemas en continues
         if(this%N_ad .gt. 0) then
            allocate(heli_ad(numNodesRankPar))
            !$acc enter data create(heli_ad(:))
            
            vol_T(:) = 0.0d0

            !$acc update host(coordPar(:,:),Ml(:))
            
            !!$acc parallel loop reduction(+:vol_T)
            do iNodeL = 1,numNodesRankPar
               x_ad =  coordPar(iNodeL,1)
               y_ad =  coordPar(iNodeL,2)
               z_ad =  coordPar(iNodeL,3)
               
               heli_ad(iNodeL) = 0

               !$acc loop seq
               do iAD=1,this%N_ad
                  if((z_ad .gt. (pos_z(iAD)-(0.025_rp*D(iAD)))) .and. (z_ad .le. (pos_z(iAD)+(0.025_rp*D(iAD))))) then
                     radius = sqrt( (y_ad-pos_y(iAD))**2 + (x_ad-pos_x(iAD))**2 )
                     if(radius .le. (D(iAD)*0.5_rp)) then 
                        vol_T(iAD) = vol_T(iAD) + real(Ml(iNodeL),8)
                        heli_ad(iNodeL) = iAD
                     end if
                  end if
               end do
            end do
            !!$acc end parallel loop
                     
            do iAD=1,this%N_ad            
               call MPI_Allreduce(vol_T(iAD),vol_T2(iAD),1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
               vol_correct(iAD) = real(vol_T2(iAD),rp)/(v_pi*((D(iAD)*0.5_rp)**2)*(D(iAD)*0.05_rp))
               if(mpi_rank.eq.0)  write(111,*) " volT",vol_T2(iAD)," corr ",vol_correct(iAD)
            end do

            !$acc update device(vol_correct(:),heli_ad(:))

            !$acc parallel loop  
            do iNodeL = 1,numNodesRankPar
               if(heli_ad(iNodeL) .ne. 0) then
                  iAD = heli_ad(iNodeL)
                  source_term(iNodeL,1) = 0.0_rp
                  source_term(iNodeL,2) = 0.0_rp
                  source_term(iNodeL,3) = 0.5_rp*(CT(iAD)/(0.1_rp*D(iAD)))*(this%vo**2)*vol_correct(iAD)
               end if
            end do
            !$acc end parallel loop
         end if
      end if


   end subroutine HelicopterSolverIncomp_afterDt

   subroutine  HelicopterSolverIncomp_initialBuffer(this)
      class(HelicopterSolverIncomp), intent(inout) :: this
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

   end subroutine  HelicopterSolverIncomp_initialBuffer

   subroutine HelicopterSolverIncomp_fill_BC_Types(this)
      class(HelicopterSolverIncomp), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine HelicopterSolverIncomp_fill_BC_Types

   subroutine HelicopterSolverIncomp_initializeSourceTerms(this)
      class(HelicopterSolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL


      allocate(source_term(numNodesRankPar,ndime))
      !$acc enter data create(source_term(:,:))


      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) =  0.00_rp
         source_term(iNodeL,2) =  0.00_rp
         source_term(iNodeL,3) =  0.00_rp
      end do
      !$acc end parallel loop

   end subroutine HelicopterSolverIncomp_initializeSourceTerms

   subroutine HelicopterSolverIncomp_initializeParameters(this)
      use json_module
      implicit none
      class(HelicopterSolverIncomp), intent(inout) :: this
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

      call json%get("doTimerAnalysis",this%doTimerAnalysis, found,.false.)
      !----------------------------------------------

      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("period_walave",period_walave, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,1.0_rp); call this%checkFound(found,found_aux)

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
      call this%readJSONAD()

   end subroutine HelicopterSolverIncomp_initializeParameters

   subroutine HelicopterSolverIncomp_evalInitialConditions(this)
      class(HelicopterSolverIncomp), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL

      call nvtxStartRange("BluffBody3D_incomp Init")
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

   end subroutine HelicopterSolverIncomp_evalInitialConditions

end module HelicopterSolverIncomp_mod
