module mod_saveFields

   use mod_constants
   use mod_numerical_params , only : nspecies
   use mod_mpi
   use mod_custom_types
   !use json_module

   !------------------------------------------------------------------------------------------------------------------------
   integer(4), parameter :: indNS_rho = 1, indNS_mu = 2, indNS_pr = 3, indNS_ener = 4, indNS_eta = 5, &
                            indNS_csound = 6, indNS_machno = 7, indNS_divU = 8, indNS_qcrit = 9, indNS_temp = 10
	integer(4), parameter :: numNodeScalarFields = 10

   integer(4), parameter :: indES_mut = 1, indES_mue = 2 !, indES_qua = 3
   integer(4), parameter :: numElGPScalarFields = 2 !3

   integer(4), parameter :: indNV_vel = 1, indNV_gradRho = 2, indNV_curlU = 3
   integer(4), parameter :: numNodeVectorFields = 3

   character(128) :: nodeScalarNameFields(numNodeScalarFields),nodeVectorNameFields(numNodeVectorFields),elGPScalarNameFields(numElGPScalarFields)
   !------------------------------------------------------------------------------------------------------------------------
   integer(4), parameter :: indANS_rho = 1, indANS_pr = 2, indANS_pr2 = 3, indANS_mueff = 4
	integer(4), parameter :: numAvgNodeScalarFields = 4

   integer(4), parameter :: indANV_vel = 1, indANV_ve2 = 2, indANV_vex = 3, indANV_tw = 4
   integer(4), parameter :: numAvgNodeVectorFields = 4

   integer(4), parameter :: numAvgElGPScalarFields = 0

   character(128) :: avgNodeScalarNameFields(numAvgNodeScalarFields),avgNodeVectorNameFields(numAvgNodeVectorFields)
   !------------------------------------------------------------------------------------------------------------------------
   integer(4), parameter :: indRF_rho = 1, indRF_ux = 2, indRF_uy = 3, indRF_uz = 4, indRF_pr = 5, &
                            indRF_ener = 6, indRF_mut = 7, indRF_mue = 8, &
                            indRF_walavex = 9, indRF_walavey = 10, indRF_walavez = 11
	integer(4), parameter :: numRestartFields = 11

   character(128) :: restartNameFields(numRestartFields)
   !------------------------------------------------------------------------------------------------------------------------

   logical :: saveFieldsInitialized = .false.

   integer(4) :: numNodeScalarFields2save,numNodeVectorFields2save,numElGPScalarFields2save
   integer(4) :: numAvgNodeScalarFields2save,numAvgNodeVectorFields2save,numAvgElGPScalarFields2save
   type(ptr_array1d_rp_save) :: nodeScalarFields2save(numNodeScalarFields),avgNodeScalarFields2save(numAvgNodeScalarFields)
   type(ptr_array2d_rp_save) :: nodeVectorFields2save(numNodeVectorFields),avgNodeVectorFields2save(numAvgNodeVectorFields)
   type(ptr_array2d_rp_save) :: elGPScalarFields2save(numElGPScalarFields),avgElGPScalarFields2save(numAvgELGPScalarFields)
   logical :: save_nodeScalarField_rho,     save_nodeScalarField_muFluid,  save_nodeScalarField_pr,       save_nodeScalarField_energy, &
              save_nodeScalarField_entropy, save_nodeScalarField_csound,   save_nodeScalarField_machno,   save_nodeScalarField_divU,   &
              save_nodeScalarField_qcrit,   save_nodeScalarField_temp
   logical :: save_elGPScalarField_muSgs, save_elGPScalarField_muEnvit!, save_elGPScalarField_quality
   logical :: save_nodeVectorField_vel, save_nodeVectorField_gradRho, save_nodeVectorField_curlU
   
   logical :: save_avgNodeScalarField_rho, save_avgNodeScalarField_pr, save_avgNodeScalarField_pr2, save_avgNodeScalarField_mueff
   logical :: save_avgNodeVectorField_vel, save_avgNodeVectorField_ve2, save_avgNodeVectorField_vex, save_avgNodeVectorField_tw

contains

   subroutine init_saveFields()
      implicit none

      call initializeNameFields()
      call initializeSaveFieldsCnt()
      call initializeDefaultSaveFields()

      saveFieldsInitialized = .true.

   end subroutine init_saveFields

   subroutine initializeSaveFieldsCnt()
      implicit none

      ! Instantaneous fields
      numNodeScalarFields2save  = 0
      numNodeVectorFields2save  = 0
      numElGPScalarFields2save  = 0

      ! Average fields
      numAvgNodeScalarFields2save = 0
      numAvgNodeVectorFields2save = 0
      numAvgElGPScalarFields2save = 0

   end subroutine initializeSaveFieldsCnt

   subroutine initializeDefaultSaveFields()
      implicit none

      call initializeSaveFieldsCnt()

      save_nodeScalarField_rho      = .false.
      save_nodeScalarField_muFluid  = .false.
      save_nodeScalarField_pr       = .false.
      save_nodeScalarField_energy   = .false.
      save_nodeScalarField_entropy  = .false.
      save_nodeScalarField_csound   = .false.
      save_nodeScalarField_machno   = .false.
      save_nodeScalarField_divU     = .false.
      save_nodeScalarField_qcrit    = .false.
      save_nodeScalarField_temp     = .false.
      save_elGPScalarField_muSgs    = .false.
      save_elGPScalarField_muEnvit  = .false.
!       save_elGPScalarField_quality  = .false.
      save_nodeVectorField_vel      = .false.
      save_nodeVectorField_gradRho  = .false.
      save_nodeVectorField_curlU    = .false.

      save_avgNodeScalarField_rho   = .false.
      save_avgNodeScalarField_pr    = .false.
      save_avgNodeScalarField_pr2   = .false.
      save_avgNodeScalarField_mueff = .false.
      save_avgNodeVectorField_vel   = .false.
      save_avgNodeVectorField_ve2   = .false.
      save_avgNodeVectorField_vex   = .false.
      save_avgNodeVectorField_tw    = .false.

      
   end subroutine initializeDefaultSaveFields

   subroutine initializeDefaultSaveFieldsCompressible()
      implicit none

      call initializeSaveFieldsCnt()

      save_nodeScalarField_rho      = .true.
      save_nodeScalarField_muFluid  = .true.
      save_nodeScalarField_pr       = .true.
      save_nodeScalarField_energy   = .true.
      save_nodeScalarField_entropy  = .true.
      save_nodeScalarField_csound   = .true.
      save_nodeScalarField_machno   = .true.
      save_nodeScalarField_divU     = .true.
      save_nodeScalarField_qcrit    = .true.
      save_nodeScalarField_temp     = .false.
      save_elGPScalarField_muSgs    = .true.
      save_elGPScalarField_muEnvit  = .true.
!       save_elGPScalarField_quality  = .false.
      save_nodeVectorField_vel      = .true.
      save_nodeVectorField_gradRho  = .true.
      save_nodeVectorField_curlU    = .true.

      save_avgNodeScalarField_rho   = .true.
      save_avgNodeScalarField_pr    = .true.
      save_avgNodeScalarField_pr2   = .true.
      save_avgNodeScalarField_mueff = .true.
      save_avgNodeVectorField_vel   = .true.
      save_avgNodeVectorField_ve2   = .true.
      save_avgNodeVectorField_vex   = .true.
      save_avgNodeVectorField_tw    = .true.

   end subroutine initializeDefaultSaveFieldsCompressible
   
   subroutine initializeDefaultSaveFieldsIncompressible()
      implicit none

      call initializeSaveFieldsCnt()

      save_nodeScalarField_rho      = .false.
      save_nodeScalarField_muFluid  = .true.
      save_nodeScalarField_pr       = .true.
      save_nodeScalarField_energy   = .false.
      save_nodeScalarField_entropy  = .true.
      save_nodeScalarField_csound   = .false.
      save_nodeScalarField_machno   = .false.
      save_nodeScalarField_divU     = .false.
      save_nodeScalarField_qcrit    = .true.
      save_nodeScalarField_temp     = .false.
      save_elGPScalarField_muSgs    = .true.
      save_elGPScalarField_muEnvit  = .true.
!       save_elGPScalarField_quality  = .false.
      save_nodeVectorField_vel      = .true.
      save_nodeVectorField_gradRho  = .false.
      save_nodeVectorField_curlU    = .true.

      save_avgNodeScalarField_rho   = .false.
      save_avgNodeScalarField_pr    = .true.
      save_avgNodeScalarField_pr2   = .true.
      save_avgNodeScalarField_mueff = .true.
      save_avgNodeVectorField_vel   = .true.
      save_avgNodeVectorField_ve2   = .true.
      save_avgNodeVectorField_vex   = .true.
      save_avgNodeVectorField_tw    = .true.

   end subroutine initializeDefaultSaveFieldsIncompressible

   subroutine initializeDefaultSaveFieldsElasticity()
      implicit none

      call initializeSaveFieldsCnt()

      save_nodeScalarField_rho      = .false.
      save_nodeScalarField_muFluid  = .false.
      save_nodeScalarField_pr       = .false.
      save_nodeScalarField_energy   = .false.
      save_nodeScalarField_entropy  = .false.
      save_nodeScalarField_csound   = .false.
      save_nodeScalarField_machno   = .false.
      save_nodeScalarField_divU     = .false.
      save_nodeScalarField_qcrit    = .false.
      save_nodeScalarField_temp     = .false.
      save_elGPScalarField_muSgs    = .false.
      save_elGPScalarField_muEnvit  = .false.
!       save_elGPScalarField_quality  = .true.
      save_nodeVectorField_vel      = .true.
      save_nodeVectorField_gradRho  = .false.
      save_nodeVectorField_curlU    = .false.

      save_avgNodeScalarField_rho   = .false.
      save_avgNodeScalarField_pr    = .false.
      save_avgNodeScalarField_pr2   = .false.
      save_avgNodeScalarField_mueff = .false.
      save_avgNodeVectorField_vel   = .false.
      save_avgNodeVectorField_ve2   = .false.
      save_avgNodeVectorField_vex   = .false.
      save_avgNodeVectorField_tw    = .false.


   end subroutine initializeDefaultSaveFieldsElasticity
   
   subroutine initializeNameFields()
      implicit none

      if(mpi_rank.eq.0) write(*,*) 'Initializing default available var names'
      !----------------------------------------------------------------------

      !---------   nodeScalars  ---------------------------------------------
      !----------------------------------------------------------------------
      nodeScalarNameFields(indNS_rho)     = 'rho'
      nodeScalarNameFields(indNS_mu)      = 'mu_fluid'
      nodeScalarNameFields(indNS_pr)      = 'pr'
      nodeScalarNameFields(indNS_ener)    = 'E'
      nodeScalarNameFields(indNS_eta)     = 'eta'
      nodeScalarNameFields(indNS_csound)  = 'csound'
      nodeScalarNameFields(indNS_machno)  = 'machno'
      nodeScalarNameFields(indNS_divU)    = 'divU'
      nodeScalarNameFields(indNS_qcrit)   = 'qcrit'
      nodeScalarNameFields(indNS_temp)    = 'temp'
      !----------  vectorScalars   ------------------------------------------
      !----------------------------------------------------------------------
      nodeVectorNameFields(indNV_vel)     = 'u'
      nodeVectorNameFields(indNV_gradRho) = 'gradRho'
      nodeVectorNameFields(indNV_curlU)   = 'curlU'
      !---------    elGPScalars   -------------------------------------------
      !----------------------------------------------------------------------
      elGPScalarNameFields(indES_mut)     = 'mut'
      elGPScalarNameFields(indES_mue)     = 'mue'
!       elGPScalarNameFields(indES_qua)     = 'quality'
      !----------------------------------------------------------------------

      !----------------------  AVERAGE FIELDS  ------------------------------
      !---------   nodeScalars  ---------------------------------------------
      !----------------------------------------------------------------------
      avgNodeScalarNameFields(indANS_rho)   = 'avrho'
      avgNodeScalarNameFields(indANS_pr)    = 'avpre'
      avgNodeScalarNameFields(indANS_pr2)   = 'avpr2'
      avgNodeScalarNameFields(indANS_mueff) = 'avmueff'
      !---------  vectorScalars   -------------------------------------------
      !----------------------------------------------------------------------
      avgNodeVectorNameFields(indANV_vel) = 'avvel'
      avgNodeVectorNameFields(indANV_ve2) = 'avve2'
      avgNodeVectorNameFields(indANV_vex) = 'avvex'
      avgNodeVectorNameFields(indANV_tw)  = 'avtw'

      !----------------------  RESTART FIELDS  ------------------------------
      !----------------------------------------------------------------------

      restartNameFields(indRF_rho)     = 'rho'
      restartNameFields(indRF_ux)      = 'u_x'
      restartNameFields(indRF_uy)      = 'u_y'
      restartNameFields(indRF_uz)      = 'u_z'
      restartNameFields(indRF_pr)      = 'pr'
      restartNameFields(indRF_ener)    = 'E'
      restartNameFields(indRF_mut)     = 'mue'
      restartNameFields(indRF_mue)     = 'mut'
      restartNameFields(indRF_walavex) = 'walave_u_x'
      restartNameFields(indRF_walavey) = 'walave_u_y'
      restartNameFields(indRF_walavez) = 'walave_u_z'

   end subroutine initializeNameFields

   subroutine read_json_saveFields(json_filename)
      use json_module
      implicit none
      character(len=*),intent(in) :: json_filename
      type(json_file) :: json_f
      logical :: isFound

      call json_f%initialize()
      call json_f%load_file(json_filename)

      !if not found in the json file, the default value is preserved!
      call json_f%get("save_nodeScalarField_rho",    save_nodeScalarField_rho    , isFound, save_nodeScalarField_rho)    
      call json_f%get("save_nodeScalarField_muFluid",save_nodeScalarField_muFluid, isFound, save_nodeScalarField_muFluid)
      call json_f%get("save_nodeScalarField_pr",     save_nodeScalarField_pr     , isFound, save_nodeScalarField_pr)
      call json_f%get("save_nodeScalarField_energy", save_nodeScalarField_energy , isFound, save_nodeScalarField_energy)
      call json_f%get("save_nodeScalarField_entropy",save_nodeScalarField_entropy, isFound, save_nodeScalarField_entropy)
      call json_f%get("save_nodeScalarField_csound", save_nodeScalarField_csound , isFound, save_nodeScalarField_csound)
      call json_f%get("save_nodeScalarField_machno", save_nodeScalarField_machno , isFound, save_nodeScalarField_machno)
      call json_f%get("save_nodeScalarField_divU",   save_nodeScalarField_divU   , isFound, save_nodeScalarField_divU)
      call json_f%get("save_nodeScalarField_qcrit",  save_nodeScalarField_qcrit  , isFound, save_nodeScalarField_qcrit)
      call json_f%get("save_nodeScalarField_temp",   save_nodeScalarField_temp   , isFound, save_nodeScalarField_temp)
      call json_f%get("save_elGPScalarField_muSgs",  save_elGPScalarField_muSgs  , isFound, save_elGPScalarField_muSgs)
      call json_f%get("save_elGPScalarField_muEnvit",save_elGPScalarField_muEnvit, isFound, save_elGPScalarField_muEnvit)
!       call json_f%get("save_elGPScalarField_quality",save_elGPScalarField_quality, isFound, save_elGPScalarField_quality)
      call json_f%get("save_nodeVectorField_vel",    save_nodeVectorField_vel    , isFound, save_nodeVectorField_vel)
      call json_f%get("save_nodeVectorField_gradRho",save_nodeVectorField_gradRho, isFound, save_nodeVectorField_gradRho)
      call json_f%get("save_nodeVectorField_curlU",  save_nodeVectorField_curlU  , isFound, save_nodeVectorField_curlU)

      call json_f%get("save_avgNodeScalarField_rho",  save_avgNodeScalarField_rho  , isFound, save_avgNodeScalarField_rho)
      call json_f%get("save_avgNodeScalarField_pr",   save_avgNodeScalarField_pr   , isFound, save_avgNodeScalarField_pr)
      call json_f%get("save_avgNodeScalarField_pr2",  save_avgNodeScalarField_pr2  , isFound, save_avgNodeScalarField_pr2)
      call json_f%get("save_avgNodeScalarField_mueff",save_avgNodeScalarField_mueff, isFound, save_avgNodeScalarField_mueff)
      call json_f%get("save_avgNodeVectorField_vel",  save_avgNodeVectorField_vel  , isFound, save_avgNodeVectorField_vel)
      call json_f%get("save_avgNodeVectorField_ve2",  save_avgNodeVectorField_ve2  , isFound, save_avgNodeVectorField_ve2)
      call json_f%get("save_avgNodeVectorField_vex",  save_avgNodeVectorField_vex  , isFound, save_avgNodeVectorField_vex)
      call json_f%get("save_avgNodeVectorField_tw",   save_avgNodeVectorField_tw   , isFound, save_avgNodeVectorField_tw)

      call json_f%destroy()
   end subroutine

   subroutine setFields2Save(rho,mu_fluid,pr,E,eta,csound,machno,divU,qcrit,Tem,u,gradRho,curlU,mu_sgs,mu_e,&
                             avrho,avpre,avpre2,avmueff,avvel,avve2,avvex,avtw,Yk,mu_e_Yk)!,quality_e)
      implicit none
      real(rp),intent(in),dimension(:) :: rho,mu_fluid,pr,E,eta,csound,machno,divU,qcrit,Tem
      real(rp),intent(in),dimension(:,:) :: u,gradRho,curlU
      real(rp),intent(in),dimension(:,:) :: mu_sgs,mu_e
      real(rp_avg),intent(in),dimension(:) :: avrho,avpre,avpre2,avmueff
      real(rp_avg),intent(in),dimension(:,:) :: avvel,avve2,avvex,avtw
      real(rp),intent(in),dimension(:,:),optional :: Yk
      real(rp),intent(in),dimension(:,:,:),optional :: mu_e_Yk
      !real(rp),intent(in),dimension(:),optional :: quality_e ! constant on element, not gauss point!
      integer(4) :: ispc
      character(128) :: dsetname


		if(saveFieldsInitialized .eqv. .false.) then
         write(*,*) "FATAL ERROR in setFields2Save! Save Fields NOT properly initialized! Aborting!"
       	call MPI_Abort(app_comm,-1,mpi_err)
      end if

      if(mpi_rank.eq.0) write(*,*) 'Setting default fields to be saved'
      !-----------------------------------------------------------------------

      !---------   nodeScalars  -----------------------------
      !------------------------------------------------------
      if(save_nodeScalarField_rho) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_rho),rho(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_muFluid) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_mu),mu_fluid(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_pr) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_pr),pr(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_energy) then
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_ener),E(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_entropy) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_eta),eta(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_csound) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_csound),csound(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_machno) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_machno),machno(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_divU) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_divU),divU(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_qcrit) then 
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_qcrit),qcrit(:))
      end if
      !------------------------------------------------------
      if(save_nodeScalarField_temp) then
         call add_nodeScalarField2save(nodeScalarNameFields(indNS_temp),Tem(:))
      end if

      !---------------  vectorScalars   -------------------------------------
      !----------------------------------------------------------------------
      if(save_nodeVectorField_vel) then 
         call add_nodeVectorField2save(nodeVectorNameFields(indNV_vel),u(:,:))
      end if
      !----------------------------------------------------------------------
      if(save_nodeVectorField_gradRho) then 
         call add_nodeVectorField2save(nodeVectorNameFields(indNV_gradRho),gradRho(:,:))
      end if
      !----------------------------------------------------------------------
      if(save_nodeVectorField_curlU) then 
         call add_nodeVectorField2save(nodeVectorNameFields(indNV_curlU),curlU(:,:))
      end if
      !----------------------------------------------------------------------

      !-------------    elGPScalars   ------------------------------------- 
      if(save_elGPScalarField_muSgs) then 
         call add_elGPScalarField2save(elGPScalarNameFields(indES_mut),mu_sgs(:,:))
      end if
      !----------------------------------------------------------------------
      if(save_elGPScalarField_muEnvit) then 
         call add_elGPScalarField2save(elGPScalarNameFields(indES_mue),mu_e(:,:))
      end if
      !----------------------------------------------------------------------
!       if(save_elGPScalarField_quality) then
!          call add_elGPScalarField2save(elGPScalarNameFields(indES_qua),quality_e(:,:))
!       end if
      !----------------------------------------------------------------------
      
      !---------------  Species   -------------------------------------
      if(present(Yk)) then
         do ispc = 1, nspecies
            dsetname=trim(adjustl('Yk'))   !//char(ispc)
            call add_nodeScalarField2save(dsetname,Yk(:,ispc))
            dsetname=trim(adjustl('mue_Yk'))   !//char(ispc)
            call add_elGPScalarField2save(dsetname,mu_e_Yk(:,:,ispc))
         end do
      end if

      !----------------------  AVERAGE FIELDS  ------------------------------

      !---------   nodeScalars  -----------------------------------
      !------------------------------------------------------------
      if(save_avgNodeScalarField_rho) then
         call add_avgNodeScalarField2save(avgNodeScalarNameFields(indANS_rho),avrho(:))
      end if
       !------------------------------------------------------------
      if(save_avgNodeScalarField_pr) then
         call add_avgNodeScalarField2save(avgNodeScalarNameFields(indANS_pr),avpre(:))
      end if
      !------------------------------------------------------------
      if(save_avgNodeScalarField_pr2) then
         call add_avgNodeScalarField2save(avgNodeScalarNameFields(indANS_pr2),avpre2(:))
      end if
      !------------------------------------------------------------
      if(save_avgNodeScalarField_mueff) then
         call add_avgNodeScalarField2save(avgNodeScalarNameFields(indANS_mueff),avmueff(:))
      end if

      !---------------  vectorScalars   -------------------------------------
      !----------------------------------------------------------------------
      if(save_avgNodeVectorField_vel) then
         call add_avgNodeVectorField2save(avgNodeVectorNameFields(indANV_vel),avvel(:,:))
      end if
      if(save_avgNodeVectorField_ve2) then
         call add_avgNodeVectorField2save(avgNodeVectorNameFields(indANV_ve2),avve2(:,:))
      end if
      if(save_avgNodeVectorField_vex) then
         call add_avgNodeVectorField2save(avgNodeVectorNameFields(indANV_vex),avvex(:,:))
      end if
      if(save_avgNodeVectorField_tw) then
         call add_avgNodeVectorField2save(avgNodeVectorNameFields(indANV_tw),avtw(:,:))
      end if
      !----------------------------------------------------------------------
      
      
   
   end subroutine setFields2Save


!--------------------------------------------------------------------------------------------------------------------------

   subroutine add_nodeScalarField2save(fieldSaveName,array2save)
      implicit none
      character(*),intent(in) :: fieldSaveName
      real(rp),target,intent(in) :: array2save(:)

      numNodeScalarFields2save = numNodeScalarFields2save + 1

      if(numNodeScalarFields2save .gt. max_num_saved_fields) then
         if(mpi_rank.eq.0) then
            write(111,*) 'WARNING! Trying to add node scalarfield ',fieldSaveName,' num ',numNodeScalarFields2save,' but max_num_saved_fields ',max_num_saved_fields
            write(111,*) 'node scalarfield NOT ADDED! If required, modify max_num_saved_fields in mod_constants.f90'
         end if
         return
      end if

      nodeScalarFields2save(numNodeScalarFields2save)%nameField = trim(adjustl(fieldSaveName))
      nodeScalarFields2save(numNodeScalarFields2save)%ptr_rp => array2save
      nodeScalarFields2save(numNodeScalarFields2save)%ptr_avg => null()

   end subroutine add_nodeScalarField2save

!--------------------------------------------------------------------------------------------------------------------------

   subroutine add_nodeVectorField2save(fieldSaveName,array2save)
      implicit none
      character(*),intent(in) :: fieldSaveName
      real(rp),target,intent(in) :: array2save(:,:)

      numNodeVectorFields2save = numNodeVectorFields2save + 1

      if(numNodeVectorFields2save .gt. max_num_saved_fields) then
         if(mpi_rank.eq.0) then
            write(111,*) 'WARNING! Trying to add node vectorfield ',fieldSaveName,' num ',numNodeVectorFields2save,' but max_num_saved_fields ',max_num_saved_fields
            write(111,*) 'node vectorfield NOT ADDED! If required, modify max_num_saved_fields in mod_constants.f90'
         end if
         return
      end if

      nodeVectorFields2save(numNodeVectorFields2save)%nameField = trim(adjustl(fieldSaveName))
      nodeVectorFields2save(numNodeVectorFields2save)%ptr_rp => array2save
      nodeVectorFields2save(numNodeVectorFields2save)%ptr_avg => null()

   end subroutine add_nodeVectorField2save

!--------------------------------------------------------------------------------------------------------------------------

   subroutine add_elGPScalarField2save(fieldSaveName,array2save)
      implicit none
      character(*),intent(in) :: fieldSaveName
      real(rp),target,intent(in) :: array2save(:,:)

      numElGPScalarFields2save = numElGPScalarFields2save + 1

      if(numElGPScalarFields2save .gt. max_num_saved_fields) then
         if(mpi_rank.eq.0) then
            write(111,*) 'WARNING! Trying to add elGP scalarfield ',fieldSaveName,' num ',numElGPScalarFields2save,' but max_num_saved_fields ',max_num_saved_fields
            write(111,*) 'elGP scalarfield NOT ADDED! If required, modify max_num_saved_fields in mod_constants.f90'
         end if
         return
      end if

      elGPScalarFields2save(numElGPScalarFields2save)%nameField = trim(adjustl(fieldSaveName))
      elGPScalarFields2save(numElGPScalarFields2save)%ptr_rp => array2save
      elGPScalarFields2save(numElGPScalarFields2save)%ptr_avg => null()

   end subroutine add_elGPScalarField2save

   subroutine add_avgNodeScalarField2save(fieldSaveName,array2save)
      implicit none
      character(*),intent(in) :: fieldSaveName
      real(rp_avg),target,intent(in) :: array2save(:)

      numAvgNodeScalarFields2save = numAvgNodeScalarFields2save + 1

      if(numAvgNodeScalarFields2save .gt. max_num_saved_fields) then
         if(mpi_rank.eq.0) then
            write(111,*) 'WARNING! Trying to add node scalarfield ',fieldSaveName,' num ',numAvgNodeScalarFields2save,' but max_num_saved_fields ',max_num_saved_fields
            write(111,*) 'node scalarfield NOT ADDED! If required, modify max_num_saved_fields in mod_constants.f90'
         end if
         return
      end if

      avgNodeScalarFields2save(numAvgNodeScalarFields2save)%nameField = trim(adjustl(fieldSaveName))
      avgNodeScalarFields2save(numAvgNodeScalarFields2save)%ptr_avg => array2save
      avgNodeScalarFields2save(numAvgNodeScalarFields2save)%ptr_rp => null()

   end subroutine add_avgNodeScalarField2save

!--------------------------------------------------------------------------------------------------------------------------

   subroutine add_avgNodeVectorField2save(fieldSaveName,array2save)
      implicit none
      character(*),intent(in) :: fieldSaveName
      real(rp_avg),target,intent(in) :: array2save(:,:)

      numAvgNodeVectorFields2save = numAvgNodeVectorFields2save + 1

      if(numAvgNodeVectorFields2save .gt. max_num_saved_fields) then
         if(mpi_rank.eq.0) then
            write(111,*) 'WARNING! Trying to add node vectorfield ',fieldSaveName,' num ',numAvgNodeVectorFields2save,' but max_num_saved_fields ',max_num_saved_fields
            write(111,*) 'node vectorfield NOT ADDED! If required, modify max_num_saved_fields in mod_constants.f90'
         end if
         return
      end if

      avgNodeVectorFields2save(numAvgNodeVectorFields2save)%nameField = trim(adjustl(fieldSaveName))
      avgNodeVectorFields2save(numAvgNodeVectorFields2save)%ptr_avg => array2save
      avgNodeVectorFields2save(numAvgNodeVectorFields2save)%ptr_rp => null()

   end subroutine add_avgNodeVectorField2save

!--------------------------------------------------------------------------------------------------------------------------

   subroutine add_avgElGPScalarField2save(fieldSaveName,array2save)
      implicit none
      character(*),intent(in) :: fieldSaveName
      real(rp_avg),target,intent(in) :: array2save(:,:)

      numAvgElGPScalarFields2save = numAvgElGPScalarFields2save + 1

      if(numAvgElGPScalarFields2save .gt. max_num_saved_fields) then
         if(mpi_rank.eq.0) then
            write(111,*) 'WARNING! Trying to add elGP scalarfield ',fieldSaveName,' num ',numAvgElGPScalarFields2save,' but max_num_saved_fields ',max_num_saved_fields
            write(111,*) 'elGP scalarfield NOT ADDED! If required, modify max_num_saved_fields in mod_constants.f90'
         end if
         return
      end if

      avgElGPScalarFields2save(numAvgElGPScalarFields2save)%nameField = trim(adjustl(fieldSaveName))
      avgElGPScalarFields2save(numAvgElGPScalarFields2save)%ptr_avg => array2save
      avgElGPScalarFields2save(numAvgElGPScalarFields2save)%ptr_rp => null()

   end subroutine add_avgElGPScalarField2save

end module mod_saveFields