module mod_arrays
      use mod_constants

      implicit none

      ! main allocatable arrays
      ! integer ---------------------------------------------------
      integer(4), allocatable :: lelpn(:),point2elem(:),bouCodes2BCType(:)
      integer(4), allocatable :: atoIJ(:),atoIJK(:),invAtoIJK(:,:,:),gmshAtoI(:),gmshAtoJ(:),gmshAtoK(:),lnbnNodes(:)
      integer(4), allocatable :: witel(:), buffstep(:),maskMapped(:),ad(:)
      integer(4), allocatable :: witGlob(:)

      ! real ------------------------------------------------------
      real(rp), allocatable :: normalsAtNodes(:,:)
      real(rp), allocatable :: helem(:),helem_l(:,:)
      real(rp), allocatable :: xgp(:,:), wgp(:), xgp_b(:,:), wgp_b(:)
      real(rp), allocatable :: Ngp(:,:), dNgp(:,:,:), Ngp_b(:,:), dNgp_b(:,:,:)
      real(rp), allocatable :: Ngp_l(:,:), dNgp_l(:,:,:),dlxigp_ip(:,:,:)
      real(rp), allocatable :: Ngp_equi(:,:)
      real(rp), allocatable :: gpvol(:,:,:),Je(:,:), He(:,:,:,:), bou_norm(:,:),Ml(:),mu_factor(:),source_term(:,:)

      real(rp), target,allocatable :: gradRho(:,:), curlU(:,:), divU(:), Qcrit(:)
      real(rp), target,allocatable :: u(:,:,:),q(:,:,:),rho(:,:),pr(:,:),E(:,:),Tem(:,:),e_int(:,:),csound(:),eta(:,:),machno(:),tauw(:,:)
      real(rp), target,allocatable :: mu_e(:,:),mu_fluid(:),mu_sgs(:,:)

      real(rp_avg), allocatable :: avrho(:), avpre(:), avpre2(:), avvel(:,:), avve2(:,:), avmueff(:),avvex(:,:),avtw(:,:)
      real(rp), allocatable :: kres(:),etot(:),au(:,:),ax1(:),ax2(:),ax3(:)
      real(rp), allocatable :: Fpr(:,:), Ftau(:,:)
      real(rp), allocatable :: witxi(:,:), Nwit(:,:), buffwit(:,:,:), bufftime(:)

      real(rp), allocatable :: u_buffer(:,:), u_mapped(:,:)

      ! exponential average for wall law
      real(rp), allocatable :: walave_u(:,:),walave_t(:)

      ! roughness for wall law
      real(rp), allocatable :: zo(:)

      ! for entropy and sgs visc.
      real(rp),  allocatable :: mue_l(:,:),al_weights(:),am_weights(:),an_weights(:)
      
      ! for scalars 
      real(rp), allocatable :: mu_e_Yk(:,:,:), eta_Yk(:,:,:), Yk(:,:,:), Yk_buffer(:,:)

end module mod_arrays

module CFDSolverBase_mod
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
      use time_integ_imex
      use time_integ_ls
      use time_integ_species_imex
      use time_integ_species
      use mod_analysis
      use mod_numerical_params
      use mod_time_ops
      use mod_fluid_viscosity
      use mod_postpro
      use mod_aver
      use mod_mpi
      use mod_mpi_mesh
      use mod_hdf5
      use mod_comms
      use mod_comms_boundaries
      use mod_custom_types
      use mod_saveFields
      use mod_witness_points
      use mod_filters
      use mod_InSitu
   implicit none
   private

   type, public :: CFDSolverBase

      ! main integer parameters
      integer(4), public :: save_logFile_first,save_logFile_step,save_logFile_next
      integer(4), public :: save_restartFile_first,save_restartFile_step,save_restartFile_next
      integer(4), public :: save_resultsFile_first,save_resultsFile_step,save_resultsFile_next
      integer(4), public :: restartFileCnt,restartFile_to_load
      integer(4), public :: initial_istep,final_istep,load_step,local_step

      integer(4), public :: currentNonLinealIter
      integer(4), public :: nwit,nwitPar,leapwit
      integer(4), public :: leapwitsave=100 !Default value, found to be robustly fast
      integer(4), public :: load_stepwit = 0
      integer(4), public :: nvarwit=5 !Default value, only to be substituted if function update_witness is modified to

      ! main logical parameters
      logical, public :: loadRestartFile=.false.,saveAvgFile=.false.,loadAvgFile=.false.,saveInitialField=.false.,saveSurfaceResults=.false.,continue_oldLogs=.false.
      logical, public :: doGlobalAnalysis=.false.,isFreshStart=.true.,doTimerAnalysis=.false.,isWallModelOn=.false.,isSymmetryOn=.false.
      logical, public :: useIntInComms=.false.,useRealInComms=.false.
      logical, public :: have_witness=.false.,wit_save_u_i=.false.,wit_save_pr=.false.,wit_save_rho=.false., continue_witness=.false.,wit_save=.true.

      ! main char variables
      character(512) :: log_file_name
      character(512) :: mesh_h5_file_path,mesh_h5_file_name
      character(512) :: results_h5_file_path,results_h5_file_name
      character(512) :: io_prepend_path,io_append_info
      character(512) :: witness_inp_file_name,witness_h5_file_name
      character(512) :: meshFile_h5_full_name,surface_meshFile_h5_full_name
      character(512) :: base_resultsFile_h5_full_name,base_avgResultsFile_h5_full_name,base_restartFile_h5_full_name

      ! main real parameters
      real(rp) , public                   :: cfl_conv,cfl_diff
      real(rp) , public                   :: leviCivi(3,3,3), surfArea, EK, VolTot, eps_D, eps_S, eps_T, maxmachno
      real(rp) , public                   :: dt, Cp, Rgas, gamma_gas,Prt
      real(rp) , public                   :: time, maxPhysTime, initial_avgTime, elapsed_avgTime
      real(rp) , public                   :: loadtimewit=0.0_rp
      logical  , public                   :: noBoundaries = .false.

   contains
      procedure, public :: printDt => CFDSolverBase_printDt
      procedure, public :: printAll => CFDSolverBase_printAll
      procedure, public :: run => CFDSolverBase_run
      !-------------------------------------------------------------------------------------------------------
      procedure, public :: initializeDefaultFileNames => CFDSolverBase_initializeDefaultFileNames
      procedure, public :: initializeDefaultVariables => CFDSolverBase_initializeDefaultVariables
      procedure, public :: initializeDefaultFlags => CFDSolverBase_initializeDefaultFlags
      procedure, public :: initializeDefaultSaveFields => CFDSolverBase_initializeDefaultSaveFields
      procedure, public :: initializeDefaultCustomSettings => CFDSolverBase_initializeDefaultCustomSettings
      procedure, public :: initializeDefaultParameters => CFDSolverBase_initializeDefaultParameters
      procedure, public :: initializeParameters => CFDSolverBase_initializeParameters
      !-------------------------------------------------------------------------------------------------------
      procedure, public :: optimizeParameters => CFDSolverBase_optimizeParameters
      procedure, public :: initializeSourceTerms => CFDSolverBase_initializeSourceTerms
      procedure, public :: openMesh => CFDSolverBase_openMesh
      procedure, public :: evalCharLength => CFDSolverBase_evalCharLength
      procedure, public :: boundaryFacesToNodes => CFDSolverBase_boundaryFacesToNodes
      procedure, public :: normalFacesToNodes => CFDSolverBase_normalFacesToNodes
      procedure, public :: fillBCTypes => CFDSolverBase_fill_BC_Types
      procedure, public :: allocateVariables => CFDSolverBase_allocateVariables
      procedure, public :: deallocateVariables => CFDSolverBase_deallocateVariables
      procedure, public :: evalOrLoadInitialConditions => CFDSolverBase_evalOrLoadInitialConditions
      procedure, public :: evalInitialConditions => CFDSolverBase_evalInitialConditions
      procedure, public :: evalInitialViscosity =>CFDSolverBase_evalInitialViscosity
      procedure, public :: evalViscosityFactor=>CFDSolverBase_evalViscosityFactor
      procedure, public :: evalInitialDt =>CFDSolverBase_evalInitialDt
      procedure, public :: evalDt =>CFDSolverBase_evalDt
      procedure, public :: evalShapeFunctions =>CFDSolverBase_evalShapeFunctions
      procedure, public :: evalBoundaryNormals =>CFDSolverBase_evalBoundaryNormals
      procedure, public :: evalJacobians =>CFDSolverBase_evalJacobians
      procedure, public :: evalAtoIJKInverse =>CFDSolverBase_evalAtoIJKInverse
      procedure, public :: eval_elemPerNode_and_nearBoundaryNode =>CFDSolverBase_eval_elemPerNode_and_nearBoundaryNode
      procedure, public :: set_mappedFaces_linkingNodes =>CFDSolverBase_set_mappedFaces_linkingNodes
      procedure, public :: evalMass=>CFDSolverBase_evalMass
      procedure, public :: evalFirstOutput =>CFDSolverBase_evalFirstOutput
      procedure, public :: evalTimeIteration =>CFDSolverBase_evalTimeIteration
      procedure, public :: callTimeIntegration =>CFDSolverBase_callTimeIntegration
      procedure, public :: saveRestartFile =>CFDSolverBase_saveRestartFile
      procedure, public :: saveAvgResultsFiles =>CFDSolverBase_saveAvgResultsFiles
      procedure, public :: saveInstResultsFiles =>CFDSolverBase_saveInstResultsFiles
      procedure, public :: afterDt =>CFDSolverBase_afterDt
      procedure, public :: beforeTimeIteration =>CFDSolverBase_beforeTimeIteration
      procedure, public :: afterTimeIteration =>CFDSolverBase_afterTimeIteration
      procedure, public :: update_witness =>CFDSolverBase_update_witness
      procedure, public :: preprocWitnessPoints =>CFDSolverBase_preprocWitnessPoints
      procedure, public :: loadWitnessPoints =>CFDSolverBase_loadWitnessPoints
      procedure, public :: save_witness =>CFDSolverBase_save_witness

      procedure, public :: initialBuffer =>CFDSolverBase_initialBuffer

      procedure, public :: initNSSolver => CFDSolverBase_initNSSolver
      procedure, public :: endNSSolver => CFDSolverBase_endNSSolver

      procedure, public :: checkFound => CFDSolverBase_checkFound
      procedure, public :: readJSONBCTypes => CFDSolverBase_readJSONBCTypes
      procedure, public :: readJSONBuffer => CFDSolverBase_readJSONBuffer
      procedure, public :: eval_vars_after_load_hdf5_resultsFile => CFDSolverBase_eval_vars_after_load_hdf5_resultsFile 
      procedure, public :: findFixPressure => CFDSolverBase_findFixPressure

      procedure :: open_log_file
      procedure :: close_log_file
      procedure :: open_analysis_files
      procedure :: close_analysis_files
      procedure :: flush_log_file
      procedure :: eval_initial_mu_sgs
      procedure :: checkIfWallModelOn
      procedure :: checkIfSymmetryOn
   end type CFDSolverBase
contains

subroutine CFDSolverBase_findFixPressure(this)
   implicit none
   class(CFDSolverBase), intent(inout) :: this
   integer(4) :: ielem, iRankl, iNodeL,iRank

   ielem = numElemsRankPar*0.5
   iRankl = mpi_size + 1
   !$acc parallel loop  
   do iNodeL = 1,numNodesRankPar
      if(maskMapped(iNodeL) .eq. 0) then
         iRankl = mpi_rank
      end if
   end do
   !$acc end parallel loop
   call MPI_Allreduce(iRankl,iRank,1,mpi_datatype_int,MPI_MIN,app_comm,mpi_err)

   if(iRank == (mpi_size+1)) iRank = 0
   if(mpi_rank.eq.iRank) then
       inode_fix_press =  connecParWork(ielem,atoIJK(nnode))
       if(maskMapped(inode_fix_press) .eq. 1) then
               do iNodeL = 1,numNodesRankPar
                     if(maskMapped(iNodeL) .eq. 0) then
                       inode_fix_press = iNodeL
                       exit
                     end if
               end do
       end if
       !inode_fix_press = 4084081  ! 11500000 mid cedval mesh
       write(111,*) '--| Node to fix pressure',inode_fix_press, " mpi_rank ",iRank,"  coord ",coordPar(inode_fix_press,1:3)," dir  ", perMapFaceDir
   else
      flag_fs_fix_pressure = .false.
   end if

end subroutine CFDSolverBase_findFixPressure

   subroutine CFDSolverBase_readJSONBuffer(this)
      use json_module
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      logical :: found, found_aux = .false.
      type(json_file) :: json
      integer :: json_nbuff,iBuff,id
      TYPE(json_core) :: jCore
      TYPE(json_value), pointer :: buffPointer, testPointer, p
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      call json%info("buffer",n_children=json_nbuff)
      call json%get_core(jCore)
      call json%get('buffer', buffPointer, found_aux)

      if(found_aux .eqv. .true.) flag_buffer_on = .true.


      do iBuff=1, json_nbuff
            call jCore%get_child(buffPointer, iBuff, testPointer, found)
            call jCore%get_child(testPointer, 'type', p, found)
            if(found) then
               call jCore%get(p,value)
               if(value .eq. "east")  then
                  flag_buffer_on_east = .true.
                  call jCore%get_child(testPointer, 'min', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_e_min)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define east min position'
                        stop 1      
                     end if
                  end if
                  call jCore%get_child(testPointer, 'size', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_e_size)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define east size'
                        stop 1      
                     end if
                  end if
               else if (value .eq. "west") then
                  flag_buffer_on_west = .true.
                  call jCore%get_child(testPointer, 'min', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_w_min)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define west min position'
                        stop 1      
                     end if
                  end if
                  call jCore%get_child(testPointer, 'size', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_w_size)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define west size'
                        stop 1      
                     end if
                  end if
               else if (value .eq. "north") then
                  flag_buffer_on_north = .true.
                  call jCore%get_child(testPointer, 'min', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_n_min)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define north min position'
                        stop 1      
                     end if
                  end if
                  call jCore%get_child(testPointer, 'size', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_n_size)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define north size'
                        stop 1      
                     end if
                  end if
               else if (value .eq. "south") then
                  flag_buffer_on_south = .true.
                  call jCore%get_child(testPointer, 'min', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_s_min)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define south min position'
                        stop 1      
                     end if
                  end if
                  call jCore%get_child(testPointer, 'size', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_s_size)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define south size'
                        stop 1      
                     end if
                  end if
               else if (value .eq. "top") then
                  flag_buffer_on_top = .true.
                  call jCore%get_child(testPointer, 'min', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_t_min)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define top min position'
                        stop 1      
                     end if
                  end if
                  call jCore%get_child(testPointer, 'size', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_t_size)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define top size'
                        stop 1      
                     end if
                  end if
               else if (value .eq. "bottom") then
                  flag_buffer_on_bottom = .true.
                  call jCore%get_child(testPointer, 'min', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_b_min)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define bottom min position'
                        stop 1      
                     end if
                  end if
                  call jCore%get_child(testPointer, 'size', p, found)
                  if(found) then
                     call jCore%get(p,flag_buffer_b_size)
                  else
                     if(mpi_rank .eq. 0) then
                        write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define bottom size'
                        stop 1      
                     end if
                  end if
               end if
            else
               if(mpi_rank .eq. 0) then
                  write(111,*) 'ERROR! JSON file error on the buffer definition, you need to define at least one buffer type'
                  stop 1      
               end if
            end if
      end do

      call json%destroy()

   end subroutine CFDSolverBase_readJSONBuffer

   subroutine CFDSolverBase_readJSONBCTypes(this)
      use json_module
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      logical :: found, found_aux = .false.
      type(json_file) :: json
      integer :: json_nbouCodes,iBouCodes,id
      TYPE(json_core) :: jCore
      TYPE(json_value), pointer :: bouCodesPointer, testPointer, p
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      call json%info("bouCodes",n_children=json_nbouCodes)
      call json%get_core(jCore)
      call json%get('bouCodes', bouCodesPointer, found_aux)

      do iBouCodes=1, json_nbouCodes
            call jCore%get_child(bouCodesPointer, iBouCodes, testPointer, found)
            call jCore%get_child(testPointer, 'id', p, found)
            if(found) then
               call jCore%get(p,id)
            else
               if(mpi_rank .eq. 0) then
                  write(111,*) 'ERROR! JSON file error on the bouCodes definition, you need to define the id'
                  stop 1      
               end if
            end if
            call jCore%get_child(testPointer, 'bc_type', p, found)
            if(found) then
               call jCore%get(p,value)
               if(value .eq. "bc_type_far_field") then
                  bouCodes2BCType(id) = bc_type_far_field
               else if(value .eq. "bc_type_outlet_incomp") then
                  bouCodes2BCType(id) = bc_type_outlet_incomp
               else if(value .eq. "bc_type_recirculation_inlet") then
                  bouCodes2BCType(id) = bc_type_recirculation_inlet
               else if(value .eq. "bc_type_non_slip_adiabatic") then
                  bouCodes2BCType(id) = bc_type_non_slip_adiabatic
               else if(value .eq. "bc_type_non_slip_hot") then
                  bouCodes2BCType(id) = bc_type_non_slip_hot
               else if(value .eq. "bc_type_non_slip_cold") then
                  bouCodes2BCType(id) = bc_type_non_slip_cold
               else if(value .eq. "bc_type_slip_adiabatic") then
                  bouCodes2BCType(id) = bc_type_slip_adiabatic
               else if(value .eq. "bc_type_slip_wall_model") then
                  bouCodes2BCType(id) = bc_type_slip_wall_model
               else if(value .eq. "bc_type_unsteady_inlet") then
                  bouCodes2BCType(id) = bc_type_unsteady_inlet
               else if(value .eq. "bc_type_far_field_SB") then
                  bouCodes2BCType(id) = bc_type_far_field_SB 
               else if(value .eq. "bc_type_slip_wall_model_iso") then
                  bouCodes2BCType(id) = bc_type_slip_wall_model_iso                  
               else if(value .eq. "bc_type_slip_atmosphere") then
                  bouCodes2BCType(id) = bc_type_slip_atmosphere
               else if(value .eq. "bc_type_symmetry") then
                  bouCodes2BCType(id) = bc_type_symmetry
               else if(value .eq. "bc_type_slip_isothermal") then
                  bouCodes2BCType(id) = bc_type_slip_isothermal
               else if(value .eq. "bc_type_far_field_supersonic") then
                  bouCodes2BCType(id) = bc_type_far_field_supersonic
               else if(value .eq. "bc_type_outlet_supersonic") then
                  bouCodes2BCType(id) = bc_type_outlet_supersonic
               else if(value .eq. "bc_type_non_slip_unsteady") then
                  bouCodes2BCType(id) = bc_type_non_slip_unsteady
               end if
            else
               if(mpi_rank .eq. 0) then
                  write(111,*) 'ERROR! JSON file error on the bouCodes definition, you need to define the bc_type'
                  stop 1      
               end if
            end if
      end do

      !$acc update device(bouCodes2BCType(:))

      call json%destroy()

      if((found_aux .eqv. .false.) .and.mpi_rank .eq. 0) then
         write(111,*) 'ERROR! JSON file missing bouCodes definition'
         stop 1      
      end if

   end subroutine CFDSolverBase_readJSONBCTypes

   subroutine CFDSolverBase_checkFound(this,found,found_aux)
      class(CFDSolverBase), intent(inout) :: this
      logical , intent(in)    :: found
      logical , intent(inout) :: found_aux

      if(found .eqv. .false.) found_aux = .true.

   end subroutine CFDSolverBase_checkFound

   subroutine CFDSolverBase_printDt(this)
      class(CFDSolverBase), intent(inout) :: this

      write(*,*) " Dt ",this%dt
   end subroutine CFDSolverBase_printDt

   subroutine CFDSolverBase_printAll(this)
      class(CFDSolverBase), intent(inout) :: this

      call this%printDt()
   end subroutine CFDSolverBase_printAll

   subroutine CFDSolverBase_initializeSourceTerms(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iNodeL

      allocate(source_term(numNodesRankPar,ndime+2))
      !$acc enter data create(source_term(:,:))
      !$acc kernels
      source_term(:,:) = 0.00_rp
      !$acc end kernels

   end subroutine CFDSolverBase_initializeSourceTerms

   subroutine CFDSolverBase_initNSSolver(this)
      class(CFDSolverBase), intent(inout) :: this
      if(flag_implicit == 1) then
         if (implicit_solver == implicit_solver_imex) then
            call init_imex_solver(numNodesRankPar,numElemsRankPar)
         end if
      else 
         if(flag_rk_ls .eqv. .false.) then
            call init_rk4_solver(numNodesRankPar)
         else
            call init_rk4_ls_solver(numNodesRankPar)
         end if
      end if
      if(flag_use_species .eqv. .true.) then
         call init_imex_species_solver(numNodesRankPar,numElemsRankPar)
         !call init_rk4_ls_species_solver(numNodesRankPar,numElemsRankPar)
      end if

   end subroutine CFDSolverBase_initNSSolver

   subroutine CFDSolverBase_endNSSolver(this)
      class(CFDSolverBase), intent(inout) :: this

      if(flag_implicit == 1) then
         if (implicit_solver == implicit_solver_imex) then
            call end_imex_solver()
         end if
      else 
         if(flag_rk_ls .eqv. .false.) then
            call end_rk4_solver()
         else
            call end_rk4_ls_solver()
         end if
      end if

   end subroutine CFDSolverBase_endNSSolver

   subroutine CFDSolverBase_initializeDefaultFileNames(this)
      class(CFDSolverBase), intent(inout) :: this

      write(this%mesh_h5_file_path,*) "./"
      write(this%mesh_h5_file_name,*) "meshFile"

      write(this%results_h5_file_path,*) "./"
      write(this%results_h5_file_name,*) "resultsFile"

      write(this%io_prepend_path,*) "./"
      write(this%io_append_info,*) ""

   end subroutine CFDSolverBase_initializeDefaultFileNames

   subroutine CFDSolverBase_initializeDefaultFlags(this)
      class(CFDSolverBase), intent(inout) :: this

      this%loadRestartFile    = .false.
      this%saveAvgFile        = .false.
      this%loadAvgFile        = .false.
      this%continue_oldLogs   = .false.
      this%doGlobalAnalysis   = .false.
      this%doTimerAnalysis    = .false.
      this%isFreshStart       = .true.
      this%saveInitialField   = .false.
      this%saveSurfaceResults = .false.
      this%isWallModelOn      = .false.
      this%isSymmetryOn       = .false.

   end subroutine CFDSolverBase_initializeDefaultFlags

   subroutine CFDSolverBase_initializeDefaultVariables(this)
      class(CFDSolverBase), intent(inout) :: this

      this%time = 0.0_rp
      this%initial_istep = 1
      this%maxPhysTime = 1.0e6_rp

      !--------------------------------------------------------------------------
      this%save_logFile_first = 1
      this%save_logFile_step = 100
      !--------------------------------------------------------------------------
      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%restartFile_to_load = 1
      this%restartFileCnt = 1
      !--------------------------------------------------------------------------
      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 100000
      !--------------------------------------------------------------------------
      this%initial_avgTime = 0.0_rp
      this%elapsed_avgTime = 0.0_rp
      !--------------------------------------------------------------------------

   end subroutine CFDSolverBase_initializeDefaultVariables

   subroutine CFDSolverBase_initializeDefaultSaveFields(this)
      class(CFDSolverBase), intent(inout) :: this

      call initializeDefaultSaveFieldsCompressible()

   end subroutine CFDSolverBase_initializeDefaultSaveFields

   subroutine CFDSolverBase_initializeDefaultCustomSettings(this)
      class(CFDSolverBase), intent(inout) :: this

      !empty !to be overwritten

   end subroutine CFDSolverBase_initializeDefaultCustomSettings

   subroutine CFDSolverBase_initializeDefaultParameters(this)
      class(CFDSolverBase), intent(inout) :: this

      call this%initializeDefaultFileNames()

      call this%initializeDefaultVariables()

      call this%initializeDefaultFlags()

      call this%initializeDefaultSaveFields()

      call this%initializeDefaultCustomSettings()


   end subroutine CFDSolverBase_initializeDefaultParameters

!--------------------------------------------------------------------------------------------------------------------------

   subroutine CFDSolverBase_initializeParameters(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_initializeParameters

   subroutine CFDSolverBase_optimizeParameters(this)
      class(CFDSolverBase), intent(inout) :: this

      !if(flag_high_mach) then !too less disipation in complex scenarious
      !   factor_comp = 1.0_rp
      !end if

   end subroutine CFDSolverBase_optimizeParameters


   subroutine CFDSolverBase_openMesh(this)
      class(CFDSolverBase), intent(inout) :: this

      call nvtxStartRange("Open mesh")

      call set_hdf5_meshFile_name(this%mesh_h5_file_path,this%mesh_h5_file_name,mpi_size,this%meshFile_h5_full_name)
      call set_hdf5_all_resultsFile_baseName(this%results_h5_file_path,this%results_h5_file_name,this%mesh_h5_file_name,mpi_size,&
                        this%base_resultsFile_h5_full_name,this%base_avgResultsFile_h5_full_name,this%base_restartFile_h5_full_name)

      this%useIntInComms=.true.
      this%useRealInComms=.true.

      call load_hdf5_meshfile(this%meshFile_h5_full_name)

      ! init comms
      call init_comms(this%useIntInComms,this%useRealInComms,1,5)
      ! init comms boundaries
      call init_comms_bnd(this%useIntInComms,this%useRealInComms,1,5)

      if (isMeshBoundaries .and. this%saveSurfaceResults) then
         call set_hdf5_surface_meshFile_name(this%mesh_h5_file_path,this%mesh_h5_file_name,mpi_size,this%surface_meshFile_h5_full_name)
         call save_surface_mesh_hdf5_file(this%meshFile_h5_full_name,this%surface_meshFile_h5_full_name,npbou,mesh_gmsh2ij,mesh_vtk2ij)
      end if

      call nvtxEndRange

      call MPI_Barrier(app_comm,mpi_err)

      ! Initialize InSitu & correct (substract 1)  connecVTK
      call init_InSitu()

   end subroutine CFDSolverBase_openMesh

   subroutine CFDSolverBase_fill_BC_Types(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_fill_BC_Types

   subroutine checkIfWallModelOn(this)
      class(CFDSolverBase), intent(inout) :: this
      integer :: iBound,bcCode,auxBoundCnt

      this%isWallModelOn = .false.
      do iBound = 1,numBoundCodes
         if(bouCodes2BCType(iBound) .eq. bc_type_slip_wall_model) then
            this%isWallModelOn = .true.
            if(mpi_rank.eq.0) write(111,*) "--| Wall-Model activated in Boundary id",iBound
         else if(bouCodes2BCType(iBound) .eq. bc_type_slip_wall_model_iso) then
            this%isWallModelOn = .true.
            if(mpi_rank.eq.0) write(111,*) "--| Wall-Model activated in Boundary id",iBound
         end if
      end do

      numBoundsWMRankPar = 0
      do iBound = 1,numBoundsRankPar
         bcCode = bouCodes2BCType(bouCodesPar(iBound))
         if((bcCode .eq. bc_type_slip_wall_model) .or. (bcCode .eq. bc_type_slip_wall_model_iso)) then
            numBoundsWMRankPar = numBoundsWMRankPar + 1
         end if
      end do
      !write(*,*) '[',mpi_rank,'] numBoundsWMRankPar',numBoundsWMRankPar

      allocate(listBoundsWallModel(numBoundsWMRankPar))
      !$acc enter data create(listBoundsWallModel(:))
      auxBoundCnt = 0
      do iBound = 1,numBoundsRankPar
         bcCode = bouCodes2BCType(bouCodesPar(iBound))
         if((bcCode .eq. bc_type_slip_wall_model) .or. (bcCode .eq. bc_type_slip_wall_model_iso)) then
            auxBoundCnt = auxBoundCnt + 1
            listBoundsWallModel(auxBoundCnt) = iBound
         end if
      end do
      !$acc update device(listBoundsWallModel(:))

   end subroutine checkIfWallModelOn

      subroutine checkIfSymmetryOn(this)
      class(CFDSolverBase), intent(inout) :: this
      integer :: iBound,bcCode,auxBoundCnt

      this%isSymmetryOn = .false.
      do iBound = 1,numBoundCodes
         if(bouCodes2BCType(iBound) .eq. bc_type_slip_adiabatic) then
            this%isSymmetryOn = .true.
            if(mpi_rank.eq.0) write(111,*) "--| Symmetry activated in Boundary id",iBound
         end if
      end do

   end subroutine checkIfSymmetryOn

   subroutine CFDSolverBase_normalFacesToNodes(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), allocatable    :: aux1(:)
      integer(4) :: iNodeL,iBound,ipbou,iElem,innerNodeL,bndNodeL,idime,iAux
      real(rp) :: aux(3), normaux,sig

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating Normals at Nodes for Wall-Model"

      allocate(normalsAtNodes(numNodesRankPar,ndime))
      !$acc enter data create(normalsAtNodes(:,:))

      !$acc kernels
      normalsAtNodes(:,:) = 0.0_rp
      !$acc end kernels

      !$acc parallel loop gang
      !do iAux = 1,numBoundsWMRankPar
      do iAux = 1,numBoundsRankPar
         !iBound = listBoundsWallModel(iAux)
         iBound = iAux
         iElem = point2elem(boundPar(iBound,npbou)) ! I use an internal face node to be sure is the correct element
         innerNodeL = connecParWork(iElem,nnode)         ! internal node
         !$acc loop vector private(aux)
         do ipbou = 1,npbou
            bndNodeL = boundPar(iBound,ipbou) ! node at the boundary
            sig=1.0_rp
            aux(1) = boundNormalPar(iBound,(ipbou-1)*ndime+1)
            aux(2) = boundNormalPar(iBound,(ipbou-1)*ndime+2)
            aux(3) = boundNormalPar(iBound,(ipbou-1)*ndime+3)
            normaux = sqrt(dot_product(aux,aux))
            if(dot_product(coordPar(innerNodeL,:)-coordPar(bndNodeL,:), aux(:)) .lt. 0.0_rp ) then
               sig=-1.0_rp
            end if
            !$acc loop seq
            do idime = 1,ndime
               aux(idime) = aux(idime)*sig/normaux
               !$acc atomic update
               normalsAtNodes(bndNodeL,idime) = normalsAtNodes(bndNodeL,idime) + aux(idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop

      if(mpi_size.ge.2) then
         call mpi_halo_boundary_atomic_update_real(normalsAtNodes(:,1))
         call mpi_halo_boundary_atomic_update_real(normalsAtNodes(:,2))
         call mpi_halo_boundary_atomic_update_real(normalsAtNodes(:,3))
      end if

      !$acc parallel loop  private(aux)
      do iNodeL = 1,numNodesRankPar
         aux(1) = normalsAtNodes(iNodeL,1)
         aux(2) = normalsAtNodes(iNodeL,2)
         aux(3) = normalsAtNodes(iNodeL,3)
         normaux = sqrt(dot_product(aux,aux))

         if(normaux .gt. 1e-10) then
            normalsAtNodes(iNodeL,1) = aux(1)/normaux
            normalsAtNodes(iNodeL,2) = aux(2)/normaux
            normalsAtNodes(iNodeL,3) = aux(3)/normaux
         end if
      end do
      !$acc end parallel loop

      !$acc update host(normalsAtNodes(:,:)) !!just in case

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_normalFacesToNodes

   subroutine CFDSolverBase_boundaryFacesToNodes(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), allocatable    :: aux1(:)
      integer(4) :: iNodeL,iBound,ipbou

      allocate(bouCodesNodesPar(numNodesRankPar))
      allocate(aux1(numNodesRankPar))
      allocate(bouCodes2BCType(numBoundCodes))
      !$acc enter data create(bouCodesNodesPar(:))
      !$acc enter data create(aux1(:))
      !$acc enter data create(bouCodes2BCType(:))

      bouCodes2BCType(:) = 0
      !$acc update device(bouCodes2BCType(:))

      call this%fillBCTypes()

      !$acc kernels
      aux1(:) = max_num_bou_codes
      bouCodesNodesPar(:) =  max_num_bou_codes
      !$acc end kernels

      !$acc parallel loop gang
      do iBound = 1,numBoundsRankPar
         !$acc loop vector
         do ipbou = 1,npbou
            !$acc atomic update
            aux1(boundPar(iBound,ipbou)) = min(aux1(boundPar(iBound,ipbou)),bouCodes2BCType(bouCodesPar(iBound)))
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop

      if((isMeshBoundaries).and.(mpi_size.ge.2)) then
         call mpi_halo_min_boundary_update_int_iSendiRcv(aux1)
      end if

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(aux1(iNodeL) .lt. max_num_bou_codes) then
            bouCodesNodesPar(iNodeL) = aux1(iNodeL)
         end if
      end do
      !$acc end parallel loop

      call this%checkIfWallModelOn()
      call this%checkIfSymmetryOn()

      !$acc exit data delete(aux1(:))
      deallocate(aux1)

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_boundaryFacesToNodes

   subroutine CFDSolverBase_evalCharLength(this)
      class(CFDSolverBase), intent(inout) :: this
      real(rp) :: he_aux
      integer(4) :: iElem

      call nvtxStartRange("Elem size compute")
      allocate(helem(numElemsRankPar))
      !$acc enter data create(helem(:))
      do iElem = 1,numElemsRankPar
         call char_length(nnode,iElem,numElemsRankPar,numNodesRankPar,connecParOrig,coordPar,he_aux)
         helem(iElem) = he_aux
      end do
      !$acc update device(helem(:))
      call nvtxEndRange

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalCharLength

   subroutine CFDSolverBase_allocateVariables(this)
      class(CFDSolverBase), intent(inout) :: this

      if(mpi_rank.eq.0) write(111,*) "--| ALLOCATING MAIN VARIABLES"
      call nvtxStartRange("Allocate main vars")
      !
      ! Last rank is for prediction-advance related to entropy viscosity,
      ! where 1 is prediction, 2 is final value
      !
      allocate(u(numNodesRankPar,ndime,4))  ! Velocity
      allocate(q(numNodesRankPar,ndime,4))  ! momentum
      allocate(rho(numNodesRankPar,4))      ! Density
      allocate(pr(numNodesRankPar,4))       ! Pressure
      allocate(E(numNodesRankPar,4))        ! Total Energy
      allocate(Tem(numNodesRankPar,2))      ! Temperature
      allocate(e_int(numNodesRankPar,2))    ! Internal Energy
      allocate(eta(numNodesRankPar,4))      ! entropy
      allocate(csound(numNodesRankPar))     ! Speed of sound
      allocate(machno(numNodesRankPar))     ! Speed of sound
      allocate(mu_fluid(numNodesRankPar))   ! Fluid viscosity
      allocate(mu_factor(numNodesRankPar))  ! Fluid viscosity
      allocate(mu_e(numElemsRankPar,ngaus))  ! Elemental viscosity
      allocate(mu_sgs(numElemsRankPar,ngaus))! SGS viscosity
      allocate(u_buffer(numNodesRankPar,ndime))  ! momentum at the buffer
      allocate(u_mapped(numNodesRankPar,ndime))  ! velocity correction in mapped nodes
      allocate(mue_l(numElemsRankPar,nnode))
      !$acc enter data create(u(:,:,:))
      !$acc enter data create(q(:,:,:))
      !$acc enter data create(rho(:,:))
      !$acc enter data create(pr(:,:))
      !$acc enter data create(E(:,:))
      !$acc enter data create(Tem(:,:))
      !$acc enter data create(e_int(:,:))
      !$acc enter data create(eta(:,:))
      !$acc enter data create(csound(:))
      !$acc enter data create(machno(:))
      !$acc enter data create(mu_fluid(:))
      !$acc enter data create(mu_factor(:))
      !$acc enter data create(mu_e(:,:))
      !$acc enter data create(mu_sgs(:,:))
      !$acc enter data create(u_buffer(:,:))
      !$acc enter data create(u_mapped(:,:))
      !$acc enter data create(mue_l(:,:))

      allocate(tauw(numNodesRankPar,ndime))  ! momentum at the buffer
      !$acc enter data create(tauw(:,:))

      !$acc kernels
      u(:,:,:) = 0.0_rp
      q(:,:,:) = 0.0_rp
      rho(:,:) = 0.0_rp
      pr(:,:) = 0.0_rp
      E(:,:) = 0.0_rp
      Tem(:,:) = 0.0_rp
      e_int(:,:) = 0.0_rp
      eta(:,:) = 0.0_rp
      csound(:) = 0.0_rp
      machno(:) = 0.0_rp
      mu_fluid(:) = 0.0_rp
      mu_factor(:) = 1.0_rp
      mu_e(:,:) = 0.0_rp
      mu_sgs(:,:) = 0.0_rp

      u_buffer(:,:) = 0.0_rp
      u_mapped(:,:) = 0.0_rp
      tauw(:,:) = 0.0_rp
      !$acc end kernels

      !ilsa
      allocate(kres(numNodesRankPar))
      allocate(etot(numNodesRankPar))
      allocate(au(numNodesRankPar,ndime))
      allocate(ax1(numNodesRankPar))
      allocate(ax2(numNodesRankPar))
      allocate(ax3(numNodesRankPar))
      !$acc enter data create(au(:,:))
      !$acc enter data create(kres(:))
      !$acc enter data create(etot(:))
      !$acc enter data create(ax1(:))
      !$acc enter data create(ax2(:))
      !$acc enter data create(ax3(:))
      !$acc kernels
      kres(:) = 0.0_rp
      etot(:) = 0.0_rp
      au(:,:) = 0.0_rp
      ax1(:) = 0.0_rp
      ax2(:) = 0.0_rp
      ax3(:) = 0.0_rp
      !$acc end kernels

      !boundary
      if(numBoundCodes .ge. 1) then
         allocate(Fpr(ndime,numBoundCodes))
         allocate(Ftau(ndime,numBoundCodes))
         !$acc enter data create(Fpr(:,:))
         !$acc enter data create(Ftau(:,:))
         !$acc kernels
         Fpr(:,:) = 0.0_rp
         Ftau(:,:) = 0.0_rp
         !$acc end kernels
      end if

      !*********************************************************************!
      ! Derivative-related fields                                           !
      !*********************************************************************!
      allocate(gradRho(numNodesRankPar,ndime))
      allocate(curlU(numNodesRankPar,ndime))
      allocate(divU(numNodesRankPar))
      allocate(Qcrit(numNodesRankPar))
      !$acc enter data create(gradRho(:,:))
      !$acc enter data create(curlU(:,:))
      !$acc enter data create(divU(:))
      !$acc enter data create(Qcrit(:))

      !$acc kernels
      gradRho(:,:) = 0.0_rp
      curlU(:,:) = 0.0_rp
      divU(:) = 0.0_rp
      Qcrit(:) = 0.0_rp
      !$acc end kernels

      allocate(avrho(numNodesRankPar))
      allocate(avpre(numNodesRankPar))
      allocate(avpre2(numNodesRankPar))
      allocate(avmueff(numNodesRankPar))
      allocate(avvel(numNodesRankPar,ndime))
      allocate(avve2(numNodesRankPar,ndime))
      allocate(avvex(numNodesRankPar,ndime))
      allocate(avtw(numNodesRankPar,ndime))
      !$acc enter data create(avrho(:))
      !$acc enter data create(avpre(:))
      !$acc enter data create(avpre2(:))
      !$acc enter data create(avmueff(:))
      !$acc enter data create(avvel(:,:))
      !$acc enter data create(avve2(:,:))
      !$acc enter data create(avvex(:,:))
      !$acc enter data create(avtw(:,:))

      !$acc kernels
      avrho(:) = 0.0_rp
      avpre(:) = 0.0_rp
      avpre2(:) = 0.0_rp
      avmueff(:) = 0.0_rp
      avvel(:,:) = 0.0_rp
      avve2(:,:) = 0.0_rp
      avvex(:,:) = 0.0_rp
      avtw(:,:) = 0.0_rp
      !$acc end kernels

      if (this%have_witness) then
         allocate(witel(this%nwit))
         allocate(witxi(this%nwit,ndime))
         allocate(Nwit(this%nwit,nnode))
      end if

      ! Exponential average velocity for wall law
      if(flag_walave) then
         allocate(walave_u(numNodesRankPar,ndime))
         !$acc enter data create(walave_u(:,:))
         !$acc kernels
         walave_u(:,:) = 0.0_rp
         !$acc end kernels

         if(flag_use_species) then
            allocate(walave_t(numNodesRankPar))
            !$acc enter data create(walave_t(:))
            !$acc kernels
            walave_t(:) = 0.0_rp
            !$acc end kernels
         end if
      else
         allocate(walave_u(0,0)) !dummy allocation
         allocate(walave_t(0)) !dummy allocation
      end if

      if(flag_type_wmles==wmles_type_abl) then
         allocate(zo(numNodesRankPar))
         allocate(ad(numNodesRankPar))
         !$acc enter data create(zo(:))
         !$acc enter data create(ad(:))
         !$acc kernels
         zo(:) = 0.0_rp
         ad(:) = 0
         !$acc end kernels
      end if

      allocate(maskMapped(numNodesRankPar))
      !$acc enter data create(maskMapped(:))
      !$acc kernels
      maskMapped(:) = 1 !1 for calculation domain / 0 for recirculating domain where buffer is not applied
      !$acc end kernels


      if(flag_use_species .eqv. .true.) then
         allocate(eta_Yk(numNodesRankPar,nspecies,4))      ! entropy
         allocate(mu_e_Yk(numElemsRankPar,ngaus,nspecies))  ! Elemental viscosity
         allocate(Yk(numNodesRankPar,nspecies,4))! SGS viscosity
         allocate(Yk_buffer(numNodesRankPar,nspecies))

         !$acc enter data create(eta_Yk(:,:,:))
         !$acc enter data create(mu_e_Yk(:,:,:))
         !$acc enter data create(Yk(:,:,:))
         !$acc enter data create(Yk_buffer(:,:))
         
         !$acc kernels
         eta_Yk(:,:,:) = 0.0_rp
         mu_e_Yk(:,:,:) = 0.0_rp
         Yk(:,:,:) = 0.0_rp
         Yk_buffer(:,:) = 0.0_rp
         !$acc end kernels
      end if


      call nvtxEndRange

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_allocateVariables

   subroutine CFDSolverBase_deallocateVariables(this)
      class(CFDSolverBase), intent(inout) :: this

      if(mpi_rank.eq.0) write(111,*) "--| DEALLOCATING MAIN VARIABLES"
      call nvtxStartRange("Deallocate main vars")

      !TO BE COMPLETED! NOT STRICTLY NECESSARY BUT IS GOOD TO DO IT AS GOOD PROGRAMMING PRACTICE :)

      call nvtxEndRange

   end subroutine CFDSolverBase_deallocateVariables

   subroutine CFDSolverBase_evalOrLoadInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

      this%save_logFile_next = this%save_logFile_first
      this%save_restartFile_next = this%save_restartFile_first
      this%save_resultsFile_next = this%save_resultsFile_first

      if(this%loadRestartFile) then
         if(mpi_rank.eq.0) write(111,*) "--| Loading restart file ",this%restartFile_to_load
         call load_hdf5_restartFile(this%base_restartFile_h5_full_name,nnode,ngaus,this%restartFile_to_load,this%load_step,flag_walave,this%time,&
                                    rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,mu_sgs,walave_u,&
                                    numRestartFields,restartNameFields,indRF_rho,indRF_ux,indRF_uy,indRF_uz,indRF_pr,&
                                    indRF_ener,indRF_mut,indRF_mue,indRF_walavex,indRF_walavey,indRF_walavez)

         if((flag_les.eq.0).and.(flag_les_ilsa.eq.0)) then
            !$acc kernels
            mu_sgs(:,:) = 0._rp
            !$acc end kernels
         end if

         if(mpi_rank.eq.0) write(111,*) "   --| Loaded results for iStep",this%load_step,"time",this%time
         call this%eval_vars_after_load_hdf5_resultsFile()

         if(this%continue_oldLogs) then
            this%initial_istep = this%load_step+1

            do while(this%save_logFile_next .le. this%load_step)
               this%save_logFile_next = this%save_logFile_next + this%save_logFile_step
            end do

            do while(this%save_restartFile_next .le. this%load_step)
               this%save_restartFile_next = this%save_restartFile_next + this%save_restartFile_step
            end do

            do while(this%save_resultsFile_next .le. this%load_step)
               this%save_resultsFile_next = this%save_resultsFile_next + this%save_resultsFile_step
            end do

            if(this%loadAvgFile) then
               if(mpi_rank.eq.0) write(111,*) "--| Loading Avg Results File (TO IMPLEMENT)",this%restartFile_to_load
               call load_avgResults_hdf5_file(this%base_avgResultsFile_h5_full_name,nnode,ngaus,Ngp_l,this%restartFile_to_load,&
                                       this%initial_avgTime,this%elapsed_avgTime,&
                                       numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
                                       numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
                                       numAvgElGPScalarFields2save,avgElGPScalarFields2save)

               if(mpi_rank.eq.0) write(111,*) "   --| Loaded Avg results! Setting initial_avgTime",this%initial_avgTime,"elapsed_avgTime",this%elapsed_avgTime
            end if

            if(mpi_rank.eq.0) then
               write(111,*) "   --| Continuing old logs..."
               write(111,*) "   --| initial_istep",this%initial_istep
               write(111,*) "   --| save_logFile_next",this%save_logFile_next
               write(111,*) "   --| save_restartFile_next",this%save_restartFile_next
               write(111,*) "   --| save_resultsFile_next",this%save_resultsFile_next
            end if



            this%isFreshStart = .false.
         else
            this%time = 0.0_rp
         end if
      else
         if(mpi_rank.eq.0) write(111,*) "--| Evaluating Initial Conditions..."
         call this%evalInitialConditions()
      end if

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalOrLoadInitialConditions

   subroutine CFDSolverBase_evalInitialConditions(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_evalInitialConditions

   subroutine CFDSolverBase_evalInitialViscosity(this)
      class(CFDSolverBase), intent(inout) :: this

        if (flag_real_diff == 1) then
           if (flag_diff_suth == 0) then
              call constant_viscosity(numNodesRankPar,incomp_viscosity,mu_fluid)
           else
              call sutherland_viscosity(numNodesRankPar,Tem(:,2),mu_factor,mu_fluid)
           end if
        else if (flag_real_diff == 0) then
           !$acc kernels
           mu_fluid(:) = 0.0_rp
           !$acc end kernels
        else
           if(mpi_rank.eq.0) write(111,*) "--| DIFFUSION FLAG MUST BE EITHER 0 OR 1, NOT: ",flag_real_diff
           stop 1
        end if

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalInitialViscosity

   subroutine CFDSolverBase_evalViscosityFactor(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iNodeL

      ! set out of the buffer zone
      ! remember that the mu_factor field has to we filled at least with the
      ! flag_mu_factor

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

   end subroutine CFDSolverBase_evalViscosityFactor

   subroutine eval_initial_mu_sgs(this)
      class(CFDSolverBase), intent(inout) :: this

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating initial mu_sgs..."

      call nvtxStartRange("MU_SGS")
      if(flag_les_ilsa == 1) then
         this%dt = 1.0_rp !To avoid 0.0 division inside sgs_ilsa_visc calc
         call sgs_ilsa_visc(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,this%dt,rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3,mue_l)
      else
         call sgs_visc(numElemsRankPar,numNodesRankPar,connecParWork,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),Ml,mu_sgs,mue_l)
      end if
      call nvtxEndRange

   end subroutine eval_initial_mu_sgs

   subroutine CFDSolverBase_evalInitialDt(this)
      class(CFDSolverBase), intent(inout) :: this

      !*********************************************************************!
      ! Compute initial time-step size                                      !
      !*********************************************************************!

      if(mpi_rank.eq.0) write(111,*) "--| Evaluating initial dt..."
      !if (flag_real_diff == 1) then
         call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,mu_e,rho(:,2),this%Cp,this%Prt)
      !else
      !   call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
      !end if
      if(mpi_rank.eq.0) write(111,*) "--| Initial time-step dt := ",this%dt,"s"

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalInitialDt

   subroutine CFDSolverBase_evalDt(this)
      class(CFDSolverBase), intent(inout) :: this

      !if (flag_real_diff == 1) then
         call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt,this%cfl_diff,mu_fluid,mu_sgs,mu_e,rho(:,2),this%Cp,this%Prt)
      !else
      !   call adapt_dt_cfl(numElemsRankPar,numNodesRankPar,connecParWork,helem,u(:,:,2),csound,this%cfl_conv,this%dt)
      !end if

   end subroutine CFDSolverBase_evalDt

   subroutine CFDSolverBase_evalShapeFunctions(this)
      class(CFDSolverBase), intent(inout) :: this
      real(8)   :: s,t,z,xi_gll(porder+1),xgp_equi(ngaus,ndime)
      integer(4) :: igaus,ii
      real(8) :: xgp_r8(ngaus,ndime),wgp_r8(ngaus),xgp_b_r8(npbou,ndime-1), wgp_b_r8(npbou)
      real(8) :: Ngp_equi_r8(ngaus,nnode), dNgp_equi_r8(ndime,nnode,ngaus)
      real(8) :: Ngp_r8(ngaus,nnode),dNgp_r8(ndime,nnode,ngaus),Ngp_b_r8(npbou,npbou),dNgp_b_r8(ndime-1,npbou,npbou)
      real(8) :: Ngp_l_r8(ngaus,nnode),dNgp_l_r8(ndime,nnode,ngaus),dlxigp_ip_r8(ngaus,ndime,porder+1)


      !*********************************************************************!
      ! Generate GLL table                                                  !
      !*********************************************************************!
      if(mpi_rank.eq.0) write(111,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."
      call nvtxStartRange("Gaussian Quadrature")

      !*********************************************************
      !           Allocating required arrays!
      allocate(atoIJ(npbou))
      allocate(atoIJK(nnode))

      allocate(xgp(ngaus,ndime),wgp(ngaus))
      !$acc enter data create(wgp(:))
      allocate(xgp_b(npbou,ndime-1),wgp_b(npbou))
      !$acc enter data create(wgp_b(:))

      allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))
      allocate(Ngp_l(ngaus,nnode),dNgp_l(ndime,nnode,ngaus))
      allocate(Ngp_b(npbou,npbou),dNgp_b(ndime-1,npbou,npbou))
      allocate(dlxigp_ip(ngaus,ndime,porder+1))
      !$acc enter data create(Ngp(:,:))
      !$acc enter data create(dNgp(:,:,:))
      !$acc enter data create(Ngp_l(:,:))
      !$acc enter data create(dNgp_l(:,:,:))
      !$acc enter data create(Ngp_b(:,:))
      !$acc enter data create(dNgp_b(:,:,:))
      !$acc enter data create(dlxigp_ip(:,:,:))
      !*********************************************************

      xi_gll(:) = 0.0d0 
      xgp_equi(:,:) = 0.0d0 
      xgp_r8(:,:) = 0.0d0 
      wgp_r8(:) = 0.0d0 
      xgp_b_r8(:,:) = 0.0d0 
      wgp_b_r8(:) = 0.0d0 
      Ngp_equi_r8(:,:) = 0.0d0 
      dNgp_equi_r8(:,:,:) = 0.0d0 
      Ngp_r8(:,:) = 0.0d0 
      dNgp_r8(:,:,:) = 0.0d0 
      Ngp_b_r8(:,:) = 0.0d0 
      dNgp_b_r8(:,:,:) = 0.0d0 
      Ngp_l_r8(:,:) = 0.0d0 
      dNgp_l_r8(:,:,:) = 0.0d0 
      dlxigp_ip_r8(:,:,:) = 0.0d0 

      !*********************************************************

      atoIJK(:) = mesh_a2ijk(:)
      atoIJ(:)  = mesh_a2ij(:)

      if(mpi_rank.eq.0) write(111,*) "  --| Generating Gauss-Lobatto-Legendre table..."

      call GaussLobattoLegendre_hex(porder,ngaus,atoIJK,xgp_r8,wgp_r8)
      xgp(:,:) = real(xgp_r8(:,:),rp)
      wgp(:) = real(wgp_r8,rp)
      !$acc update device(wgp(:))

      call GaussLobattoLegendre_qua(porder,npbou,atoIJ,xgp_b_r8,wgp_b_r8)
      xgp_b(:,:) = real(xgp_b_r8(:,:),rp)
      wgp_b(:) = real(wgp_b_r8(:),rp)
      !$acc update device(wgp_b(:))

      !-------------------------------------------------------------------------------
      ! Generating Ngp_equi to interpolate from GLL nodes mesh to Equispace nodes mesh
      allocate(Ngp_equi(ngaus,nnode))
      !$acc enter data create(Ngp_equi(:,:))

      call getGaussPoints_equispaced_hex(porder,ngaus,atoIJK,xgp_equi)
      call getGaussLobattoLegendre_roots(porder,xi_gll)

      do igaus = 1,ngaus
         s = xgp_equi(igaus,1)
         t = xgp_equi(igaus,2)
         z = xgp_equi(igaus,3)

         call TripleTensorProduct(porder,nnode,xi_gll,s,t,z,atoIJK,Ngp_equi_r8(igaus,:),dNgp_equi_r8(:,:,igaus))
      end do
      Ngp_equi(:,:) = real(Ngp_equi_r8(:,:),rp)
      !$acc update device(Ngp_equi(:,:))

      !-------------------------------------------------------------------------------

      call nvtxEndRange

      !*********************************************************************!
      ! Generate N and dN for all GP                                        !
      !*********************************************************************!

      ! TODO: Allow for more element types

      if(mpi_rank.eq.0) write(111,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."
      call nvtxStartRange("N and dN")

      do igaus = 1,ngaus
         s = xgp(igaus,1)
         t = xgp(igaus,2)
         z = xgp(igaus,3)
         call hex_highorder(porder,nnode,s,t,z,atoIJK,Ngp_r8(igaus,:),dNgp_r8(:,:,igaus),Ngp_l_r8(igaus,:),dNgp_l_r8(:,:,igaus),dlxigp_ip_r8(igaus,:,:))
      end do
      Ngp(:,:)         = real(Ngp_r8(:,:),rp)
      dNgp(:,:,:)      = real(dNgp_r8(:,:,:),rp)
      Ngp_l(:,:)       = real(Ngp_l_r8(:,:),rp)
      dNgp_l(:,:,:)    = real(dNgp_l_r8(:,:,:),rp)
      dlxigp_ip(:,:,:) = real(dlxigp_ip_r8(:,:,:),rp)
      !$acc update device(Ngp(:,:))
      !$acc update device(dNgp(:,:,:))
      !$acc update device(Ngp_l(:,:))
      !$acc update device(dNgp_l(:,:,:))
      !$acc update device(dlxigp_ip(:,:,:))
      !
      ! Compute N andd dN for boundary elements
      !
      do igaus = 1,npbou
         s = xgp_b(igaus,1)
         t = xgp_b(igaus,2)
         call quad_highorder(porder,npbou,s,t,atoIJ,Ngp_b_r8(igaus,:),dNgp_b_r8(:,:,igaus))
      end do
      Ngp_b(:,:)    = Ngp_b_r8(:,:)
      dNgp_b(:,:,:) = dNgp_b_r8(:,:,:)
      !$acc update device(Ngp_b(:,:))
      !$acc update device(dNgp_b(:,:,:))

      call nvtxEndRange

      !
      ! Compute Levi-Civita tensor
      !
      this%leviCivi = 0.0_rp
      this%leviCivi(2,3,1) =  1.0_rp
      this%leviCivi(3,2,1) = -1.0_rp
      this%leviCivi(1,3,2) = -1.0_rp
      this%leviCivi(3,1,2) =  1.0_rp
      this%leviCivi(1,2,3) =  1.0_rp
      this%leviCivi(2,1,3) = -1.0_rp

      !
      ! Compute al,am,an weights and convertIJK
      !

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalShapeFunctions

   subroutine CFDSolverBase_evalBoundaryNormals(this)
      class(CFDSolverBase), intent(inout) :: this

      if (isMeshBoundaries) then
         if(mpi_rank.eq.0) write(111,*) "--| COMPUTING BOUNDARY ELEMENT NORMALS"
         allocate(boundNormalPar(numBoundsRankPar,ndime*npbou))
         !$acc enter data create(boundNormalPar(:,:))
         call nvtxStartRange("Bou normals")
         call boundary_normals(npbou,numNodesRankPar,numBoundsRankPar,boundParOrig,this%leviCivi,coordPar,dNgp_b,boundNormalPar)
         call nvtxEndRange
         !$acc update device(boundNormalPar(:,:))
      end if

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalBoundaryNormals

   subroutine CFDSolverBase_evalJacobians(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: ielem, igaus,iElemG
      real(8) :: vol_rank, vol_tot_d

      !*********************************************************************!
      ! Generate Jacobian related information                               !
      !*********************************************************************!

      if(mpi_rank.eq.0) write(111,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

      call nvtxStartRange("Jacobian info")
      allocate(He(ndime,ndime,ngaus,numElemsRankPar))
      allocate(gpvol(1,ngaus,numElemsRankPar))
      !$acc enter data create(He(:,:,:,:))
      !$acc enter data create(gpvol(:,:,:))

      call elem_jacobian(numElemsRankPar,numNodesRankPar,connecParOrig,coordPar,dNgp,wgp,gpvol,He)
      call  nvtxEndRange
      vol_rank  = 0.0d8
      vol_tot_d = 0.0d8
      !$acc parallel loop reduction(+:vol_rank)
      do ielem = 1,numElemsRankPar
         !$acc loop vector
         do igaus = 1,ngaus
            vol_rank = vol_rank+real(gpvol(1,igaus,ielem),8)
         end do
      end do
      !$acc end parallel loop

      call MPI_Allreduce(vol_rank,vol_tot_d,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

      this%VolTot = real(vol_tot_d,rp)

      call MPI_Barrier(app_comm,mpi_err)
      if(mpi_rank.eq.0) write(111,*) '--| DOMAIN VOLUME := ',this%VolTot

   end subroutine CFDSolverBase_evalJacobians

   subroutine CFDSolverBase_evalAtoIJKInverse(this)
      class(CFDSolverBase), intent(inout) :: this

      allocate(invAtoIJK(porder+1,porder+1,porder+1))
      allocate(gmshAtoI(nnode))
      allocate(gmshAtoJ(nnode))
      allocate(gmshAtoK(nnode))
      !$acc enter data create(invAtoIJK(:,:,:))
      !$acc enter data create(gmshAtoI(:))
      !$acc enter data create(gmshAtoJ(:))
      !$acc enter data create(gmshAtoK(:))

      call atoIJKInverse(porder,nnode,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)
      !$acc update device(invAtoIJK(:,:,:))
      !$acc update device(gmshAtoI(:))
      !$acc update device(gmshAtoJ(:))
      !$acc update device(gmshAtoK(:))

   end subroutine CFDSolverBase_evalAtoIJKInverse

   subroutine CFDSolverBase_eval_elemPerNode_and_nearBoundaryNode(this)
      class(CFDSolverBase), intent(inout) :: this

      !*********************************************************************!
      ! Compute list of elements per node (connectivity index)              !
      !*********************************************************************!
      ! evaluate near boundaries for the inlets and outlets
      !not the best place Oriol!
      if(mpi_rank.eq.0) write(111,*) "--| Doing near boundary calculations..."
      allocate(lelpn(numNodesRankPar))
      !$acc enter data create(lelpn(:))
      allocate(point2elem(numNodesRankPar))
      !$acc enter data create(point2elem(:))
      if(mpi_rank.eq.0) write(111,*) '  --| Evaluating point2elem array...'
      call elemPerNode(nnode,numElemsRankPar,numNodesRankPar,connecParWork,lelpn,point2elem)

      if(mpi_rank.eq.0) write(111,*) '  --| Evaluating lnbnNodes arrays...'
      allocate(lnbnNodes(numNodesRankPar))
      !$acc enter data create(lnbnNodes(:))
      call nearBoundaryNode(porder,nnode,npbou,numElemsRankPar,numNodesRankPar,numBoundsRankPar,connecParWork,coordPar,boundPar,bouCodesNodesPar,point2elem,atoIJK,lnbnNodes)

      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_eval_elemPerNode_and_nearBoundaryNode

   subroutine CFDSolverBase_set_mappedFaces_linkingNodes(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iNode,iNodePer,iNodeMap,perMapFaceAbsDir

      if(.not.isMappedFaces) return

      if(mpi_rank.eq.0) write(111,*) '  --| Setting mapped faces linking nodes in lbnbNodes array...'

		!$acc parallel loop  
		do iNode = 1,numPerMapLinkedNodesRankPar
         iNodePer = perMapLinkedNodesRankPar(iNode,1)
         iNodeMap = perMapLinkedNodesRankPar(iNode,2)
         !if(mpi_rank.eq.0) write(111,*) 'i',iNode,'iNodeMap',iNodeMap,'xyz',coordPar(iNodeMap,:)
         !if(mpi_rank.eq.0) write(111,*) 'i',iNode,'iNodePer',iNodePer,'xyz',coordPar(iNodePer,:)
         lnbnNodes(iNodeMap) = iNodePer                  
		end do
		!$acc end parallel loop

      !---------------------------------------------------------------------------------------------------------------------

      if(mpi_rank.eq.0) write(111,*) '  --| Applying Mask Mapped from pos',perMapFaceGapCoord,'at dir',perMapFaceDir
      !maskMapped by default is 1, and here we set it to 0 in the precursor domain so the buffer is not applied
      !i.e.: if maskMapped=1: buffer IS applied // if maskMapped=0: buffer IS NOT applied

      perMapFaceAbsDir = abs(perMapFaceDir)
      !$acc parallel loop
      do iNode = 1,numNodesRankPar
         if(perMapFaceDir.ge.0) then !positive x,y,z
            if (coordPar(iNode,perMapFaceAbsDir) .le. perMapFaceGapCoord) maskMapped(iNode) = 0
         else !negative x,y,z
            if (coordPar(iNode,perMapFaceAbsDir) .ge. perMapFaceGapCoord) maskMapped(iNode) = 0
         end if
      end do
      !$acc end parallel loop

      !$acc update host(maskMapped(:))

   end subroutine CFDSolverBase_set_mappedFaces_linkingNodes

   subroutine CFDSolverBase_evalMass(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: iElem

      !*********************************************************************!
      ! Compute mass matrix (Lumped and Consistent) and set solver type     !
      !*********************************************************************!

      if(mpi_rank.eq.0) write(111,*) '--| COMPUTING LUMPED MASS MATRIX...'
      call nvtxStartRange("Lumped mass compute")
      allocate(Ml(numNodesRankPar))
      !$acc enter data create(Ml(:))
      call lumped_mass_spectral(numElemsRankPar,numNodesRankPar,connecParWork,gpvol,Ml)
      call nvtxEndRange

      !charecteristic length for spectral elements for the entropy
      !stablisation
      allocate(helem_l(numElemsRankPar,nnode))
      !$acc enter data create(helem_l(:,:))
      do iElem = 1,numElemsRankPar
         call char_length_spectral(nnode,iElem,numElemsRankPar,numNodesRankPar,connecParOrig,coordPar,Ml,helem_l)
      end do
      !$acc update device(helem_l(:,:))
      call MPI_Barrier(app_comm,mpi_err)

   end subroutine CFDSolverBase_evalMass

   subroutine CFDSolverBase_evalFirstOutput(this)
      class(CFDSolverBase), intent(inout) :: this
      character(500) :: tmpname
      integer(4) :: iCode

      if(this%saveInitialField) then
         call this%saveInstResultsFiles(0)
      end if
      !*********************************************************************!
      ! Compute surface forces and area                                                                !
      !*********************************************************************!
      if (isMeshBoundaries) then
         do iCode = 1,numBoundCodes
            call nvtxStartRange("Surface info")
            call surfInfo(0,0.0_rp,numElemsRankPar,numNodesRankPar,numBoundsRankPar,iCode,connecParWork,boundPar,point2elem,&
               bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
               mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(:,iCode),Ftau(:,iCode),.TRUE.)
            call nvtxEndRange
         end do
      end if

      call compute_fieldDerivs(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)
      call volAvg_EK(numElemsRankPar,numNodesRankPar,connecParWork,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
      call visc_dissipationRate(numElemsRankPar,numNodesRankPar,connecParWork,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
      call maxMach(numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,machno,this%maxmachno)
      call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
      if(mpi_rank.eq.0) then
         write(111,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
         write(111,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
         20 format(6(F16.8,2X))
      end if

   end subroutine CFDSolverBase_evalFirstOutput

   subroutine CFDSolverBase_callTimeIntegration(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in) :: istep

      if(mpi_rank.eq.0) write(111,*) " Time integration should be overwritted"
      stop 1

   end subroutine CFDSolverBase_callTimeIntegration

   subroutine CFDSolverBase_saveRestartFile(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in) :: istep
      
      call save_hdf5_restartFile(this%base_restartFile_h5_full_name,nnode,ngaus,this%restartFileCnt,istep,flag_walave,this%time,&
                                 rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,mu_sgs,walave_u,&
                                 numRestartFields,restartNameFields,indRF_rho,indRF_ux,indRF_uy,indRF_uz,indRF_pr,&
                                 indRF_ener,indRF_mut,indRF_mue,indRF_walavex,indRF_walavey,indRF_walavez)

      if(this%restartFileCnt .eq. 1) then
         this%restartFileCnt = 2
      else if(this%restartFileCnt .eq. 2) then
         this%restartFileCnt = 1
      else
         if(mpi_rank.eq.0) write(111,*) 'Wrong value in restartFileCnt! Setting it to default value 1'
         this%restartFileCnt = 1
      end if

   end subroutine CFDSolverBase_saveRestartFile

   subroutine CFDSolverBase_saveAvgResultsFiles(this)
      class(CFDSolverBase), intent(inout) :: this

      call save_avgResults_hdf5_file(this%base_avgResultsFile_h5_full_name,this%meshFile_h5_full_name,nnode,ngaus,Ngp_equi,&
               this%restartFileCnt,this%initial_avgTime,this%elapsed_avgTime,&
               numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
               numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
               numAvgElGPScalarFields2save,avgElGPScalarFields2save)

      if (isMeshBoundaries .and. this%saveSurfaceResults) then
         call save_surface_avgResults_hdf5_file(this%base_avgResultsFile_h5_full_name,this%meshFile_h5_full_name,&
                  this%surface_meshFile_h5_full_name,this%restartFileCnt,&
                  numAvgNodeScalarFields2save,avgNodeScalarFields2save,&
                  numAvgNodeVectorFields2save,avgNodeVectorFields2save,&
                  numAvgElGPScalarFields2save,avgElGPScalarFields2save)
      end if

   end subroutine CFDSolverBase_saveAvgResultsFiles

   subroutine CFDSolverBase_beforeTimeIteration(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_beforeTimeIteration

   subroutine CFDSolverBase_afterTimeIteration(this)
      class(CFDSolverBase), intent(inout) :: this

   end subroutine CFDSolverBase_afterTimeIteration

   subroutine CFDSolverBase_saveInstResultsFiles(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in) :: istep

      call save_instResults_hdf5_file(this%base_resultsFile_h5_full_name,this%meshFile_h5_full_name,&
               nnode,ngaus,Ngp_equi,iStep,this%time,&
               numNodeScalarFields2save,nodeScalarFields2save,&
               numNodeVectorFields2save,nodeVectorFields2save,&
               numElGPScalarFields2save,elGPScalarFields2save)

      if (isMeshBoundaries .and. this%saveSurfaceResults) then
         call save_surface_instResults_hdf5_file(this%base_resultsFile_h5_full_name,this%meshFile_h5_full_name,this%surface_meshFile_h5_full_name,istep,&
               numNodeScalarFields2save,nodeScalarFields2save,&
               numNodeVectorFields2save,nodeVectorFields2save,&
               numElGPScalarFields2save,elGPScalarFields2save)
      end if

   end subroutine CFDSolverBase_saveInstResultsFiles

   subroutine CFDSolverBase_afterDt(this,istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in) :: istep

   end subroutine CFDSolverBase_afterDt

   subroutine CFDSolverBase_initialBuffer(this)
      class(CFDSolverBase), intent(inout) :: this

      !$acc kernels
      u_buffer(:,1) = nscbc_u_inf
      !$acc end kernels

      if(flag_use_species .eqv. .true.) then
            !$acc kernels
            Yk_buffer(:,:) = 0.0_rp
            !$acc end kernels
      end if

   end subroutine CFDSolverBase_initialBuffer

   subroutine CFDSolverBase_evalTimeIteration(this)
      class(CFDSolverBase), intent(inout) :: this
      integer(4) :: icode,istep,inonLineal,iwitstep=0
      character(4) :: timeStep
      real(8) :: iStepTimeRank,iStepTimeMax,iStepEndTime,iStepStartTime,iStepAvgTime
      real(rp) :: inv_iStep
      real(rp) :: dtfact,avwei

      call MPI_Barrier(app_comm,mpi_err)

      call nvtxStartRange("Start time-integration")
      if(mpi_rank.eq.0) then
         write(*,*) 'Strarting evalTimeItarion! All info will be written in the log file: ',this%log_file_name
         write(111,*) 'Doing evalTimeIteration. Ini step:',this%initial_istep,'| End step:',this%final_istep
      end if

      call this%initNSSolver()

      this%local_step=1
      do istep = this%initial_istep,this%final_istep
         !if (istep==this%nsave.and.mpi_rank.eq.0) write(111,*) '   --| STEP: ', istep
         call nvtxStartRange("Init pred")
         !$acc kernels
         rho(:,1) = rho(:,2)
         u(:,:,1) = u(:,:,2)
         q(:,:,1) = q(:,:,2)
         pr(:,1) = pr(:,2)
         E(:,1) = E(:,2)
         Tem(:,1) = Tem(:,2)
         e_int(:,1) = e_int(:,2)
         eta(:,1) = eta(:,2)
         !$acc end kernels
         if(flag_use_species .eqv. .true.) then
            !$acc kernels      
               eta_Yk(:,:,1) = eta_Yk(:,:,2)
               Yk(:,:,1) = Yk(:,:,2)
            !$acc end kernels
         end if
         call nvtxEndRange

         !
         ! Exponential averaging for wall law
         !
         call nvtxStartRange("Wall Average")
         if(flag_walave) then
            !
            ! outside acc kernels following pseudo_cfl in next loop
            !
            dtfact = this%dt/(this%dt+period_walave)
            avwei  = 1.0_rp - dtfact
            !$acc kernels
            walave_u(:,:) = dtfact*u(:,:,2) + avwei*walave_u(:,:)
            !$acc end kernels

            if(flag_use_species) then
               !$acc kernels
               walave_t(:) = dtfact*(Yk(:,1,2) - Yk_buffer(:,1)) + avwei*walave_t(:)
               !$acc end kernels
            end if

         end if
         call nvtxEndRange

         call nvtxStartRange("Time-step")

         if(this%doTimerAnalysis) iStepStartTime = MPI_Wtime()
         call this%callTimeIntegration(istep)

         !$acc kernels
         rho(:,4) = rho(:,3)
         E(:,4) = E(:,3)
         q(:,:,4) = q(:,:,3)
         eta(:,4) = eta(:,3)
         u(:,:,4) = u(:,:,3)
         pr(:,4) = pr(:,3)

         rho(:,3) = rho(:,1)
         E(:,3) = E(:,1)
         q(:,:,3) = q(:,:,1)
         eta(:,3) = eta(:,1)
         u(:,:,3) = u(:,:,1)
         pr(:,3) = pr(:,1)
         !$acc end kernels

         if(flag_use_species .eqv. .true.) then
            !$acc kernels      
               eta_Yk(:,:,4) = eta_Yk(:,:,3)
               Yk(:,:,4) = Yk(:,:,3)

               eta_Yk(:,:,3) = eta_Yk(:,:,1)
               Yk(:,:,3) = Yk(:,:,1)
            !$acc end kernels
         end if

         if(this%doTimerAnalysis) then
            iStepEndTime = MPI_Wtime()
            iStepTimeRank = iStepEndTime - iStepStartTime
            call MPI_Allreduce(iStepTimeRank,iStepTimeMax,1,mpi_datatype_real8,MPI_MAX,app_comm,mpi_err)
            inv_iStep = 1.0_rp/real(istep)
            iStepAvgTime = (iStepAvgTime*(istep-1)+iStepTimeMax)*inv_iStep

            if((mpi_rank.eq.0).and.(this%save_logFile_next==istep)) then
               write(123,*) istep,iStepTimeMax,iStepAvgTime
               call flush(123)
            end if
         end if

         this%time = this%time+this%dt
         global_time = this%time

         if ((this%save_logFile_next==istep) .and. (this%doGlobalAnalysis)) then
            call volAvg_EK(numElemsRankPar,numNodesRankPar,connecParWork,gpvol,Ngp,nscbc_rho_inf,rho(:,2),u(:,:,2),this%EK)
            call visc_dissipationRate(numElemsRankPar,numNodesRankPar,connecParWork,this%leviCivi,nscbc_rho_inf,mu_fluid,mu_e,u(:,:,2),this%VolTot,gpvol,He,dNgp,this%eps_S,this%eps_D,this%eps_T)
            call maxMach(numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,machno,this%maxmachno)
            call write_EK(this%time,this%EK,this%eps_S,this%eps_D,this%eps_T,this%maxmachno)
            if(mpi_rank.eq.0) then
               write(111,*) "--| time     EK     eps_S     eps_D     eps_T     max(Ma)"
               write(111,20) this%time, this%EK, this%eps_S, this%eps_D, this%eps_T, this%maxmachno
               20 format(6(F16.8,2X))
               call flush(666)
            end if
         end if

         call this%evalDt()

         call this%afterDt(istep)
         call nvtxEndRange

         if(this%saveAvgFile) then
            ! Update the accumulators for averaging
            if(this%time .ge. this%initial_avgTime) then
               call nvtxStartRange("Accumulate")

               call eval_average_iter(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,this%dt,this%elapsed_avgTime,&
                                 rho,u,pr,mu_fluid,mu_e,mu_sgs,tauw,avrho,avpre,avpre2,avvel,avve2,avvex,avmueff,avtw)
               call nvtxEndRange
            end if
         end if

         if(this%save_logFile_next==istep) then
            if(mpi_rank.eq.0) write(111,*) "step",istep,"time:",this%time,"s (dt",this%dt,"s)"

            if (isMeshBoundaries) then
               do icode = 1,numBoundCodes
                  call nvtxStartRange("Surface info")
                  call surfInfo(istep,this%time,numElemsRankPar,numNodesRankPar,numBoundsRankPar,icode,connecParWork,boundPar,point2elem, &
                     bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
                     mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(:,iCode),Ftau(:,iCode),.TRUE.)

                  call nvtxEndRange
                  if(mpi_rank.eq.0) call flush(888+icode)
               end do
            end if
         end if

         ! ---- SAVING RESTART FILE ---------------------------------------------------------------------------
         if (this%save_restartFile_next == istep) then
            this%save_restartFile_next = this%save_restartFile_next + this%save_restartFile_step
            if (mpi_rank.eq.0) write(111,*) ' - Saving restart file',this%restartFileCnt,'in step',istep,'(next to save',this%save_restartFile_next,')'

            if(this%saveAvgFile) then
               if (mpi_rank.eq.0) write(111,*) '   - Saving avgResults file step:',this%restartFileCnt
               call nvtxStartRange("Output AVG")
               call this%saveAvgResultsFiles
               call nvtxEndRange
            end if

            call nvtxStartRange("Saving_restart_file")
            call this%saveRestartFile(istep)
            call nvtxEndRange

         end if

         ! ---- SAVING INST RESULTS FILE -----------------------------------------------------------------------
         if (this%save_resultsFile_next == istep) then
            this%save_resultsFile_next = this%save_resultsFile_next + this%save_resultsFile_step
            if (mpi_rank.eq.0) write(111,*) ' - Saving results file step:',istep,'(next to save',this%save_resultsFile_next,')'
            call nvtxStartRange("Output")
            call compute_fieldDerivs(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)
            call this%saveInstResultsFiles(istep)
            call nvtxEndRange
         end if

         ! InSitu - pass data
         call run_InSitu(u(:,:,2),istep)

!         call this%afterDt(istep)

         if(this%save_logFile_next==istep) then
            this%save_logFile_next = this%save_logFile_next + this%save_logFile_step
            if(mpi_rank.eq.0) call flush(111)
         end if

         !!! Witness points interpolation !!!
         if(this%have_witness) then
            if (this%continue_oldLogs) then
               if (mod(istep-this%load_step,this%leapwit)==0) then
                  iwitstep = iwitstep+1
                  call this%update_witness(istep, iwitstep)
               end if
               if ((istep-this%load_step > 0) .and. (mod((istep-this%load_step),this%leapwitsave*this%leapwit)==0)) then
                  if (this%wit_save) call this%save_witness(istep)
                  iwitstep = 0
               end if
            else
               if (mod(istep,this%leapwit)==0) then
                  iwitstep = iwitstep+1
                  call this%update_witness(istep, iwitstep)
               end if
               if ((istep > 0) .and. (mod((istep),this%leapwitsave*this%leapwit)==0)) then
                  if (this%wit_save) call this%save_witness(istep)
                  iwitstep = 0
               end if
            end if
         end if


         ! End simulation when physical time is reached (user defined)
         if (this%time .ge. this%maxPhysTime) then
            write(111,*) "--| Time integration finished at step: ",istep,"| time: ",this%time
            ! TODO: check if we want to save the last step
            exit
         end if
         this%local_step = this%local_step + 1
      end do
      call nvtxEndRange

      call this%endNSSolver()

   end subroutine CFDSolverBase_evalTimeIteration

   subroutine CFDSolverBase_update_witness(this, istep, iwitstep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in)              :: istep, iwitstep
      integer(4)                          :: iwit, iwitglobal, itewit, inode
      real(rp)                            :: start, finish, auxux, auxuy, auxuz, auxpr, auxrho

      !$acc parallel loop gang
      do iwit = 1,this%nwitPar
         auxux  = 0.0_rp
         auxuy  = 0.0_rp
         auxuz  = 0.0_rp
         auxpr  = 0.0_rp
         auxrho = 0.0_rp
         !$acc loop vector reduction(+:auxux,auxuy,auxuz,auxpr,auxrho)
         do inode = 1,nnode
            auxux  = auxux + Nwit(iwit,inode)*u(connecParOrig(witel(iwit),inode),1,2)
            auxuy  = auxuy + Nwit(iwit,inode)*u(connecParOrig(witel(iwit),inode),2,2)
            auxuz  = auxuz + Nwit(iwit,inode)*u(connecParOrig(witel(iwit),inode),3,2)
            auxpr  = auxpr + Nwit(iwit,inode)*pr(connecParOrig(witel(iwit),inode),2)
            auxrho = auxrho + Nwit(iwit,inode)*rho(connecParOrig(witel(iwit),inode),2)
         end do
         buffwit(iwit,iwitstep,1) = auxux
         buffwit(iwit,iwitstep,2) = auxuz
         buffwit(iwit,iwitstep,3) = auxuy
         buffwit(iwit,iwitstep,4) = auxpr
         buffwit(iwit,iwitstep,5) = auxrho
      end do
      !$acc end loop
      bufftime(iwitstep) = this%time
      buffstep(iwitstep) = istep
   end subroutine CFDSolverBase_update_witness

   subroutine CFDSolverBase_save_witness(this, istep)
      class(CFDSolverBase), intent(inout) :: this
      integer(4), intent(in)              :: istep
      integer(4)                          :: iwit, iwitglobal, itewit
      real(rp)                            :: start, finish
      character(512)                      :: fname

      if ((this%continue_witness .eqv. .false.) .AND. (this%continue_oldLogs .eqv. .false.)) then
         itewit = istep/(this%leapwit)
      end if
      if ((this%continue_witness .eqv. .false.) .AND. (this%continue_oldLogs .eqv. .true.)) then
         itewit = (istep - this%load_step)/(this%leapwit)
      end if
      if ((this%continue_witness .eqv. .true.) .AND. (this%continue_oldLogs .eqv. .true.)) then
         itewit = this%load_stepwit + (istep - this%load_step)/(this%leapwit)
      end if
      fname = trim(adjustl(this%witness_h5_file_name))
      call update_witness_hdf5(itewit, this%leapwitsave, buffwit, this%nwit, this%nwitPar, this%nvarwit, fname, bufftime, buffstep, this%wit_save_u_i, this%wit_save_pr, this%wit_save_rho)
   end subroutine CFDSolverBase_save_witness

   subroutine CFDSolverBase_preprocWitnessPoints(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer(4)                          :: iwit, jwit, ielem, inode, ifound, nwitParCand, icand, nwitFound, nwit2find, icount=0, imiss=0, myrank, mpi_rank_found
      integer(4)                          :: witGlobCand(this%nwit), witGlobFound(this%nwit*mpi_size)
      integer(4), allocatable             :: witGlobFound2(:), witGlobMiss(:)
      real(rp)                            :: xi(ndime), radwit(numElemsRankPar), maxL, center(numElemsRankPar,ndime), aux1, aux2, aux3, auxvol, helemmax(numElemsRankPar), Niwit(nnode), dist(numElemsRankPar), xyzwit(ndime), mindist
      real(rp), parameter                 :: wittol=1e-7
      real(rp)                            :: witxyz(this%nwit,ndime), witxyzPar(this%nwit,ndime), witxyzParCand(this%nwit,ndime)
      real(rp)                            :: xmin, ymin, zmin, xmax, ymax, zmax
      real(rp)                            :: xminloc, yminloc, zminloc, xmaxloc, ymaxloc, zmaxloc
      logical                             :: isinside, found
      character(512)                      :: input, output

      type two_real
	      real(4) :: real_1
	      real(4) :: real_2
      end type
      
      type(two_real)                      :: locdist, globdist

      if(mpi_rank.eq.0) then
         write(*,*) "--| Preprocessing witness points"
      end if

      allocate(witGlob(this%nwit))

      !$acc parallel loop gang
      do ielem = 1, numElemsRankPar
         aux1   = 0.0_rp
         aux2   = 0.0_rp
         aux3   = 0.0_rp
         auxvol = 0.0_rp
         !$acc loop vector reduction(+:aux1, aux2, aux3, auxvol)
         do inode = 1, nnode
            aux1   = aux1 + coordPar(connecParOrig(ielem,inode),1)
            aux2   = aux2 + coordPar(connecParOrig(ielem,inode),2)
            aux3   = aux3 + coordPar(connecParOrig(ielem,inode),3)
            auxvol = auxvol+gpvol(1,inode,ielem) !nnode = ngaus
         end do
         center(ielem,1) = aux1/nnode
         center(ielem,2) = aux2/nnode
         center(ielem,3) = aux3/nnode
         helemmax(ielem) = auxvol**(1.0/3.0)
      end do
      !$acc end loop
      maxL = maxval(abs(helemmax))

      xminloc = minval(coordPar(:,1)) - wittol
      yminloc = minval(coordPar(:,2)) - wittol
      zminloc = minval(coordPar(:,3)) - wittol
      xmaxloc = maxval(coordPar(:,1)) + wittol
      ymaxloc = maxval(coordPar(:,2)) + wittol
      zmaxloc = maxval(coordPar(:,3)) + wittol

      call MPI_Allreduce(xminloc, xmin, 1, mpi_datatype_real, MPI_MIN, app_comm, mpi_err)
      call MPI_Allreduce(yminloc, ymin, 1, mpi_datatype_real, MPI_MIN, app_comm, mpi_err)
      call MPI_Allreduce(zminloc, zmin, 1, mpi_datatype_real, MPI_MIN, app_comm, mpi_err)
      call MPI_Allreduce(xmaxloc, xmax, 1, mpi_datatype_real, MPI_MAX, app_comm, mpi_err)
      call MPI_Allreduce(ymaxloc, ymax, 1, mpi_datatype_real, MPI_MAX, app_comm, mpi_err)
      call MPI_Allreduce(zmaxloc, zmax, 1, mpi_datatype_real, MPI_MAX, app_comm, mpi_err)

      !$acc kernels
      witGlobCand(:) = 0
      witGlob(:) = 0
      witxyzPar(:,:) = 0.0_rp
      !$acc end kernels
      ifound  = 0
      icand   = 0
      input = trim(adjustl(this%witness_inp_file_name))
      call read_points(input, this%nwit, witxyz)
      do iwit = 1, this%nwit
	 if (witxyz(iwit,1) < xmin .OR. witxyz(iwit,2) < ymin .OR. witxyz(iwit,3) < zmin .OR. witxyz(iwit,1) > xmax .OR. witxyz(iwit,2) > ymax .OR. witxyz(iwit,3) > zmax) then
		write(*,*) "FATAL ERROR!! Witness point out of bounds", witxyz(iwit,:)
         	call MPI_Abort(app_comm,-1,mpi_err)
	 end if
	 if (witxyz(iwit,1) > xminloc .AND. witxyz(iwit,2) > yminloc .AND. witxyz(iwit,3) > zminloc .AND. witxyz(iwit,1) < xmaxloc .AND. witxyz(iwit,2) < ymaxloc .AND. witxyz(iwit,3) < zmaxloc) then
            icand = icand + 1
            witGlobCand(icand) = iwit
            witxyzParCand(icand,:) = witxyz(iwit,:)
         end if
      end do
      nwitParCand = icand

      do iwit = 1, nwitParCand
	 !$acc kernels
         radwit(:) = ((witxyzParCand(iwit, 1)-center(:,1))*(witxyzParCand(iwit, 1)-center(:,1))+(witxyzParCand(iwit, 2)-center(:,2))*(witxyzParCand(iwit, 2)-center(:,2))+(witxyzParCand(iwit, 3)-center(:,3))*(witxyzParCand(iwit, 3)-center(:,3)))-maxL*maxL
         !$acc end kernels
         do ielem = 1, numElemsRankPar
            if (radwit(ielem) < 0) then
               call isocoords(coordPar(connecParOrig(ielem,:),:), witxyzParCand(iwit,:), atoIJK, xi, isinside, Niwit)
               if (isinside .AND. (abs(xi(1)) < 1.0_rp+wittol) .AND. (abs(xi(2)) < 1.0_rp+wittol) .AND. (abs(xi(3)) < 1.0_rp+wittol)) then
                  ifound = ifound+1
                  witel(ifound)   = ielem
                  witxi(ifound,:) = xi(:)
                  witxyzPar(ifound,:)  = witxyzParCand(iwit, :)
                  witGlob(ifound) = witGlobCand(iwit)
                  Nwit(ifound,:) = Niwit(:)
                  exit
               end if
            end if
         end do
      end do
      this%nwitPar = ifound
      !Check that all witness points have been found
      call MPI_Allreduce(this%nwitPar, nwitFound, 1, MPI_INTEGER, MPI_SUM, app_comm,mpi_err)
      if (nwitFound < this%nwit) then
         nwit2find = this%nwit - nwitFound
         if (mpi_rank .eq. 0) then
            write(*,*) "WARNING!!!! The following witness points were not found inside any element, taking the element with the closest centroid as the one they belong to. Make sure they are inside the domain"
         endif
         call MPI_Allgather(witGlob, this%nwit, MPI_INTEGER, witGlobFound, this%nwit, MPI_INTEGER, app_comm,mpi_err)
         allocate(witGlobFound2(nwitFound))
         allocate(witGlobMiss(nwit2Find))
         do iwit = 1, this%nwit*mpi_size
            if (witGlobFound(iwit) > 0) then
               icount = icount + 1
		      witGlobFound2(icount) = witGlobFound(iwit)
	         end if
	      end do
	      do iwit = 1, this%nwit
	 	      found = .false.
		      do jwit = 1, nwitFound
	 		      if (witGlobFound2(jwit) == iwit) then
			      	found = .true.
			      	exit
			      end if
		      end do
		      if (found .eqv. .false.) then
		      	imiss = imiss + 1
		      	witGlobMiss(imiss) = iwit
		      end if
	      end do
         do iwit = 1, nwit2find
            xyzwit(:)  = witxyz(witGlobMiss(iwit),:)
            !$acc kernels
            dist(:)    = (center(:,1)-xyzwit(1))*(center(:,1)-xyzwit(1))+(center(:,2)-xyzwit(2))*(center(:,2)-xyzwit(2))+(center(:,3)-xyzwit(3))*(center(:,3)-xyzwit(3))
            !$acc end kernels
            ielem      = minloc(dist(:),1)
            locdist % real_1 = dist(ielem)
            locdist % real_2 = mpi_rank
            call MPI_Allreduce(locdist, globdist, 1, MPI_2REAL, MPI_MINLOC, app_comm, mpi_err)
            mpi_rank_found = int(globdist % real_2, kind=4)
            if (mpi_rank .eq. mpi_rank_found) then
	            write(*,*) "[NOT FOUND WITNESS] ", xyzwit(:)
               this%nwitPar              = this%nwitPar+1
               witGlob(this%nwitPar)     = witGlobMiss(iwit)
               witel(this%nwitPar)       = ielem
	            witxyzPar(this%nwitPar,:) = witxyz(witGlobMiss(iwit),:)
	            call isocoords(coordPar(connecParOrig(ielem,:),:), witxyzPar(this%nwitPar,:), atoIJK, witxi(this%nwitPar,:), isinside, Nwit(this%nwitPar,:)) 
            end if
         end do
         deallocate(witGlobFound2)
         deallocate(witGlobMiss)
      end if

      allocate(buffwit(this%nwitPar,this%leapwitsave,this%nvarwit))
      allocate(bufftime(this%leapwitsave))
      allocate(buffstep(this%leapwitsave))
      output = trim(adjustl(this%witness_h5_file_name))
      if (this%wit_save) call create_witness_hdf5(output, nnode, witxyzPar, witel, witxi, Nwit, this%nwit, this%nwitPar, witGlob, this%wit_save_u_i, this%wit_save_pr, this%wit_save_rho)
      if(mpi_rank.eq.0) then
         write(*,*) "--| End of preprocessing witness points"
      end if
   end subroutine CFDSolverBase_preprocWitnessPoints

   subroutine CFDSolverBase_loadWitnessPoints(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      character(512)                      :: output
      
      output = trim(adjustl(this%witness_h5_file_name))
      call load_witness_hdf5(output, nnode, this%nwit, this%load_step, this%load_stepwit, this%nwitPar, witel, witxi, Nwit)
      allocate(buffwit(this%nwitPar,this%leapwitsave,this%nvarwit))
      allocate(bufftime(this%leapwitsave))
      allocate(buffstep(this%leapwitsave))
   end subroutine CFDSolverBase_loadWitnessPoints

   subroutine open_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      character(len=1024) :: aux_string_mpisize
      integer :: iCode

      if(mpi_rank.eq.0) then
         write(aux_string_mpisize,'(I0)') mpi_size

         if(len_trim(adjustl(this%io_append_info)) == 0) then
            this%io_append_info = trim(adjustl(this%io_append_info))
        else
            this%io_append_info = "."//trim(adjustl(this%io_append_info))
         end if

         this%log_file_name = trim(adjustl(this%io_prepend_path))//'sod2d_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//trim(this%io_append_info)//'.log'
         if(this%continue_oldLogs) then
            open(unit=111,file=this%log_file_name,status='old',position='append')
         else
            open(unit=111,file=this%log_file_name,status='replace')
         end if
      end if

      if(mpi_rank.eq.0) then
         write(111,*) "--| Current case settings:"
         write(111,*) "--------------------------------------------"
         write(111,*) "  # Constants:"
         write(111,*) "    rp: ",               rp
         write(111,*) "    rp_vtk: ",           rp_vtk
         write(111,*) "    rp_avg: ",           rp_avg
         write(111,*) "    porder: ",           porder
         write(111,*) "    flag_real_diff: ",   flag_real_diff
         write(111,*) "    flag_diff_suth: ",   flag_diff_suth
         write(111,*) "--------------------------------------------"
         write(111,*) "  # Numerical parameters:"
         write(111,*) "    flag_implicit: ",            flag_implicit
         write(111,*) "    flag_rk_order: ",            flag_rk_order
         write(111,*) "    flag_les: ",                 flag_les
         write(111,*) "    flag_les_ilsa: ",            flag_les_ilsa
         write(111,*) "    flag_solver_type: ",         flag_solver_type
         write(111,*) "    flag_spectralElem: ",        flag_spectralElem
         write(111,*) "    flag_normalise_entropy: ",   flag_normalise_entropy
         write(111,*) "    flag_walave: ",              flag_walave
         write(111,*) "    flag_buffer_on: ",           flag_buffer_on
         write(111,*) "--------------------------------------------"
         write(111,*) "    ce: ",      ce
         write(111,*) "    cmax: ",    cmax
         write(111,*) "    c_sgs: ",   c_sgs
         write(111,*) "    cglob: ",   cglob
         write(111,*) "    c_rho: ",   c_rho
         write(111,*) "    c_ener: ",  c_ener
         write(111,*) "    stau: ",    stau
         write(111,*) "    T_ilsa: ",  T_ilsa
         write(111,*) "--------------------------------------------"
         write(111,*) "    cfl_conv: ",         this%cfl_conv
         write(111,*) "    cfl_diff: ",         this%cfl_diff
         write(111,*) "    maxIter: ",          maxIter
         write(111,*) "    tol: ",              tol
         write(111,*) "--------------------------------------------"
      end if

   end subroutine open_log_file

   subroutine open_analysis_files(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      character(len=1024) :: filenameAnalysis,fileNameTimer,filenameBound,aux_string_mpisize,aux_string_code
      integer :: iCode

      if(mpi_rank.eq.0) then
         write(aux_string_mpisize,'(I0)') mpi_size

         if(this%doGlobalAnalysis) then
            filenameAnalysis = trim(adjustl(this%io_prepend_path))//'analysis_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//trim(this%io_append_info)//'.dat'
            if(this%continue_oldLogs) then
               open(unit=666,file=filenameAnalysis,status='old',position='append')
            else
               open(unit=666,file=filenameAnalysis,status='replace')
            end if
         end if

         if(this%doTimerAnalysis) then
            fileNameTimer = trim(adjustl(this%io_prepend_path))//'timer_'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//trim(this%io_append_info)//'.dat'
            open(unit=123,file=fileNameTimer,status='replace')
            write(123,*) "iter iteTime iteTimeAvg"
         end if

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               write(aux_string_code,'(I0)') iCode
               filenameBound = trim(adjustl(this%io_prepend_path))//'surf_code_'//trim(aux_string_code)//'-'//trim(adjustl(this%mesh_h5_file_name))//'-'//trim(aux_string_mpisize)//trim(this%io_append_info)//'.dat'
               if(this%continue_oldLogs) then
                  open(unit=888+iCode,form='formatted',file=filenameBound,status='old',position='append')
               else
                  open(unit=888+iCode,form='formatted',file=filenameBound,status='replace')
                  write(888+iCode,60) "ITER", "TIME", "AREA", "FPRE_X", "FPRE_Y", "FPRE_Z", "FTAU_X", "FTAU_Y", "FTAU_Z"
                  60 format(9(3X,A,5X))
               end if
            end do
         end if
      end if

   end subroutine open_analysis_files

   subroutine close_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         close(unit=111)
      end if
   end subroutine close_log_file

   subroutine close_analysis_files(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         if(this%doGlobalAnalysis) close(unit=666)
         if(this%doTimerAnalysis)  close(unit=123)

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               close(unit=888+iCode)
            end do
         end if
      end if
   end subroutine close_analysis_files

   subroutine flush_log_file(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iCode

      if(mpi_rank.eq.0) then
         flush(unit=111)

         if(this%doGlobalAnalysis) then
            flush(unit=666)
         end if

         if (isMeshBoundaries) then
            do iCode = 1,numBoundCodes
               flush(unit=888+iCode)
            end do
         end if
      end if

   end subroutine flush_log_file

   subroutine CFDSolverBase_eval_vars_after_load_hdf5_resultsFile(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this
      integer :: iNodeL,idime
      real(rp) :: umag

      !values loaded -> rho,u,pr,E,mu_e,mu_sgs

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         e_int(iNodeL,2) = E(iNodeL,2)/rho(iNodeL,2) - 0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))
         !e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) =  eta(iNodeL,2)
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         umag = 0.0_rp
         !$acc loop seq
         do idime = 1,ndime
            umag = umag + u(iNodeL,idime,2)**2
         end do
         umag = sqrt(umag)
         machno(iNodeL) = umag/csound(iNodeL)
         !machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)
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
      !$acc end kernels

   end subroutine CFDSolverBase_eval_vars_after_load_hdf5_resultsFile

   subroutine CFDSolverBase_run(this)
      implicit none
      class(CFDSolverBase), intent(inout) :: this

      ! Init MPI
      call init_mpi()

      ! Init HDF5 interface
      call init_hdf5_interface()

      ! Init Save Fields vars and arrays
      call init_saveFields()

      ! Main simulation parameters
      call this%initializeDefaultParameters()
      call this%initializeParameters()

      call this%optimizeParameters()

      call read_json_saveFields(json_filename)

      ! Open log file
      call this%open_log_file()

      ! read the mesh
      call this%openMesh()

      ! Init hdf5 auxiliar saving arrays
      call init_hdf5_auxiliar_saving_arrays()

      ! Open analysis files
      call this%open_analysis_files

      ! Eval shape Functions
      call this%evalShapeFunctions()

      ! Allocate variables
      call this%allocateVariables()

      ! Allocating and defining filters
      call init_filters()

      ! Setting fields to be saved
      if(flag_use_species .eqv. .true.) then
         call setFields2Save(rho(:,2),mu_fluid,pr(:,2),E(:,2),eta(:,2),csound,machno,divU,qcrit,Tem(:,2),&
                             u(:,:,2),gradRho,curlU,mu_sgs,mu_e,&
                             avrho,avpre,avpre2,avmueff,avvel,avve2,avvex,avtw,Yk(:,:,2),mu_e_Yk)
      else
         call setFields2Save(rho(:,2),mu_fluid,pr(:,2),E(:,2),eta(:,2),csound,machno,divU,qcrit,Tem(:,2),&
                             u(:,:,2),gradRho,curlU,mu_sgs,mu_e,&
                             avrho,avpre,avpre2,avmueff,avvel,avve2,avvex,avtw)
      end if

      ! Eval or load initial conditions
      call this%evalOrLoadInitialConditions()

      ! Eval  viscosty factor
      call this%evalViscosityFactor()

      ! Eval initial viscosty
      call this%evalInitialViscosity()

      ! Compute characteristic size of the elements
      call this%evalCharLength()

      ! Eval boundary element normal
      call this%evalBoundaryNormals()

      ! Eval Jacobian information
      call this%evalJacobians()

      ! Eval AtoIJK inverse
      call this%evalAtoIJKInverse()

      ! Eval BoundaryFacesToNodes
      call  this%boundaryFacesToNodes()

      ! Eval list Elems per Node and Near Boundary Node
      call this%eval_elemPerNode_and_nearBoundaryNode()

      call this%set_mappedFaces_linkingNodes()   

      ! Init of the source terms
      call this%initializeSourceTerms()

      ! Eval mass
      call this%evalMass()

      ! Find Pressure Node to fix
      if(flag_fs_fix_pressure) call this%findFixPressure()

      ! Preprocess witness points
      if (this%have_witness) then
         if (this%continue_witness) then
            call this%loadWitnessPoints() ! Load witness points and continue them
         else
            call this%preprocWitnessPoints()
         end if
      end if

      ! Eval first output
      if(this%isFreshStart) call this%evalFirstOutput()
      call this%flush_log_file()

      call  this%normalFacesToNodes()

      ! Eval initial time step
      call this%evalInitialDt()

      call this%initialBuffer()

      call this%beforeTimeIteration()

      call this%flush_log_file()

      ! Do the time iteration
      call this%evalTimeIteration()

      ! After the time iteration
      call this%afterTimeIteration()

      call this%close_log_file()
      call this%close_analysis_files()

      ! Deallocate the variables
      call this%deallocateVariables()

      ! Deallocate the filters
      call deallocate_filters()

      ! End hdf5 auxiliar saving arrays
      call end_hdf5_auxiliar_saving_arrays()

      ! End hdf5 interface
      call end_hdf5_interface()

      ! Finalize InSitu   - bettre done before end_comms
      call end_InSitu()

      ! End comms
      call end_comms()
      call end_comms_bnd()

      ! End MPI
      call end_mpi()

   end subroutine CFDSolverBase_run

end module CFDSolverBase_mod
