#define ACTUATION 1
#define IMPLICIT 0
#define NEUMANN 1
#define SB 1


module BLTSBDRLFlowSolver_mod
   use mod_arrays
   use mod_nvtx
#ifndef NOACC
   use cudafor
#endif
   use mod_veclen

   use elem_qua
   use elem_hex
   use jacobian_oper
   use quadrature_rules
   use mesh_reader
   use inicond_reader
   use mass_matrix
   use mod_geom
   use mod_output
   use mod_period
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
   use mod_constants
   use mod_smartredis
   implicit none
   private

   real(rp), allocatable, dimension(:) :: eta_b,f,f_prim ! auxiliary data arrays needs to be here because how cuda works
   real(rp), allocatable, dimension(:,:) :: rectangleControl ! two point coordinates that define each rectangle
   integer(4), allocatable, dimension(:) :: actionMask ! mask that contains whether a point is in a rectangle control or not

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BLTSBDRLFlowSolver
      real(rp) , public  ::  M, d0, U0, rho0, Red0, Re, to, po, mu, amp_tbs, x_start, x_rise, x_end, x_fall, x_rerise, x_restart, coeff_tbs
      character(len=8), public :: tag="0"
      logical :: db_clustered = .false.
      character(512), public :: fileControlName ! file that contains the points defining the rectangle controls
      integer(4), public :: nRectangleControl
      real(rp), public :: previousActuationTime, periodActuation, frequencyActuation, timeBeginActuation ! parameters of the actuation

   contains
      procedure, public :: fillBCTypes => BLTSBDRLFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BLTSBDRLFlowSolver_initializeParameters
      procedure, public :: evalInitialConditions => BLTSBDRLFlowSolver_evalInitialConditions
      procedure, public :: initialBuffer => BLTSBDRLFlowSolver_initialBuffer
      procedure, public :: fillBlasius => BLTSBDRLFlowSolver_fillBlasius
      procedure, public :: smoothStep => BLTSBDRLFlowSolver_smoothStep
      procedure, public :: afterDt => BLTSBDRLFlowSolver_afterDt
      procedure, public :: getControlNodes => BLTSBDRLFlowSolver_getControlNodes
      procedure, public :: readControlRectangles => BLTSBDRLFlowSolver_readControlRectangles
#ifdef SMARTREDIS
      procedure, public :: initSmartRedis  => BLTSBDRLFlowSolver_initSmartRedis
      procedure, public :: afterTimeIteration => BLTSBDRLFlowSolver_afterTimeIteration
      procedure, public :: smoothControlFunction => BLTSBDRLFlowSolver_smoothControlFunction
#endif
   end type BLTSBDRLFlowSolver

contains

#ifdef SMARTREDIS
   subroutine BLTSBDRLFlowSolver_initSmartRedis(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this

      this%previousActuationTime = this%time
      open(unit=443,file="control_fortran_raw.txt",status='replace')
      open(unit=444,file="control_fortran_smooth.txt",status='replace')
      call init_smartredis(client, this%nwitPar, this%nRectangleControl, this%db_clustered)
      call write_step_type(client, 1, "step_type")
   end subroutine BLTSBDRLFlowSolver_initSmartRedis

   subroutine BLTSBDRLFlowSolver_afterTimeIteration(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this

      call write_step_type(client, 0, "step_type")
      close(443)
      close(444)
   end subroutine BLTSBDRLFlowSolver_afterTimeIteration
#endif

   subroutine BLTSBDRLFlowSolver_fill_BC_Types(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
#if (ACTUATION)
      bouCodes2BCType(1) = bc_type_unsteady_inlet ! wall + actuation
#else
      bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! wall
#endif
#if (NEUMANN)
      bouCodes2BCType(2) = bc_type_far_field_SB
#else
      bouCodes2BCType(2) = bc_type_far_field         ! Upper part of the domain
#endif
      bouCodes2BCType(3) = bc_type_far_field          ! inlet part of the domain
      bouCodes2BCType(4) = bc_type_far_field          ! outlet part of the domain
      !$acc update device(bouCodes2BCType(:))

   end subroutine BLTSBDRLFlowSolver_fill_BC_Types

   subroutine BLTSBDRLFlowSolver_initializeParameters(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      real(rp) :: mur
      integer :: num_args, equal_pos, iarg, ierr
      character(len=64) :: arg, output_dir
      character(len=8) :: restart_step_str="", db_clustered_str="0", frequencyActuation_str="1.0"
      logical :: output_dir_exists

      ! get command line args, ie: mpirun -n 4 sod2d --tag=12 --restart_step=2500
      num_args = command_argument_count()
      do iarg = 1, num_args
         call get_command_argument(iarg, arg)
         equal_pos = scan(adjustl(trim(arg)), "=")
         if (adjustl(trim(arg(:equal_pos-1))) .eq. "--tag") then
            this%tag = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--restart_step") then
            restart_step_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--db_clustered") then
            db_clustered_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--f_action") then
            frequencyActuation_str = trim(adjustl(arg(equal_pos+1:)))
         else
            stop "Unknown command line argument"
         end if
      end do

      if (this%tag == "") this%tag = "0"
      if (db_clustered_str == "" .or. db_clustered_str == "0") this%db_clustered = .false.
      read(frequencyActuation_str,*,iostat=ierr) this%frequencyActuation
      this%periodActuation = 1.0_rp / this%frequencyActuation

      ! create output dir if not existing
      output_dir = "./output_"//trim(adjustl(this%tag))//"/"
      inquire(file=trim(adjustl(output_dir)), exist=output_dir_exists)
      if (.not. output_dir_exists .and. mpi_rank .eq. 0) call execute_command_line("mkdir -p "//trim(adjustl(output_dir)))

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "bl_les"
      write(this%results_h5_file_path,*) "./output_"//trim(adjustl(this%tag))//"/"
      write(this%results_h5_file_name,*) "results_"//trim(adjustl(this%tag))
      write(this%io_prepend_path,*) output_dir

      if (mpi_rank .eq. 0) then
         write(*,*) "Received arguments are:"
         write(*,*) "--tag:  ", this%tag
         write(*,*) "--restart_step:  ", restart_step_str
         write(*,*) "--db_clustered:  ", db_clustered_str
         write(*,*) "--freq_action:  ", this%frequencyActuation
      end if

      !----------------------------------------------
      ! I/O params
      this%final_istep = 1000 ! 10000001

      this%save_logFile_first = 1
      this%save_logFile_step  = 50
#if (IMPLICIT)
      this%save_restartFile_first = 1
      this%save_restartFile_step = 200
      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 200
#else
      this%save_restartFile_first = 1
      this%save_restartFile_step = 1000
      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 1000
#endif

      if (restart_step_str == "" .or. restart_step_str == "0") then
         this%loadRestartFile = .false.
         this%restartFile_to_load = 1 ! not used
         this%continue_oldLogs = .false.
      else
         this%loadRestartFile = .true.
         read(restart_step_str, *) this%restartFile_to_load ! 1 (start) or 2 (last time step)
         this%continue_oldLogs = .false.
      end if

      this%initial_avgTime = 0.0 ! 3000.0_rp
      this%saveAvgFile = .true.
      this%loadAvgFile = .false.

      !----------------------------------------------
      ! Witness points parameters
      this%have_witness          = .true.
      this%witness_inp_file_name = "witness.txt"
      this%witness_h5_file_name  = "resultwit.h5"
      this%leapwit               = 1
      this%leapwitsave           = 1
      this%wit_save              = .false.
      this%wit_save_u_i          = .false.
      this%wit_save_pr           = .false.
      this%wit_save_rho          = .false.
      this%continue_witness      = .false.
      ! this%load_step             = 1

      !----------------------------------------------
      ! Control parameters
#if (ACTUATION)
      write(this%fileControlName ,*) "rectangleControl.txt"
      this%timeBeginActuation = 0.0 ! 2000.0_rp
#endif

      !----------------------------------------------
      ! numerical params
      flag_les = 1
#if (IMPLICIT)
      flag_implicit = 1
      maxIterNonLineal=200
      tol=1e-3
      pseudo_cfl =1.95_rp
      flag_implicit_repeat_dt_if_not_converged = 0
      flag_use_constant_dt = 1
      this%dt = 0.5_rp
#else
      flag_implicit = 0
      this%cfl_conv = 1.0_rp !M0.1 1.5, M0.3 0.95
      this%cfl_diff = 1.0_rp
#endif

      this%Cp   = 1004.0_rp
      this%Prt  = 0.71_rp
      this%M    = 0.1_rp
      this%d0   = 1.0_rp
      this%U0   = 1.0_rp
      this%rho0 = 1.0_rp

      ! Coefficients for the wall-normal boundary condition on the top
      !this%amp_tbs   = 1.0
      this%amp_tbs   = 0.5
      this%x_start   = 200.0_rp   ! x coordinate where 1st step function on the free streamwise velocity starts
      this%x_rise    = 20.0_rp   ! x length of the smoothing of the 1st step function for the free streamwise velocity
      this%x_end     = 435.0_rp  ! x coordinate where 2nd step function on the free streamwise velocity ends
      this%x_fall    = 115.0_rp   ! x half length of the smoothing of the 2nd step function for the free streamwise velocity
      this%x_rerise  = 100.0_rp   ! x length of the smoothing of the 3rd step function for the free streamwise velocity
      this%x_restart = 3.0_rp*(this%x_end-this%x_fall) - 2.0_rp*this%x_start - this%x_rise - 0.5_rp*this%x_rerise  ! x coordinate
      !where 3rd step function on the free streamwise velocity starts. The location is calculated to have a zero vertical mass flow
      !rate
      this%coeff_tbs = 0.5_rp

      this%Red0  = 450.0_rp
      this%gamma_gas = 1.40_rp

      this%mu    = this%rho0*this%d0*this%U0/this%Red0

      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%U0*this%U0/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = this%mu/mur

      nscbc_u_inf = this%U0
      nscbc_p_inf = this%po
      nscbc_rho_inf = this%rho0
      nscbc_gamma_inf = this%gamma_gas
      nscbc_c_inf = sqrt(this%gamma_gas*this%po/this%rho0)
      nscbc_Rgas_inf = this%Rgas


      flag_buffer_on = .true.
      ! x outlet
      flag_buffer_on_east = .true.
      flag_buffer_e_min = 950.0_rp
      flag_buffer_e_size = 50.0_rp

      ! x inlet
      flag_buffer_on_west = .true.
      flag_buffer_w_min = -50.0_rp
      flag_buffer_w_size = 50.0_rp

#if (NEUMANN)
      flag_include_neumann_flux = 1
#else
      ! y top
      flag_buffer_on_north = .true.
      flag_buffer_n_min =  100.0_rp
      flag_buffer_n_size = 20.0_rp
#endif

      ! -------- Instantaneous results file -------------
      this%save_scalarField_rho        = .true.
      this%save_scalarField_muFluid    = .false.
      this%save_scalarField_pressure   = .true.
      this%save_scalarField_energy     = .false.
      this%save_scalarField_entropy    = .false.
      this%save_scalarField_csound     = .false.
      this%save_scalarField_machno     = .false.
      this%save_scalarField_divU       = .false.
      this%save_scalarField_qcrit      = .true.
      this%save_scalarField_muSgs      = .false.
      this%save_scalarField_muEnvit    = .false.
      this%save_vectorField_vel        = .true.
      this%save_vectorField_gradRho    = .false.
      this%save_vectorField_curlU      = .true.

      ! -------- Average results file -------------
      this%save_avgScalarField_rho     = .true.
      this%save_avgScalarField_pr      = .true.
      this%save_avgScalarField_mueff   = .false.
      this%save_avgVectorField_vel     = .true.
      this%save_avgVectorField_ve2     = .false.
      this%save_avgVectorField_vex     = .true.
      this%save_avgVectorField_vtw     = .true.

      !Blasius analytical function
      call this%fillBlasius()
   end subroutine BLTSBDRLFlowSolver_initializeParameters

   subroutine BLTSBDRLFlowSolver_getControlNodes(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this

      integer(4) :: iNodeL, iRectangleControl, bcode
      real(rp)   :: xPoint, zPoint, x1RectangleControl, x2RectangleControl, z1RectangleControl, z2RectangleControl

      call this%readControlRectangles()
      allocate(actionMask(numNodesRankPar))
      !$acc enter data create(actionMask(:))

      !this%countPar = 0
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         actionMask(iNodeL) = 0
         if (bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL)
            if (bcode == bc_type_unsteady_inlet) then ! we are on the wall and we need to check if node is in the control rectangle
               do iRectangleControl = 1,this%nRectangleControl
                  xPoint = coordPar(iNodeL,1)
                  zPoint = coordPar(iNodeL,3)
                  x1RectangleControl = rectangleControl(1,2*iRectangleControl-1)
                  z1RectangleControl = rectangleControl(2,2*iRectangleControl-1)
                  x2RectangleControl = rectangleControl(1,2*iRectangleControl  )
                  z2RectangleControl = rectangleControl(2,2*iRectangleControl  )
                  if (xPoint >= x1RectangleControl .and. xPoint <= x2RectangleControl .and. zPoint >= z1RectangleControl .and. zPoint <= z2RectangleControl) then
                     actionMask(iNodeL) = iRectangleControl
                     !this%countPar = this%countPar + 1
                     exit
                  endif
               end do
            end if
         end if
      end do
      !$acc end parallel loop
      !$acc update device(actionMask(:))
   end subroutine BLTSBDRLFlowSolver_getControlNodes

   subroutine BLTSBDRLFlowSolver_readControlRectangles(this)
      ! This subroutine reads the file that contains the two points defining a rectangle parallel to the X-Z axis. Several rectangles
      ! can be introduced. In this rectangles is where control will be applied
      class(BLTSBDRLFlowSolver), intent(inout) :: this

      integer(rp)                           :: ii

      open(unit=99, file=this%fileControlName, status='old', action='read')
      read(99,*) this%nRectangleControl
      allocate(rectangleControl(2,2*this%nRectangleControl))
      !$acc enter data create(rectangleControl(:,:))
      do ii = 1, this%nRectangleControl
         read(99, *) rectangleControl(:,2*ii-1)  ! First point
         read(99, *) rectangleControl(:,2*ii  )  ! Second point
         read(99, *)
      end do
      close(99)

      !$acc update device(rectangleControl(:,:))
   end subroutine BLTSBDRLFlowSolver_readControlRectangles

   subroutine BLTSBDRLFlowSolver_afterDt(this,istep)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      integer(4), intent(in) :: istep

      integer(4) :: iNodeL, bcode, rectangleId
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3
      integer(4) :: k, j, ielem, iCen
      real(rp) :: dy, fx1, fx2, x
#if ACTUATION
#ifdef SMARTREDIS
      real(rp) :: action_global_instant(action_global_size)
#endif
#endif
      cd = 1.0_rp
      lx = this%d0*3.5_rp
      ly = this%d0*3.5_rp
      xmin = -15.0_rp*this%d0

      xmax = xmin+lx

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(coordPar(iNodeL,2) < ly) then
            if((coordPar(iNodeL,1)<xmax)  .and. (coordPar(iNodeL,1)>xmin)) then
               source_term(iNodeL,1) = -0.5_rp*rho(iNodeL,2)*cd*u(iNodeL,1,2)*abs(u(iNodeL,1,2))/lx
               source_term(iNodeL,2) = 0.00_rp
               source_term(iNodeL,3) = 0.00_rp
            end if
         end if
      end do
      !$acc end parallel loop
#if ACTUATION
#ifdef SMARTREDIS
      if (this%time .gt. this%timeBeginActuation) then
         write(*,*) this%time, this%time - this%previousActuationTime, this%periodActuation
            ! check if a new action is needed
         if (this%time - this%previousActuationTime .gt. this%periodActuation) then
            ! save old action values and time - useful for interpolating to new action_global values
            !$acc kernels
            action_global_previous(:) = action_global(:)
            this%previousActuationTime = this%previousActuationTime + this%periodActuation
            !$acc end kernels
            write(*,*) "Sod2D to write time: ", this%time
            call write_time(client, this%time, "time") ! writes current time into database
            write(*,*) "Sod2D wrote time: ", this%time
            call read_action(client, "action") ! modifies action_global (the target control values)
            write(*,*) "Sod2D read action: ", action_global
            write(443,*) action_global(1), this%time+this%periodActuation
         end if

         call this%smoothControlFunction(action_global_instant)
         write(444,*) action_global_instant(1), this%time


         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if ((bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) .and. (actionMask(iNodeL) .gt. 0)) then
               bcode = bouCodesNodesPar(iNodeL)
               if (bcode == bc_type_unsteady_inlet) then
                  u_buffer(iNodeL,1) = 0.0_rp
                  u_buffer(iNodeL,2) = action_global(actionMask(iNodeL))
                  u_buffer(iNodeL,3) = 0.0_rp
               end if
            end if
         end do
         !$acc end parallel loop
      end if

      if (step_type_mod .eq. 1) then
         call write_step_type(client, 2, "step_type")
         write(*,*) "Sod2D wrote step: 2"
      end if
#endif
#endif
      if(flag_include_neumann_flux == 1) then
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if(coordPar(iNodeL,2)  .gt. 100.0_rp) then
#if (SB)
               u_buffer_flux(iNodeL,1) = this%mu*((0.158531_rp-0.00110576_rp*coordPar(iNodeL,1)+1.8030232059983043_rp*10.0_rp**(-6.0_rp)*coordPar(iNodeL,1)**2.0_rp)*exp(-0.00008192_rp*(306.641_rp- coordPar(iNodeL,1))**2.0_rp))
#else
               x=coordPar(iNodeL,1)
               fx1 = 0.9_rp*exp(-((x-171.9_rp)/(0.3375_rp*100.0_rp))**2)
               x=coordPar(iNodeL,1)+11.0_rp
               fx2 = 0.9_rp*exp(-((x-171.9_rp)/(0.3375_rp*100.0_rp))**2)

               u_buffer_flux(iNodeL,1) = this%mu*((fx2-fx1)/11.0_rp)
#endif
            end if
         end do
         !$acc end parallel loop
      end if
   end subroutine BLTSBDRLFlowSolver_afterDt

#if ACTUATION
#ifdef SMARTREDIS
   subroutine BLTSBDRLFlowSolver_smoothControlFunction(this, action_global_instant)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      real(rp), intent(inout) :: action_global_instant(action_global_size)

      real(rp) :: f1, f2, f3

      f1 = exp(-1.0_rp / ((this%time - this%previousActuationTime) / this%periodActuation))
      f2 = exp(-1.0_rp / (1.0_rp - (this%time - this%previousActuationTime) / this%periodActuation))
      f3 = f1 / (f1 + f2)

      !$acc kernels
      action_global_instant(:) = action_global_previous(:) + f3 * (action_global(:) - action_global_previous(:))
      !$acc end kernels
   end subroutine BLTSBDRLFlowSolver_smoothControlFunction
#endif
#endif

   subroutine BLTSBDRLFlowSolver_initialBuffer(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      integer(4) :: iNodeL,k,j,bcode,ielem,iCen
      real(rp) :: yp,eta_y,f_y,f_prim_y, f1, f2, f3,x,dy

#if (ACTUATION)
      call this%getControlNodes()
#endif

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2)
         eta_y = yp !with our normalisation is sqrt(U/(nu x ) is actually 1 for the inlet)
         j = 45
         !$acc loop seq
         label1:do k=1,45
            if(eta_y<eta_b(k)) then
               j = k
               exit label1
            end if
         end do label1
         if(j == 1) then
            u_buffer(iNodeL,1) = 0.0_rp
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp
         else if(j==45) then
            u_buffer(iNodeL,1) = this%U0
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp
         else
            f_y      = f(j-1)      + (f(j)-f(j-1))*(eta_y-eta_b(j-1))/(eta_b(j)-eta_b(j-1))
            f_prim_y = f_prim(j-1) + (f_prim(j)-f_prim(j-1))*(eta_y-eta_b(j-1))/(eta_b(j)-eta_b(j-1))

            u_buffer(iNodeL,1) = f_prim_y
            u_buffer(iNodeL,2) = 0.5_rp*sqrt(1.0/(450.0_rp*450.0_rp))*(eta_y*f_prim_y-f_y)
            u_buffer(iNodeL,3) = 0.0_rp
         end if
         if(yp .gt. 100.0_rp) then
#if (SB)
            u_buffer(iNodeL,2) =  0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                              -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
#else
            u_buffer(iNodeL,2) =  0.9_rp*exp(-((coordPar(iNodeL,1)-171.9_rp)/(0.3375_rp*100.0_rp))**2)
#endif
         end if
      end do
      !$acc end parallel loop

   end subroutine BLTSBDRLFlowSolver_initialBuffer

   subroutine BLTSBDRLFlowSolver_evalInitialConditions(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      integer(4) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime, j,k,bcode
      real(rp) :: yp,eta_y,f_y,f_prim_y, f1, f2, f3, x
      integer(4)   :: iLine,iNodeGSrl,auxCnt

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2)
         eta_y = yp !with our normalisation is sqrt(U/(nu x ) is actually 1 for the inlet)
         j = 45
         !$acc loop seq
         label1:do k=1,45
            if(eta_y<eta_b(k)) then
               j = k
               exit label1
            end if
         end do label1
         if(j == 1) then
            u(iNodeL,1,2) = 0.0_rp
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = 0.0_rp
         else if(j==45) then
            u(iNodeL,1,2) = this%U0
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = 0.0_rp
         else
            f_y      = f(j-1)      + (f(j)-f(j-1))*(eta_y-eta_b(j-1))/(eta_b(j)-eta_b(j-1))
            f_prim_y = f_prim(j-1) + (f_prim(j)-f_prim(j-1))*(eta_y-eta_b(j-1))/(eta_b(j)-eta_b(j-1))

            u(iNodeL,1,2) = f_prim_y
            u(iNodeL,2,2) = 0.5_rp*sqrt(1.0/(450.0_rp*450.0_rp))*(eta_y*f_prim_y-f_y)
            u(iNodeL,3,2) = 0.0_rp
         end if
#if(!NEUMANN)
         if(yp .gt. 100.0_rp) then
#if (SB)
            u(iNodeL,2,2) =  0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                              -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
#else
            u(iNodeL,2,2) =  0.9_rp*exp(-((coordPar(iNodeL,1)-171.9_rp)/(0.3375_rp*100.0_rp))**2)
#endif
         end if
#endif
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%rho0
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)

         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3) = rho(iNodeL,2)
         E(iNodeL,3) =  E(iNodeL,2)
         eta(iNodeL,3) = eta(iNodeL,2)
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
   end subroutine BLTSBDRLFlowSolver_evalInitialConditions

   subroutine BLTSBDRLFlowSolver_smoothStep(this,x,y)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      real(rp), intent(in)  :: x
      real(rp), intent(out) :: y

      if(x<=0.0_rp) then
         y = 0.0_rp
      elseif(x<1.0_rp) then
         y = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
      else
         y = 1
      end if

   end subroutine BLTSBDRLFlowSolver_smoothStep

   subroutine BLTSBDRLFlowSolver_fillBlasius(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this

      allocate(eta_b(45), f(45), f_prim(45))
      !$acc enter data create(eta_b(:))
      !$acc enter data create(f(:))
      !$acc enter data create(f_prim(:))

      eta_b(1) = real(0.0E+00,rp)
      eta_b(2) = real(2.0E-01,rp)
      eta_b(3) = real(4.0E-01,rp)
      eta_b(4) = real(6.0E-01,rp)
      eta_b(5) = real(8.0E-01,rp)
      eta_b(6) = real(1.0E+00,rp)
      eta_b(7) = real(1.2E+00,rp)
      eta_b(8) = real(1.4E+00,rp)
      eta_b(9) = real(1.6E+00,rp)
      eta_b(10) = real(1.8E+00,rp)
      eta_b(11) = real(2.0E+00,rp)
      eta_b(12) = real(2.2E+00,rp)
      eta_b(13) = real(2.4E+00,rp)
      eta_b(14) = real(2.6E+00,rp)
      eta_b(15) = real(2.8E+00,rp)
      eta_b(16) = real(3.0E+00,rp)
      eta_b(17) = real(3.2E+00,rp)
      eta_b(18) = real(3.4E+00,rp)
      eta_b(19) = real(3.6E+00,rp)
      eta_b(20) = real(3.8E+00,rp)
      eta_b(21) = real(4.0E+00,rp)
      eta_b(22) = real(4.2E+00,rp)
      eta_b(23) = real(4.4E+00,rp)
      eta_b(24) = real(4.6E+00,rp)
      eta_b(25) = real(4.8E+00,rp)
      eta_b(26) = real(5.0E+00,rp)
      eta_b(27) = real(5.2E+00,rp)
      eta_b(28) = real(5.4E+00,rp)
      eta_b(29) = real(5.6E+00,rp)
      eta_b(30) = real(5.8E+00,rp)
      eta_b(31) = real(6.0E+00,rp)
      eta_b(32) = real(6.2E+00,rp)
      eta_b(33) = real(6.4E+00,rp)
      eta_b(34) = real(6.6E+00,rp)
      eta_b(35) = real(6.8E+00,rp)
      eta_b(36) = real(7.0E+00,rp)
      eta_b(37) = real(7.2E+00,rp)
      eta_b(38) = real(7.4E+00,rp)
      eta_b(39) = real(7.6E+00,rp)
      eta_b(40) = real(7.8E+00,rp)
      eta_b(41) = real(8.0E+00,rp)
      eta_b(42) = real(8.2E+00,rp)
      eta_b(43) = real(8.4E+00,rp)
      eta_b(44) = real(8.6E+00,rp)
      eta_b(45) = real(8.8E+00,rp)

      f(1)  = real(0.000000000E+00,rp)
      f(2)  = real(6.640999715E-03,rp)
      f(3)  = real(2.655988402E-02,rp)
      f(4)  = real(5.973463750E-02,rp)
      f(5)  = real(1.061082208E-01,rp)
      f(6)  = real(1.655717258E-01,rp)
      f(7)  = real(2.379487173E-01,rp)
      f(8)  = real(3.229815738E-01,rp)
      f(9)  = real(4.203207655E-01,rp)
      f(10) = real(5.295180377E-01,rp)
      f(11) = real(6.500243699E-01,rp)
      f(12) = real(7.811933370E-01,rp)
      f(13) = real(9.222901256E-01,rp)
      f(14) = real(1.072505977E+00,rp)
      f(15) = real(1.230977302E+00,rp)
      f(16) = real(1.396808231E+00,rp)
      f(17) = real(1.569094960E+00,rp)
      f(18) = real(1.746950094E+00,rp)
      f(19) = real(1.929525170E+00,rp)
      f(20) = real(2.116029817E+00,rp)
      f(21) = real(2.305746418E+00,rp)
      f(22) = real(2.498039663E+00,rp)
      f(23) = real(2.692360938E+00,rp)
      f(24) = real(2.888247990E+00,rp)
      f(25) = real(3.085320655E+00,rp)
      f(26) = real(3.283273665E+00,rp)
      f(27) = real(3.481867612E+00,rp)
      f(28) = real(3.680919063E+00,rp)
      f(29) = real(3.880290678E+00,rp)
      f(30) = real(4.079881939E+00,rp)
      f(31) = real(4.279620923E+00,rp)
      f(32) = real(4.479457297E+00,rp)
      f(33) = real(4.679356615E+00,rp)
      f(34) = real(4.879295811E+00,rp)
      f(35) = real(5.079259772E+00,rp)
      f(36) = real(5.279238811E+00,rp)
      f(37) = real(5.479226847E+00,rp)
      f(38) = real(5.679220147E+00,rp)
      f(39) = real(5.879216466E+00,rp)
      f(40) = real(6.079214481E+00,rp)
      f(41) = real(6.279213431E+00,rp)
      f(42) = real(6.479212887E+00,rp)
      f(43) = real(6.679212609E+00,rp)
      f(44) = real(6.879212471E+00,rp)
      f(45) = real(7.079212403E+00,rp)

      f_prim(1)  = real(0.000000000E+00,rp)
      f_prim(2)  = real(6.640779210E-02,rp)
      f_prim(3)  = real(1.327641608E-01,rp)
      f_prim(4)  = real(1.989372524E-01,rp)
      f_prim(5)  = real(2.647091387E-01,rp)
      f_prim(6)  = real(3.297800312E-01,rp)
      f_prim(7)  = real(3.937761044E-01,rp)
      f_prim(8)  = real(4.562617647E-01,rp)
      f_prim(9)  = real(5.167567844E-01,rp)
      f_prim(10) = real(5.747581439E-01,rp)
      f_prim(11) = real(6.297657365E-01,rp)
      f_prim(12) = real(6.813103772E-01,rp)
      f_prim(13) = real(7.289819351E-01,rp)
      f_prim(14) = real(7.724550211E-01,rp)
      f_prim(15) = real(8.115096232E-01,rp)
      f_prim(16) = real(8.460444437E-01,rp)
      f_prim(17) = real(8.760814552E-01,rp)
      f_prim(18) = real(9.017612214E-01,rp)
      f_prim(19) = real(9.233296659E-01,rp)
      f_prim(20) = real(9.411179967E-01,rp)
      f_prim(21) = real(9.555182298E-01,rp)
      f_prim(22) = real(9.669570738E-01,rp)
      f_prim(23) = real(9.758708321E-01,rp)
      f_prim(24) = real(9.826835008E-01,rp)
      f_prim(25) = real(9.877895262E-01,rp)
      f_prim(26) = real(9.915419002E-01,rp)
      f_prim(27) = real(9.942455354E-01,rp)
      f_prim(28) = real(9.961553040E-01,rp)
      f_prim(29) = real(9.974777682E-01,rp)
      f_prim(30) = real(9.983754937E-01,rp)
      f_prim(31) = real(9.989728724E-01,rp)
      f_prim(32) = real(9.993625417E-01,rp)
      f_prim(33) = real(9.996117017E-01,rp)
      f_prim(34) = real(9.997678702E-01,rp)
      f_prim(35) = real(9.998638190E-01,rp)
      f_prim(36) = real(9.999216041E-01,rp)
      f_prim(37) = real(9.999557173E-01,rp)
      f_prim(38) = real(9.999754577E-01,rp)
      f_prim(39) = real(9.999866551E-01,rp)
      f_prim(40) = real(9.999928812E-01,rp)
      f_prim(41) = real(9.999962745E-01,rp)
      f_prim(42) = real(9.999980875E-01,rp)
      f_prim(43) = real(9.999990369E-01,rp)
      f_prim(44) = real(9.999995242E-01,rp)
      f_prim(45) = real(9.999997695E-01,rp)

   !$acc update device(eta_b(:))
   !$acc update device(f(:))
   !$acc update device(f_prim(:))

  end subroutine BLTSBDRLFlowSolver_fillBlasius


end module BLTSBDRLFlowSolver_mod
