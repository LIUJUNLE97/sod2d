#define SMOOTH 0
#define ACTUATION 0

module BLTSBFlowSolver_mod
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
   implicit none
   private

   real(rp), allocatable, dimension(:)   :: eta_b,f,f_prim     !auxiliary data arrays needs to be here because how cuda works
   real(rp), allocatable    :: rectangleControl(:,:)           ! Two point coordinates that define each rectangle
   integer(4),  allocatable :: actionMask(:)                   ! Mask that contains whether a point is in a rectangle control or not
   real(rp), allocatable    :: actionPar(:), actionNodePar(:)  ! Variables to define the control nodes

                             
   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BLTSBFlowSolver

      real(rp) , public  ::  M, d0, U0, rho0, Red0, Re, to, po, mu, amp_tbs, x_start, x_rise, x_end, x_fall, x_rerise, x_restart, coeff_tbs
      integer(4), public       :: countPar                                   ! Number of points in a rectangle of control per partition
      character(512), public   :: fileControlName                            ! Input: path to the file that contains the points defining the rectangle controls
      integer(4), public         :: nRectangleControl                          ! Number of rectangle control
      real(rp), public         :: amplitudeActuation, frequencyActuation     ! Parameters of the actuation
   contains
      procedure, public :: fillBCTypes           => BLTSBFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BLTSBFlowSolver_initializeParameters
      procedure, public :: evalInitialConditions => BLTSBFlowSolver_evalInitialConditions
      procedure, public :: initialBuffer => BLTSBFlowSolver_initialBuffer
      procedure, public :: fillBlasius => BLTSBFlowSolver_fillBlasius
      procedure, public :: smoothStep => BLTSBFlowSolver_smoothStep
      procedure, public :: afterDt => BLTSBFlowSolver_afterDt
      procedure, public :: getControlNodes => BLTSBFlowSolver_getControlNodes
      procedure, public :: readControlRectangles => BLTSBFlowSolver_readControlRectangles
   end type BLTSBFlowSolver
contains

   subroutine BLTSBFlowSolver_afterDt(this,istep)
      class(BLTSBFlowSolver), intent(inout) :: this
      integer(4)             , intent(in)   :: istep
      integer(4) :: iNodeL, bcode
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3,x

      cd = 1.0_rp
#if (SMOOTH)      
      lx = this%d0*3.5_rp
      ly = this%d0*1.5_rp
      xmin = -15.0_rp*this%d0
#else
      !lx = this%d0*1.5_rp
      !ly = this%d0*1.5_rp
      !xmin = 0.0_rp*this%d0
      lx = this%d0*3.5_rp
      ly = this%d0*1.5_rp
      xmin = -15.0_rp*this%d0
#endif
      xmax = xmin+lx

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(coordPar(iNodeL,2) < ly) then
            if((coordPar(iNodeL,1)<xmax)  .and. (coordPar(iNodeL,1)>xmin)) then
               source_term(iNodeL,1) = -0.5_rp*rho(iNodeL,2)*cd*u(iNodeL,1,2)*abs(u(iNodeL,1,2))/lx
               source_term(iNodeL,2) = 0.00_rp
               source_term(iNodeL,3) = 0.00_rp
#if (ACTUATION)               
            elseif (actionMask(iNodeL)>0) then  ! Control
               u_buffer(iNodeL,1) = 0
               u_buffer(iNodeL,2) = this%amplitudeActuation*sin(2.0_rp*3.14159265_rp*this%frequencyActuation*this%time)
               u_buffer(iNodeL,3) = 0
#endif               
            end if
         end if
      end do
      !$acc end parallel loop
   end subroutine BLTSBFlowSolver_afterDt

   subroutine BLTSBFlowSolver_readControlRectangles(this)
      ! This subroutine reads the file that contains the two points defining a rectanle paralel to the X-Z axis. Several rectangles
      ! can be introduced. In this rectangles is where control will be applied
      class(BLTSBFlowSolver), intent(inout) :: this

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
      
   end subroutine BLTSBFlowSolver_readControlRectangles

   subroutine BLTSBFlowSolver_getControlNodes(this)
      class(BLTSBFlowSolver), intent(inout) :: this

      integer(4) :: iNodeL, iRectangleControl, bcode
      real(rp)   :: xPoint, zPoint, x1RectangleControl, x2RectangleControl, z1RectangleControl, z2RectangleControl
      
      call this%readControlRectangles()
      
      allocate(actionMask (numNodesRankPar))
      !$acc enter data create(actionMask(:))

      !this%countPar = 0
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if (bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL)
            if (bcode == bc_type_unsteady_inlet) then ! we are on the wall and we need to check if node is in the control rectanle
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
       !!$acc end parallel loop

   end subroutine BLTSBFlowSolver_getControlNodes

   subroutine BLTSBFlowSolver_fill_BC_Types(this)
      class(BLTSBFlowSolver), intent(inout) :: this
#if (ACTUATION)
      bouCodes2BCType(1) = bc_type_unsteady_inlet ! wall + actuation
#else
      bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! wall
#endif      
      bouCodes2BCType(2) = bc_type_far_field          ! Upper part of the domain
      bouCodes2BCType(3) = bc_type_far_field          ! inlet part of the domain
      bouCodes2BCType(4) = bc_type_far_field          ! outlet part of the domain
      !$acc update device(bouCodes2BCType(:))

   end subroutine BLTSBFlowSolver_fill_BC_Types

   subroutine BLTSBFlowSolver_fillBlasius(this)
      class(BLTSBFlowSolver), intent(inout) :: this

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

  end subroutine BLTSBFlowSolver_fillBlasius

   subroutine BLTSBFlowSolver_initializeParameters(this)
      class(BLTSBFlowSolver), intent(inout) :: this
      real(rp) :: mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "bl_les"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 10000001

      this%save_logFile_first = 1
      this%save_logFile_step  = 50
      this%save_restartFile_first = 1
      this%save_restartFile_step = 20000
      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 20000

      this%loadRestartFile = .false.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%initial_avgTime = 3000.0_rp
      this%saveAvgFile = .true.
      this%loadAvgFile = .false.
      !----------------------------------------------

      ! numerical params
      flag_les = 1
      flag_implicit = 0
      maxIterNonLineal=200
      tol=1e-3
      pseudo_cfl =1.95_rp 
      flag_implicit_repeat_dt_if_not_converged = 0

      this%cfl_conv = 1.5_rp !M0.1 1.5, M0.3 0.95
      this%cfl_diff = 1.5_rp 
      !flag_use_constant_dt = 1
      !this%dt = 1e-1

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
      
      ! y top
      flag_buffer_on_north = .true.
      flag_buffer_n_min =  160.0_rp
      flag_buffer_n_size = 20.0_rp

      ! -------- Instantaneous results file -------------
      this%save_scalarField_rho        = .true.
      this%save_scalarField_muFluid    = .false.
      this%save_scalarField_pressure   = .true.
      this%save_scalarField_energy     = .false.
      this%save_scalarField_entropy    = .false.
      this%save_scalarField_csound     = .false.
      this%save_scalarField_machno     = .true.
      this%save_scalarField_divU       = .true.
      this%save_scalarField_qcrit      = .true.
      this%save_scalarField_muSgs      = .false.
      this%save_scalarField_muEnvit    = .false.
      this%save_vectorField_vel        = .true.
      this%save_vectorField_gradRho    = .false.
      this%save_vectorField_curlU      = .true.

      ! -------- Average results file -------------
      this%save_avgScalarField_rho     = .true.
      this%save_avgScalarField_pr      = .true.
      this%save_avgScalarField_mueff   = .true.
      this%save_avgVectorField_vel     = .true.
      this%save_avgVectorField_ve2     = .true.
      this%save_avgVectorField_vex     = .true.
      this%save_avgVectorField_vtw     = .true.

      ! control parameters
#if (ACTUATION)      
      this%fileControlName = 'rectangleControl.dat'
      this%amplitudeActuation = 0.05
      this%frequencyActuation = 0.0025 
#endif      

   end subroutine BLTSBFlowSolver_initializeParameters

   subroutine BLTSBFlowSolver_initialBuffer(this)
      class(BLTSBFlowSolver), intent(inout) :: this
      integer(4) :: iNodeL,k,j,bcode
      real(rp) :: yp,eta_y,f_y,f_prim_y, f1, f2, f3,x


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
#if (SMOOTH) 
         u_buffer(iNodeL,2) = 0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                              -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
#else
         !if(yp .gt. 160.0_rp) then
         if(yp .gt. 100.0_rp) then
            x = (coordPar(iNodeL,1)-this%x_start  )/        this%x_rise
            if(x<=0.0_rp) then
               f1 = 0.0_rp
            elseif(x<1.0_rp) then
               f1 = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
            else
               f1 = 1
            end if

            x = (coordPar(iNodeL,1)-this%x_end    )/(2.0_rp*this%x_fall  )+1.0_rp
            if(x<=0.0_rp) then
               f2 = 0.0_rp
            elseif(x<1.0_rp) then
               f2 = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
            else
               f2 = 1
            end if

            x = (coordPar(iNodeL,1)-this%x_restart)/        this%x_rerise
            if(x<=0.0_rp) then
               f3 = 0.0_rp
            elseif(x<1.0_rp) then
               f3 = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
            else
               f3 = 1
            end if               

            u_buffer(iNodeL,2) = (f1 - (1.0_rp + this%coeff_tbs)*f2 + this%coeff_tbs*f3)*this%amp_tbs
         end if
#endif         
      end do
      !$acc end parallel loop

      !u_buffer(inodel,1) = (-1.0_rp/160000.0_rp)*yp*yp + yp*(1.0_rp/200.0_rp)

   end subroutine BLTSBFlowSolver_initialBuffer

   subroutine BLTSBFlowSolver_evalInitialConditions(this)
      class(BLTSBFlowSolver), intent(inout) :: this
      integer(4) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL, idime, j,k,bcode
      real(rp) :: yp,eta_y,f_y,f_prim_y, f1, f2, f3, x
      integer(4)   :: iLine,iNodeGSrl,auxCnt

      call this%fillBlasius()
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
#if (SMOOTH) 
           u(iNodeL,2,2) = 0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                              -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
#else       
         !if(yp .gt. 160.0_rp) then
         if(yp .gt. 100.0_rp) then
            x = (coordPar(iNodeL,1)-this%x_start  )/        this%x_rise
            if(x<=0.0_rp) then
               f1 = 0.0_rp
            elseif(x<1.0_rp) then
               f1 = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
            else
               f1 = 1
            end if

            x = (coordPar(iNodeL,1)-this%x_end    )/(2.0_rp*this%x_fall  )+1.0_rp
            if(x<=0.0_rp) then
               f2 = 0.0_rp
            elseif(x<1.0_rp) then
               f2 = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
            else
               f2 = 1
            end if

            x = (coordPar(iNodeL,1)-this%x_restart)/        this%x_rerise
            if(x<=0.0_rp) then
               f3 = 0.0_rp
            elseif(x<1.0_rp) then
               f3 = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
            else
               f3 = 1
            end if               

            u(iNodeL,2,2) = (f1 - (1.0_rp + this%coeff_tbs)*f2 + this%coeff_tbs*f3)*this%amp_tbs
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
   end subroutine BLTSBFlowSolver_evalInitialConditions

   subroutine BLTSBFlowSolver_smoothStep(this,x,y)
      class(BLTSBFlowSolver), intent(inout) :: this
      real(rp), intent(in)  :: x
      real(rp), intent(out) :: y

      if(x<=0.0_rp) then
         y = 0.0_rp
      elseif(x<1.0_rp) then
         y = 1.0_rp/(1.0_rp+exp(1.0_rp/(x-1.0_rp)+1.0_rp/x))
      else
         y = 1
      end if

   end subroutine BLTSBFlowSolver_smoothStep


end module BLTSBFlowSolver_mod
