module BluffBodySolverIncompDRL_mod
#ifdef SMARTREDIS
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
   use mod_smartredis
   use CFDSolverPeriodicWithBoundariesIncomp_mod
   implicit none
   private

   real(rp), allocatable, dimension(:,:)  :: rectangleControl
   integer(4),  allocatable, dimension(:) :: actionMask

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BluffBodySolverIncompDRL

      real(rp) , public  :: vo, delta, rho0, Re, aoa

      character(len=64), public :: tag, output_dir
      logical :: db_clustered = .false.
      integer(4), public :: nRectangleControl, n_pseudo_envs
      character(512), public :: fileControlName
      real(rp), public :: timeBeginActuation, frequencyActuation, alphaActuation
      real(rp), public :: periodEpisode, previousActuationTime, nextActuationTime, periodActuation, time_coeff
      real(rp), public :: spanLength, alphaReward, Cd_base
      real(rp), allocatable, public :: Cd_avg(:), Cl_avg(:), reward(:)

   contains
      procedure, public :: readControlRectangles => BluffBodySolverIncompDRL_readControlRectangles
      procedure, public :: getControlNodes       => BluffBodySolverIncompDRL_getControlNodes
      procedure, public :: beforeTimeIteration => BluffBodySolverIncompDRL_beforeTimeIteration
      procedure, public :: afterTimeIteration => BluffBodySolverIncompDRL_afterTimeIteration
      procedure, public :: afterDt => BluffBodySolverIncompDRL_afterDt
      procedure, public :: fillBCTypes           => BluffBodySolverIncompDRL_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolverIncompDRL_initializeParameters
      procedure, public :: evalInitialConditions => BluffBodySolverIncompDRL_evalInitialConditions
      procedure, public :: initialBuffer => BluffBodySolverIncompDRL_initialBuffer
      procedure, public :: initSmartRedis  => BluffBodySolverIncompDRL_initSmartRedis
      procedure, public :: computeClCd => BluffBodySolverIncompDRL_computeClCd
      procedure, public :: smoothControlFunction => BluffBodySolverIncompDRL_smoothControlFunction
   end type BluffBodySolverIncompDRL
contains

   subroutine BluffBodySolverIncompDRL_readControlRectangles(this)

      class(BluffBodySolverIncompDRL), intent(inout) :: this

      integer(rp) :: ii

      open(unit=99, file=trim(adjustl(this%fileControlName)), status='old', action='read')

      read(99,*) this%nRectangleControl

      allocate(rectangleControl(3,2*this%nRectangleControl))
      !$acc enter data create(rectangleControl(:,:))
      do ii = 1, this%nRectangleControl
         read(99, *) rectangleControl(:,2*ii-1)  ! First point [xMin, yMin, zMin]
         read(99, *) rectangleControl(:,2*ii  )  ! Second point [xMax, yMax, zMax]
         read(99, *)
      end do
      close(99)
      !$acc update device(rectangleControl(:,:))

   end subroutine BluffBodySolverIncompDRL_readControlRectangles

   subroutine BluffBodySolverIncompDRL_getControlNodes(this)

      class(BluffBodySolverIncompDRL), intent(inout) :: this

      integer(4) :: iNodeL, iRectangleControl, bcode
      real(rp)   :: xPoint, yPoint, zPoint, x1RectangleControl, x2RectangleControl, y1RectangleControl, y2RectangleControl, z1RectangleControl, z2RectangleControl

      call this%readControlRectangles()

      allocate(actionMask(numNodesRankPar))
      !$acc enter data create(actionMask(:))

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         actionMask(iNodeL) = 0

         if (bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL)
         
            if (bcode == bc_type_unsteady_inlet) then

               do iRectangleControl = 1,this%nRectangleControl

                  ! Mesh coordinates
                  xPoint = coordPar(iNodeL,1)
                  yPoint = coordPar(iNodeL,2)
                  zPoint = coordPar(iNodeL,3)

                  ! Control surfaces xMin, yMin, zMin (1) and xMax, yMax, zMax (2) coordinatess
                  x1RectangleControl = rectangleControl(1,2*iRectangleControl-1)
                  y1RectangleControl = rectangleControl(2,2*iRectangleControl-1)
                  z1RectangleControl = rectangleControl(3,2*iRectangleControl-1)

                  x2RectangleControl = rectangleControl(1,2*iRectangleControl)
                  y2RectangleControl = rectangleControl(2,2*iRectangleControl)
                  z2RectangleControl = rectangleControl(3,2*iRectangleControl)

                  ! Check if the mesh point is within the control surface limits and create mask
                  if (xPoint >= x1RectangleControl .and. xPoint <= x2RectangleControl .and. &
                     yPoint >= y1RectangleControl .and. yPoint <= y2RectangleControl .and. &
                     zPoint >= z1RectangleControl .and. zPoint <= z2RectangleControl) then

                     actionMask(iNodeL) = iRectangleControl

                     exit

                  endif
               end do
            end if
         end if
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolverIncompDRL_getControlNodes

   subroutine BluffBodySolverIncompDRL_initSmartRedis(this)
      class(BluffBodySolverIncompDRL), intent(inout) :: this

      ! Open files to save the action and reward
      if (mpi_rank .eq. 0) open(unit=443,file=trim(adjustl(this%output_dir))//"control_action.txt",status='replace')
      if (mpi_rank .eq. 0) open(unit=445,file=trim(adjustl(this%output_dir))//"control_reward.txt",status='replace')

      ! Initialise smartredis database 
      call init_smartredis(client, this%nwitPar, this%nRectangleControl, this%n_pseudo_envs, trim(adjustl(this%tag)), this%db_clustered)

      ! Indicate in the log that Smartredis has been initialised
      if (mpi_rank .eq. 0) write(111, *) "SmartRedis initialised"

      ! Write the current step type in the database: 1 == Initialising time step
      call write_step_type(client, 1, "ensemble_"//trim(adjustl(this%tag))//".step_type")

   end subroutine BluffBodySolverIncompDRL_initSmartRedis

   subroutine BluffBodySolverIncompDRL_beforeTimeIteration(this)
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      integer(4) :: bcode, iNodeL

      ! Open file to save the instantaneous action
      if (mpi_rank .eq. 0) open(unit=444,file=trim(adjustl(this%output_dir))//"smooth_control_action.txt",status='replace')

      ! Open file to save the lift and drag coefficients
      if (mpi_rank .eq. 0) open(unit=446,file=trim(adjustl(this%output_dir))//"ClCd.txt",status='replace')
      if (mpi_rank .eq. 0) open(unit=447,file=trim(adjustl(this%output_dir))//"ClCd_avg.txt",status='replace')

      ! Obtain the mask of the control nodes
      call this%getControlNodes()

      ! Initially force U=0 at the walls
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if (bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL)
            if (bcode == bc_type_unsteady_inlet) then
         
               u_buffer(iNodeL,1) = 0.0_rp
               u_buffer(iNodeL,2) = 0.0_rp
               u_buffer(iNodeL,3) = 0.0_rp 

            end if
         end if 
      end do  
      !$acc end parallel loop

      ! Initialise Smartredis
      call this%initSmartRedis()

      ! Initialise averaged coefficients
      allocate(this%Cd_avg(this%n_pseudo_envs))
      !$acc enter data create(this%Cd_avg(:))
      allocate(this%Cl_avg(this%n_pseudo_envs))
      !$acc enter data create(this%Cl_avg(:))
      allocate(this%reward(this%n_pseudo_envs))
      !$acc enter data create(this%reward(:))

      this%Cd_avg(:) = 0.0d0
      !$acc update device(this%Cd_avg(:))
      this%Cl_avg(:) = 0.0d0 
      !$acc update device(this%Cl_avg(:))

      this%time_coeff = 0.0

      this%nextActuationTime = this%timeBeginActuation ! The first "nextActuationTime" is the initial actuation time

   end subroutine BluffBodySolverIncompDRL_beforeTimeIteration

   subroutine BluffBodySolverIncompDRL_afterTimeIteration(this)
      class(BluffBodySolverIncompDRL), intent(inout) :: this

      if (mpi_rank .eq. 0) close(444)
      if (mpi_rank .eq. 0) close(446)
      if (mpi_rank .eq. 0) close(443)
      if (mpi_rank .eq. 0) close(445)
      if (mpi_rank .eq. 0) close(447)
      call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
      call end_smartredis(client)

   end subroutine BluffBodySolverIncompDRL_afterTimeIteration

   subroutine BluffBodySolverIncompDRL_afterDt(this,istep)
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      integer(4), intent(in) :: istep
      integer(4) :: icode
      real(rp) :: CL(this%n_pseudo_envs), CD(this%n_pseudo_envs)

      integer(4) :: iNodeL, invert
      real(rp)   :: zPoint
      real(rp) :: action_classic

      real(rp) :: action_global_instant(action_global_size)
      real(rp) :: eliti, ave  ! Weights used to perform the averaging

      real(rp)   :: xPoint, yPoint, theta ! BORRAR


      ! Update the surface forces (Fpr, Ftau)
      if (isMeshBoundaries) then
         do icode = 1,numBoundCodes
            call surfInfo(istep,this%time,numElemsRankPar,numNodesRankPar,numBoundsRankPar,icode,connecParWork,boundPar,point2elem, &
               bouCodesPar,boundNormalPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coordPar, &
               mu_fluid,mu_e,mu_sgs,rho(:,2),u(:,:,2),pr(:,2),this%surfArea,Fpr(:,iCode),Ftau(:,iCode),.FALSE.)
         end do
      end if

      ! Compute drag and lift in the current time (at each pseudo-environment)
      call this%computeClCd(CL, CD)
      if(this%save_logFile_next==istep) then
         if (mpi_rank .eq. 0) then
            write(446,'(*(ES16.6,:,","))') this%time, CL, CD
         end if
         call flush(446)
      end if

      ! Update the step_type to 2 == Running
      if (step_type_mod .eq. 1) then
         call write_step_type(client, 2, "ensemble_"//trim(adjustl(this%tag))//".step_type")
      end if

      !  Start actuating when t > timeBeginActuation
      if (this%time .gt. this%timeBeginActuation) then

         ! Compute average Cl and Cd
         ave = this%dt / (this%time_coeff + this%dt)
         eliti = this%time_coeff / (this%time_coeff + this%dt)

         this%time_coeff = this%time_coeff + this%dt
         
         this%Cd_avg(:) = ave * CD(:) + eliti * this%Cd_avg(:)
         this%Cl_avg(:) = ave * CL(:) + eliti * this%Cl_avg(:)

         ! Compute reward
         this%reward(:) = this%Cd_base - this%Cd_avg(:) - this%alphaReward * abs(this%Cl_avg(:))

         ! Check if a new action is needed
         if (this%time .ge. this%nextActuationTime) then
            if (mpi_rank .eq. 0) write(111, *) "Performing SmartRedis comms"

            ! reset coefficient averaging time
            this%time_coeff = 0.0_rp

            ! save old action values and time - useful for interpolating to new action_global values
            !$acc kernels
            action_global_previous(:) = action_global(:)
            this%previousActuationTime = this%nextActuationTime
            this%nextActuationTime = this%previousActuationTime + this%periodActuation
            !$acc end kernels

            ! Write current state and reward into the smartedis database and then read action
            call this%update_witness(istep, 1) ! manual update of witness points
            !call write_state(client, buffwit(:, 1, 1), "ensemble_"//trim(adjustl(this%tag))//".state") ! the streamwise velocity u
            call write_state(client, buffwit(:, 1, 4), "ensemble_"//trim(adjustl(this%tag))//".state") ! pressure

            call write_reward(client, this%reward, "ensemble_"//trim(adjustl(this%tag))//".reward")
            if (mpi_rank .eq. 0) write(445,'(*(ES16.6,:,","))') this%time, this%reward
            if (mpi_rank .eq. 0) call flush(445)
            if (mpi_rank .eq. 0) write(447,'(*(ES16.6,:,","))') this%time, this%Cl_avg, this%Cd_avg
            if (mpi_rank .eq. 0) call flush(447)

            call read_action(client, "ensemble_"//trim(adjustl(this%tag))//".action") ! modifies action_global
            if (mpi_rank .eq. 0) write(443,'(*(ES16.6,:,","))') this%time+this%periodActuation, action_global
            if (mpi_rank .eq. 0) call flush(443)

            ! if the next time that we require actuation value is the last one, write now step_type=0 into database
            if (this%time + 1.0 * this%periodActuation .gt. this%maxPhysTime) then
               call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
            end if
         end if

         ! apply actions to control surfaces
         call this%smoothControlFunction(action_global_instant)
         if(this%save_logFile_next==istep) then
            if (mpi_rank .eq. 0) then
               write(444,'(*(ES16.6,:,","))') this%time, action_global_instant
            end if
            call flush(444)
         end if

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then

               !u_buffer(iNodeL,1) = action_global_instant(actionMask(iNodeL))*sin(this%alphaActuation*v_pi/180_rp)
               !u_buffer(iNodeL,2) = action_global_instant(actionMask(iNodeL))*cos(this%alphaActuation*v_pi/180_rp)
               !u_buffer(iNodeL,3) = 0.0_rp

               ! With a specific jet velocity function (action_global_instant here is Q==mass flow rate)
               xPoint = coordPar(iNodeL,1)
               yPoint = coordPar(iNodeL,2)
               if (yPoint < 0.0_rp) then ! For the lower surface, compute the jet velocity as for the upper surface (take the symmetry point)
                  yPoint = -yPoint
               end if
               theta = atan2((yPoint - 0.0_rp), (xPoint - 0.0_rp)) ! Angle of the current point

               u_buffer(iNodeL,1) = (action_global_instant(actionMask(iNodeL)) * ((v_pi)/(this%rho0*(10 *v_pi/180)*this%delta)) * cos((v_pi/(10 *v_pi/180))*(theta-(90 *v_pi/180)))) * cos(theta)
               u_buffer(iNodeL,2) = (action_global_instant(actionMask(iNodeL)) * ((v_pi)/(this%rho0*(10 *v_pi/180)*this%delta)) * cos((v_pi/(10 *v_pi/180))*(theta-(90 *v_pi/180)))) * sin(theta)
               u_buffer(iNodeL,3) = 0.0_rp

            end if
         end do
         !$acc end parallel loop

      end if

   end subroutine BluffBodySolverIncompDRL_afterDt

   subroutine BluffBodySolverIncompDRL_computeClCd(this, CL, CD)
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      integer(4) :: surfCode
      real(rp), intent(out)  :: CL(this%n_pseudo_envs), CD(this%n_pseudo_envs)
      real(rp) :: Fx(this%n_pseudo_envs), Fy(this%n_pseudo_envs)

      ! Compute the lift and drag coefficients (CL and CD) at each pseudo-environment
      !$acc loop seq
      do surfCode = 1, this%n_pseudo_envs

         ! Summ the pressure and viscous forces in the x and y direction
         !Fy(surfCode) = 2.0_rp * (Fpr(2,surfCode) + Ftau(2,surfCode)) / (this%rho0 * (this%delta*this%spanLength/this%n_pseudo_envs) * (this%vo*this%vo))
         !Fx(surfCode) = 2.0_rp * (Fpr(1,surfCode) + Ftau(1,surfCode)) / (this%rho0 * (this%delta*this%spanLength/this%n_pseudo_envs) * (this%vo*this%vo))

         ! Now consider the AoA
         !CL(surfCode) = -Fx(surfCode)*sin(this%AoA) + Fy(surfCode)*cos(this%AoA)
         !CD(surfCode) = Fx(surfCode)*cos(this%AoA) + Fy(surfCode)*sin(this%AoA)

         ! ONLY CONSIDER THE PRESSURE FORCES 
         CL(surfCode) = 2.0_rp * Fpr(2,surfCode) / (this%rho0 * (this%delta*this%spanLength/this%n_pseudo_envs) * (this%vo*this%vo))
         CD(surfCode) = 2.0_rp * Fpr(1,surfCode) / (this%rho0 * (this%delta*this%spanLength/this%n_pseudo_envs) * (this%vo*this%vo))  
         
      end do

   end subroutine BluffBodySolverIncompDRL_computeClCd

   subroutine BluffBodySolverIncompDRL_smoothControlFunction(this, action_global_instant)
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      real(rp), intent(inout) :: action_global_instant(action_global_size)
      
      real(rp) :: f1, f2, f3
      
      f1 = exp(-1.0_rp / ((this%time - this%previousActuationTime) / this%periodActuation))
      f2 = exp(-1.0_rp / (1.0_rp - (this%time - this%previousActuationTime) / this%periodActuation))
      f3 = f1 / (f1 + f2)
      !$acc kernels
      action_global_instant(:) = action_global_previous(:) + f3 * (action_global(:) - action_global_previous(:))
      !$acc end kernels

   end subroutine BluffBodySolverIncompDRL_smoothControlFunction

   subroutine BluffBodySolverIncompDRL_fill_BC_Types(this)
      class(BluffBodySolverIncompDRL), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BluffBodySolverIncompDRL_fill_BC_Types

   subroutine BluffBodySolverIncompDRL_initializeParameters(this)
      use json_module
      implicit none
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      real(rp) :: mul, mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      integer :: num_args, equal_pos, iarg, ierr
      character(len=64) :: arg
      character(len=8) :: restart_step_str="", periodEpisode_str="1.0", frequencyActuation_str="1.0", timeBeginActuation_str="0.0", db_clustered_str="0"
      logical :: output_dir_exists

      ! Get command line args, i.e.: mpirun -n 4 sod2d --restart_step=1 --t_episode=100.0
      num_args = command_argument_count()
      do iarg = 2, num_args                        ! ARGUMENT 1 IS THE JSON FILE
         call get_command_argument(iarg, arg)
         equal_pos = scan(adjustl(trim(arg)), "=") ! Find the position in the string where "=" is located
         if (adjustl(trim(arg(:equal_pos-1))) .eq. "--restart_step") then
            restart_step_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--t_episode") then
            periodEpisode_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--f_action") then
            frequencyActuation_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--t_begin_control") then
            timeBeginActuation_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--db_clustered") then
            db_clustered_str = trim(adjustl(arg(equal_pos+1:)))
         else
            stop "Unknown command line argument"
         end if
      end do

      write(this%tag, *) app_color

      ! db_clustered?
      if (db_clustered_str == "" .or. db_clustered_str == "0") this%db_clustered = .false.

      ! Which restart we will load?
      if (restart_step_str == "" .or. restart_step_str == "0") then
         this%loadRestartFile = .false.
         this%restartFile_to_load = 0
      else
         this%loadRestartFile = .true.
         read(restart_step_str, *) this%restartFile_to_load ! 1: baseline restart 1, 2: baseline restart 2
      end if
      
      ! Period of an episode
      read(periodEpisode_str,*,iostat=ierr) this%periodEpisode

      ! Frequency of an actuation
      read(frequencyActuation_str,*,iostat=ierr) this%frequencyActuation
      this%periodActuation = 1.0_rp / this%frequencyActuation

      ! Time to start actuation
      read(timeBeginActuation_str,*,iostat=ierr) this%timeBeginActuation

      ! Set the final simulation time as the episode duration
      this%maxPhysTime =  this%timeBeginActuation + this%periodEpisode

      if (mpi_rank .eq. 0) then
         write(*,*) "RL simulation params:"
         write(*,*) "--tag: ", adjustl(trim(this%tag))
         write(*,*) "--restart_step: ", this%restartFile_to_load
         write(*,*) "--f_action: ", this%frequencyActuation
         write(*,*) "--t_episode: ", this%periodEpisode
         write(*,*) "--t_begin_control: ", this%timeBeginActuation
         write(*,*) "--db_clustered: ", db_clustered_str
      end if


      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)

      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"cylin"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value

      call json%get("results_path",value, found,""); call this%checkFound(found,found_aux)
      this%output_dir = trim(adjustl(value))//"output_"//trim(adjustl(this%tag))//"/"

      write(this%results_h5_file_path,*) this%output_dir ! path to save h5
      write(this%io_prepend_path,*) this%output_dir  ! path to save surface and log files

      inquire(file=trim(adjustl(this%output_dir)), exist=output_dir_exists)
      if (mpi_rank .eq. 0) then
         if (.not. output_dir_exists) call system("mkdir -p "//trim(adjustl(this%output_dir)))             ! Create results directory
         if (this%restartFile_to_load .gt. 0) call system("cp restart/* "//trim(adjustl(this%output_dir))) ! Copy baseline restarts to start from those
      end if


      call json%get("results_h5_file_name",value, found,"results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value//"_"//trim(adjustl(this%tag))


      !  --------------  I/O params -------------
      call json%get("final_istep",this%final_istep, found,1000001); call this%checkFound(found,found_aux)

      call json%get("save_logFile_first",this%save_logFile_first, found, 1); call this%checkFound(found,found_aux)
      call json%get("save_logFile_step",this%save_logFile_step, found, 10); call this%checkFound(found,found_aux)

      call json%get("save_resultsFile_first",this%save_resultsFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_resultsFile_step" ,this%save_resultsFile_step, found,10000); call this%checkFound(found,found_aux)

      call json%get("save_restartFile_first",this%save_restartFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_restartFile_step" ,this%save_restartFile_step, found,10000); call this%checkFound(found,found_aux)

      call json%get("continue_oldLogs" ,this%continue_oldLogs, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveAvgFile" ,this%saveAvgFile, found, .true.); call this%checkFound(found,found_aux)
      call json%get("loadAvgFile" ,this%loadAvgFile, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)
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
      call json%get("aoa",this%aoa, found,0.0_rp); call this%checkFound(found,found_aux)

      call json%get("c_sgs",c_sgs, found,0.025_rp); 

      ! fixed by the type of base class parameters
      incomp_viscosity = (this%rho0*this%delta*this%vo)/this%Re
      flag_mu_factor = 1.0_rp

      nscbc_u_inf = this%vo
      nscbc_p_inf = 0.0_rp
      nscbc_rho_inf = this%rho0

      ! Actuation Parameters
      call json%get("fileControlName",value, found,"rectangleControl.txt"); call this%checkFound(found,found_aux)       ! Jet surface limits (file name)
      write(this%fileControlName,*) value
      call json%get("spanLength",this%spanLength, found,1.0_rp); call this%checkFound(found,found_aux)                  ! Spanwise length
      call json%get("n_pseudo_envs",this%n_pseudo_envs, found,1); call this%checkFound(found,found_aux)            ! Number of pseudo-environments
      call json%get("alphaActuation",this%alphaActuation, found,0.0_rp); call this%checkFound(found,found_aux)          ! Jet angle with respect the y-axis
      call json%get("Cd_base",this%Cd_base, found,0.0_rp); call this%checkFound(found,found_aux)                        ! Baseline drag (compute reward)
      call json%get("alphaReward",this%alphaReward, found,0.0_rp); call this%checkFound(found,found_aux)                ! Lift contr. weight (compute reward)

      !Witness points parameters
      call json%get("have_witness",this%have_witness, found,.false.)
      if(this%have_witness .eqv. .true.) then
         call json%get("witness_inp_file_name",value, found,"witness.txt"); call this%checkFound(found,found_aux)
         write(this%witness_inp_file_name,*) value
         call json%get("witness_h5_file_name",value, found,"resultwit.h5"); call this%checkFound(found,found_aux)
         write(this%witness_h5_file_name,*) trim(adjustl(this%output_dir))//value

         call json%get("leapwit",this%leapwit, found,1); call this%checkFound(found,found_aux)
         call json%get("leapwitsave",this%leapwitsave, found,100); call this%checkFound(found,found_aux)
         call json%get("nwit",this%nwit, found,17986); call this%checkFound(found,found_aux)
         call json%get("wit_save",this%wit_save, found,.true.); call this%checkFound(found,found_aux)
         call json%get("wit_save_u_i",this%wit_save_u_i, found,.true.); call this%checkFound(found,found_aux)
         call json%get("wit_save_pr",this%wit_save_pr, found,.true.); call this%checkFound(found,found_aux)
         call json%get("wit_save_rho",this%wit_save_rho, found,.true.); call this%checkFound(found,found_aux)
         call json%get("continue_witness",this%continue_witness, found,.false.); call this%checkFound(found,found_aux)
      end if  

      call this%readJSONBuffer()

      call json%destroy()

   end subroutine BluffBodySolverIncompDRL_initializeParameters

   subroutine BluffBodySolverIncompDRL_evalInitialConditions(this)
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo*cos(this%aoa*v_pi/180.0_rp)
         u(iNodeL,2,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)
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

   end subroutine BluffBodySolverIncompDRL_evalInitialConditions


   subroutine BluffBodySolverIncompDRL_initialBuffer(this)
      class(BluffBodySolverIncompDRL), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo*cos(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolverIncompDRL_initialBuffer

#endif
end module BluffBodySolverIncompDRL_mod
