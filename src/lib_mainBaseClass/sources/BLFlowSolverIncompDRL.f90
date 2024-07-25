module BLFlowSolverIncompDRL_mod
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
   real(rp), public :: eta_b(45), f(45), f_prim(45)

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BLFlowSolverIncompDRL

      real(rp) , public  ::   d0, U0, rho0, Red0, mu, vmax_SB, xc_SB, eps_SB, psi_SB, xmin_tripping, lx_tripping, ly_tripping

      character(len=64), public :: tag, output_dir
      logical :: db_clustered = .false.
      character(512), public :: fileControlName
      integer(4), public :: nRectangleControl, n_pseudo_envs
      real(rp), public :: periodEpisode, previousActuationTime, nextActuationTime, frequencyActuation, periodActuation, timeBeginActuation, spanLength, time_coeff
      real(rp), allocatable, public ::  reward(:), lx_recirculation_avg(:)

   contains
      procedure, public :: fillBCTypes           => BLFlowSolverIncompDRL_fill_BC_Types
      procedure, public :: initializeParameters  => BLFlowSolverIncompDRL_initializeParameters
      procedure, public :: evalInitialConditions => BLFlowSolverIncompDRL_evalInitialConditions
      procedure, public :: initialBuffer => BLFlowSolverIncompDRL_initialBuffer
      procedure, public :: fillBlasius => BLFlowSolverIncompDRL_fillBlasius
      procedure, public :: afterDt => BLFlowSolverIncompDRL_afterDt
      procedure, public :: beforeTimeIteration => BLFlowSolverIncompDRL_beforeTimeIteration
      procedure, public :: afterTimeIteration => BLFlowSolverIncompDRL_afterTimeIteration
      procedure, public :: computeTauW => BLFlowSolverIncompDRL_computeTauW
      procedure, public :: getControlNodes => BLFlowSolverIncompDRL_getControlNodes
      procedure, public :: readControlRectangles => BLFlowSolverIncompDRL_readControlRectangles
      procedure, public :: initSmartRedis  => BLFlowSolverIncompDRL_initSmartRedis
      procedure, public :: smoothControlFunction => BLFlowSolverIncompDRL_smoothControlFunction

   end type BLFlowSolverIncompDRL
contains

   subroutine BLFlowSolverIncompDRL_readControlRectangles(this)

      class(BLFlowSolverIncompDRL), intent(inout) :: this

      integer(rp) :: ii

      open(unit=99, file=trim(adjustl(this%fileControlName)), status='old', action='read')

      read(99,*) this%nRectangleControl
      allocate(rectangleControl(2,2*this%nRectangleControl))
      !$acc enter data create(rectangleControl(:,:))
      do ii = 1, this%nRectangleControl
         read(99, *) rectangleControl(:,2*ii-1)  ! First point [xMin, zMin]
         read(99, *) rectangleControl(:,2*ii  )  ! Second point [xMax, zMax]
         read(99, *)
      end do
      close(99)
      !$acc update device(rectangleControl(:,:))

   end subroutine BLFlowSolverIncompDRL_readControlRectangles

   subroutine BLFlowSolverIncompDRL_getControlNodes(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this

      integer(4) :: iNodeL, iRectangleControl, bcode
      real(rp)   :: xPoint, zPoint, x1RectangleControl, x2RectangleControl, z1RectangleControl, z2RectangleControl

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

                  ! Mesh node coordinates
                  xPoint = coordPar(iNodeL,1)
                  zPoint = coordPar(iNodeL,3)

                  ! Control surfaces xMin, zMin, xMax, zMax coordinatess
                  x1RectangleControl = rectangleControl(1,2*iRectangleControl-1)
                  z1RectangleControl = rectangleControl(2,2*iRectangleControl-1)
                  x2RectangleControl = rectangleControl(1,2*iRectangleControl  )
                  z2RectangleControl = rectangleControl(2,2*iRectangleControl  )

                  ! Check if the mesh point is within the control surface limits
                  if (xPoint >= x1RectangleControl .and. xPoint <= x2RectangleControl .and. zPoint >= z1RectangleControl .and. zPoint <= z2RectangleControl) then
                     actionMask(iNodeL) = iRectangleControl
                     exit
                  endif
               end do
            end if
         end if
      end do
      !$acc end parallel loop
   end subroutine BLFlowSolverIncompDRL_getControlNodes

   subroutine BLFlowSolverIncompDRL_initSmartRedis(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this

      ! Open files to save the action and reward
      open(unit=443,file=trim(adjustl(this%output_dir))//"control_action.txt",status='replace')
      open(unit=445,file=trim(adjustl(this%output_dir))//"control_reward.txt",status='replace')

      ! Initialise smartredis database
      call init_smartredis(client, this%nwitPar, this%nRectangleControl, this%n_pseudo_envs, trim(adjustl(this%tag)), this%db_clustered)
      if (mpi_rank .eq. 0) write(111, *) "SmartRedis initialised"

      ! Write the current step type in the database: 1 == Initialising time step
      call write_step_type(client, 1, "ensemble_"//trim(adjustl(this%tag))//".step_type")

   end subroutine BLFlowSolverIncompDRL_initSmartRedis


   subroutine BLFlowSolverIncompDRL_beforeTimeIteration(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this
      integer(4) :: iboun,bcode,ipbou,iBoundNode,iNodeL

     !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL) ! Boundary element code
            if (bcode == bc_type_far_field_SB) then ! top far field for BL
               u_buffer(iNodeL,2) =  this%vmax_SB*this%U0*(2.0_rp**0.5_rp)*(this%xc_SB-coordPar(iNodeL,1))/this%eps_SB*exp(this%psi_SB-((this%xc_SB &
                                 -coordPar(iNodeL,1))/this%eps_SB)**2_rp) ! v_top to force separation bubble
               u(iNodeL,2,2) = u_buffer(iNodeL,2)
            end if
         end if
      end do
      !$acc end parallel loop

      ! Open file to save the instantaneous action
      if (mpi_rank .eq. 0) open(unit=444,file=trim(adjustl(this%output_dir))//"control_action_smooth.txt",status='replace')

      ! Open file to save the instantaneous separation length
      if (mpi_rank .eq. 0) open(unit=446,file=trim(adjustl(this%output_dir))//"lx.txt",status='replace')

      ! Open file to save the averaged the separation length (over an action period)
      if (mpi_rank .eq. 0) open(unit=447,file=trim(adjustl(this%output_dir))//"lx_avg.txt",status='replace')

      ! Obtain the mask of the control nodes (>0)
      call this%getControlNodes()

      ! Initialise Smartredis
      call this%initSmartRedis()

      allocate(this%lx_recirculation_avg(this%n_pseudo_envs))
      !$acc enter data create(this%lx_recirculation_avg(:))
      allocate(this%reward(this%n_pseudo_envs))
      !$acc enter data create(this%reward(:))
      
      this%lx_recirculation_avg(:) = 0.0d0
      !$acc update device(this%lx_recirculation_avg(:))

      this%time_coeff = 0.0

      ! The first "nextActuationTime" is the initial actuation time
      this%nextActuationTime = this%timeBeginActuation 

   end subroutine BLFlowSolverIncompDRL_beforeTimeIteration

   subroutine BLFlowSolverIncompDRL_afterTimeIteration(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this

      if (mpi_rank .eq. 0) close(444)
      if (mpi_rank .eq. 0) close(446)
      if (mpi_rank .eq. 0) close(443)
      if (mpi_rank .eq. 0) close(445)
      if (mpi_rank .eq. 0) close(447)
      call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
      call end_smartredis(client)

   end subroutine BLFlowSolverIncompDRL_afterTimeIteration

   subroutine BLFlowSolverIncompDRL_afterDt(this,istep)
      class(BLFlowSolverIncompDRL), intent(inout) :: this
      integer(4) , intent(in)   :: istep
      real(rp) :: cd, lx, ly, xmin, xmax
      integer(4) :: ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime,iNodeL,bcode,isoII, isoJJ, isoKK,type_ijk
      integer(4) :: invert
      real(rp) :: xp, yp, yc
      real(rp)  :: gradIsoV(ndime)
      real(rp)  :: gradV(ndime),vl(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
      real(rp) :: lx_recirculation(this%n_pseudo_envs)

      real(rp) :: action_global_instant(action_global_size)
      real(rp) :: eliti, ave  ! Weights used to perform the averaging

      ! Tripping transition
      cd = 1.0_rp

      lx = this%lx_tripping*this%d0
      ly = this%ly_tripping*this%d0

      xmin = this%xmin_tripping*this%d0
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

      !$acc parallel loop private(vl,dlxi_ip, dleta_ip, dlzeta_ip,gradIsoV,gradV)
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2)
         xp = coordPar(iNodeL,1)

         if(bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL) ! Boundary element code
            if (bcode == bc_type_far_field_SB) then ! top far field for BL
               u_buffer(iNodeL,1) = 0.0_rp
               ielem = point2elem(iNodeL)
               !$acc loop seq
               do inode = 1,nnode
                  vl(inode) = u(connecParWork(ielem,inode),2,2)
               end do
               igaus = minloc(abs(connecParWork(ielem,:)-iNodeL),1)
               !$acc loop seq
               do ii=1,porder+1
                  dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                  dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                  dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
               end do
               isoI = gmshAtoI(igaus)
               isoJ = gmshAtoJ(igaus)
               isoK = gmshAtoK(igaus)
               type_ijk = 0
               if((isoI .eq. 1) .or. (isoI .eq. 2)) then
                  type_ijk = 1
                  if(isoI .eq. 1) then
                     isoII = 2
                  else
                     isoII = 1
                  end if
               else if ((isoJ .eq. 1) .or. (isoJ .eq. 2)) then
                  type_ijk = 2
                  if(isoJ .eq. 1) then
                     isoJJ = 2
                  else
                     isoJJ = 1
                  end if
               else
                  type_ijk = 3
                  if(isoK .eq. 1) then
                     isoKK = 2
                  else
                     isoKK = 1
                  end if
               end if

               gradIsoV(:) = 0.0_rp
               !$acc loop seq
               do ii=1,porder+1
                  gradIsoV(1) = gradIsoV(1) + dlxi_ip(ii)*vl(invAtoIJK(ii,isoJ,isoK))
                  gradIsoV(2) = gradIsoV(2) + dleta_ip(ii)*vl(invAtoIJK(isoI,ii,isoK))
                  gradIsoV(3) = gradIsoV(3) + dlzeta_ip(ii)*vl(invAtoIJK(isoI,isoJ,ii))
               end do

               gradV(:) = 0.0_rp
               !$acc loop seq
               do idime=1, ndime
                  !$acc loop seq
                  do jdime=1, ndime
                     gradV(idime) = gradV(idime) + He(idime,jdime,igaus,ielem) * gradIsoV(jdime)
                  end do
               end do

               if(type_ijk .eq. 1) then
                  iCen = connecParWork(ielem,invAtoIJK(isoII,isoJ,isoK))
               else if(type_ijk .eq. 2) then
                  iCen = connecParWork(ielem,invAtoIJK(isoI,isoJJ,isoK))
               else
                  iCen = connecParWork(ielem,invAtoIJK(isoI,isoJ,isoKK))
               end if

               yc =  coordPar(iCen,2)
               u_buffer(iNodeL,1)= u(iCen,1,2)-gradV(1)*(yc-yp)
            end if
         end if
      end do
      !$acc end parallel loop
      !Filtering u buffer at the top of the domain
      if(mpi_size.ge.2) then
         call nvtxStartRange("MPI_comms_tI")
         call mpi_halo_max_boundary_update_real_iSendiRcv(u_buffer(:,1))
         call nvtxEndRange
      end if

      ! Compute wall shear stress
      call this%computeTauW(lx_recirculation)
      if (this%save_logFile_next .eq. istep .and. mpi_rank .eq. 0) then
         write(446,'(*(ES16.6,:,","))') this%time, lx_recirculation
         call flush(446)
      end if

      !  Start actuating when t > timeBeginActuation
      if (this%time .gt. this%timeBeginActuation) then

         ! Compute average separation length
         ave = this%dt/(this%time_coeff+this%dt)
         eliti = this%time_coeff/(this%time_coeff+this%dt)

         this%time_coeff = this%time_coeff + this%dt

         this%lx_recirculation_avg(:) = ave * lx_recirculation(:) + eliti * this%lx_recirculation_avg(:)

         ! Compute reward
         this%reward(:) = - this%lx_recirculation_avg(:)
   
         ! Check if a new action is needed
         if (this%time .ge. this%nextActuationTime) then
            if (mpi_rank .eq. 0) write(111, *) "Performing SmartRedis comms"

            ! reset averaging time coefficient
            this%time_coeff = 0.0_rp

            ! save old action values and time - useful for interpolating to new action_global values
            !$acc kernels
            action_global_previous(:) = action_global(:)
            this%previousActuationTime = this%nextActuationTime
            this%nextActuationTime = this%previousActuationTime + this%periodActuation
            !$acc end kernels

            ! Write current state and reward into the smartedis database and then read action
            call this%update_witness(istep, 1) ! manual update of witness points
            call write_state(client, buffwit(:, 1, 1), "ensemble_"//trim(adjustl(this%tag))//".state") ! the streamwise velocity u

            call write_reward(client, this%reward, "ensemble_"//trim(adjustl(this%tag))//".reward") ! the streamwise component tw_x
            if (mpi_rank .eq. 0) write(445,'(*(ES16.6,:,","))') this%time, this%reward
            if (mpi_rank .eq. 0) call flush(445)
            if (mpi_rank .eq. 0) write(447,'(*(ES16.6,:,","))') this%time, this%lx_recirculation_avg
            if (mpi_rank .eq. 0) call flush(447)

            call read_action(client, "ensemble_"//trim(adjustl(this%tag))//".action") ! modifies action_global (the target control values)
            if (mpi_rank .eq. 0) write(443,'(*(ES16.6,:,","))') this%time+this%periodActuation, action_global
            if (mpi_rank .eq. 0) call flush(443)

            ! if the next time that we require actuation value is the last one, write now step_type=0 into database
            if (this%time + 1.0 * this%periodActuation .gt. this%maxPhysTime) then
               call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
            end if
         end if

         ! apply actions to control surfaces
         call this%smoothControlFunction(action_global_instant)
         if (this%save_logFile_next .eq. istep .and. mpi_rank .eq. 0) then
            write(444,'(*(ES16.6,:,","))') this%time, action_global_instant
            call flush(444)
         end if

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then
                  u_buffer(iNodeL,1) = 0.0_rp
                  u_buffer(iNodeL,2) = action_global_instant(actionMask(iNodeL))
                  u_buffer(iNodeL,3) = 0.0_rp
            end if
         end do
         !$acc end parallel loop

      end if

   end subroutine BLFlowSolverIncompDRL_afterDt

   subroutine BLFlowSolverIncompDRL_computeTauW(this, lx_recirculation)
      class(BLFlowSolverIncompDRL), intent(inout) :: this
      real(rp), intent(out) :: lx_recirculation(this%n_pseudo_envs)
      real(rp) :: lx_r
      integer(4) :: surfCode

      !$acc loop seq
      do surfCode=1, this%n_pseudo_envs
         call twInfo(numElemsRankPar, numNodesRankPar, numBoundsRankPar, surfCode, 1, connecParWork, boundPar, &
            point2elem, bouCodesPar, boundNormalPar, invAtoIJK, gmshAtoI, gmshAtoJ, gmshAtoK, wgp_b, dlxigp_ip, He, coordPar, &
            mu_fluid, mu_e, mu_sgs, rho(:,2), u(:,:,2), lx_r)
         lx_recirculation(surfCode) = lx_r / (this%spanLength / dble(this%n_pseudo_envs))
      end do
   end subroutine BLFlowSolverIncompDRL_computeTauW

   subroutine BLFlowSolverIncompDRL_smoothControlFunction(this, action_global_instant)
      class(BLFlowSolverIncompDRL), intent(inout) :: this
      real(rp), intent(inout) :: action_global_instant(action_global_size)
      
      real(rp) :: f1, f2, f3
      
      f1 = exp(-1.0_rp / ((this%time - this%previousActuationTime) / this%periodActuation))
      f2 = exp(-1.0_rp / (1.0_rp - (this%time - this%previousActuationTime) / this%periodActuation))
      f3 = f1 / (f1 + f2)
      !$acc kernels
      action_global_instant(:) = action_global_previous(:) + f3 * (action_global(:) - action_global_previous(:))
      !$acc end kernels

   end subroutine BLFlowSolverIncompDRL_smoothControlFunction

   subroutine BLFlowSolverIncompDRL_fill_BC_Types(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BLFlowSolverIncompDRL_fill_BC_Types

   subroutine BLFlowSolverIncompDRL_fillBlasius(this)
     class(BLFlowSolverIncompDRL), intent(inout) :: this

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

  end subroutine BLFlowSolverIncompDRL_fillBlasius

   subroutine BLFlowSolverIncompDRL_initializeParameters(this)
      use json_module
      implicit none
      class(BLFlowSolverIncompDRL), intent(inout) :: this
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
      call json%get("mesh_h5_file_name",value, found,"bl"); call this%checkFound(found,found_aux)
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

      call json%get("initial_avgTime",this%initial_avgTime, found, 0.0); call this%checkFound(found,found_aux)
      call json%get("saveAvgFile" ,this%saveAvgFile, found, .false.); call this%checkFound(found,found_aux)
      call json%get("loadAvgFile" ,this%loadAvgFile, found, .false.); call this%checkFound(found,found_aux)

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)    

      !----------------------------------------------

      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("U0",this%U0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("d0",this%d0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Red0",this%Red0, found,100.0_rp); call this%checkFound(found,found_aux)

      this%mu = this%rho0*this%d0*this%U0/this%Red0

      incomp_viscosity = this%mu
      flag_mu_factor = 1.0_rp

      nscbc_u_inf = this%U0
      nscbc_p_inf = 0.0_rp
      nscbc_rho_inf = this%rho0

      ! Separation bubble top velocity parameters
      call json%get("vmax_SB",this%vmax_SB, found,0.4); call this%checkFound(found,found_aux)
      call json%get("xc_SB",this%xc_SB, found,306.640625); call this%checkFound(found,found_aux)
      call json%get("eps_SB",this%eps_SB, found,110.485435); call this%checkFound(found,found_aux)
      call json%get("psi_SB",this%psi_SB, found,0.95); call this%checkFound(found,found_aux)

      ! Tripping region
      call json%get("xmin_tripping",this%xmin_tripping, found,0.0); call this%checkFound(found,found_aux)
      call json%get("ymin_tripping",this%xmin_tripping, found,0.0); call this%checkFound(found,found_aux)
      call json%get("lx_tripping",this%lx_tripping, found,1.0); call this%checkFound(found,found_aux)
      call json%get("ly_tripping",this%ly_tripping, found,1.0); call this%checkFound(found,found_aux)

      ! Actuation Parameters
      call json%get("fileControlName",value, found,"rectangleControl.txt"); call this%checkFound(found,found_aux)       ! Jet surface limits (file name)
      write(this%fileControlName,*) value

      ! Parameters to compute mean separation length
      call json%get("spanLength",this%spanLength, found,1.0_rp); call this%checkFound(found,found_aux)  ! Spanwise length
      call json%get("n_pseudo_envs",this%n_pseudo_envs, found,1); call this%checkFound(found,found_aux)  ! Number of pseudo-environemnts (extrusions)

      ! witness points
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

   end subroutine BLFlowSolverIncompDRL_initializeParameters

   subroutine BLFlowSolverIncompDRL_initialBuffer(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this
      integer(4) :: iNodeL,k,j
      real(rp) :: yp,eta_y,f_y,f_prim_y


      call this%fillBlasius()

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
      end do
      !$acc end parallel loop

   end subroutine BLFlowSolverIncompDRL_initialBuffer

   subroutine BLFlowSolverIncompDRL_evalInitialConditions(this)
      class(BLFlowSolverIncompDRL), intent(inout) :: this
      integer(4) :: iNodeL, idime, j,k
      real(rp) :: yp,eta_y,f_y,f_prim_y

      call this%fillBlasius()

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

         pr(iNodeL,2) = 0.0_rp
         rho(iNodeL,2) = this%rho0
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
   end subroutine BLFlowSolverIncompDRL_evalInitialConditions

#endif
end module BLFlowSolverIncompDRL_mod
