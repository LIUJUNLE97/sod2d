#define ACTUATION 1

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

   real(rp), allocatable, dimension(:)    :: eta_b,f,f_prim         !auxiliary data arrays needs to be here because how cuda works
   real(rp), allocatable, dimension(:,:)  :: rectangleControl       ! Two point coordinates that define each rectangle
   integer(4),  allocatable, dimension(:) :: actionMask             ! Mask that contains whether a point is in a rectangle control or not


   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BLTSBDRLFlowSolver

      real(rp) , public  ::  M, d0, U0, rho0, Red0, Re, to, po, mu, amp_tbs, x_start, x_rise, x_end, x_fall, x_rerise, x_restart, coeff_tbs
      character(len=64), public :: tag
      logical :: db_clustered = .false.
      integer(4), public       :: countPar                                   ! Number of points in a rectangle of control per partition
      character(512), public   :: fileControlName                            ! Input: path to the file that contains the points defining the rectangle controls
      integer(4), public         :: nRectangleControl                          ! Number of rectangle control
      real(rp), public :: periodEpisode, previousActuationTime, periodActuation, frequencyActuation, timeBeginActuation ! parameters of the actuation
      integer(4), public :: tw_write_interval
   contains
      procedure, public :: fillBCTypes => BLTSBDRLFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters => BLTSBDRLFlowSolver_initializeParameters
      procedure, public :: evalInitialConditions => BLTSBDRLFlowSolver_evalInitialConditions
      procedure, public :: initialBuffer => BLTSBDRLFlowSolver_initialBuffer
      procedure, public :: fillBlasius => BLTSBDRLFlowSolver_fillBlasius
      procedure, public :: smoothStep => BLTSBDRLFlowSolver_smoothStep
      procedure, public :: afterDt => BLTSBDRLFlowSolver_afterDt
      procedure, public :: getControlNodes => BLTSBDRLFlowSolver_getControlNodes
      procedure, public :: readControlRectangles => BLTSBDRLFlowSolver_readControlRectangles
      procedure, public :: computeReward => BLTSBDRLFlowSolver_computeReward
      procedure, public :: beforeTimeIteration => BLTSBDRLFlowSolver_beforeTimeIteration
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

      open(unit=443,file="./output_"//trim(adjustl(this%tag))//"/"//"control_fortran_raw.txt",status='replace')
      open(unit=444,file="./output_"//trim(adjustl(this%tag))//"/"//"control_fortran_smooth.txt",status='replace')
      open(unit=445,file="./output_"//trim(adjustl(this%tag))//"/"//"control_tw.txt",status='replace')
      call init_smartredis(client, this%nwitPar, this%nRectangleControl, trim(adjustl(this%tag)), this%db_clustered)
      if (mpi_rank .eq. 0) write(111, *) "SmartRedis initialised"
      call write_step_type(client, 1, "ensemble_"//trim(adjustl(this%tag))//".step_type")
   end subroutine BLTSBDRLFlowSolver_initSmartRedis
#endif

   subroutine BLTSBDRLFlowSolver_beforeTimeIteration(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this

#if (ACTUATION)
      call this%getControlNodes()
#endif

   end subroutine BLTSBDRLFlowSolver_beforeTimeIteration

   subroutine BLTSBDRLFlowSolver_afterTimeIteration(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
#ifdef SMARTREDIS
      call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
      if (mpi_rank .eq. 0) close(443)
      if (mpi_rank .eq. 0) close(444)
      if (mpi_rank .eq. 0) close(445)
      call end_smartredis(client)
#endif
      if (mpi_rank .eq. 0) close(446)
   end subroutine BLTSBDRLFlowSolver_afterTimeIteration

   subroutine BLTSBDRLFlowSolver_afterDt(this,istep)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      integer(4)             , intent(in)   :: istep
      integer(4) :: iNodeL, bcode,iNodeL2
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3
      integer(4) :: k,j,ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime
      real(rp) :: dy,fx1,fx2,xp
      real(rp) :: mul , yp,yc
      real(rp)  :: gradIsoV(ndime),gradIsoU(ndime)
      real(rp)  :: gradV(ndime),vl(nnode),fact,targ,gradU(ndime),ul(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
      real(8) :: Ftau_neg(ndime), Ftau_pos(ndime)
#if ACTUATION
#ifdef SMARTREDIS
      real(rp) :: action_global_instant(action_global_size)
#endif
#endif
      ! test with teh new condition

      cd = 1.0_rp
      lx = this%d0*2.5_rp
      ly = this%d0*2.5_rp
      xmin = -40.0_rp*this%d0

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
         if (step_type_mod .eq. 1) then
            call write_step_type(client, 2, "ensemble_"//trim(adjustl(this%tag))//".step_type")
            ! write(*,*) "Sod2D wrote step: 2"
         end if

         ! check if a new action is needed
         if (this%time - this%previousActuationTime .ge. this%periodActuation) then
            if (mpi_rank .eq. 0) write(111, *) "Performing SmartRedis comms"
            ! save old action values and time - useful for interpolating to new action_global values
            !$acc kernels
            action_global_previous(:) = action_global(:)
            this%previousActuationTime = this%previousActuationTime + this%periodActuation
            !$acc end kernels

            call this%update_witness(istep, 1) ! manual update of witness points
            call write_state(client, buffwit(:, 1, 1), "ensemble_"//trim(adjustl(this%tag))//".state") ! the streamwise velocity u
            ! write(*,*) "Sod2D wrote state(1:5): ", buffwit(1:5, 1, 1)
            call this%computeReward(bc_type_unsteady_inlet, Ftau_neg, Ftau_pos)
            call write_reward(client, Ftau_neg(1), Ftau_pos(1), "ensemble_"//trim(adjustl(this%tag))//".reward") ! the streamwise component tw_x
            ! write(*,*) "Sod2D wrote reward: ", Ftau_neg(1), Ftau_pos(1)
            if (mpi_rank .eq. 0) write(445,'(*(ES12.4,:,","))') this%time, Ftau_neg(1), Ftau_pos(1)
            if (mpi_rank .eq. 0) call flush(445)

            call read_action(client, "ensemble_"//trim(adjustl(this%tag))//".action") ! modifies action_global (the target control values)
            ! write(*,*) "Sod2D read action: ", action_global
            if (mpi_rank .eq. 0) write(443,'(*(ES12.4,:,","))') this%time+this%periodActuation, action_global(1)
            if (mpi_rank .eq. 0) call flush(443)

            ! if the next time that we require actuation value is the last one, write now step_type=0 into database
            if (this%time + 2.0 * this%periodActuation .gt. this%maxPhysTime) then
               call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
               ! write(*,*) "Sod2D wrote step: 0"
            end if
         end if

         call this%smoothControlFunction(action_global_instant)
         if (mpi_rank .eq. 0) write(444,'(*(ES12.4,:,","))') this%time, action_global_instant(1)
         call flush(444)
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then
                  u_buffer(iNodeL,1) = 0.0_rp
                  u_buffer(iNodeL,2) = action_global(actionMask(iNodeL))
                  u_buffer(iNodeL,3) = 0.0_rp
            end if
         end do
         !$acc end parallel loop
      end if
#endif
#endif
      if (mod(istep, this%tw_write_interval) .eq. 0) then
         call this%computeReward(bc_type_unsteady_inlet, Ftau_neg, Ftau_pos)
         if (mpi_rank .eq. 0) write(446,'(*(ES12.4,:,","))') this%time, Ftau_neg(1), Ftau_pos(1)
         if (mpi_rank .eq. 0) call flush(446)
      end if

      !$acc parallel loop private(vl,dlxi_ip, dleta_ip, dlzeta_ip,gradIsoV,gradV,gradIsoU,gradU,ul)
      do iNodeL2 = 1,numWorkingNodesRankPar
         iNodeL = workingNodesPar(iNodeL2)
         yp = coordPar(iNodeL,2)
         xp = coordPar(iNodeL,1)
         if((xp .gt.-50.0_rp) .and. (xp .lt. 950.0_rp) .and. (yp .gt. 100.0_rp) .and. (yp .lt. 110.0_rp)) then
               ielem = point2elem(iNodeL)
               !$acc loop seq
               do inode = 1,nnode
                  vl(inode) = u(connecParWork(ielem,inode),2,2)
                  !ul(inode) = u(connecParWork(ielem,inode),1,2)
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

               gradIsoV(:) = 0.0_rp
              !gradIsoU(:) = 0.0_rp
               !$acc loop seq
               do ii=1,porder+1
                  gradIsoV(1) = gradIsoV(1) + dlxi_ip(ii)*vl(invAtoIJK(ii,isoJ,isoK))
                  gradIsoV(2) = gradIsoV(2) + dleta_ip(ii)*vl(invAtoIJK(isoI,ii,isoK))
                  gradIsoV(3) = gradIsoV(3) + dlzeta_ip(ii)*vl(invAtoIJK(isoI,isoJ,ii))
                  !gradIsoU(1) = gradIsoU(1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK))
                  !gradIsoU(2) = gradIsoU(2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK))
                  !gradIsoU(3) = gradIsoU(3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii))
               end do

               gradV(:) = 0.0_rp
               !gradU(:) = 0.0_rp
               !$acc loop seq
               do idime=1, ndime
                  !$acc loop seq
                  do jdime=1, ndime
                     gradV(idime) = gradV(idime) + He(idime,jdime,igaus,ielem) * gradIsoV(jdime)
                     !gradU(idime) = gradU(idime) + He(idime,jdime,igaus,ielem) * gradIsoU(jdime)
                  end do
                end do

               !if(gradU(2) .lt. gradV(1)) then
               !   u_buffer(iNodeL,1) = u_buffer(iNodeL,1)*1.01_rp
               !else
               !   u_buffer(iNodeL,1) = u_buffer(iNodeL,1)*0.99_rp
               !endif
               !u_buffer(iNodeL,1) = max(u_buffer(iNodeL,1),-this%U0)
               !u_buffer(iNodeL,1) = min(u_buffer(iNodeL,1),this%U0)
               iCen = connecParWork(ielem,atoIJK(nnode))
               yc =  coordPar(iCen,2)
               u_buffer(iNodeL,1) = (gradV(1))*(yp-yc) + u(iCen,1,2)
         end if
      end do
      !$acc end parallel loop
      !Filtering u buffer at the top of the domain
      if(mpi_size.ge.2) then
         call nvtxStartRange("MPI_comms_tI")
         call mpi_halo_max_boundary_update_real_iSendiRcv(u_buffer(:,1))
         call nvtxEndRange
      end if
      !$acc parallel loop gang
      do ielem = 1, numElemsRankPar
         iCen = connecParWork(ielem,atoIJK(nnode))
         yp = coordPar(iCen,2)
         xp = coordPar(iCen,1)

         !if((xp .gt.-50.0_rp) .and. (xp .lt. 950.0_rp) .and. (yp .gt. 80.0_rp)) then
         if((yp .gt. 100.0_rp)) then
            fact = 0.0_rp
            !$acc loop vector reduction(+:fact)
            do inode = 1, nnode
               fact = fact + u_buffer(connecParWork(ielem,inode),1)
            end do
            fact = fact/real(nnode,rp)
            !$acc loop vector
            do inode = 1, nnode
               !$acc atomic write
               u_buffer(connecParWork(ielem,inode),1) = fact
               !$acc end atomic
            end do
         end if
      end do
      !$acc end loop
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

   subroutine BLTSBDRLFlowSolver_computeReward(this, surfCode, Ftau_neg, Ftau_pos)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      integer(4), intent(in) :: surfCode
      real(8), intent(out) :: Ftau_neg(ndime), Ftau_pos(ndime)

      call twInfo(numElemsRankPar, numNodesRankPar, numBoundsRankPar, surfCode, connecParWork, boundPar, &
         point2elem, bouCodesPar, boundNormalPar, invAtoIJK, gmshAtoI, gmshAtoJ, gmshAtoK, wgp_b, dlxigp_ip, He, coordPar, &
         mu_fluid, mu_e, mu_sgs, rho(:,2), u(:,:,2), Ftau_neg, Ftau_pos)

   end subroutine BLTSBDRLFlowSolver_computeReward

   subroutine BLTSBDRLFlowSolver_readControlRectangles(this)
      ! This subroutine reads the file that contains the two points defining a rectanle paralel to the X-Z axis. Several rectangles
      ! can be introduced. In this rectangles is where control will be applied
      class(BLTSBDRLFlowSolver), intent(inout) :: this

      integer(rp)                           :: ii

      open(unit=99, file=trim(adjustl(this%fileControlName)), status='old', action='read')

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
      !$acc end parallel loop
      !$acc update device(actionMask(:))

   end subroutine BLTSBDRLFlowSolver_getControlNodes

   subroutine BLTSBDRLFlowSolver_fill_BC_Types(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
#if (ACTUATION)
      bouCodes2BCType(1) = bc_type_unsteady_inlet ! wall + actuation
#else
      bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! wall
#endif
      bouCodes2BCType(2) = bc_type_far_field         ! Upper part of the domain
      bouCodes2BCType(3) = bc_type_far_field          ! inlet part of the domain
      bouCodes2BCType(4) = bc_type_far_field          ! outlet part of the domain
      !$acc update device(bouCodes2BCType(:))

   end subroutine BLTSBDRLFlowSolver_fill_BC_Types

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

   subroutine BLTSBDRLFlowSolver_initializeParameters(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      real(rp) :: mur
      integer :: num_args, equal_pos, iarg, ierr
      character(len=64) :: arg, output_dir
      character(len=8) :: restart_step_str="", db_clustered_str="0", frequencyActuation_str="1.0", periodEpisode_str="1.0", timeBeginActuation_str="0.0"
      logical :: output_dir_exists

      ! get command line args, ie: mpirun -n 4 sod2d --restart_step=2500
      num_args = command_argument_count()
      do iarg = 1, num_args
         call get_command_argument(iarg, arg)
         equal_pos = scan(adjustl(trim(arg)), "=")
         if (adjustl(trim(arg(:equal_pos-1))) .eq. "--restart_step") then
            restart_step_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--db_clustered") then
            db_clustered_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--f_action") then
            frequencyActuation_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--t_episode") then
            periodEpisode_str = trim(adjustl(arg(equal_pos+1:)))
         else if (adjustl(trim(arg(:equal_pos-1))) .eq. "--t_begin_control") then
            timeBeginActuation_str = trim(adjustl(arg(equal_pos+1:)))
         else
            stop "Unknown command line argument"
         end if
      end do

      write(this%tag, *) app_color
      if (db_clustered_str == "" .or. db_clustered_str == "0") this%db_clustered = .false.
      if (restart_step_str == "" .or. restart_step_str == "0") then
         this%loadRestartFile = .false.
         this%restartFile_to_load = 0
      else
         this%loadRestartFile = .true.
         read(restart_step_str, *) this%restartFile_to_load ! 1: baseline, 2: last episode
      end if

      read(frequencyActuation_str,*,iostat=ierr) this%frequencyActuation
      this%periodActuation = 1.0_rp / this%frequencyActuation
      read(periodEpisode_str,*,iostat=ierr) this%periodEpisode
      read(timeBeginActuation_str,*,iostat=ierr) this%timeBeginActuation

      ! create output dir if not existing and copy the baseline restarts.
      ! when a random restart is selected and it is not the first episode, it will only create a
      ! copy of the *_1.h5 file, so that the random selection from python is either 1 (baseline)
      ! or 2 (last episode result) saved as *_2.h5
      output_dir = "./output_"//trim(adjustl(this%tag))//"/"
      inquire(file=trim(adjustl(output_dir)), exist=output_dir_exists)
      if (mpi_rank .eq. 0) then
         if (.not. output_dir_exists) call system("mkdir -p "//trim(adjustl(output_dir)))
         if (this%restartFile_to_load .gt. 0) call system("cp restart/* "//trim(adjustl(output_dir)))
      end if

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "bl_les"
      write(this%results_h5_file_path,*) "./output_"//trim(adjustl(this%tag))//"/"
      write(this%results_h5_file_name,*) "results_"//trim(adjustl(this%tag))
      write(this%io_prepend_path,*) output_dir


      if (mpi_rank .eq. 0) then
         write(*,*) "Received arguments are:"
         write(*,*) "--tag: ", adjustl(trim(this%tag))
         write(*,*) "--restart_step: ", this%restartFile_to_load
         write(*,*) "--db_clustered: ", db_clustered_str
         write(*,*) "--f_action: ", this%frequencyActuation
         write(*,*) "--t_episode: ", this%periodEpisode
         write(*,*) "--t_begin_control: ", this%timeBeginActuation
      end if

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 10000001
      this%maxPhysTime =  this%timeBeginActuation + this%periodEpisode

      this%save_logFile_first = 1
      this%save_logFile_step = 250

      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 10000

      this%continue_oldLogs = .false.

      this%initial_avgTime = 0.0 ! 3000.0_rp
      this%saveAvgFile = .true.
      this%loadAvgFile = .false. ! .true.

      ! wall shear stress output
      if (mpi_rank .eq. 0) open(unit=446,file="./output_"//trim(adjustl(this%tag))//"/"//"tw.txt",status='replace')
      this%tw_write_interval = 10
      !----------------------------------------------

      ! numerical params
      flag_les = 1
      flag_implicit = 0
      this%cfl_conv = 1.0_rp !M0.1 1.5, M0.3 0.95
      this%cfl_diff = 1.0_rp

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

      !! x inlet
      flag_buffer_on_west = .true.
      flag_buffer_w_min = -50.0_rp
      flag_buffer_w_size = 50.0_rp

      ! y top
      flag_buffer_on_north = .true.
      flag_buffer_n_min =  100.0_rp
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
      write(this%fileControlName ,*) "rectangleControl.txt"
      this%previousActuationTime = this%timeBeginActuation

      !Blasius analytical function
      call this%fillBlasius()

      this%have_witness          = .true.
      this%witness_inp_file_name = "witness.txt"
      this%witness_h5_file_name  = "resultwit.h5"
      this%leapwit               = 1 ! (update witness ever n dts) | in this class, we update the witness points manually
      this%leapwitsave           = 10 ! how many dts are stored in buffer
      this%wit_save              = .true. ! save witness or not
      this%wit_save_u_i          = .true.
      this%wit_save_pr           = .false.
      this%wit_save_rho          = .false.
      this%continue_witness      = .false.

   end subroutine BLTSBDRLFlowSolver_initializeParameters

   subroutine BLTSBDRLFlowSolver_initialBuffer(this)
      class(BLTSBDRLFlowSolver), intent(inout) :: this
      integer(4) :: iNodeL,k,j,bcode,ielem,iCen
      real(rp) :: yp,eta_y,f_y,f_prim_y, f1, f2, f3,x,dy

!#if (ACTUATION)
!      call this%getControlNodes()
!#endif

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
            u_buffer(iNodeL,2) =  0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                              -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
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
         if(yp .gt. 100.0_rp) then
            u(iNodeL,2,2) =  0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                              -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
         end if
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


end module BLTSBDRLFlowSolver_mod
