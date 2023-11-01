#define ACTUATION 1

module BLMARLFlowSolverIncomp_mod
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

   real(rp), allocatable, dimension(:,:)  :: rectangleControl       ! Two point coordinates that define each rectangle
   integer(4),  allocatable, dimension(:) :: actionMask             ! Mask that contains whether a point is in a rectangle control or not

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BLMARLFlowSolverIncomp

      real(rp) , public  ::   d0, U0, rho0, Red0, Re, mu, Lz
      real(rp), public   :: eta_b(45), f(45), f_prim(45)

      character(len=64), public :: tag
      logical :: db_clustered = .false.
      character(512), public   :: fileControlName                            ! Input: path to the file that contains the points defining the rectangle controls
      integer(4), public         :: nRectangleControl                          ! Number of rectangle control
      real(rp), public :: amplitudeActuation, periodEpisode, previousActuationTime, periodActuation, frequencyActuation, timeBeginActuation ! parameters of the actuation
      integer(4), public :: tw_write_interval, n_pseudo_envs
      real(rp), public :: time_tw
      real(8), allocatable, public :: lx_recirculation_avg(:)

   contains
      procedure, public :: fillBCTypes           => BLMARLFlowSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => BLMARLFlowSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => BLMARLFlowSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer => BLMARLFlowSolverIncomp_initialBuffer
      procedure, public :: fillBlasius => BLMARLFlowSolverIncomp_fillBlasius
      procedure, public :: afterDt => BLMARLFlowSolverIncomp_afterDt
      procedure, public :: beforeTimeIteration => BLMARLFlowSolverIncomp_beforeTimeIteration
      procedure, public :: afterTimeIteration => BLMARLFlowSolverIncomp_afterTimeIteration
      procedure, public :: computeTauW => BLMARLFlowSolverIncomp_computeTauW
#if ACTUATION
      procedure, public :: getControlNodes => BLMARLFlowSolverIncomp_getControlNodes
      procedure, public :: readControlRectangles => BLMARLFlowSolverIncomp_readControlRectangles
#ifdef SMARTREDIS
      procedure, public :: initSmartRedis  => BLMARLFlowSolverIncomp_initSmartRedis
      procedure, public :: smoothControlFunction => BLMARLFlowSolverIncomp_smoothControlFunction
#endif
#endif
   end type BLMARLFlowSolverIncomp
contains

#if ACTUATION
   subroutine BLMARLFlowSolverIncomp_readControlRectangles(this)
      ! This subroutine reads the file that contains the two points defining a rectanle paralel to the X-Z axis. Several rectangles
      ! can be introduced. In this rectangles is where control will be applied
      class(BLMARLFlowSolverIncomp), intent(inout) :: this

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

   end subroutine BLMARLFlowSolverIncomp_readControlRectangles

   subroutine BLMARLFlowSolverIncomp_getControlNodes(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this

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
                     exit
                  endif
               end do
            end if
         end if
      end do
      !$acc end parallel loop
   end subroutine BLMARLFlowSolverIncomp_getControlNodes

#ifdef SMARTREDIS
   subroutine BLMARLFlowSolverIncomp_initSmartRedis(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this

      open(unit=443,file="./output_"//trim(adjustl(this%tag))//"/"//"control_fortran_raw.txt",status='replace')
      open(unit=445,file="./output_"//trim(adjustl(this%tag))//"/"//"control_tw.txt",status='replace')
      call init_smartredis(client, this%nwitPar, this%nRectangleControl, this%n_pseudo_envs, trim(adjustl(this%tag)), this%db_clustered)
      if (mpi_rank .eq. 0) write(111, *) "SmartRedis initialised"
      call write_step_type(client, 1, "ensemble_"//trim(adjustl(this%tag))//".step_type")
   end subroutine BLMARLFlowSolverIncomp_initSmartRedis
#endif
#endif

   subroutine BLMARLFlowSolverIncomp_beforeTimeIteration(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
      integer(4)                 :: iboun,bcode,ipbou,iBoundNode,iNodeL

     !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(bouCodesNodesPar(iNodeL) .lt. max_num_bou_codes) then
            bcode = bouCodesNodesPar(iNodeL) ! Boundary element code
            if (bcode == bc_type_far_field_SB) then ! top far field for BL
               u_buffer(iNodeL,2) =  0.470226_rp*(306.640625_rp-coordPar(iNodeL,1))/110.485435_rp*exp(0.95_rp-((306.640625_rp &
                                 -coordPar(iNodeL,1))/110.485435_rp)**2_rp)
               u(iNodeL,2,2) = u_buffer(iNodeL,2)
            end if
         end if
      end do
      !$acc end parallel loop

#if ACTUATION
      open(unit=444,file="./output_"//trim(adjustl(this%tag))//"/"//"control_fortran_smooth.txt",status='replace')
      call this%getControlNodes()
#ifdef SMARTREDIS
      if (mpi_rank .eq. 0) open(unit=447,file="./output_"//trim(adjustl(this%tag))//"/"//"lx_avg.txt",status='replace')
      call this%initSmartRedis()
#endif
#endif

      if (mpi_rank .eq. 0) open(unit=446,file="./output_"//trim(adjustl(this%tag))//"/"//"lx.txt",status='replace')

      allocate(this%lx_recirculation_avg(this%n_pseudo_envs))
      !$acc enter data create(this%lx_recirculation_avg(:))
      this%lx_recirculation_avg(:) = 0.0d0
      !$acc update device(this%lx_recirculation_avg(:))
      this%time_tw = 0.0
   end subroutine BLMARLFlowSolverIncomp_beforeTimeIteration

   subroutine BLMARLFlowSolverIncomp_afterTimeIteration(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
#if ACTUATION
      if (mpi_rank .eq. 0) close(444)
#ifdef SMARTREDIS
      call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
      if (mpi_rank .eq. 0) close(443)
      if (mpi_rank .eq. 0) close(445)
      if (mpi_rank .eq. 0) close(447)
      call end_smartredis(client)
#endif
#endif
      if (mpi_rank .eq. 0) close(446)
   end subroutine BLMARLFlowSolverIncomp_afterTimeIteration

   subroutine BLMARLFlowSolverIncomp_afterDt(this,istep)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3
      integer(4) :: k,j,ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime,iNodeL,iNodeL2,bcode,isoII, isoJJ, isoKK,type_ijk
      real(rp) :: dy,fx1,fx2,xp
      real(rp) :: mul , yc
      real(rp)  :: gradIsoV(ndime),gradIsoU(ndime)
      real(rp)  :: gradV(ndime),vl(nnode),fact,targ,gradU(ndime),ul(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
      real(rp) :: yp,eta_y,f_y,f_prim_y
      real(rp) :: eliti, ave
      real(8) :: lx_recirculation(this%n_pseudo_envs)
      integer(4) :: i
#if ACTUATION
#ifdef SMARTREDIS
      real(rp) :: action_global_instant(action_global_size)
#else
      real(rp) :: action_classic
#endif
#endif

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

      ! wall shear stress
      call this%computeTauW(lx_recirculation)
      if (mod(istep, this%tw_write_interval) .eq. 0 .and. mpi_rank .eq. 0) then
         write(446,'(*(ES12.4,:,","))') this%time, lx_recirculation
         call flush(446)
      end if

#if ACTUATION
#ifdef SMARTREDIS
      if (step_type_mod .eq. 1) then
         call write_step_type(client, 2, "ensemble_"//trim(adjustl(this%tag))//".step_type")
      end if
#endif
#endif

#if ACTUATION
      if (this%time .gt. this%timeBeginActuation) then
#ifdef SMARTREDIS
         ! average wall shear stress
         ave = this%dt/(this%time_tw+this%dt)
         eliti = this%time_tw/(this%time_tw+this%dt)
         this%time_tw = this%time_tw+this%dt
         ! !$acc kernels
         this%lx_recirculation_avg(:) = ave * lx_recirculation(:) + eliti * this%lx_recirculation_avg(:)
         ! !$acc end kernels
         if (mod(istep, this%tw_write_interval) .eq. 0 .and. mpi_rank .eq. 0) then
            write(447,'(*(ES12.4,:,","))') this%time, this%lx_recirculation_avg
            call flush(447)
         end if

         ! check if a new action is needed
         if (this%time - this%previousActuationTime .ge. this%periodActuation) then
            if (mpi_rank .eq. 0) write(111, *) "Performing SmartRedis comms"

            ! reset wall shear stress averaging time
            this%time_tw = 0.0_rp

            ! save old action values and time - useful for interpolating to new action_global values
            !$acc kernels
            action_global_previous(:) = action_global(:)
            this%previousActuationTime = this%previousActuationTime + this%periodActuation
            !$acc end kernels

            call this%update_witness(istep, 1) ! manual update of witness points
            call write_state(client, buffwit(:, 1, 1), "ensemble_"//trim(adjustl(this%tag))//".state") ! the streamwise velocity u
            call write_reward(client, this%lx_recirculation_avg, "ensemble_"//trim(adjustl(this%tag))//".reward") ! the streamwise component tw_x
            if (mpi_rank .eq. 0) write(445,'(*(ES12.4,:,","))') this%time, this%lx_recirculation_avg
            if (mpi_rank .eq. 0) call flush(445)

            call read_action(client, "ensemble_"//trim(adjustl(this%tag))//".action") ! modifies action_global (the target control values)
            if (mpi_rank .eq. 0) write(443,'(*(ES12.4,:,","))') this%time+this%periodActuation, action_global
            if (mpi_rank .eq. 0) call flush(443)

            ! if the next time that we require actuation value is the last one, write now step_type=0 into database
            if (this%time + 2.0 * this%periodActuation .gt. this%maxPhysTime) then
               call write_step_type(client, 0, "ensemble_"//trim(adjustl(this%tag))//".step_type")
            end if
         end if

         ! apply actions to control surfaces
         call this%smoothControlFunction(action_global_instant)
         if (mpi_rank .eq. 0) write(444,'(*(ES12.4,:,","))') this%time, action_global_instant
         call flush(444)

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then
                  u_buffer(iNodeL,1) = 0.0_rp
                  u_buffer(iNodeL,2) = action_global_instant(actionMask(iNodeL))
                  u_buffer(iNodeL,3) = 0.0_rp
            end if
         end do
         !$acc end parallel loop
#else
         action_classic = this%amplitudeActuation*sin(2.0_rp*v_pi*this%frequencyActuation*this%time)
         if (mpi_rank .eq. 0) write(444,'(*(ES12.4,:,","))') this%time, action_classic
         call flush(444)
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then
               u_buffer(iNodeL,1) = 0.0_rp
               u_buffer(iNodeL,2) = action_classic
               u_buffer(iNodeL,3) = 0.0_rp
            end if
         end do
         !$acc end parallel loop
#endif
      end if
#endif
   end subroutine BLMARLFlowSolverIncomp_afterDt

#if ACTUATION
#ifdef SMARTREDIS
      subroutine BLMARLFlowSolverIncomp_smoothControlFunction(this, action_global_instant)
         class(BLMARLFlowSolverIncomp), intent(inout) :: this
         real(rp), intent(inout) :: action_global_instant(action_global_size)

         real(rp) :: f1, f2, f3

         f1 = exp(-1.0_rp / ((this%time - this%previousActuationTime) / this%periodActuation))
         f2 = exp(-1.0_rp / (1.0_rp - (this%time - this%previousActuationTime) / this%periodActuation))
         f3 = f1 / (f1 + f2)
         !$acc kernels
         action_global_instant(:) = action_global_previous(:) + f3 * (action_global(:) - action_global_previous(:))
         !$acc end kernels
      end subroutine BLMARLFlowSolverIncomp_smoothControlFunction
#endif
#endif

   subroutine BLMARLFlowSolverIncomp_computeTauW(this, lx_recirculation)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
      real(8), intent(out) :: lx_recirculation(this%n_pseudo_envs)
      real(8) :: lx_r
      integer(4) :: surfCode

      !$acc loop seq
      do surfCode=1, this%n_pseudo_envs
         call twInfo(numElemsRankPar, numNodesRankPar, numBoundsRankPar, surfCode, 1, connecParWork, boundPar, &
            point2elem, bouCodesPar, boundNormalPar, invAtoIJK, gmshAtoI, gmshAtoJ, gmshAtoK, wgp_b, dlxigp_ip, He, coordPar, &
            mu_fluid, mu_e, mu_sgs, rho(:,2), u(:,:,2), lx_r)
         lx_recirculation(surfCode) = lx_r / this%Lz
      end do
   end subroutine BLMARLFlowSolverIncomp_computeTauW

   subroutine BLMARLFlowSolverIncomp_fill_BC_Types(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this

#if ACTUATION
      ! wall + actuation
      bouCodes2BCType(1) = bc_type_unsteady_inlet ! wall 1
      bouCodes2BCType(2) = bc_type_unsteady_inlet ! wall 2
      bouCodes2BCType(3) = bc_type_unsteady_inlet ! wall 3
      bouCodes2BCType(4) = bc_type_unsteady_inlet ! wall 4
      bouCodes2BCType(5) = bc_type_unsteady_inlet ! wall 5
#else
      ! wall
      bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! wall 1
      bouCodes2BCType(2) = bc_type_non_slip_adiabatic ! wall 2
      bouCodes2BCType(3) = bc_type_non_slip_adiabatic ! wall 3
      bouCodes2BCType(4) = bc_type_non_slip_adiabatic ! wall 4
      bouCodes2BCType(5) = bc_type_non_slip_adiabatic ! wall 5
#endif
      bouCodes2BCType(6) = bc_type_far_field_SB ! upper part of the domain
      bouCodes2BCType(7) = bc_type_far_field ! inlet part of the domain
      bouCodes2BCType(8) = bc_type_outlet_incomp ! outlet part of the domain
      !$acc update device(bouCodes2BCType(:))

   end subroutine BLMARLFlowSolverIncomp_fill_BC_Types

   subroutine BLMARLFlowSolverIncomp_fillBlasius(this)
     class(BLMARLFlowSolverIncomp), intent(inout) :: this

     this%eta_b(1) = real(0.0E+00,rp)
     this%eta_b(2) = real(2.0E-01,rp)
     this%eta_b(3) = real(4.0E-01,rp)
     this%eta_b(4) = real(6.0E-01,rp)
     this%eta_b(5) = real(8.0E-01,rp)
     this%eta_b(6) = real(1.0E+00,rp)
     this%eta_b(7) = real(1.2E+00,rp)
     this%eta_b(8) = real(1.4E+00,rp)
     this%eta_b(9) = real(1.6E+00,rp)
     this%eta_b(10) = real(1.8E+00,rp)
     this%eta_b(11) = real(2.0E+00,rp)
     this%eta_b(12) = real(2.2E+00,rp)
     this%eta_b(13) = real(2.4E+00,rp)
     this%eta_b(14) = real(2.6E+00,rp)
     this%eta_b(15) = real(2.8E+00,rp)
     this%eta_b(16) = real(3.0E+00,rp)
     this%eta_b(17) = real(3.2E+00,rp)
     this%eta_b(18) = real(3.4E+00,rp)
     this%eta_b(19) = real(3.6E+00,rp)
     this%eta_b(20) = real(3.8E+00,rp)
     this%eta_b(21) = real(4.0E+00,rp)
     this%eta_b(22) = real(4.2E+00,rp)
     this%eta_b(23) = real(4.4E+00,rp)
     this%eta_b(24) = real(4.6E+00,rp)
     this%eta_b(25) = real(4.8E+00,rp)
     this%eta_b(26) = real(5.0E+00,rp)
     this%eta_b(27) = real(5.2E+00,rp)
     this%eta_b(28) = real(5.4E+00,rp)
     this%eta_b(29) = real(5.6E+00,rp)
     this%eta_b(30) = real(5.8E+00,rp)
     this%eta_b(31) = real(6.0E+00,rp)
     this%eta_b(32) = real(6.2E+00,rp)
     this%eta_b(33) = real(6.4E+00,rp)
     this%eta_b(34) = real(6.6E+00,rp)
     this%eta_b(35) = real(6.8E+00,rp)
     this%eta_b(36) = real(7.0E+00,rp)
     this%eta_b(37) = real(7.2E+00,rp)
     this%eta_b(38) = real(7.4E+00,rp)
     this%eta_b(39) = real(7.6E+00,rp)
     this%eta_b(40) = real(7.8E+00,rp)
     this%eta_b(41) = real(8.0E+00,rp)
     this%eta_b(42) = real(8.2E+00,rp)
     this%eta_b(43) = real(8.4E+00,rp)
     this%eta_b(44) = real(8.6E+00,rp)
     this%eta_b(45) = real(8.8E+00,rp)

     this%f(1)  = real(0.000000000E+00,rp)
     this%f(2)  = real(6.640999715E-03,rp)
     this%f(3)  = real(2.655988402E-02,rp)
     this%f(4)  = real(5.973463750E-02,rp)
     this%f(5)  = real(1.061082208E-01,rp)
     this%f(6)  = real(1.655717258E-01,rp)
     this%f(7)  = real(2.379487173E-01,rp)
     this%f(8)  = real(3.229815738E-01,rp)
     this%f(9)  = real(4.203207655E-01,rp)
     this%f(10) = real(5.295180377E-01,rp)
     this%f(11) = real(6.500243699E-01,rp)
     this%f(12) = real(7.811933370E-01,rp)
     this%f(13) = real(9.222901256E-01,rp)
     this%f(14) = real(1.072505977E+00,rp)
     this%f(15) = real(1.230977302E+00,rp)
     this%f(16) = real(1.396808231E+00,rp)
     this%f(17) = real(1.569094960E+00,rp)
     this%f(18) = real(1.746950094E+00,rp)
     this%f(19) = real(1.929525170E+00,rp)
     this%f(20) = real(2.116029817E+00,rp)
     this%f(21) = real(2.305746418E+00,rp)
     this%f(22) = real(2.498039663E+00,rp)
     this%f(23) = real(2.692360938E+00,rp)
     this%f(24) = real(2.888247990E+00,rp)
     this%f(25) = real(3.085320655E+00,rp)
     this%f(26) = real(3.283273665E+00,rp)
     this%f(27) = real(3.481867612E+00,rp)
     this%f(28) = real(3.680919063E+00,rp)
     this%f(29) = real(3.880290678E+00,rp)
     this%f(30) = real(4.079881939E+00,rp)
     this%f(31) = real(4.279620923E+00,rp)
     this%f(32) = real(4.479457297E+00,rp)
     this%f(33) = real(4.679356615E+00,rp)
     this%f(34) = real(4.879295811E+00,rp)
     this%f(35) = real(5.079259772E+00,rp)
     this%f(36) = real(5.279238811E+00,rp)
     this%f(37) = real(5.479226847E+00,rp)
     this%f(38) = real(5.679220147E+00,rp)
     this%f(39) = real(5.879216466E+00,rp)
     this%f(40) = real(6.079214481E+00,rp)
     this%f(41) = real(6.279213431E+00,rp)
     this%f(42) = real(6.479212887E+00,rp)
     this%f(43) = real(6.679212609E+00,rp)
     this%f(44) = real(6.879212471E+00,rp)
     this%f(45) = real(7.079212403E+00,rp)

     this%f_prim(1)  = real(0.000000000E+00,rp)
     this%f_prim(2)  = real(6.640779210E-02,rp)
     this%f_prim(3)  = real(1.327641608E-01,rp)
     this%f_prim(4)  = real(1.989372524E-01,rp)
     this%f_prim(5)  = real(2.647091387E-01,rp)
     this%f_prim(6)  = real(3.297800312E-01,rp)
     this%f_prim(7)  = real(3.937761044E-01,rp)
     this%f_prim(8)  = real(4.562617647E-01,rp)
     this%f_prim(9)  = real(5.167567844E-01,rp)
     this%f_prim(10) = real(5.747581439E-01,rp)
     this%f_prim(11) = real(6.297657365E-01,rp)
     this%f_prim(12) = real(6.813103772E-01,rp)
     this%f_prim(13) = real(7.289819351E-01,rp)
     this%f_prim(14) = real(7.724550211E-01,rp)
     this%f_prim(15) = real(8.115096232E-01,rp)
     this%f_prim(16) = real(8.460444437E-01,rp)
     this%f_prim(17) = real(8.760814552E-01,rp)
     this%f_prim(18) = real(9.017612214E-01,rp)
     this%f_prim(19) = real(9.233296659E-01,rp)
     this%f_prim(20) = real(9.411179967E-01,rp)
     this%f_prim(21) = real(9.555182298E-01,rp)
     this%f_prim(22) = real(9.669570738E-01,rp)
     this%f_prim(23) = real(9.758708321E-01,rp)
     this%f_prim(24) = real(9.826835008E-01,rp)
     this%f_prim(25) = real(9.877895262E-01,rp)
     this%f_prim(26) = real(9.915419002E-01,rp)
     this%f_prim(27) = real(9.942455354E-01,rp)
     this%f_prim(28) = real(9.961553040E-01,rp)
     this%f_prim(29) = real(9.974777682E-01,rp)
     this%f_prim(30) = real(9.983754937E-01,rp)
     this%f_prim(31) = real(9.989728724E-01,rp)
     this%f_prim(32) = real(9.993625417E-01,rp)
     this%f_prim(33) = real(9.996117017E-01,rp)
     this%f_prim(34) = real(9.997678702E-01,rp)
     this%f_prim(35) = real(9.998638190E-01,rp)
     this%f_prim(36) = real(9.999216041E-01,rp)
     this%f_prim(37) = real(9.999557173E-01,rp)
     this%f_prim(38) = real(9.999754577E-01,rp)
     this%f_prim(39) = real(9.999866551E-01,rp)
     this%f_prim(40) = real(9.999928812E-01,rp)
     this%f_prim(41) = real(9.999962745E-01,rp)
     this%f_prim(42) = real(9.999980875E-01,rp)
     this%f_prim(43) = real(9.999990369E-01,rp)
     this%f_prim(44) = real(9.999995242E-01,rp)
     this%f_prim(45) = real(9.999997695E-01,rp)

  end subroutine BLMARLFlowSolverIncomp_fillBlasius

   subroutine BLMARLFlowSolverIncomp_initializeParameters(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
      real(rp) :: mur
      integer :: num_args, equal_pos, iarg, ierr
      character(len=64) :: arg, output_dir
      character(len=8) :: restart_step_str="", periodEpisode_str="1.0", frequencyActuation_str="1.0", timeBeginActuation_str="0.0", db_clustered_str="0"
      logical :: output_dir_exists

      ! get command line args, ie: mpirun -n 4 sod2d --restart_step=1 --t_episode=100.0
      num_args = command_argument_count()
      do iarg = 1, num_args
         call get_command_argument(iarg, arg)
         equal_pos = scan(adjustl(trim(arg)), "=")
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
      if (db_clustered_str == "" .or. db_clustered_str == "0") this%db_clustered = .false.
      if (restart_step_str == "" .or. restart_step_str == "0") then
         this%loadRestartFile = .false.
         this%restartFile_to_load = 0
      else
         this%loadRestartFile = .true.
         read(restart_step_str, *) this%restartFile_to_load ! 1: baseline, 2: last episode
      end if
      if (periodEpisode_str == "") then
         this%periodEpisode = 10000000.0_rp
      else
         read(periodEpisode_str,*,iostat=ierr) this%periodEpisode
      end if
      read(frequencyActuation_str,*,iostat=ierr) this%frequencyActuation
      this%periodActuation = 1.0_rp / this%frequencyActuation
      read(periodEpisode_str,*,iostat=ierr) this%periodEpisode
      read(timeBeginActuation_str,*,iostat=ierr) this%timeBeginActuation

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
         write(*,*) "RL simulation params:"
         write(*,*) "--tag: ", adjustl(trim(this%tag))
         write(*,*) "--restart_step: ", this%restartFile_to_load
         write(*,*) "--f_action: ", this%frequencyActuation
         write(*,*) "--t_episode: ", this%periodEpisode
         write(*,*) "--t_begin_control: ", this%timeBeginActuation
         write(*,*) "--db_clustered: ", db_clustered_str
      end if

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 9000000
      this%maxPhysTime =  this%timeBeginActuation + this%periodEpisode

      this%save_logFile_first = 1
      this%save_logFile_step  = 250

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 15000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 15000
      this%continue_oldLogs = .false.

      this%initial_avgTime = this%timeBeginActuation
      this%saveAvgFile = .true.
      this%loadAvgFile = .false.
      ! -------- Average results file -------------
      this%save_avgScalarField_rho     = .false.
      this%save_avgScalarField_pr      = .true.
      this%save_avgScalarField_mueff   = .true.
      this%save_avgVectorField_vel     = .true.
      this%save_avgVectorField_ve2     = .false.
      this%save_avgVectorField_vex     = .false.
      this%save_avgVectorField_vtw     = .false.
      !----------------------------------------------

      ! wall shear stress output
      this%n_pseudo_envs = 3
      this%Lz = 125.0_rp / this%n_pseudo_envs
      this%tw_write_interval = 1

      ! numerical params
      flag_les = 1

      maxIter = 30
      tol = 1e-2

      this%cfl_conv = 0.95_rp
      this%cfl_diff = 0.95_rp
      ! flag_use_constant_dt = 1
      ! this%dt = 0.065_rp

      this%d0   = 1.0_rp
      this%U0   = 1.0_rp
      this%rho0 = 1.0_rp
      this%Red0  = 450.0_rp

      this%mu    = this%rho0*this%d0*this%U0/this%Red0

      incomp_viscosity = this%mu
      flag_mu_factor = 1.0_rp

      nscbc_u_inf = this%U0
      nscbc_p_inf = 0.0_rp
      nscbc_rho_inf = this%rho0

      flag_fs_fix_pressure = .false.

      flag_buffer_on = .true.
      ! x outlet
      flag_buffer_on_east = .true.
      flag_buffer_e_min = 950.0_rp
      flag_buffer_e_size = 50.0_rp
      !! x inlet
      flag_buffer_on_west = .true.
      flag_buffer_w_min = -50.0_rp
      flag_buffer_w_size = 50.0_rp

      ! witness points
      this%have_witness          = .true.
      this%nwit                  = 240
      this%witness_inp_file_name = "witness.txt"
      this%witness_h5_file_name  = "./output_"//trim(adjustl(this%tag))//"/resultwit.h5"
      this%leapwit               = 5 ! (update witness ever n dts) | in this class, we update the witness points manually
      this%leapwitsave           = 100 ! how many dts are stored in buffer
      this%wit_save              = .true. ! save witness or not
      this%wit_save_u_i          = .true.
      this%wit_save_pr           = .true.
      this%wit_save_rho          = .false.
      this%continue_witness      = .false.

#if ACTUATION
      write(this%fileControlName,*) "rectangleControl.txt"
#ifndef SMARTREDIS
      this%amplitudeActuation = 0.05
      this%frequencyActuation = 0.0030_rp
#endif
#endif
   end subroutine BLMARLFlowSolverIncomp_initializeParameters

   subroutine BLMARLFlowSolverIncomp_initialBuffer(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
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
            if(eta_y<this%eta_b(k)) then
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
            f_y      = this%f(j-1)      + (this%f(j)-this%f(j-1))*(eta_y-this%eta_b(j-1))/(this%eta_b(j)-this%eta_b(j-1))
            f_prim_y = this%f_prim(j-1) + (this%f_prim(j)-this%f_prim(j-1))*(eta_y-this%eta_b(j-1))/(this%eta_b(j)-this%eta_b(j-1))

            u_buffer(iNodeL,1) = f_prim_y
            u_buffer(iNodeL,2) = 0.5_rp*sqrt(1.0/(450.0_rp*450.0_rp))*(eta_y*f_prim_y-f_y)
            u_buffer(iNodeL,3) = 0.0_rp
         end if
      end do
      !$acc end parallel loop

   end subroutine BLMARLFlowSolverIncomp_initialBuffer

   subroutine BLMARLFlowSolverIncomp_evalInitialConditions(this)
      class(BLMARLFlowSolverIncomp), intent(inout) :: this
      integer(4) :: iNodeL, idime, j,k
      real(rp) :: yp,eta_y,f_y,f_prim_y
      integer(4)   :: iLine,iNodeGSrl,auxCnt
      integer(4)                 :: iboun,bcode,ipbou,iBoundNode

      call this%fillBlasius()

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2)
         eta_y = yp !with our normalisation is sqrt(U/(nu x ) is actually 1 for the inlet)
         j = 45
         !$acc loop seq
         label1:do k=1,45
            if(eta_y<this%eta_b(k)) then
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
            f_y      = this%f(j-1)      + (this%f(j)-this%f(j-1))*(eta_y-this%eta_b(j-1))/(this%eta_b(j)-this%eta_b(j-1))
            f_prim_y = this%f_prim(j-1) + (this%f_prim(j)-this%f_prim(j-1))*(eta_y-this%eta_b(j-1))/(this%eta_b(j)-this%eta_b(j-1))

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
   end subroutine BLMARLFlowSolverIncomp_evalInitialConditions

end module BLMARLFlowSolverIncomp_mod
