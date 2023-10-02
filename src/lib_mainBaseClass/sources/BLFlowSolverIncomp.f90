#define ACTUATION 1

module BLFlowSolverIncomp_mod
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
   use CFDSolverPeriodicWithBoundariesIncomp_mod
   implicit none
   private

   real(rp), allocatable, dimension(:,:)  :: rectangleControl       ! Two point coordinates that define each rectangle
   integer(4),  allocatable, dimension(:) :: actionMask             ! Mask that contains whether a point is in a rectangle control or not 

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BLFlowSolverIncomp

      real(rp) , public  ::   d0, U0, rho0, Red0, Re, mu
      real(rp), public   :: eta_b(45), f(45), f_prim(45)
      integer(4), public       :: countPar                                   ! Number of points in a rectangle of control per partition
      character(512), public   :: fileControlName                            ! Input: path to the file that contains the points defining the      
      integer(4), public         :: nRectangleControl                          ! Number of rectangle control
      real(rp), public         :: amplitudeActuation, frequencyActuation , timeBeginActuation    ! Parameters of the actuation
   contains
      procedure, public :: fillBCTypes           => BLFlowSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => BLFlowSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => BLFlowSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer => BLFlowSolverIncomp_initialBuffer
      procedure, public :: fillBlasius => BLFlowSolverIncomp_fillBlasius
      procedure, public :: afterDt => BLFlowSolverIncomp_afterDt
      procedure, public :: beforeTimeIteration => BLFlowSolverIncomp_beforeTimeIteration
      procedure, public :: getControlNodes => BLFlowSolverIncomp_getControlNodes
      procedure, public :: readControlRectangles => BLFlowSolverIncomp_readControlRectangles      
   end type BLFlowSolverIncomp
contains

   subroutine BLFlowSolverIncomp_readControlRectangles(this)
      ! This subroutine reads the file that contains the two points defining a rectanle paralel to the X-Z axis. Several rectangles
      ! can be introduced. In this rectangles is where control will be applied
      class(BLFlowSolverIncomp), intent(inout) :: this

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

   end subroutine BLFlowSolverIncomp_readControlRectangles

   subroutine BLFlowSolverIncomp_getControlNodes(this)
      class(BLFlowSolverIncomp), intent(inout) :: this

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

   end subroutine BLFlowSolverIncomp_getControlNodes

   subroutine BLFlowSolverIncomp_beforeTimeIteration(this)
      class(BLFlowSolverIncomp), intent(inout) :: this
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

#if (ACTUATION)
      call this%getControlNodes()
#endif
   end subroutine BLFlowSolverIncomp_beforeTimeIteration

   subroutine BLFlowSolverIncomp_afterDt(this,istep)
      class(BLFlowSolverIncomp), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3
      integer(4) :: k,j,ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime,iNodeL,iNodeL2,bcode,isoII, isoJJ, isoKK,type_ijk
      real(rp) :: dy,fx1,fx2,xp
      real(rp) :: mul , yc
      real(rp)  :: gradIsoV(ndime),gradIsoU(ndime)
      real(rp)  :: gradV(ndime),vl(nnode),fact,targ,gradU(ndime),ul(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
      real(rp) :: yp,eta_y,f_y,f_prim_y


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

#if (ACTUATION)
      if (this%time .gt. this%timeBeginActuation) then
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then
                  u_buffer(iNodeL,1) = 0.0_rp
                  u_buffer(iNodeL,2) = this%amplitudeActuation*sin(2.0_rp*v_pi*this%frequencyActuation*this%time)
                  u_buffer(iNodeL,3) = 0.0_rp
            end if
         end do
         !$acc end parallel loop
      end if
#endif

   end subroutine BLFlowSolverIncomp_afterDt

   subroutine BLFlowSolverIncomp_fill_BC_Types(this)
      class(BLFlowSolverIncomp), intent(inout) :: this

#if (ACTUATION)
      bouCodes2BCType(1) = bc_type_unsteady_inlet ! wall + actuation
#else
      bouCodes2BCType(1) = bc_type_non_slip_adiabatic ! wall
#endif      
      bouCodes2BCType(2) = bc_type_far_field_SB       ! Upper part of the domain
      bouCodes2BCType(3) = bc_type_far_field      ! inlet part of the domain
      bouCodes2BCType(4) = bc_type_outlet_incomp  ! outlet part of the domain
      !$acc update device(bouCodes2BCType(:))

   end subroutine BLFlowSolverIncomp_fill_BC_Types

   subroutine BLFlowSolverIncomp_fillBlasius(this)
      class(BLFlowSolverIncomp), intent(inout) :: this

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

  end subroutine BLFlowSolverIncomp_fillBlasius

   subroutine BLFlowSolverIncomp_initializeParameters(this)
      class(BLFlowSolverIncomp), intent(inout) :: this
      real(rp) :: mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "bl"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 9000000 

      this%save_logFile_first = 1 
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 10000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%loadRestartFile = .true.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .false.

      this%initial_avgTime = 2000.0_rp
      this%saveAvgFile = .true.
      this%loadAvgFile = .false.
      !----------------------------------------------

      ! numerical params
      flag_les = 1
      
      maxIter = 20
      tol = 1e-2

      this%cfl_conv = 0.95_rp
      this%cfl_diff = 0.95_rp
      !flag_use_constant_dt = 1
      !this%dt = 0.1_rp

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

#if (ACTUATION)
      write(this%fileControlName ,*) "rectangleControl.dat"
      this%amplitudeActuation = 0.2
      this%frequencyActuation = 0.0025_rp
      this%timeBeginActuation = 0.0_rp
#endif

   end subroutine BLFlowSolverIncomp_initializeParameters

   subroutine BLFlowSolverIncomp_initialBuffer(this)
      class(BLFlowSolverIncomp), intent(inout) :: this
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

   end subroutine BLFlowSolverIncomp_initialBuffer

   subroutine BLFlowSolverIncomp_evalInitialConditions(this)
      class(BLFlowSolverIncomp), intent(inout) :: this
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
   end subroutine BLFlowSolverIncomp_evalInitialConditions

end module BLFlowSolverIncomp_mod
