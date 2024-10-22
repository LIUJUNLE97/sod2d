module BluffBodySolverIncompAFC_mod
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

   real(rp), allocatable, dimension(:,:)  :: rectangleControl
   integer(4),  allocatable, dimension(:) :: actionMask

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BluffBodySolverIncompAFC

      real(rp) , public  :: vo, delta, rho0, Re, aoa

      integer(4), public :: nRectangleControl, flag_MassConservationStrategy
      character(512), public :: fileControlName
      real(rp), public :: timeBeginActuation, amplitudeActuation, frequencyActuation, alphaActuation, spanFreqActuation

   contains
      procedure, public :: readControlRectangles => BluffBodySolverIncompAFC_readControlRectangles
      procedure, public :: getControlNodes       => BluffBodySolverIncompAFC_getControlNodes
      procedure, public :: beforeTimeIteration => BluffBodySolverIncompAFC_beforeTimeIteration
      procedure, public :: afterTimeIteration => BluffBodySolverIncompAFC_afterTimeIteration
      procedure, public :: afterDt => BluffBodySolverIncompAFC_afterDt
      procedure, public :: fillBCTypes           => BluffBodySolverIncompAFC_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolverIncompAFC_initializeParameters
      procedure, public :: evalInitialConditions => BluffBodySolverIncompAFC_evalInitialConditions
      procedure, public :: initialBuffer => BluffBodySolverIncompAFC_initialBuffer
   end type BluffBodySolverIncompAFC
contains

   subroutine BluffBodySolverIncompAFC_readControlRectangles(this)

      class(BluffBodySolverIncompAFC), intent(inout) :: this

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

   end subroutine BluffBodySolverIncompAFC_readControlRectangles

   subroutine BluffBodySolverIncompAFC_getControlNodes(this)

      class(BluffBodySolverIncompAFC), intent(inout) :: this

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

                  ! Mesh node coordinates
                  xPoint = coordPar(iNodeL,1)
                  yPoint = coordPar(iNodeL,2)
                  zPoint = coordPar(iNodeL,3)

                  ! Control surfaces xMin, yMin, zMin, xMax, yMax, zMax coordinatess
                  x1RectangleControl = rectangleControl(1,2*iRectangleControl-1)
                  y1RectangleControl = rectangleControl(2,2*iRectangleControl-1)
                  z1RectangleControl = rectangleControl(3,2*iRectangleControl-1)

                  x2RectangleControl = rectangleControl(1,2*iRectangleControl)
                  y2RectangleControl = rectangleControl(2,2*iRectangleControl)
                  z2RectangleControl = rectangleControl(3,2*iRectangleControl)

                  ! Check if the mesh point is within the control surface limits
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

   end subroutine BluffBodySolverIncompAFC_getControlNodes

   subroutine BluffBodySolverIncompAFC_beforeTimeIteration(this)
      class(BluffBodySolverIncompAFC), intent(inout) :: this
      integer(4) :: bcode, iNodeL

      ! Open file to save the instantaneous action
      if (mpi_rank .eq. 0) open(unit=444,file="instantaneous_action.txt",status='replace')

      ! Obtain the mask of the control nodes (>0)
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

   end subroutine BluffBodySolverIncompAFC_beforeTimeIteration

   subroutine BluffBodySolverIncompAFC_afterTimeIteration(this)
      class(BluffBodySolverIncompAFC), intent(inout) :: this

      if (mpi_rank .eq. 0) close(444)

   end subroutine BluffBodySolverIncompAFC_afterTimeIteration

   subroutine BluffBodySolverIncompAFC_afterDt(this,istep)
      class(BluffBodySolverIncompAFC), intent(inout) :: this
      integer(4), intent(in) :: istep

      integer(4) :: iNodeL, invert
      real(rp)   :: zPoint
      real(rp) :: action_classic

      real(rp)   :: xPoint, yPoint, theta ! BORRAR

      !  Start actuating when t > timeBeginActuation
      if (this%time .gt. this%timeBeginActuation) then

         ! Apply CLASSIC ACTUATION 
         action_classic = this%amplitudeActuation*sin(2.0_rp*v_pi*this%frequencyActuation*(this%time - this%timeBeginActuation))

         if(this%save_logFile_next==istep) then
            !if (mpi_rank .eq. 0) write(444,'(*(ES16.6,:,","))') this%time, action_classic
            if (mpi_rank .eq. 0) then
               write(444,"(I8,1X,A)", ADVANCE="NO") istep, ","
               write(444,50) this%time, ",", action_classic, ""
               50 format(16(1X,F16.8,1X,A))
            end if
            call flush(444)
         end if

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            if (actionMask(iNodeL) .gt. 0) then

               ! Apply mass conservation strategy: Consecutive actuators are paired so that the opposite mass flow rate is forced (when the flag is activated)
               ! EVEN NUMBER OF ACTUATORS REQUIRED !
               if (MOD(actionMask(iNodeL),2) .eq. 0 .and. this%flag_MassConservationStrategy .eq. 1) then
                  invert = -1
               else 
                  invert = 1
               end if

               if  (this%spanFreqActuation .eq. 0) then

                  u_buffer(iNodeL,1) = action_classic*sin(this%alphaActuation*v_pi/180_rp) * invert
                  u_buffer(iNodeL,2) = action_classic*cos(this%alphaActuation*v_pi/180_rp) * invert
                  u_buffer(iNodeL,3) = 0.0_rp * invert

               else

                  zPoint = coordPar(iNodeL,3)

                  u_buffer(iNodeL,1) = action_classic*sin(2.0_rp*v_pi*this%spanFreqActuation*zPoint)*sin(this%alphaActuation*v_pi/180_rp)
                  u_buffer(iNodeL,2) = action_classic*sin(2.0_rp*v_pi*this%spanFreqActuation*zPoint)*cos(this%alphaActuation*v_pi/180_rp)
                  u_buffer(iNodeL,3) = 0.0_rp

               end if
               
            end if
         end do
         !$acc end parallel loop

      end if

   end subroutine BluffBodySolverIncompAFC_afterDt

   subroutine BluffBodySolverIncompAFC_fill_BC_Types(this)
      class(BluffBodySolverIncompAFC), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BluffBodySolverIncompAFC_fill_BC_Types

   subroutine BluffBodySolverIncompAFC_initializeParameters(this)
      use json_module
      implicit none
      class(BluffBodySolverIncompAFC), intent(inout) :: this
      real(rp) :: mul, mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)

      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"cylin"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value

      call json%get("results_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%results_h5_file_path,*) value
      call json%get("results_h5_file_name",value, found,"results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value

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

      call json%get("saveSurfaceResults",this%saveSurfaceResults, found,.false.); call this%checkFound(found,found_aux)
      !----------------------------------------------
      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,1.0_rp); call this%checkFound(found,found_aux)

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
      call json%get("amplitudeActuation",this%amplitudeActuation, found,1.0_rp); call this%checkFound(found,found_aux)  ! Actuation amplitude
      call json%get("frequencyActuation",this%frequencyActuation, found,1.0_rp); call this%checkFound(found,found_aux)  ! Actuation frequency
      call json%get("alphaActuation",this%alphaActuation, found,0.0_rp); call this%checkFound(found,found_aux)          ! Jet angle with respect the y-axis
      call json%get("spanFreqActuation",this%spanFreqActuation, found,0.0_rp); call this%checkFound(found,found_aux)    ! Spanwise frequency
      call json%get("timeBeginActuation",this%timeBeginActuation, found,0.0_rp); call this%checkFound(found,found_aux)  ! Time start actuation
      call json%get("flag_MassConservationStrategy",this%flag_MassConservationStrategy, found,0); call this%checkFound(found,found_aux)   ! Apply mass conservation strategy

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

      call this%readJSONBuffer()

      call json%destroy()

   end subroutine BluffBodySolverIncompAFC_initializeParameters

   subroutine BluffBodySolverIncompAFC_evalInitialConditions(this)
      class(BluffBodySolverIncompAFC), intent(inout) :: this
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

   end subroutine BluffBodySolverIncompAFC_evalInitialConditions


   subroutine BluffBodySolverIncompAFC_initialBuffer(this)
      class(BluffBodySolverIncompAFC), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo*cos(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,2) = this%vo*sin(this%aoa*v_pi/180.0_rp)
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolverIncompAFC_initialBuffer

end module BluffBodySolverIncompAFC_mod
