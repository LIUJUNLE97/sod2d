module BluffBodySolver_mod
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
   use CFDSolverPeriodicWithBoundaries_mod
   implicit none
   private

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BluffBodySolver

      real(rp) , public  :: vo, M, delta, rho0, Re, to, po, userDefinedDt = -1
      logical, public    :: useFixDt = .false.

      character(100) :: initialConditionType    ! FromZero -> evalInitialConditionsFromZero
                                                ! FromFile -> evalInitialConditionsFromFile
                                                ! FromPrevRun -> evalInitialConditionsFromPreviousRun

   contains
      procedure, public :: fillBCTypes           =>BluffBodySolver_fill_BC_Types
      procedure, public :: initializeParameters  => BluffBodySolver_initializeParameters
!      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditionsFromZero
!      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditionsFromFile
!      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditionsFromPreviousRun
      procedure, public :: evalInitialConditions => BluffBodySolver_chooseInitialConditionsType

      procedure, private:: readConfigFile       =>  BluffBodySolver_readConfigFile
      procedure, private:: evalInitialConditionsFromZero   =>      BluffBodySolver_evalInitialConditionsFromZero
      procedure, private:: evalInitialConditionsFromFile   =>      BluffBodySolver_evalInitialConditionsFromFile
      procedure, private:: evalInitialConditionsFromPreviousRun => BluffBodySolver_evalInitialConditionsFromPreviousRun
!      procedure, public :: evalInitialConditions => BluffBodySolver_evalInitialConditions
      procedure, public :: initialBuffer => BluffBodySolver_initialBuffer
      procedure, public :: afterDt => BluffBodySolver_afterDt
   end type BluffBodySolver
contains

   subroutine BluffBodySolver_fill_BC_Types(this)
      class(BluffBodySolver), intent(inout) :: this

      bouCodes2BCType(1) = bc_type_far_field
      bouCodes2BCType(2) = bc_type_far_field
      bouCodes2BCType(3) = bc_type_far_field
      bouCodes2BCType(4) = bc_type_far_field
      bouCodes2BCType(5) = bc_type_non_slip_adiabatic
      !$acc update device(bouCodes2BCType(:))

   end subroutine BluffBodySolver_fill_BC_Types

   subroutine BluffBodySolver_initializeParameters(this)
      class(BluffBodySolver), intent(inout) :: this
      real(rp) :: mul, mur

      write(this%mesh_h5_file_path,*) ""
      write(this%mesh_h5_file_name,*) "mesh"
      !write(this%mesh_h5_file_name,*) "naca"

      write(this%results_h5_file_path,*) ""
      write(this%results_h5_file_name,*) "results"

      ! numerical params
      flag_les = 1
      flag_implicit = 0
      flag_les_ilsa=0
      flag_rk_order=4
      tol=1e-3

      this%cfl_conv = 0.5_rp 
      this%cfl_diff = 0.5_rp 

      !----------------------------------------------
      !  --------------  I/O params -------------
      this%final_istep = 1000000 

      this%save_logFile_first = 1
      this%save_logFile_step  = 10

      this%save_resultsFile_first = 1
      this%save_resultsFile_step = 10000

      this%save_restartFile_first = 1
      this%save_restartFile_step = 10000
      this%loadRestartFile = .true.
      this%restartFile_to_load = 1 !1 or 2
      this%continue_oldLogs = .true.

      this%saveAvgFile = .true.
      this%loadAvgFile = .false.
      !----------------------------------------------

      this%Cp = 1004.0_rp
      this%Prt = 0.71_rp
      this%vo = 2.0_rp
      this%M  = 1.0_rp
      this%delta  = 1.0_rp
      this%rho0   = 1.0_rp
      this%gamma_gas = 1.40_rp
      this%Re     =  10000.0_rp


      call this%readConfigFile("sod2d_config.txt")

      if(mpi_rank.eq.0) then
        write(*,*) "Delta ", this%delta
        write(*,*) "vo ", this%vo
        write(*,*) "Re ", this%Re

        write(*,*) "Reference time ", this%delta/this%vo
      endif

      mul    = (this%rho0*this%delta*this%vo)/this%Re
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = mul/mur

      if(mpi_rank.eq.0) then
        write(*,*) "To ", this%to
        write(*,*) "Po ", this%po
        write(*,*) "mur ", mur
        write(*,*) "flag_mu_factor ", flag_mu_factor


      endif
      nscbc_u_inf = this%vo
      nscbc_p_inf = this%po
      nscbc_rho_inf = this%rho0
      nscbc_gamma_inf = this%gamma_gas
      nscbc_c_inf = sqrt(this%gamma_gas*this%po/this%rho0)
      nscbc_Rgas_inf = this%Rgas

      !Witness points parameters
!      this%have_witness          = .false.
!      this%witness_inp_file_name = "witness.txt"
!      this%witness_h5_file_name  = "resultwit.h5"
!      this%leapwit               = 1
!      this%nwit                  = 17986
!      this%wit_save_u_i          = .true.
!      this%wit_save_pr           = .true.
!      this%wit_save_rho          = .true.
!      this%continue_witness      = .false.     

 
!      flag_buffer_on = .true.
!     !cylinder
!     flag_buffer_on_east = .true.
!     flag_buffer_e_min = 40.0_rp
!     flag_buffer_e_size = 10.0_rp

!     flag_buffer_on_west = .true.
!     flag_buffer_w_min = -20.0_rp
!     flag_buffer_w_size = 10.0_rp

!     flag_buffer_on_north = .true.
!     flag_buffer_n_min = 20.0_rp
!     flag_buffer_n_size = 10.0_rp

!     flag_buffer_on_south = .true.
!     flag_buffer_s_min = -20.0_rp
!     flag_buffer_s_size = 10.0_rp

      !naca
     !flag_buffer_on_east = .true.
     !flag_buffer_e_min = 10.0_rp
     !flag_buffer_e_size = 5.0_rp 

     !flag_buffer_on_west = .true.
     !flag_buffer_w_min = -10.0_rp
     !flag_buffer_w_size = 2.5_rp 
     
     !flag_buffer_on_north = .true.
     !flag_buffer_n_min = 10.0_rp
     !flag_buffer_n_size = 2.5_rp 
     
     !flag_buffer_on_south = .true.
     !flag_buffer_s_min = -10.0_rp
     !flag_buffer_s_size = 2.5_rp 

     if (this%useFixDt)then
         if(this%userDefinedDt.lt.0.)then
             write(*,*) "A condition to fix the timestep has been set, but no dt value has been provided to userDefinedDt"
             stop
         endif
     endif
   end subroutine BluffBodySolver_initializeParameters

   subroutine BluffBodySolver_chooseInitialConditionsType(this)
        class(BluffBodySolver), intent(inout) :: this

       select case (this%initialConditionType)
        case('FromZero') 
            call this%evalInitialConditionsFromZero()
        case('FromAlya') 
            call this%evalInitialConditionsFromFile()
        case('FromPrevRun') 
            call this%evalInitialConditionsFromPreviousRun()
       case default
         ! Ignore unknown variable
         if(mpi_rank.eq.0) write(111,*) "Type of initialisation (", this%initialConditionType ,") not recognised. Aborting."
         stop
       end select


   end subroutine BluffBodySolver_chooseInitialConditionsType

   subroutine BluffBodySolver_evalInitialConditionsFromPreviousRun(this)
        class(BluffBodySolver), intent(inout) :: this
        integer(4) :: iNodeL
!        real(rp) :: vo_old

!       if(mpi_rank.eq.0) write(111,*) "Initialising the case based on a previous SOD2D case. Scaling the old case velocities and clipping." 
       if(mpi_rank.eq.0) write(111,*) "Initialising the case based on a previous SOD2D case. Loading only velocities." 

!        vo_old = 10._rp
!        vo_old = 20._rp

        if(mpi_rank.eq.0) write(111,*) "--| Loading results load_step ",this%load_step
!         call load_hdf5_restartFile(            this%restartFile_to_load,this%load_step,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,mu_sgs)
        call load_hdf5_restartFile(nnode,ngaus,this%restartFile_to_load,this%load_step,flag_walave,this%time,rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,mu_sgs,walave_u)
        if(mpi_rank.eq.0) write(111,*) "   --| Loaded results for load_step",this%load_step,"time",this%time

!      !$acc parallel loop
!      do iNodeL = 1,numNodesRankPar
!        ! Scale read values from input file based on the reference velocity 
!         u(iNodeL,1,2) = u(iNodeL,1,2)*this%vo/vo_old
!         u(iNodeL,2,2) = u(iNodeL,2,2)*this%vo/vo_old
!         u(iNodeL,3,2) = u(iNodeL,3,2)*this%vo/vo_old
!
!         ! Clip values if simulation is behaving oddly
!         if( abs(u(iNodeL,1,2)) .gt. 2.*this%vo)  u(iNodeL,1,2) = sign(2.*this%vo,u(iNodeL,1,2))
!
!         if( abs(u(iNodeL,2,2)) .gt. 2.*this%vo)  u(iNodeL,2,2) = sign(2.*this%vo,u(iNodeL,2,2))
!
!         if( abs(u(iNodeL,3,2)) .gt. 2.*this%vo)  u(iNodeL,3,2) = sign(2.*this%vo,u(iNodeL,3,2))
!
!      end do
!      !$acc end parallel loop

!      call this%eval_vars_after_load_hdf5_resultsFile()
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%po/this%Rgas/this%to
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)
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

      if(this%continue_oldLogs) then
!          this%initial_istep = this%load_step
!          this%nsave = this%load_step+this%nleap
!          this%nsaveAVG = this%load_step+this%nleapAVG
!          this%nsave2 = this%load_step+this%nleap2
         this%isFreshStart = .false.
      else
         this%time = 0.0_rp
      end if


   end subroutine BluffBodySolver_evalInitialConditionsFromPreviousRun

   subroutine BluffBodySolver_evalInitialConditionsFromFile(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles
      character(512) :: initialField_filePath

      if(mpi_rank.eq.0) write(111,*) "Initialising the case using a VELOC.alya file"

      call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
      write(initialField_filePath,*) ""
      call read_veloc_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,u(:,:,2),connecParOrig,Ngp_l,matGidSrlOrdered)
!      call read_densi_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,rho(:,2),connecParOrig,Ngp_l,matGidSrlOrdered)
!      call read_press_from_file_Par(numElemsRankPar,numNodesRankPar,totalNumNodesSrl,initialField_filePath,pr(:,2),connecParOrig,Ngp_l,matGidSrlOrdered)

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
        ! Scale read values from input file based on the reference velocity 
         u(iNodeL,1,2) = u(iNodeL,1,2)*this%vo
         u(iNodeL,2,2) = u(iNodeL,2,2)*this%vo
         u(iNodeL,3,2) = u(iNodeL,3,2)*this%vo
!         u(iNodeL,1,2) = this%vo
!         u(iNodeL,2,2) = 0.0_rp
!         u(iNodeL,3,2) = 0.0_rp
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%po/this%Rgas/this%to
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*this%Rgas)
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)
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

   end subroutine BluffBodySolver_evalInitialConditionsFromFile

   subroutine BluffBodySolver_evalInitialConditionsFromZero(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(rp) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4) :: iNodeL
      logical :: readFiles

      if(mpi_rank.eq.0) write(111,*) "Initialising the case from 0"
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         u(iNodeL,1,2) = this%vo
         u(iNodeL,2,2) = 0.0_rp
         u(iNodeL,3,2) = 0.0_rp
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%po/this%Rgas/this%to
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

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolver_evalInitialConditionsFromZero


   subroutine BluffBodySolver_initialBuffer(this)
      class(BluffBodySolver), intent(inout) :: this
      integer(4) :: iNodeL

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
            u_buffer(iNodeL,1) = this%vo
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
      end do
      !$acc end parallel loop

   end subroutine BluffBodySolver_initialBuffer
   subroutine BluffBodySolver_afterDt(this,istep)
      class(BluffBodySolver), intent(inout) :: this
      integer(4), intent(in) :: istep
      character(4) :: timeStep

!         if (mpi_rank.eq.0) write(111,*) "Call BluffBodySolver_afterDt"
!         if (mpi_rank.eq.0) write(111,*) "this%dt ", this%dt, " this%userDefinedDt ", this%userDefinedDt, " ", this%useFixDt, " ",  (this%dt .gt. this%userDefinedDt)

      if(this%useFixDt) then
          if (this%dt .gt. this%userDefinedDt) then
              ! Only oversubscribe the dt if dt defined by the user is smaller 
              ! than the one required by stability
              this%dt = this%userDefinedDt
          else
              if (mpi_rank.eq.0) write(111,*) "Warning at iteration", istep ,", user defined dt too large, using the one dictated by stability: ", this%dt
          endif
      endif

      if(this%dt .gt. 2e-4) then
         if (mpi_rank.eq.0) write(111,*) "At iteration ", istep, " next time step ", this%dt, " . Capping it at 2e-4 ." 


         if (mpi_rank.eq.0) write(111,*) ' - Extra saving results file step:',istep
         call nvtxStartRange("Output "//timeStep,istep)
         call compute_fieldDerivs(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,workingNodesPar,connecParWork,lelpn,He,dNgp,this%leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),gradRho,curlU,divU,Qcrit)
         call this%saveInstResultsFiles(istep)
         call nvtxEndRange

         this%dt = 2e-4

         stop
      endif


   end subroutine BluffBodySolver_afterDt


   subroutine BluffBodySolver_readConfigFile(this, filename)
     class(BluffBodySolver), intent(inout) :: this
     character(len=*), intent(in) :: filename
     character(len=256) :: line
     integer :: status
     integer :: ieq, iexcl
     character(len=100) :: varname, varvalue
!     real :: varvalue

     if(mpi_rank.eq.0) write(*,*) "--| Reading configuration file ", trim(filename)

     open(10, file=filename, status='old', action='read')
     do
       read(10, '(A)', iostat=status) line
       if (status /= 0) exit

!        line = trim(line) ! Remove leading/trailing spaces
        if (len(trim(line)) == 0) cycle ! Skip empty lines
!        if (index(line, "!") == 1) cycle ! Skip comment lines

       ! Ignore comments at end of line
       iexcl = index(trim(line), '!')       ! Position where comments start
                                            ! Index returns 0 is the searched char is not found
       if (iexcl > 0) line = trim(line(1:iexcl-1))  
       
       if (iexcl.eq.1) cycle                           ! If comments at the beginning of the line, skip line

       ! Parse variable name and value from line
       ieq = index(line, '=')
       varname = trim(line(1:ieq-1))
       varvalue = trim(line(ieq+1:))
       read(line(ieq+1:), *) varvalue

       ! Assign value to class variable
       select case (varname)

       case('saveInitialField') ! Save initial fla
        read(varvalue, *) this%saveInitialField
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%saveInitialField
!        case('loadResults') ! Restart from a previous run using the .h5 results file
!         read(varvalue, *) this%loadResults
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%loadResults
!        case('continue_oldLogs') ! If restarting, continue the log files from the previous run
!         read(varvalue, *) this%continue_oldLogs
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%continue_oldLogs
       case('load_step')  ! If restarting, load the results from this step
        read(varvalue, *) this%load_step
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%load_step

       case('iniCondType')  ! If not restaring, define the type of initialisation
        read(varvalue, *) this%initialConditionType
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", trim(this%initialConditionType)

       case('ce')  ! If not restaring, define the type of initialisation
        read(varvalue, *) ce
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", ce

!        case('nstep')  ! Number of steps to run
!         read(varvalue, *) this%nstep
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nstep
       case('cfl_conv') !0.5_rp !0.1_rp    
        read(varvalue, *) this%cfl_conv
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%cfl_conv
       case('cfl_diff') !0.5_rp !0.1_rp
        read(varvalue, *) this%cfl_diff
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%cfl_diff

!        case('nsave')                   ! First step to save, 
!         read(varvalue, *) this%nsave
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nsave
!        case('nsave2')                    ! First step to save,
!         read(varvalue, *) this%nsave2
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nsave2
!        case('nsaveAVG') 
!         read(varvalue, *) this%nsaveAVG
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nsaveAVG

!        case('nleap')                ! Saving interval in iterations,
!         read(varvalue, *) this%nleap
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nleap
!        case('tleap')                ! Saving interval in time units,
!         read(varvalue, *) this%tleap
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%tleap
!        case('nleap2')               ! Saving interval,
!         read(varvalue, *) this%nleap2
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nleap2
!        case('nleapAVG') 
!         read(varvalue, *) this%nleapAVG
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nleapAVG

       case('final_istep')                
        read(varvalue, *) this%final_istep
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%final_istep
       case('save_logFile_first')                
        read(varvalue, *) this%save_logFile_first
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_logFile_first
       case('save_logFile_step')               ! Saving interval,
        read(varvalue, *) this%save_logFile_step
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_logFile_step
       case('save_resultsFile_first') 
        read(varvalue, *) this%save_resultsFile_first
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_resultsFile_first
       case('save_resultsFile_step') 
        read(varvalue, *) this%save_resultsFile_step
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_resultsFile_step


       case('save_restartFile_first') 
        read(varvalue, *) this%save_restartFile_first
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_restartFile_first
       case('save_restartFile_step') 
        read(varvalue, *) this%save_restartFile_step
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_restartFile_step
       case('loadRestartFile') 
        read(varvalue, *) this%loadRestartFile
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%loadRestartFile
       case('restartFile_to_load') 
        read(varvalue, *) this%restartFile_to_load
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%restartFile_to_load

       case('continue_oldLogs') 
        read(varvalue, *) this%continue_oldLogs
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%continue_oldLogs
       case('saveAvgFile') 
        read(varvalue, *) this%saveAvgFile
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%saveAvgFile
       case('loadAvgFile') 
        read(varvalue, *) this%loadAvgFile
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%loadAvgFile


       case('Cp') 
        read(varvalue, *) this%Cp
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%Cp
       case('Prt') 
        read(varvalue, *) this%Prt
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%Prt
       case('vo') ! 10.0_rp
        read(varvalue, *) this%vo
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%vo
       case('M')   ! 0.03_rp
        read(varvalue, *) this%M
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%M
       case('delta')  
        read(varvalue, *) this%delta
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%delta
       case('rho0')   
        read(varvalue, *) this%rho0
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%rho0
       case('gamma_gas') 
        read(varvalue, *) this%gamma_gas
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%gamma_gas
       case('Re')     
        read(varvalue, *) this%Re
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%Re
       case('flag_implicit')
        read(varvalue, *) flag_implicit
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_implicit
       case('flag_les')
        read(varvalue, *) flag_les
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_les
!        case('pseudo_cfl')
!         read(varvalue, *) pseudo_cfl
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", pseudo_cfl
!        case('pseudo_ftau')
!         read(varvalue, *) pseudo_ftau
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", pseudo_ftau
!        case('maxIterNonLineal')
!         read(varvalue, *) maxIterNonLineal
!         if(mpi_rank.eq.0) write(*,*) trim(varname), " ", maxIterNonLineal
       case('tol')
        read(varvalue, *) tol
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", tol
       case('useFixDt')     
        read(varvalue, *) this%useFixDt
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%useFixDt
       case('userDefinedDt')     
        read(varvalue, *) this%userDefinedDt
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%userDefinedDt

       case('flag_buffer_on') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_on
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_on

       case('flag_buffer_on_east') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_on_east
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_on_east
       case('flag_buffer_e_min') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_e_min
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_e_min
       case('flag_buffer_e_size') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_e_size
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_e_size

       case('flag_buffer_on_west') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_on_west
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_on_west
       case('flag_buffer_w_min') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_w_min
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_w_min
       case('flag_buffer_w_size') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_w_size
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_w_size

       case('flag_buffer_on_north') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_on_north
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_on_north
       case('flag_buffer_n_min') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_n_min
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_n_min
       case('flag_buffer_n_size') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_n_size
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_n_size
        
       case('flag_buffer_on_south') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_on_south
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_on_south
       case('flag_buffer_s_min') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_s_min
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_s_min
       case('flag_buffer_s_size') ! Restart from a previous run using the .h5 results file
        read(varvalue, *) flag_buffer_s_size
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", flag_buffer_s_size


       case('save_scalarField_rho')     
        read(varvalue, *) this%save_scalarField_rho
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_rho
       case('save_scalarField_muFluid')     
        read(varvalue, *) this%save_scalarField_muFluid
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_muFluid
       case('save_scalarField_pressure')     
        read(varvalue, *) this%save_scalarField_pressure
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_pressure
       case('save_scalarField_energy')     
        read(varvalue, *) this%save_scalarField_energy
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_energy
       case('save_scalarField_entropy')     
        read(varvalue, *) this%save_scalarField_entropy
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_entropy
       case('save_scalarField_csound')     
        read(varvalue, *) this%save_scalarField_csound
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_csound
       case('save_scalarField_machno')     
        read(varvalue, *) this%save_scalarField_machno
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_machno
       case('save_scalarField_divU')     
        read(varvalue, *) this%save_scalarField_divU
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_divU
       case('save_scalarField_qcrit')     
        read(varvalue, *) this%save_scalarField_qcrit
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_qcrit
       case('save_scalarField_muSgs')     
        read(varvalue, *) this%save_scalarField_muSgs
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_muSgs
       case('save_scalarField_muEnvit')     
        read(varvalue, *) this%save_scalarField_muEnvit
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_scalarField_muEnvit
       case('save_vectorField_vel')     
        read(varvalue, *) this%save_vectorField_vel
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_vectorField_vel
       case('save_vectorField_gradRho')     
        read(varvalue, *) this%save_vectorField_gradRho
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_vectorField_gradRho
       case('save_vectorField_curlU')     
        read(varvalue, *) this%save_vectorField_curlU
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_vectorField_curlU


       case('save_avgScalarField_rho')     
        read(varvalue, *) this%save_avgScalarField_rho
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgScalarField_rho
       case('save_avgScalarField_pr')     
        read(varvalue, *) this%save_avgScalarField_pr
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgScalarField_pr
       case('save_avgScalarField_mueff')     
        read(varvalue, *) this%save_avgScalarField_mueff
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgScalarField_mueff
       case('save_avgVectorField_vel')     
        read(varvalue, *) this%save_avgVectorField_vel
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgVectorField_vel
       case('save_avgVectorField_ve2')     
        read(varvalue, *) this%save_avgVectorField_ve2
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgVectorField_ve2
       case('save_avgVectorField_vex')     
        read(varvalue, *) this%save_avgVectorField_vex
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgVectorField_vex
       case('save_avgVectorField_vtw')     
        read(varvalue, *) this%save_avgVectorField_vtw
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%save_avgVectorField_vtw


       case('have_witness')     
        read(varvalue, *) this%have_witness
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%have_witness
       case('witness_inp_file_name')     
        read(varvalue, *) this%witness_inp_file_name
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", trim(this%witness_inp_file_name)
       case('witness_h5_file_name')     
        read(varvalue, *) this%witness_h5_file_name
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", trim(this%witness_h5_file_name)
       case('leapwit')     
        read(varvalue, *) this%leapwit
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%leapwit
       case('leapwitsave')     
        read(varvalue, *) this%leapwitsave
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%leapwitsave
       case('nwit')     
        read(varvalue, *) this%nwit
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%nwit
       case('wit_save_u_i')     
        read(varvalue, *) this%wit_save_u_i
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%wit_save_u_i
       case('wit_save_pr')     
        read(varvalue, *) this%wit_save_pr
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%wit_save_pr
       case('wit_save_rho')     
        read(varvalue, *) this%wit_save_rho
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%wit_save_rho
       case('continue_witness')     
        read(varvalue, *) this%continue_witness
        if(mpi_rank.eq.0) write(*,*) trim(varname), " ", this%continue_witness
        
       case default
         ! Ignore unknown variable
         if(mpi_rank.eq.0) write(*,*) "Variable ", trim(varname), " not recognised. Skipping it."
       end select
     end do
     close(10)
   end subroutine BluffBodySolver_readConfigFile

end module BluffBodySolver_mod
