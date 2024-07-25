module BLFlowSolver_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: BLFlowSolver

      real(rp) , public  ::  M, d0, U0, rho0, Red0, Re, to, po, mu, xmin_tripping, lx_tripping, ly_tripping
      real(rp), public   :: eta_b(45), f(45), f_prim(45)

   contains
      procedure, public :: fillBCTypes           => BLFlowSolver_fill_BC_Types
      procedure, public :: initializeParameters  => BLFlowSolver_initializeParameters
      procedure, public :: evalInitialConditions => BLFlowSolver_evalInitialConditions
      procedure, public :: initialBuffer => BLFlowSolver_initialBuffer
      procedure, public :: fillBlasius => BLFlowSolver_fillBlasius
      procedure, public :: afterDt => BLFlowSolver_afterDt
   end type BLFlowSolver
contains

   subroutine BLFlowSolver_afterDt(this,istep)
      class(BLFlowSolver), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      integer(4) :: iNodeL 
      real(rp) :: cd, lx, ly, xmin, xmax

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
               source_term(iNodeL,3) = -0.5_rp*rho(iNodeL,2)*cd*u(iNodeL,1,2)*abs(u(iNodeL,1,2))/lx
               source_term(iNodeL,4) = 0.00_rp
               source_term(iNodeL,5) = 0.00_rp
            end if
         end if
      end do
      !$acc end parallel loop
   end subroutine BLFlowSolver_afterDt

   subroutine BLFlowSolver_fill_BC_Types(this)
      class(BLFlowSolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BLFlowSolver_fill_BC_Types

   subroutine BLFlowSolver_fillBlasius(this)
      class(BLFlowSolver), intent(inout) :: this

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

  end subroutine BLFlowSolver_fillBlasius

   subroutine BLFlowSolver_initializeParameters(this)
      use json_module
      implicit none
      class(BLFlowSolver), intent(inout) :: this
      real(rp) :: mul, mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)

      call json%get("mesh_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("mesh_h5_file_name",value, found,"bl"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value

      call json%get("results_h5_file_path",value, found,""); call this%checkFound(found,found_aux)
      write(this%results_h5_file_path,*) value
      call json%get("results_h5_file_name",value, found,"results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value

      !----------------------------------------------
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
      !----------------------------------------------

      ! numerical params
      call json%get("flag_les",flag_les, found,0); call this%checkFound(found,found_aux)
      call json%get("flag_implicit",flag_implicit, found,1); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("Cp",this%Cp, found,1004.0_rp); call this%checkFound(found,found_aux)
      call json%get("Prt",this%Prt, found,0.71_rp); call this%checkFound(found,found_aux)
      call json%get("M",this%M, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("d0",this%d0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("U0",this%U0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Red0",this%Red0, found,100.0_rp); call this%checkFound(found,found_aux)
      call json%get("gamma_gas",this%gamma_gas, found,1.4_rp); call this%checkFound(found,found_aux)

      call json%get("flag_rk_order",flag_rk_ls_stages, found,4);

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

      ! Tripping region
      call json%get("xmin_tripping",this%xmin_tripping, found,0.0); call this%checkFound(found,found_aux)
      call json%get("lx_tripping",this%lx_tripping, found,1.0); call this%checkFound(found,found_aux)
      call json%get("ly_tripping",this%ly_tripping, found,1.0); call this%checkFound(found,found_aux)

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

   end subroutine BLFlowSolver_initializeParameters

   subroutine BLFlowSolver_initialBuffer(this)
      class(BLFlowSolver), intent(inout) :: this
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

      !u_buffer(inodel,1) = (-1.0_rp/160000.0_rp)*yp*yp + yp*(1.0_rp/200.0_rp)

   end subroutine BLFlowSolver_initialBuffer

   subroutine BLFlowSolver_evalInitialConditions(this)
      class(BLFlowSolver), intent(inout) :: this
      integer(4) :: iNodeL, j,k
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

         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%rho0
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
   end subroutine BLFlowSolver_evalInitialConditions

end module BLFlowSolver_mod
