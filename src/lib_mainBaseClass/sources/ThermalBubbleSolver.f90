module ThermalBubbleSolver_mod
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

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: ThermalBubbleSolver

      real(rp) , public  :: Cv, rho0, mu, to, po, g0, Tc, rc, xc, yc, zc

   contains
      procedure, public :: fillBCTypes           => ThermalBubbleSolver_fill_BC_Types
      procedure, public :: initializeParameters  => ThermalBubbleSolver_initializeParameters
      procedure, public :: initializeSourceTerms => ThermalBubbleSolver_initializeSourceTerms
      procedure, public :: evalInitialConditions => ThermalBubbleSolver_evalInitialConditions
      procedure, public :: afterDt               => ThermalBubbleSolver_afterDt
      procedure, public :: setFields2Save        => ThermalBubbleSolver_setFields2Save
   end type ThermalBubbleSolver
contains

   subroutine ThermalBubbleSolver_fill_BC_Types(this)
      class(ThermalBubbleSolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine ThermalBubbleSolver_fill_BC_Types

   subroutine ThermalBubbleSolver_initializeSourceTerms(this)
      class(ThermalBubbleSolver), intent(inout) :: this
      integer(4) :: iNodeL
      real(rp)   :: source_mom, source_ener

      allocate(source_term(numNodesRankPar,ndime+1))
      !$acc enter data create(source_term(:,:))

      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_mom  = -rho(iNodeL,2)*this%g0 ! rho*g
         source_ener = -q(iNodeL,3,2)*this%g0 ! rho*w*g
         ! Three components of momentum
         source_term(iNodeL,1) = 0.0_rp
         source_term(iNodeL,2) = 0.0_rp
         source_term(iNodeL,3) = source_mom
         ! Energy
         source_term(iNodeL,4) = source_ener
      end do
      !$acc end parallel loop

   end subroutine ThermalBubbleSolver_initializeSourceTerms

   subroutine ThermalBubbleSolver_afterDt(this,istep)
      class(ThermalBubbleSolver), intent(inout) :: this
      integer(4), intent(in) :: istep
      integer(4) :: iNodeL
      real(rp)   :: source_mom, source_ener

      !$acc parallel loop  
      do iNodeL = 1,numNodesRankPar
         source_mom  = -rho(iNodeL,2)*this%g0 ! rho*g
         source_ener = -q(iNodeL,3,2)*this%g0 ! rho*w*g
         ! Three components of momentum
         source_term(iNodeL,1) = 0.0_rp
         source_term(iNodeL,2) = 0.0_rp
         source_term(iNodeL,3) = source_mom
         ! Energy
         source_term(iNodeL,4) = source_ener
      end do
      !$acc end parallel loop

   end subroutine ThermalBubbleSolver_afterDt

   subroutine ThermalBubbleSolver_initializeParameters(this)
      use json_module
      implicit none
      class(ThermalBubbleSolver), intent(inout) :: this
      real(rp) :: mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      call json%initialize()
      call json%load_file(json_filename)

      ! get(label,target,is found?, default value)
     
      !  --------------  File parameters -------------

      call json%get("FileParameters.mesh_h5_file_path", value, found, ""); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_path,*) value
      call json%get("FileParameters.mesh_h5_file_name", value, found, "bubble"); call this%checkFound(found,found_aux)
      write(this%mesh_h5_file_name,*) value

      call json%get("FileParameters.results_h5_file_path", value, found, ""); call this%checkFound(found,found_aux)
      write(this%results_h5_file_path,*) value
      call json%get("FileParameters.results_h5_file_name", value, found, "results"); call this%checkFound(found,found_aux)
      write(this%results_h5_file_name,*) value

      !  --------------  I/O parameters -------------

      call json%get("IOParameters.maxPhysTime", this%maxPhysTime, found, 1200.0_rp); call this%checkFound(found,found_aux)
      call json%get("IOParameters.final_istep", this%final_istep, found, 5000001);   call this%checkFound(found,found_aux)

      call json%get("IOParameters.save_logFile_first", this%save_logFile_first, found, 1);   call this%checkFound(found,found_aux)
      call json%get("IOParameters.save_logFile_step", this%save_logFile_step,   found, 10);  call this%checkFound(found,found_aux)

      call json%get("IOParameters.save_resultsFile_first", this%save_resultsFile_first, found, 1);     call this%checkFound(found,found_aux)
      call json%get("IOParameters.save_resultsFile_step" , this%save_resultsFile_step,  found, 10000); call this%checkFound(found,found_aux)

      call json%get("IOParameters.save_restartFile_first", this%save_restartFile_first, found, 1);     call this%checkFound(found,found_aux)
      call json%get("IOParameters.save_restartFile_step" , this%save_restartFile_step,  found, 10000); call this%checkFound(found,found_aux)

      call json%get("IOParameters.loadRestartFile",     this%loadRestartFile,     found, .true.); call this%checkFound(found,found_aux)
      call json%get("IOParameters.restartFile_to_load", this%restartFile_to_load, found, 1);      call this%checkFound(found,found_aux)

      call json%get("IOParameters.continue_oldLogs", this%continue_oldLogs, found, .false.); call this%checkFound(found,found_aux)

      this%saveAvgFile = .false. ! Do not save averages
      this%loadAvgFile = .false.

      this%saveSurfaceResults = .false. ! No surface results either

      !  --------------  Numerical parameters -------------

      call json%get("NumericalParameters.flag_les",flag_les, found, 1);        call this%checkFound(found,found_aux)
      call json%get("NumericalParameters.c_sgs",   c_sgs,    found, 0.025_rp); call this%checkFound(found,found_aux)
      
      call json%get("NumericalParameters.flag_implicit", flag_implicit, found, 1);       call this%checkFound(found,found_aux)
      call json%get("NumericalParameters.maxIter",       maxIter,       found, 20);      call this%checkFound(found,found_aux)
      call json%get("NumericalParameters.tol",           tol,           found, 0.001d0); call this%checkFound(found,found_aux)
      call json%get("NumericalParameters.envit_ce",      ce,            found, 1.0_rp);  call this%checkFound(found,found_aux)
       
      call json%get("NumericalParameters.cfl_conv", this%cfl_conv, found, 1.5_rp); call this%checkFound(found,found_aux)
      call json%get("NumericalParameters.cfl_diff", this%cfl_diff, found, 1.5_rp); call this%checkFound(found,found_aux)

      call json%get("NumericalParameters.flag_rk_ls",flag_rk_ls, found,.false.); 
      call json%get("NumericalParameters.flag_rk_ls_stages",flag_rk_ls_stages, found,5); 

      flag_walave   = .false.
      !period_walave = 200.0_rp

      !  --------------  Thermodynamic parameters -------------

      call json%get("ThermodynamicParameters.Cp",  this%Cp,   found, 1005.2_rp);   call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.Cv",  this%Cv,   found, 717.1_rp);    call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.R" ,  this%Rgas, found, 287.0_rp);    call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.Prt", this%Prt,  found, 0.71_rp);     call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.mu",  this%mu,   found, 1.0_rp);      call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.po",  this%po,   found, 101325.0_rp); call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.to",  this%to,   found, 288.15_rp);   call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.g",   this%g0,   found, 9.81_rp);     call this%checkFound(found,found_aux)  
      
      if (this%mu .eq. 0.0_rp) flag_real_diff = 0
      flag_diff_suth = 0 ! Deactivate Sutherland viscosity
      ! fixed by the type of base class parameters
      this%gamma_gas = this%Cp/this%Cv
      this%rho0      = this%po/this%Rgas/this%to

      !  --------------  Thermodynamic parameters -------------

      call json%get("BubbleParameters.Tc", this%Tc, found, 0.5_rp);   call this%checkFound(found,found_aux)  
      call json%get("BubbleParameters.rc", this%rc, found, 250.0_rp); call this%checkFound(found,found_aux)  
      call json%get("BubbleParameters.xc", this%xc, found, 500.0_rp); call this%checkFound(found,found_aux)  
      call json%get("BubbleParameters.yc", this%yc, found, 500.0_rp); call this%checkFound(found,found_aux)  
      call json%get("BubbleParameters.zc", this%zc, found, 260.0_rp); call this%checkFound(found,found_aux)  
    
      !  --------------  Witness parameters -------------

      this%have_witness = .false.

      !  --------------  Boundary parameters -------------

      nscbc_u_inf     = 0.0_rp
      nscbc_rho_inf   = this%rho0
      nscbc_p_inf     = this%po
      nscbc_Rgas_inf  = this%Rgas
      nscbc_gamma_inf = this%gamma_gas
      nscbc_T_C       = this%to
      nscbc_g         = this%g0

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'

   end subroutine ThermalBubbleSolver_initializeParameters

   subroutine ThermalBubbleSolver_evalInitialConditions(this)
      class(ThermalBubbleSolver), intent(inout) :: this
      real(rp)   :: x, y, z
      integer(4) :: iNodeL

      ! Set up thermal bubble
      ! TODO: add bubble shapes
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         ! Spatial coordinates
         x = coordPar(iNodeL,1)
         y = coordPar(iNodeL,2)
         z = coordPar(iNodeL,3)
         Tem(iNodeL,2) = this%to ! Background temperature
         ! Now set the thermal bubble shape (TODO: hardcoded for now)
         if(((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc) + (z-this%zc)*(z-this%zc)) .le. this%rc*this%rc) then
            Tem(iNodeL,2) = Tem(iNodeL,2) + this%Tc/2.0_rp*(1.0_rp + cos(v_pi*sqrt((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc) + (z-this%zc)*(z-this%zc))/this%rc))
         end if
      end do
      !$acc end parallel loop

      ! TODO: add atmospheric profiles
      ! Set up initial conditions
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         z = coordPar(iNodeL,3)
         u(iNodeL,1:ndime,2) = 0.0_rp  ! Velocity is set to zero
         ! Navas-Montilla (2023) eq. 71
         pr(iNodeL,2)        = this%po*(1 - (this%gamma_gas-1.0_rp)*this%g0*z/this%gamma_gas/this%Rgas/this%to)**(this%gamma_gas/(this%gamma_gas-1))
         rho(iNodeL,2)       = (this%po/this%Rgas/Tem(iNodeL,2))*(1 - (this%gamma_gas-1.0_rp)*this%g0*z/this%gamma_gas/this%Rgas/this%to)**(this%gamma_gas/(this%gamma_gas-1))
         
         e_int(iNodeL,2)     = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp)) ! Internal energy
         E(iNodeL,2)         = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2)) + e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL)      = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2)       = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
                 
         q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
         rho(iNodeL,3)       = rho(iNodeL,2)
         E(iNodeL,3)         = E(iNodeL,2)
         eta(iNodeL,3)       = eta(iNodeL,2) 
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         machno(iNodeL) = dot_product(u(iNodeL,:,2),u(iNodeL,:,2))/csound(iNodeL)
      end do
      !$acc end parallel loop

      !$acc kernels
      mu_e(:,:)   = 0.0_rp ! Element syabilization viscosity
      mu_sgs(:,:) = 0.0_rp
      kres(:)     = 0.0_rp
      etot(:)     = 0.0_rp
      ax1(:)      = 0.0_rp
      ax2(:)      = 0.0_rp
      ax3(:)      = 0.0_rp
      au(:,:)     = 0.0_rp
      !$acc end kernels
      call nvtxEndRange

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         mu_factor(iNodeL) = flag_mu_factor
      end do
      !$acc end parallel loop

   end subroutine ThermalBubbleSolver_evalInitialConditions

   subroutine ThermalBubbleSolver_setFields2Save(this)
      implicit none
      class(ThermalBubbleSolver), intent(inout) :: this

      if(mpi_rank.eq.0) write(*,*) 'Setting fields to be saved'
      !-----------------------------------------------------------------------

      !---------   nodeScalars  -----------------------------
      !------------------------------------------------------
      call this%add_nodeScalarField2save('Tem',Tem(:,2))
      !------------------------------------------------------
      if(this%save_scalarField_rho) then !rho(numNodesRankPar,3)
         call this%add_nodeScalarField2save('rho',rho(:,2))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_muFluid) then !mu_fluid(numNodesRankPar)
         call this%add_nodeScalarField2save('mu_fluid',mu_fluid(:))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_pressure) then !pr(numNodesRankPar,2)
         call this%add_nodeScalarField2save('pr',pr(:,2))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_energy) then !E(numNodesRankPar,3)
         call this%add_nodeScalarField2save('E',E(:,2))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_entropy) then !eta(numNodesRankPar,3)
         call this%add_nodeScalarField2save('eta',eta(:,2))
      end if
      !------------------------------------------------------
!      if(this%save_scalarField_csound) then !csound(numNodesRankPar)
!         call this%add_nodeScalarField2save('csound',csound(:))
!      end if
!      !------------------------------------------------------
!      if(this%save_scalarField_machno) then !machno(numNodesRankPar)
!         call this%add_nodeScalarField2save('machno',machno(:))
!      end if
!      !------------------------------------------------------
!      if(this%save_scalarField_divU) then !divU(numNodesRankPar)
!         call this%add_nodeScalarField2save('divU',divU(:))
!      end if
!      !------------------------------------------------------
!      if(this%save_scalarField_qcrit) then !qcrit(numNodesRankPar)
!         call this%add_nodeScalarField2save('qcrit',qcrit(:))
!      end if

      !---------------  vectorScalars   -------------------------------------
      !----------------------------------------------------------------------
      if(this%save_vectorField_vel) then !u(numNodesRankPar,ndime,2)
         call this%add_nodeVectorField2save('u',u(:,:,2))
      end if
!      !----------------------------------------------------------------------
!      if(this%save_vectorField_gradRho) then !gradRho(numNodesRankPar,ndime)
!         call this%add_nodeVectorField2save('gradRho',gradRho(:,:))
!      end if
!      !----------------------------------------------------------------------
!      if(this%save_vectorField_curlU) then !curlU(numNodesRankPar,ndime)
!         call this%add_nodeVectorField2save('curlU',curlU(:,:))
!      end if
      !----------------------------------------------------------------------

      !-------------    elemGpScalars   -------------------------------------
      !----------------------------------------------------------------------
      if(this%save_scalarField_muSgs) then !mu_sgs(numElemsRankPar,ngaus)
         call this%add_elemGpScalarField2save('mut',mu_sgs(:,:))
      end if
      !----------------------------------------------------------------------
      if(this%save_scalarField_muEnvit) then !mu_e(numElemsRankPar,ngaus)
         call this%add_elemGpScalarField2save('mue',mu_e(:,:))
      end if
      !----------------------------------------------------------------------

   end subroutine ThermalBubbleSolver_setFields2Save


end module ThermalBubbleSolver_mod
