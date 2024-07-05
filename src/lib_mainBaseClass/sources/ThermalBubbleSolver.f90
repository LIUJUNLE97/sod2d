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

   ! Types for atmospheric profile
   integer(4), parameter :: atmos_type_adiabatic = 1

   ! Types for bubble shape
   integer(4), parameter :: bubble_shape_cosine2D   = 1
   integer(4), parameter :: bubble_shape_gaussian2D = 2
   integer(4), parameter :: bubble_shape_cosine3D   = 3   

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: ThermalBubbleSolver

      real(rp) , public  :: Cv, rho0, mu, T0, p0, g0, Tc, rc, xc, yc, zc
      integer(4), public :: atmos_type, bubble_shape
      logical, public    :: save_scalarField_Tem


   contains
      procedure, public :: fillBCTypes           => ThermalBubbleSolver_fill_BC_Types
      procedure, public :: initializeParameters  => ThermalBubbleSolver_initializeParameters
      procedure, public :: evalInitialConditions => ThermalBubbleSolver_evalInitialConditions
      procedure, public :: setFields2Save        => ThermalBubbleSolver_setFields2Save
   end type ThermalBubbleSolver
contains

   subroutine ThermalBubbleSolver_fill_BC_Types(this)
      class(ThermalBubbleSolver), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine ThermalBubbleSolver_fill_BC_Types

   subroutine ThermalBubbleSolver_initializeParameters(this)
      use json_module
      implicit none
      class(ThermalBubbleSolver), intent(inout) :: this
      real(rp) :: mur
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value

      flag_high_mach       = .false.
      flag_bouyancy_effect = .true.

      call json%initialize()
      call json%load_file(json_filename)

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

      call json%get("IOParameters.save_scalarField_rho",      this%save_scalarField_rho,      found, .true.); 
      call json%get("IOParameters.save_scalarField_Tem",      this%save_scalarField_Tem,      found, .true.); 
      call json%get("IOParameters.save_scalarField_pressure", this%save_scalarField_pressure, found, .true.); 
      call json%get("IOParameters.save_scalarField_energy",   this%save_scalarField_energy,   found, .true.); 
      call json%get("IOParameters.save_scalarField_muFluid",  this%save_scalarField_muFluid,  found, .false.); 
      call json%get("IOParameters.save_scalarField_entropy",  this%save_scalarField_entropy,  found, .true.); 
      call json%get("IOParameters.save_scalarField_csound",   this%save_scalarField_csound,   found, .false.); 
      call json%get("IOParameters.save_scalarField_machno",   this%save_scalarField_machno,   found, .false.); 
      call json%get("IOParameters.save_scalarField_divU",     this%save_scalarField_divU,     found, .false.); 
      call json%get("IOParameters.save_scalarField_qcrit",    this%save_scalarField_qcrit,    found, .false.); 
      call json%get("IOParameters.save_scalarField_muSgs",    this%save_scalarField_muSgs,    found, .false.); 
      call json%get("IOParameters.save_scalarField_muEnvit",  this%save_scalarField_muEnvit,  found, .true.); 
      call json%get("IOParameters.save_vectorField_vel",      this%save_vectorField_vel,      found, .true.); 
      call json%get("IOParameters.save_vectorField_gradRho",  this%save_vectorField_gradRho,  found, .false.); 
      call json%get("IOParameters.save_vectorField_curlU",    this%save_vectorField_curlU,    found, .false.); 

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

      call json%get("ThermodynamicParameters.Cp",  this%Cp,   found, 1005.2_rp); call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.Cv",  this%Cv,   found, 717.1_rp);  call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.R" ,  this%Rgas, found, 287.0_rp);  call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.Prt", this%Prt,  found, 0.71_rp);   call this%checkFound(found,found_aux)
      call json%get("ThermodynamicParameters.mu",  this%mu,   found, 1.0_rp);    call this%checkFound(found,found_aux)
      
      if (this%mu .eq. 0.0_rp) then
         flag_real_diff = 0
         flag_mu_factor = 0.0_rp
      end if
      flag_diff_suth = 0 ! Deactivate Sutherland viscosity
      ! fixed by the type of base class parameters
      this%gamma_gas = this%Cp/this%Cv

      !  --------------  Atmosphere parameters -------------
      
      call json%get("AtmosphericParameters.type", value, found, "adiabatic"); call this%checkFound(found,found_aux)
      if (value .eq. "adiabatic") then
         this%atmos_type = atmos_type_adiabatic
      else
         write(*,*) "INVALID ATMOSPHERIC TYPE!"
         stop 1
      end if
      call json%get("AtmosphericParameters.p", this%p0, found, 100000.0_rp); call this%checkFound(found,found_aux)
      call json%get("AtmosphericParameters.T", this%T0, found, 288.15_rp);   call this%checkFound(found,found_aux)
      call json%get("AtmosphericParameters.g", this%g0, found, 9.81_rp);     call this%checkFound(found,found_aux)  

      this%rho0 = this%p0/this%Rgas/this%T0

      !  --------------  Thermodynamic parameters -------------

      call json%get("BubbleParameters.shape", value, found, "cosine2D"); call this%checkFound(found,found_aux)
      if (value .eq. "cosine2D") then
         this%bubble_shape = bubble_shape_cosine2D 
         flag_force_2D     = .true. ! Sets the z component to zero
      elseif (value .eq. "gaussian2D") then
         this%bubble_shape = bubble_shape_gaussian2D
         flag_force_2D     = .true. ! Sets the z component to zero
      elseif (value .eq. "cosine3D") then
         this%bubble_shape = bubble_shape_cosine3D
      else
         write(*,*) "INVALID BUBBLE SHAPE!"
         stop 1
      end if
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
      nscbc_p_inf     = this%p0
      nscbc_Rgas_inf  = this%Rgas
      nscbc_Cp_inf    = this%Cp
      nscbc_gamma_inf = this%gamma_gas
      nscbc_T_C       = this%T0
      nscbc_g_x       = 0.0_rp
      nscbc_g_y       = this%g0
      nscbc_g_z       = 0.0_rp

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'

   end subroutine ThermalBubbleSolver_initializeParameters

   subroutine ThermalBubbleSolver_evalInitialConditions(this)
      class(ThermalBubbleSolver), intent(inout) :: this
      real(rp)   :: x, y, z, aux, dtheta
      integer(4) :: iNodeL

      ! Set up the atmosphere
      if (this%atmos_type .eq. atmos_type_adiabatic) then ! Adiabatic atmosphere
         aux = (this%gamma_gas - 1.0_rp)/this%gamma_gas*this%g0/this%Rgas/this%T0
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            ! Spatial coordinates
            y = coordPar(iNodeL,2) 
            ! Set up the atmosphere, Navas-Montilla (2023) eq. 71 (adiabatic atmosphere)
            ! https://farside.ph.utexas.edu/teaching/sm1/lectures/node56.html
            Tem(iNodeL,2) = this%T0*(1.0_rp - aux*y)
            pr(iNodeL,2)  = this%p0*(1.0_rp - aux*y)**(this%Cp/this%Rgas) !(this%gamma_gas/(this%gamma_gas-1.0_rp))
            ! WARNING: here for numerical approximation is important to either use
            ! (gamma-1)/gamma or Cp/R as otherwise we incur in approximation issues!
         end do
         !$acc end parallel loop
      ! TODO: add more atmospheres
      endif

      ! Set up thermal bubble
      if (this%bubble_shape .eq. bubble_shape_cosine2D) then
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            ! Spatial coordinates
            x = coordPar(iNodeL,1)
            y = coordPar(iNodeL,2)
            z = coordPar(iNodeL,3)
            ! Now set the thermal bubble shape
            if(((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc)) .le. this%rc*this%rc) then
               aux = this%Tc/2.0_rp*(1.0_rp + cos(v_pi*sqrt((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc))/this%rc))
            end if
            ! We computed the perturbation in terms of potential temperature now convert it to static temperature
            ! https://en.wikipedia.org/wiki/Potential_temperature
            Tem(iNodeL,2) = Tem(iNodeL,2) + aux*(pr(iNodeL,2)/this%p0)**(this%Rgas/this%Cp)
            ! WARNING: here for numerical approximation is important to either use
            ! (gamma/gamma-1) or R/Cp as otherwise we incur in approximation issues!
         end do
         !$acc end parallel loop        
      elseif (this%bubble_shape .eq. bubble_shape_gaussian2D) then
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            ! Spatial coordinates
            x = coordPar(iNodeL,1)
            y = coordPar(iNodeL,2)
            z = coordPar(iNodeL,3)
            ! Now set the thermal bubble shape
            if(((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc)) .le. this%rc*this%rc) then
               aux = this%Tc
            else
               aux = this%Tc*exp(-sqrt(((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc)) - 50.0_rp)*sqrt(((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc)) - 50.0_rp)/100.0_rp/100.0_rp)
            endif
            ! We computed the perturbation in terms of potential temperature now convert it to static temperature
            ! https://en.wikipedia.org/wiki/Potential_temperature
            Tem(iNodeL,2) = Tem(iNodeL,2) + aux*(pr(iNodeL,2)/this%p0)**(this%Rgas/this%Cp)
            ! WARNING: here for numerical approximation is important to either use
            ! (gamma/gamma-1) or R/Cp as otherwise we incur in approximation issues!
         end do
         !$acc end parallel loop   
      elseif (this%bubble_shape .eq. bubble_shape_cosine3D) then
         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            ! Spatial coordinates
            x = coordPar(iNodeL,1)
            y = coordPar(iNodeL,2)
            z = coordPar(iNodeL,3)
            ! Now set the thermal bubble shape
            if(((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc) + (z-this%zc)*(z-this%zc)) .le. this%rc*this%rc) then
               aux = this%Tc/2.0_rp*(1.0_rp + cos(v_pi*sqrt((x-this%xc)*(x-this%xc) + (y-this%yc)*(y-this%yc) + (z-this%zc)*(z-this%zc))/this%rc))
            end if
            ! We computed the perturbation in terms of potential temperature now convert it to static temperature
            ! https://en.wikipedia.org/wiki/Potential_temperature
            Tem(iNodeL,2) = Tem(iNodeL,2) + aux*(pr(iNodeL,2)/this%p0)**(this%Rgas/this%Cp)
            ! WARNING: here for numerical approximation is important to either use
            ! (gamma/gamma-1) or R/Cp as otherwise we incur in approximation issues!
         end do
         !$acc end parallel loop
      endif

      ! Set up initial conditions
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         ! Velocity is set to zero
         u(iNodeL,1:ndime,2) = 0.0_rp
         ! Density is set with the equation of state
         rho(iNodeL,2)       = pr(iNodeL,2)/this%Rgas/Tem(iNodeL,2)

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
      mu_e(:,:)   = 0.0_rp ! Element stabilization viscosity
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
      if(this%save_scalarField_Tem) then !Tem(numNodesRankPar,3)
         call this%add_nodeScalarField2save('Tem',Tem(:,2))
      end if
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
      if(this%save_scalarField_csound) then !csound(numNodesRankPar)
         call this%add_nodeScalarField2save('csound',csound(:))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_machno) then !machno(numNodesRankPar)
         call this%add_nodeScalarField2save('machno',machno(:))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_divU) then !divU(numNodesRankPar)
         call this%add_nodeScalarField2save('divU',divU(:))
      end if
      !------------------------------------------------------
      if(this%save_scalarField_qcrit) then !qcrit(numNodesRankPar)
         call this%add_nodeScalarField2save('qcrit',qcrit(:))
      end if

      !---------------  vectorScalars   -------------------------------------
      !----------------------------------------------------------------------
      if(this%save_vectorField_vel) then !u(numNodesRankPar,ndime,2)
         call this%add_nodeVectorField2save('u',u(:,:,2))
      end if
      !----------------------------------------------------------------------
      if(this%save_vectorField_gradRho) then !gradRho(numNodesRankPar,ndime)
         call this%add_nodeVectorField2save('gradRho',gradRho(:,:))
      end if
      !----------------------------------------------------------------------
      if(this%save_vectorField_curlU) then !curlU(numNodesRankPar,ndime)
         call this%add_nodeVectorField2save('curlU',curlU(:,:))
      end if
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
