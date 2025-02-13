module SupersonicNozzle_mod
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
   use mod_operators
   use mod_hdf5
   use CFDSolverPeriodicWithBoundaries_mod
   implicit none
   private

   real(rp),  allocatable, dimension(:,:) :: gradXVel, gradYVel,gradZVel

   type, public, extends(CFDSolverPeriodicWithBoundaries) :: SupersonicNozzle

      real(rp) , public  :: rho0,to, po, mu, Re, delta, M, vo  , Ccoles, D0, utau, RetRef, Uinfp, dia, thsl, aL, yL, deltaBL

   contains
      procedure, public :: fillBCTypes           =>SupersonicNozzle_fill_BC_Types
      procedure, public :: initializeParameters  =>SupersonicNozzle_initializeParameters
      procedure, public :: evalInitialConditions =>SupersonicNozzle_evalInitialConditions
      procedure, public :: initialBuffer =>SupersonicNozzle_initialBuffer
      procedure, public :: initializeSourceTerms => SupersonicNozzle_initializeSourceTerms
      procedure, public :: afterDt => SupersonicNozzle_afterDt
      procedure, public :: beforeTimeIteration => SupersonicNozzle_beforeTimeIteration


   end type SupersonicNozzle
contains

   subroutine SupersonicNozzle_beforeTimeIteration(this)
      class(SupersonicNozzle), intent(inout) :: this
      integer(4)                 :: iboun,bcode,ipbou,iBoundNode,iNodeL

#if 1
!0 to reint very experimental, 1 for the first time step
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(maskMapped(iNodeL)==1) then
            u(iNodeL,1,2) = 0.0_rp
            u(iNodeL,2,2) = 0.0_rp
            u(iNodeL,3,2) = 0.0_rp
            if(coordPar(iNodeL,1) .lt. 3.0_rp) then ! not ideal
               u(iNodeL,1,2) = 1.0_rp
               u(iNodeL,2,2) = 0.0_rp
               u(iNodeL,3,2) = 0.0_rp
            end if
         
            pr(iNodeL,2) = this%po
            rho(iNodeL,2) = this%rho0
            e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
            Tem(iNodeL,2) = this%to
            E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
            q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
            csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
            eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
            !eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(rho(iNodeL,2)*e_int(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
            machno(iNodeL) = sqrt(dot_product(u(iNodeL,:,2),u(iNodeL,:,2)))/csound(iNodeL)

            q(iNodeL,1:ndime,3) = q(iNodeL,1:ndime,2)
            rho(iNodeL,3) = rho(iNodeL,2)
            E(iNodeL,3) =  E(iNodeL,2)
            eta(iNodeL,3) = eta(iNodeL,2) 
         end if
      end do
      !$acc end parallel loop
#endif

   end subroutine SupersonicNozzle_beforeTimeIteration


   subroutine SupersonicNozzle_initializeSourceTerms(this)
      class(SupersonicNozzle), intent(inout) :: this
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3
      integer(4) :: k,j,ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime,iNodeL,iNodeL2,bcode,isoII, isoJJ, isoKK,type_ijk
      real(rp) :: dy,fx1,fx2,xp
      real(rp) :: mul , yc
      real(rp)  :: gradIsoV(ndime),gradIsoU(ndime)
      real(rp)  :: gradV(ndime),vl(nnode),fact,targ,gradU(ndime),ul(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
      real(rp) :: yp,up,vp,wp,eta_y,f_y,f_prim_y

      
      allocate(source_term(numNodesRankPar,ndime+2),gradXVel(numNodesRankPar,3),gradYVel(numNodesRankPar,3),gradZVel(numNodesRankPar,3))
      !$acc enter data create(source_term(:,:),gradXVel(:,:),gradYVel(:,:),gradZVel(:,:))

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         source_term(iNodeL,1) = 0.0_rp
         source_term(iNodeL,2) = 0.0_rp
         source_term(iNodeL,3) = 0.0_rp
         source_term(iNodeL,4) = 0.0_rp
         source_term(iNodeL,5) = 0.0_rp
      end do
      !$acc end parallel loop
   end subroutine SupersonicNozzle_initializeSourceTerms


   subroutine SupersonicNozzle_afterDt(this,istep)
      class(SupersonicNozzle), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      real(rp) :: cd, lx, ly, xmin, xmax, f1, f2, f3
      integer(4) :: k,j,ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime,iNodeL,iNodeL2,bcode,isoII, isoJJ, isoKK,type_ijk
      real(rp) :: dy,fx1,fx2,xp
      real(rp) :: mul , yc
      real(rp)  :: gradIsoV(ndime),gradIsoU(ndime)
      real(rp)  :: gradV(ndime),vl(nnode),fact,targ,gradU(ndime),ul(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
      real(rp) :: yp,up,vp,wp,eta_y,f_y,f_prim_y,sig = 1.0_rp

      call eval_gradient(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,connecParWork,workingNodesPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,q(:,1,2),gradXVel,.true.)
      call eval_gradient(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,connecParWork,workingNodesPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,q(:,2,2),gradYVel,.true.)
      call eval_gradient(numElemsRankPar,numNodesRankPar,numWorkingNodesRankPar,connecParWork,workingNodesPar,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,q(:,3,2),gradZVel,.true.)

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         yp = coordPar(iNodeL,2)
         if(maskMapped(iNodeL)==0) then
            yp = 0.5_rp*this%delta-abs(coordPar(iNodeL,2))

            if(coordPar(iNodeL,2) .gt. 0) sig = -1.0_rp

            source_term(iNodeL,3) = (this%D0)*yp*gradXVel(iNodeL,2)*sig
            source_term(iNodeL,4) = (this%D0)*yp*gradYVel(iNodeL,2)*sig
            source_term(iNodeL,5) = (this%D0)*yp*gradZVel(iNodeL,2)*sig
         end if
      end do
      !$acc end parallel loop
   end subroutine SupersonicNozzle_afterDt


   subroutine SupersonicNozzle_fill_BC_Types(this)
      class(SupersonicNozzle), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine SupersonicNozzle_fill_BC_Types

   subroutine SupersonicNozzle_initializeParameters(this)
      use json_module
      implicit none
      class(SupersonicNozzle), intent(inout) :: this
      logical :: found, found_aux = .false.
      type(json_file) :: json
      character(len=:) , allocatable :: value
      real(rp) :: mul, mur

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

      call json%get("saveInitialField",this%saveInitialField, found,.true.); call this%checkFound(found,found_aux)

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

      call json%get("doTimerAnalysis",this%doTimerAnalysis, found,.false.)

      ! numerical params
      call json%get("flag_les",flag_les, found,1); call this%checkFound(found,found_aux)
      call json%get("flag_implicit",flag_implicit, found,0); call this%checkFound(found,found_aux)
      call json%get("flag_imex_stages",flag_imex_stages, found,4); call this%checkFound(found,found_aux)
      call json%get("maxIter",maxIter, found,20); call this%checkFound(found,found_aux)
      call json%get("tol",tol, found,0.001d0); call this%checkFound(found,found_aux)
      call json%get("flag_high_mach",flag_high_mach, found,.true.); call this%checkFound(found,found_aux)
      call json%get("flag_les_ilsa",flag_les_ilsa, found,0); call this%checkFound(found,found_aux)
      call json%get("stau",stau, found,0.022_rp); call this%checkFound(found,found_aux)
      call json%get("T_ilsa",T_ilsa, found,1.0_rp); call this%checkFound(found,found_aux)
       
      call json%get("period_walave",period_walave, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("flag_type_wmles",flag_type_wmles, found,1); call this%checkFound(found,found_aux)
      call json%get("wmles_walex",wmles_walex, found,0.1_rp); !optional depending of the model

      call json%get("cfl_conv",this%cfl_conv, found,0.95_rp); call this%checkFound(found,found_aux)
      call json%get("cfl_diff",this%cfl_diff, found,0.95_rp); call this%checkFound(found,found_aux)

      call json%get("flag_rk_ls",flag_rk_ls, found,.true.); 
      call json%get("flag_rk_ls_stages",flag_rk_ls_stages, found,5); 
      call json%get("flag_high_mach",flag_high_mach, found,.true.); call this%checkFound(found,found_aux)


      call json%get("Cp",this%Cp, found,1004.0_rp); call this%checkFound(found,found_aux)
      call json%get("Prt",this%Prt, found,0.71_rp); call this%checkFound(found,found_aux)
      call json%get("v0",this%vo, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("M",this%M, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("delta",this%delta, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("rho0",this%rho0, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("Re",this%Re, found,90000.0_rp); call this%checkFound(found,found_aux)
      call json%get("gamma_gas",this%gamma_gas, found,1.4_rp); call this%checkFound(found,found_aux)
      call json%get("RetRef",this%RetRef, found,1000.0_rp); call this%checkFound(found,found_aux)
      call json%get("Uinfp",this%Uinfp, found,22.0_rp); call this%checkFound(found,found_aux)
      call json%get("deltaBL",this%deltaBL, found,0.1_rp); call this%checkFound(found,found_aux)

      call json%get("flag_lps_stab",flag_lps_stab, found,.true.); call this%checkFound(found,found_aux)


      this%mu    = (this%rho0*this%delta*this%vo)/this%Re
      this%Rgas = this%Cp*(this%gamma_gas-1.0_rp)/this%gamma_gas
      this%to = this%vo*this%vo/(this%gamma_gas*this%Rgas*this%M*this%M)
      this%po = this%rho0*this%Rgas*this%to
      mur = 0.000001458_rp*(this%to**1.50_rp)/(this%to+110.40_rp)
      flag_mu_factor = this%mu/mur

      this%Ccoles = 4.05_rp

      this%D0     = (this%vo/this%Uinfp)/(this%deltaBL*this%Ccoles)
      this%dia    = 0.5_rp*this%deltaBL
      this%thsl   = 54.0_rp*this%mu/this%vo
      
      this%aL     = 1.0_rp
      this%yL     = 3.5_rp

      nscbc_c_inf = sqrt(this%gamma_gas*this%Rgas*this%to)
      nscbc_rho_inf = this%rho0
      nscbc_p_inf = this%po
      nscbc_Rgas_inf = this%Rgas
      nscbc_gamma_inf = this%gamma_gas
      nscbc_T_C = this%to
      nscbc_u_inf = this%vo
      nscbc_sign_ux = -1.0_rp

      call this%readJSONBuffer()

      call json%destroy()

      if(found_aux .and.mpi_rank .eq. 0) write(111,*) 'WARNING! JSON file missing a parameter, overwrtting with the default value'
   end subroutine SupersonicNozzle_initializeParameters

   subroutine SupersonicNozzle_evalInitialConditions(this)
      class(SupersonicNozzle), intent(inout) :: this
      integer(8) :: matGidSrlOrdered(numNodesRankPar,2)
      integer(4)   :: iLine,iNodeGSrl,auxCnt, iNodeL
      logical :: readFiles
      real(rp) :: velo, rti(3), yp,velo_aux1

      call nvtxStartRange("SupersonicNozzle Init")
      call order_matrix_globalIdSrl(numNodesRankPar,globalIdSrl,matGidSrlOrdered)
      auxCnt = 1          
  
      serialLoop : do iLine = 1,totalNumNodesSrl
         call random_number(rti)
         if(iLine.eq.matGidSrlOrdered(auxCnt,2)) then
            iNodeL = matGidSrlOrdered(auxCnt,1)
            auxCnt=auxCnt+1
            if(coordPar(iNodeL,1) .lt. 0.0_rp) then ! not ideal
               
               yp = 0.5_rp*this%delta-abs(coordPar(iNodeL,2))
               
               velo = 0.5_rp*this%vo*(1.0_rp + tanh(0.5_rp*this%dia/this%thsl*(1.0_rp-yp/this%dia)))
               u(iNodeL,1,2) = velo*(1.0_rp + 0.1_rp*(rti(1) -0.5_rp))
               u(iNodeL,2,2) = velo*(0.1_rp*(rti(2) -0.5_rp))
               u(iNodeL,3,2) = velo*(0.1_rp*(rti(3) -0.5_rp))
            else 
               u(iNodeL,1,2) = 0.0_rp
               u(iNodeL,2,2) = 0.0_rp
               u(iNodeL,3,2) = 0.0_rp
            end if
         end if
         if(auxCnt.gt.numNodesRankPar) then
            exit serialLoop
         end if
      end do serialLoop
      !$acc update device(u(:,:,:))

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         pr(iNodeL,2) = this%po
         rho(iNodeL,2) = this%rho0
         e_int(iNodeL,2) = pr(iNodeL,2)/(rho(iNodeL,2)*(this%gamma_gas-1.0_rp))
         Tem(iNodeL,2) = this%to
         E(iNodeL,2) = rho(iNodeL,2)*(0.5_rp*dot_product(u(iNodeL,:,2),u(iNodeL,:,2))+e_int(iNodeL,2))
         q(iNodeL,1:ndime,2) = rho(iNodeL,2)*u(iNodeL,1:ndime,2)
         csound(iNodeL) = sqrt(this%gamma_gas*pr(iNodeL,2)/rho(iNodeL,2))
         eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(pr(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         !eta(iNodeL,2) = (rho(iNodeL,2)/(this%gamma_gas-1.0_rp))*log(rho(iNodeL,2)*e_int(iNodeL,2)/(rho(iNodeL,2)**this%gamma_gas))
         machno(iNodeL) = sqrt(dot_product(u(iNodeL,:,2),u(iNodeL,:,2)))/csound(iNodeL)

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

   end subroutine SupersonicNozzle_evalInitialConditions


   subroutine SupersonicNozzle_initialBuffer(this)
      class(SupersonicNozzle), intent(inout) :: this
      integer(4) :: iNodeL, bcode

      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         if(maskMapped(iNodeL)==0) then
            u_buffer(iNodeL,1) = 1.0_rp
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
         else
            u_buffer(iNodeL,1) = 0.0_rp
            u_buffer(iNodeL,2) = 0.0_rp
            u_buffer(iNodeL,3) = 0.0_rp  
         end if
         
         u_mapped(iNodeL,1) = 1.0_rp
         u_mapped(iNodeL,2) = 0.0_rp
         u_mapped(iNodeL,3) = 0.0_rp  

      end do
      !$acc end parallel loop

   end subroutine SupersonicNozzle_initialBuffer

end module SupersonicNozzle_mod
