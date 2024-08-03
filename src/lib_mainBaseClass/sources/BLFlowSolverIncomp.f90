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

   real(rp), public :: eta_b(45), f(45), f_prim(45)

   type, public, extends(CFDSolverPeriodicWithBoundariesIncomp) :: BLFlowSolverIncomp

      real(rp) , public  ::   d0, U0, rho0, Red0, mu, vmax_SB, xc_SB, eps_SB, psi_SB, xmin_tripping, lx_tripping, ly_tripping

   contains
      procedure, public :: fillBCTypes           => BLFlowSolverIncomp_fill_BC_Types
      procedure, public :: initializeParameters  => BLFlowSolverIncomp_initializeParameters
      procedure, public :: evalInitialConditions => BLFlowSolverIncomp_evalInitialConditions
      procedure, public :: initialBuffer => BLFlowSolverIncomp_initialBuffer
      procedure, public :: fillBlasius => BLFlowSolverIncomp_fillBlasius
      procedure, public :: afterDt => BLFlowSolverIncomp_afterDt
      procedure, public :: beforeTimeIteration => BLFlowSolverIncomp_beforeTimeIteration
      end type BLFlowSolverIncomp
contains

   subroutine BLFlowSolverIncomp_beforeTimeIteration(this)
      class(BLFlowSolverIncomp), intent(inout) :: this
      integer(4) :: bcode,iNodeL

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

   end subroutine BLFlowSolverIncomp_beforeTimeIteration

   subroutine BLFlowSolverIncomp_afterDt(this,istep)
      class(BLFlowSolverIncomp), intent(inout) :: this
      integer(4)              , intent(in)   :: istep
      real(rp) :: cd, lx, ly, xmin, xmax
      integer(4) :: ielem,iCen,inode,igaus, isoI, isoJ, isoK,ii,jdime,idime,iNodeL,bcode,isoII, isoJJ, isoKK,type_ijk
      real(rp) :: xp, yp, yc
      real(rp)  :: gradIsoV(ndime)
      real(rp)  :: gradV(ndime),vl(nnode)
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

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

   end subroutine BLFlowSolverIncomp_afterDt

   subroutine BLFlowSolverIncomp_fill_BC_Types(this)
      class(BLFlowSolverIncomp), intent(inout) :: this

      call this%readJSONBCTypes()

   end subroutine BLFlowSolverIncomp_fill_BC_Types

   subroutine BLFlowSolverIncomp_fillBlasius(this)
      class(BLFlowSolverIncomp), intent(inout) :: this

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

  end subroutine BLFlowSolverIncomp_fillBlasius

   subroutine BLFlowSolverIncomp_initializeParameters(this)
      use json_module
      implicit none
      class(BLFlowSolverIncomp), intent(inout) :: this
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
      call json%get("maxPhysTime",this%maxPhysTime, found,1000.0_rp); call this%checkFound(found,found_aux)

      call json%get("save_logFile_first",this%save_logFile_first, found, 1); call this%checkFound(found,found_aux)
      call json%get("save_logFile_step",this%save_logFile_step, found, 10); call this%checkFound(found,found_aux)

      call json%get("save_resultsFile_first",this%save_resultsFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_resultsFile_step" ,this%save_resultsFile_step, found,10000); call this%checkFound(found,found_aux)

      call json%get("save_restartFile_first",this%save_restartFile_first, found,1); call this%checkFound(found,found_aux)
      call json%get("save_restartFile_step" ,this%save_restartFile_step, found,10000); call this%checkFound(found,found_aux)

      call json%get("loadRestartFile" ,this%loadRestartFile, found, .false.); call this%checkFound(found,found_aux)
      call json%get("restartFile_to_load" ,this%restartFile_to_load, found,1); call this%checkFound(found,found_aux)

      call json%get("continue_oldLogs" ,this%continue_oldLogs, found, .false.); call this%checkFound(found,found_aux)

      call json%get("initial_avgTime",this%initial_avgTime, found, 0.0_rp); call this%checkFound(found,found_aux)
      call json%get("saveAvgFile" ,this%saveAvgFile, found, .true.); call this%checkFound(found,found_aux)
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
      call json%get("vmax_SB",this%vmax_SB, found,0.4_rp); call this%checkFound(found,found_aux)
      call json%get("xc_SB",this%xc_SB, found,306.640625_rp); call this%checkFound(found,found_aux)
      call json%get("eps_SB",this%eps_SB, found,110.485435_rp); call this%checkFound(found,found_aux)
      call json%get("psi_SB",this%psi_SB, found,0.95_rp); call this%checkFound(found,found_aux)

      ! Tripping region
      call json%get("xmin_tripping",this%xmin_tripping, found,0.0_rp); call this%checkFound(found,found_aux)
      call json%get("lx_tripping",this%lx_tripping, found,1.0_rp); call this%checkFound(found,found_aux)
      call json%get("ly_tripping",this%ly_tripping, found,1.0_rp); call this%checkFound(found,found_aux)

      ! witness points
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

   end subroutine BLFlowSolverIncomp_initialBuffer

   subroutine BLFlowSolverIncomp_evalInitialConditions(this)
      class(BLFlowSolverIncomp), intent(inout) :: this
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
   end subroutine BLFlowSolverIncomp_evalInitialConditions

end module BLFlowSolverIncomp_mod
