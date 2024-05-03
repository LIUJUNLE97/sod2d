
module time_integ_ls

   use mod_nvtx
   use elem_convec
   use elem_diffu
   use elem_source
   use mod_solver
   use mod_entropy_viscosity
   use mod_numerical_params
   use mod_fluid_viscosity
   use mod_sgs_viscosity
   use mod_sgs_ilsa_viscosity
   use mod_bc_routines
   use mod_wall_model
   use time_integ, only : updateBuffer

   implicit none

   real(rp), allocatable, dimension(:)   :: aux_h
   real(rp), allocatable, dimension(:)   :: K2mass,K2ener
   real(rp), allocatable, dimension(:,:) :: K2mom,f_eta,Rwmles
   real(rp), allocatable, dimension(:)   :: Rmass,Rener
   real(rp), allocatable, dimension(:,:) :: Rmom
   real(rp), allocatable, dimension(:,:) :: Reta
   real(rp), allocatable, dimension(:) :: auxReta

   real(rp), allocatable, dimension(:)   :: a_i, b_i
   logical :: firstTimeStep = .true.

   contains

   subroutine init_rk4_ls_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin

      allocate(aux_h(npoin))
      !$acc enter data create(aux_h(:))

      allocate(K2mass(npoin),K2ener(npoin),auxReta(npoin),Rmass(npoin),Rener(npoin))
      !$acc enter data create(K2mass(:))
      !$acc enter data create(K2ener(:))
      !$acc enter data create(Rmass(:))
      !$acc enter data create(Rener(:))
      !$acc enter data create(auxReta(:))

      allocate(K2mom(npoin,ndime),f_eta(npoin,ndime),Reta(npoin,2),Rmom(npoin,ndime),Rwmles(npoin,ndime))
      !$acc enter data create(K2mom(:,:))
      !$acc enter data create(Rmom(:,:))
      !$acc enter data create(f_eta(:,:))
      !$acc enter data create(Reta(:,:))
      !$acc enter data create(Rwmles(:,:))

      allocate(a_i(flag_rk_ls_stages),b_i(flag_rk_ls_stages))
      !$acc enter data create(a_i(:))
      !$acc enter data create(b_i(:))

      if (flag_rk_ls_stages == 5) then
         a_i = [0.0_rp, real( -0.4178904745d0,rp) , real( -1.192151694643d0,rp), real( -1.697784692471d0,rp), real( -1.514183444257d0,rp)]
         b_i = [real( 0.1496590219993d0,rp), real( 0.3792103129999d0,rp), real( 0.8229550293869d0,rp), real( 0.6994504559488d0,rp), real( 0.1530572479681d0,rp)]
      else if (flag_rk_ls_stages == 14) then

         a_i(1)   = real( 0.0000000000000000d0  ,rp) 
         a_i(2)   = real(-0.7188012108672410d0  ,rp) 
         a_i(3)   = real(-0.7785331173421570d0  ,rp) 
         a_i(4)   = real(-0.0053282796654044d0  ,rp) 
         a_i(5)   = real(-0.8552979934029281d0  ,rp) 
         a_i(6)   = real(-3.9564138245774565d0  ,rp) 
         a_i(7)   = real(-1.5780575380587385d0  ,rp) 
         a_i(8)   = real(-2.0837094552574054d0  ,rp) 
         a_i(9)   = real(-0.7483334182761610d0  ,rp) 
         a_i(10)  = real(-0.7032861106563359d0  ,rp) 
         a_i(11)  = real( 0.0013917096117681d0  ,rp) 
         a_i(12)  = real(-0.0932075369637460d0  ,rp) 
         a_i(13)  = real(-0.9514200470875948d0  ,rp) 
         a_i(14)  = real(-7.1151571693922548d0  ,rp) 

         b_i(1)   = real( 0.0367762454319673d0 ,rp) 
         b_i(2)   = real( 0.3136296607553959d0 ,rp) 
         b_i(3)   = real( 0.1531848691869027d0 ,rp) 
         b_i(4)   = real( 0.0030097086818182d0 ,rp) 
         b_i(5)   = real( 0.3326293790646110d0 ,rp) 
         b_i(6)   = real( 0.2440251405350864d0 ,rp) 
         b_i(7)   = real( 0.3718879239592277d0 ,rp) 
         b_i(8)   = real( 0.6204126221582444d0 ,rp) 
         b_i(9)   = real( 0.1524043173028741d0 ,rp) 
         b_i(10)  = real( 0.0760894927419266d0 ,rp) 
         b_i(11)  = real( 0.0077604214040978d0 ,rp) 
         b_i(12)  = real( 0.0024647284755382d0 ,rp) 
         b_i(13)  = real( 0.0780348340049386d0 ,rp) 
         b_i(14)  = real( 5.5059777270269628d0 ,rp) 
      else
         write(1,*) "--| NOT CODED FOR RK_LS 4 stages > 5 YET!"
         stop 1
      end if
      !$acc update device(a_i(:))
      !$acc update device(b_i(:))

      call nvtxEndRange

   end subroutine init_rk4_ls_solver

   subroutine end_rk4_ls_solver()
      implicit none

      !$acc exit data delete(aux_h(:))
      deallocate(aux_h)
      
      !$acc exit data delete(K2mass(:))
      !$acc exit data delete(K2ener(:))
      !$acc exit data delete(Rmass(:))
      !$acc exit data delete(Rener(:))
      !$acc exit data delete(auxReta(:))
      deallocate(K2mass,K2ener,auxReta,Rmass,Rener)

      !$acc exit data delete(K2mom(:,:))
      !$acc exit data delete(Rmom(:,:))
      !$acc exit data delete(f_eta(:,:))
      !$acc exit data delete(Reta(:,:))
      !$acc exit data delete(Rwmles(:,:))
      deallocate(K2mom,f_eta,Reta,Rmom,Rwmles)

      !$acc exit data delete(a_i(:))
      !$acc exit data delete(b_i(:))
      deallocate(a_i,b_i)

   end subroutine end_rk4_ls_solver

         subroutine rk_4_ls_main(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,4)
            real(rp),             intent(inout) :: u(npoin,ndime,4)
            real(rp),             intent(inout) :: q(npoin,ndime,4)
            real(rp),             intent(inout) :: pr(npoin,4)
            real(rp),             intent(inout) :: E(npoin,4)
            real(rp),             intent(inout) :: Tem(npoin,2)
            real(rp),             intent(inout) :: e_int(npoin,2)
            real(rp),             intent(inout) :: eta(npoin,4)
            real(rp),             intent(inout) :: mu_fluid(npoin)
            real(rp),             intent(inout) :: csound(npoin)
            real(rp),             intent(inout) :: machno(npoin)
            real(rp),             intent(inout) :: mu_e(nelem,ngaus)
            real(rp),             intent(inout) :: mu_sgs(nelem,ngaus)
            real(rp),             intent(inout) :: kres(npoin)
            real(rp),             intent(inout) :: etot(npoin)
            real(rp),             intent(inout) :: au(npoin,ndime)
            real(rp),             intent(inout) :: ax1(npoin)
            real(rp),             intent(inout) :: ax2(npoin)
            real(rp),             intent(inout) :: ax3(npoin)
            real(rp),             intent(inout) :: mue_l(nelem,nnode)
            real(rp),             intent(in)    :: coord(npoin,ndime)
            real(rp),             intent(in)  ::  wgp(ngaus)
            integer(4),            intent(in)    :: numBoundsWM
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(in)      :: source_term(npoin,ndime)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode
            real(rp),    dimension(npoin)       :: Rrho
            real(rp)                            :: umag, rho_min, rho_avg

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value

            if(firstTimeStep .eqv. .true.) then
               firstTimeStep = .false.
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
                  end do
               end do
               !$acc end parallel loop
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta(:,1))

               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:,1))
               end if
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,1))   

               call nvtxStartRange("Entropy viscosity evaluation")
               call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,Reta(:,1),Ngp,coord,dNgp,gpvol,wgp, &
                  gamma_gas,rho(:,1),u(:,:,1),csound,Tem(:,1),eta(:,1),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
               call nvtxEndRange
            end if

            !$acc parallel loop
            do ipoin = 1,npoin
               K2mass(ipoin) = 0.0_rp
               K2ener(ipoin) = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  K2mom(ipoin,idime) = 0.0_rp
                  Rwmles(ipoin,idime) = 0.0_rp
               end do
            end do
            !$acc end parallel loop

                           !
               ! Evaluate wall models

            if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
               call nvtxStartRange("WALL MODEL")
               if(flag_type_wmles == wmles_type_reichardt) then
                  call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                  bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,rho(:,pos),walave_u(:,:),tauw,Rwmles,-1.0_rp)
               else if (flag_type_wmles == wmles_type_abl) then
                  call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                                       bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                                       rho(:,pos),walave_u(:,:),zo,tauw,Rwmles,-1.0_rp)
               end if   
               call nvtxEndRange                              
            end if

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               do idime = 1,ndime
                  call mpi_halo_atomic_update_real(Rwmles(:,idime))
               end do
               call nvtxEndRange
            end if

            !
            ! Call lumped mass matrix solver
            !
            call nvtxStartRange("Call solver")
            call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rwmles(:,:))
            call nvtxEndRange

            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")

            do istep = 1,flag_rk_ls_stages

               if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,pos),q(:,:,pos),E(:,pos),u_buffer)

               !
               ! Apply bcs after update
               !
               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer)
                  call nvtxEndRange
               end if


               !
               ! Update velocity and equations of state
               !
               call nvtxStartRange("Update u and EOS")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  umag = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     u(lpoin_w(ipoin),idime,pos) = q(lpoin_w(ipoin),idime,pos)/rho(lpoin_w(ipoin),pos)
                     umag = umag + (u(lpoin_w(ipoin),idime,pos)*u(lpoin_w(ipoin),idime,pos))
                  end do
                  e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos)) - 0.5_rp*umag
                  pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),pos)
                  aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*pr(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos)
                  Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)
               end do
               !$acc end parallel loop
               call nvtxEndRange
                
               ! Compute viscosities and diffusion
               !
               !
               ! Update viscosity if Sutherland's law is active
               !
               if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
                  call nvtxStartRange("MU_SUT")
                  call sutherland_viscosity(npoin,Tem(:,pos),mu_factor,mu_fluid)
                  call nvtxEndRange
               end if
               !
               ! Compute diffusion terms with values at current substep
               !
               call nvtxStartRange("CONVDIFFS")

               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),rho(:,pos),u(:,:,pos),Tem(:,pos),mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener,.true.,-1.0_rp)
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),q(:,:,pos),rho(:,pos),pr(:,pos),aux_h,Rmass,Rmom,Rener,.false.,-1.0_rp)               
               
               call nvtxEndRange
               !
               ! Call source term if applicable
               !
               if(present(source_term)) then
                  call nvtxStartRange("SOURCE TERM")
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),source_term,Rmom,-1.0_rp)
                  call nvtxEndRange
               end if

               !TESTING NEW LOCATION FOR MPICOMMS
               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(Rmass(:))
                  call mpi_halo_atomic_update_real(Rener(:))
                  do idime = 1,ndime
                     call mpi_halo_atomic_update_real(Rmom(:,idime))
                  end do
                  call nvtxEndRange
               end if

               !
               ! Call lumped mass matrix solver
               !
               call nvtxStartRange("Call solver")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass(:))
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener(:))
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom(:,:))
               call nvtxEndRange
               !
               ! Accumulate the residuals
               !
               call nvtxStartRange("Accumulate residuals")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  K2mass(lpoin_w(ipoin)) = a_i(istep)*K2mass(lpoin_w(ipoin)) + Rmass(lpoin_w(ipoin))*dt
                  rho(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos) + b_i(istep)*K2mass(lpoin_w(ipoin))
                  K2ener(lpoin_w(ipoin)) = a_i(istep)*K2ener(lpoin_w(ipoin)) + Rener(lpoin_w(ipoin))*dt
                  E(lpoin_w(ipoin),pos) = E(lpoin_w(ipoin),pos) + b_i(istep)*K2ener(lpoin_w(ipoin))
                  !$acc loop seq
                  do idime = 1,ndime
                     K2mom(lpoin_w(ipoin),idime) = a_i(istep)*K2mom(lpoin_w(ipoin),idime) + (Rmom(lpoin_w(ipoin),idime)+Rwmles(lpoin_w(ipoin),idime))*dt
                     q(lpoin_w(ipoin),idime,pos) = q(lpoin_w(ipoin),idime,pos) + b_i(istep)*K2mom(lpoin_w(ipoin),idime)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange
               !call limit_rho(nelem,npoin,connec,rho(:,pos),epsilon(umag))
            end do
            call nvtxEndRange



            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,pos),q(:,:,pos),E(:,pos),u_buffer)

            !
            ! Apply bcs after update
            !
            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer)
               call nvtxEndRange
            end if

            !
            ! Update velocity and equations of state
            !
            call nvtxStartRange("Update u and EOS")

            !$acc parallel loop
            do ipoin = 1,npoin_w
               umag = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  u(lpoin_w(ipoin),idime,pos) = q(lpoin_w(ipoin),idime,pos)/rho(lpoin_w(ipoin),pos)
                  umag = umag + u(lpoin_w(ipoin),idime,pos)**2
               end do
               e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))- &
                  0.5_rp*umag
               pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),pos)
               csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))
               umag = sqrt(umag)
               machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
               Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)

               eta(lpoin_w(ipoin),1) = eta(lpoin_w(ipoin),2)
               eta(lpoin_w(ipoin),2) = (rho(lpoin_w(ipoin),2)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)**gamma_gas))
               !$acc loop seq
               do idime = 1,ndime
                  f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,2)*eta(lpoin_w(ipoin),2)
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,2),u(:,:,2),Reta(:,2))

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(Reta(:,2))
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,2))

            call nvtxStartRange("Entropy residual")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               auxReta(lpoin_w(ipoin)) = (1.5_rp*Reta(lpoin_w(ipoin),2)-0.5_rp*Reta(lpoin_w(ipoin),1)) !+ &
                                              !(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
               Reta(lpoin_w(ipoin),1) = Reta(lpoin_w(ipoin),2)            
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,auxReta)
               call nvtxEndRange
            end if

            !
            ! If using Sutherland viscosity model:
            !
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call nvtxStartRange("Sutherland viscosity")
               call sutherland_viscosity(npoin,Tem(:,pos),mu_factor,mu_fluid)
               call nvtxEndRange
            end if

            !
            ! Compute entropy viscosity
            !
            call nvtxStartRange("Entropy viscosity evaluation")
            call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
            call nvtxEndRange

            !
            ! Compute subgrid viscosity if active
            !
            if(flag_les == 1) then
               call nvtxStartRange("MU_SGS")
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,pos),u(:,:,pos),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3,mue_l) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),u(:,:,pos),Ml,mu_sgs,mue_l)
               end if
               call nvtxEndRange
            end if

         end subroutine rk_4_ls_main

         subroutine limit_rho(nelem,npoin,connec,rho,eps)  

            implicit none

            integer(4),           intent(in)    :: nelem,npoin,connec(nelem,nnode)
            real(rp),             intent(in)    :: eps
            real(rp),             intent(inout) :: rho(npoin)
            real(rp)                            :: rho_min, rho_avg, fact
            integer(4)                          :: ielem, inode, ipoin

            !$acc parallel loop gang 
            do ielem = 1,nelem
                rho_min = real(1e6, rp)
                rho_avg = 0.0_rp
                !$acc loop seq
                do inode = 1,nnode
                        ipoin = connec(ielem,inode)
                        rho_min = min(rho_min, rho(ipoin))
                        rho_avg = rho_avg + rho(ipoin)
                end do
                rho_avg = rho_avg/real(nnode,rp)
                fact = min(1.0_rp,(rho_avg-eps)/(rho_avg-rho_min))
                !$acc loop vector
                do inode = 1,nnode
                        ipoin = connec(ielem,inode)
                        rho(ipoin) = rho_avg + fact*(rho(ipoin)-rho_avg)
                end do
             end do
            !$acc end parallel loop

            end subroutine limit_rho

      end module time_integ_ls
