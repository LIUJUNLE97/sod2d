
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

   real(rp), allocatable, dimension(:)   :: K1mass,K1ener,K1eta
   real(rp), allocatable, dimension(:,:) :: K1mom
   real(rp), allocatable, dimension(:)   :: aux_pr, aux_Tem, aux_e_int,aux_h
   real(rp), allocatable, dimension(:,:) :: aux_u
   real(rp), allocatable, dimension(:)   :: K2mass,K2ener,K2eta
   real(rp), allocatable, dimension(:,:) :: K2mom,f_eta

   real(rp), allocatable, dimension(:)   :: a_i, b_i

   contains

   subroutine init_rk4_ls_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin

      allocate(K1mass(npoin),K1ener(npoin),K1eta(npoin),K2mom(npoin,ndime))
      !$acc enter data create(K1mass(:))
      !$acc enter data create(K1ener(:))
      !$acc enter data create(K1eta(:))
      !$acc enter data create(K2mom(:,:))

      allocate(aux_pr(npoin),aux_Tem(npoin),aux_e_int(npoin),aux_h(npoin))
      !$acc enter data create(aux_pr(:))
      !$acc enter data create(aux_Tem(:))
      !$acc enter data create(aux_e_int(:))
      !$acc enter data create(aux_h(:))

      allocate(aux_u(npoin,ndime))
      !$acc enter data create(aux_u(:,:))

      allocate(K2mass(npoin),K2ener(npoin),K2eta(npoin))
      !$acc enter data create(K2mass(:))
      !$acc enter data create(K2ener(:))
      !$acc enter data create(K2eta(:))

      allocate(K2mom(npoin,ndime),f_eta(npoin,ndime))
      !$acc enter data create(K2mom(:,:))
      !$acc enter data create(f_eta(:,:))

      allocate(a_i(flag_rk_ls_stages),b_i(flag_rk_ls_stages))
      !$acc enter data create(a_i(:))
      !$acc enter data create(b_i(:))

      if (flag_rk_ls_stages == 5) then
         a_i = [0.0_rp, real( -0.7274361725534d0,rp) , real( -1.192151694643d0,rp), real( -1.697784692471d0,rp), real( -1.514183444257d0,rp)]
         b_i = [real( 0.1496590219993d0,rp), real( 0.3792103129999d0,rp), real( 0.8229550293869d0,rp), real( 0.6994504559488d0,rp), real( 0.1530572479681d0,rp)]
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

      !$acc exit data delete(K1mass(:))
      !$acc exit data delete(K1ener(:))
      !$acc exit data delete(K1eta(:))
      !$acc exit data delete(K1mom(:,:))
      deallocate(K1mass,K1ener,K1eta,K1mom)

      !$acc exit data delete(aux_pr(:))
      !$acc exit data delete(aux_Tem(:))
      !$acc exit data delete(aux_e_int(:))
      !$acc exit data delete(aux_h(:))
      deallocate(aux_pr,aux_Tem,aux_e_int,aux_h)
      
      !$acc exit data delete(aux_u(:,:))
      deallocate(aux_u)

      !$acc exit data delete(K2mass(:))
      !$acc exit data delete(K2ener(:))
      !$acc exit data delete(K2eta(:))
      deallocate(K2mass,K2ener,K2eta)

      !$acc exit data delete(K2mom(:,:))
      !$acc exit data delete(f_eta(:,:))
      deallocate(K2mom,f_eta)

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
            real(rp)                            :: umag


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value

            !
            ! Initialize variables to zero
            !
            call nvtxStartRange("Initialize variables")
            !$acc parallel loop
            do ipoin = 1,npoin            
               aux_u(ipoin,1:ndime) = 0.0_rp
               aux_pr(ipoin) = 0.0_rp
               aux_Tem(ipoin) = 0.0_rp
               aux_e_int(ipoin) = 0.0_rp
               K1mass(ipoin) = rho(ipoin,pos)
               K1mom(ipoin,1:ndime) = q(ipoin,1:ndime,pos)
               K1ener(ipoin) =  E(ipoin,pos)
               K1eta(ipoin) = eta(ipoin,pos)
               K2mass(ipoin) = 0.0_rp
               K2ener(ipoin) = 0.0_rp
               K2eta(ipoin) = 0.0_rp
               K2mom(ipoin,1:ndime) = 0.0_rp
            end do
            !$acc end parallel loop
            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")

            do istep = 1,flag_rk_ls_stages
               !
               ! Update velocity and equations of state
               !
               call nvtxStartRange("Update u and EOS")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  umag = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_u(lpoin_w(ipoin),idime) = K1mom(lpoin_w(ipoin),idime)/K1mass(lpoin_w(ipoin))
                     umag = umag + (aux_u(lpoin_w(ipoin),idime)*aux_u(lpoin_w(ipoin),idime))
                  end do
                  aux_e_int(lpoin_w(ipoin)) = (K1ener(lpoin_w(ipoin))/K1mass(lpoin_w(ipoin)))-0.5_rp*umag
                  aux_pr(lpoin_w(ipoin)) = K1mass(lpoin_w(ipoin))*(gamma_gas-1.0_rp)*aux_e_int(lpoin_w(ipoin))
                  aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*aux_pr(lpoin_w(ipoin))/K1mass(lpoin_w(ipoin))
                  aux_Tem(lpoin_w(ipoin)) = aux_pr(lpoin_w(ipoin))/(K1mass(lpoin_w(ipoin))*Rgas)
                  K1eta(lpoin_w(ipoin)) = (K1mass(lpoin_w(ipoin))/(gamma_gas-1.0_rp))* &
                  log(aux_pr(lpoin_w(ipoin))/(K1mass(lpoin_w(ipoin))**gamma_gas))
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime) = aux_u(lpoin_w(ipoin),idime)*K1eta(lpoin_w(ipoin))
                  end do

                  K2mass(lpoin_w(ipoin)) = a_i(istep)*K2mass(lpoin_w(ipoin))
                  K2ener(lpoin_w(ipoin)) = a_i(istep)*K2ener(lpoin_w(ipoin)) 
                  K2eta(lpoin_w(ipoin)) = a_i(istep)*K2eta(lpoin_w(ipoin))
                  !$acc loop seq
                  do idime = 1,ndime
                     K2mom(lpoin_w(ipoin),idime) = a_i(istep)*K2mom(lpoin_w(ipoin),idime)
                  end do

               end do
               !$acc end parallel loop
               call nvtxEndRange

               call nvtxStartRange("Entropy convection")
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,K1eta(:),aux_u(:,:),K2eta(:))
               call nvtxEndRange

               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(K2eta(:))
                  call nvtxEndRange
               end if

               if (noBoundaries .eqv. .false.) then
                  call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,K2eta)
               end if

               call nvtxStartRange("Lumped solver")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,K2eta)
               call nvtxEndRange
               !
               ! Compute viscosities and diffusion
               !
               !
               ! Update viscosity if Sutherland's law is active
               !
               if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
                  call nvtxStartRange("MU_SUT")
                  call sutherland_viscosity(npoin,aux_Tem,mu_factor,mu_fluid)
                  call nvtxEndRange
               end if
               !
               ! Compute diffusion terms with values at current substep
               !
               call nvtxStartRange("DIFFUSIONS")
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,K1mass,K1mass,aux_u,aux_Tem,mu_fluid,mu_e,mu_sgs,Ml,K2mass,K2mom,K2ener,.false.,-1.0_rp)
               call nvtxEndRange
               !
               ! Call source term if applicable
               !
               if(present(source_term)) then
                  call nvtxStartRange("SOURCE TERM")
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,K2mom)
                  call nvtxEndRange
               end if
               !
               ! Evaluate wall models

               if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  call nvtxStartRange("WALL MODEL")
                  if(flag_walave) then
                     if(flag_type_wmles == wmles_type_reichardt) then
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,K1mass(:),walave_u(:,:),tauw,K2mom)
                     else if (flag_type_wmles == wmles_type_abl) then
                        call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                                             bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                                             K1mass(:),walave_u(:,:),zo,tauw,K2mom)
                     end if   
                  else
                     if(flag_type_wmles == wmles_type_reichardt) then
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,K1mass(:),aux_u(:,:),tauw,K2mom)
                     else
                        write(1,*) "--| Only Reichardt wall model can work without time filtering!"
                        stop 1
                     end if
                  end if
                  call nvtxEndRange
               end if
               !
               !
               ! Compute convective terms
               !
               call nvtxStartRange("CONVECTIONS")
               if(flag_total_enthalpy .eqv. .true.) then
                  call full_convec_ijk_H(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,K1mom,K1mass,aux_pr,K1ener,K2mass(:),K2mom(:,:),K2ener(:),.false.,-1.0_rp)
               else 
                  call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,K1mom,K1mass,aux_pr,aux_h,K2mass(:),K2mom(:,:),K2ener(:), .false.,-1.0_rp)
               end if
               call nvtxEndRange

               !TESTING NEW LOCATION FOR MPICOMMS
               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(K2mass(:))
                  call mpi_halo_atomic_update_real(K2ener(:))
                  do idime = 1,ndime
                     call mpi_halo_atomic_update_real(K2mom(:,idime))
                  end do
                  call nvtxEndRange
               end if

               !
               ! Call lumped mass matrix solver
               !
               call nvtxStartRange("Call solver")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,K2mass(:))
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,K2ener(:))
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,K2mom(:,:))
               call nvtxEndRange
               !
               ! Accumulate the residuals
               !
               call nvtxStartRange("Accumulate residuals")
               !$acc parallel loop
               do ipoin = 1,npoin
                  K2mass(ipoin) = K2mass(ipoin)*dt
                  K2ener(ipoin) = K2ener(ipoin)*dt
                  !$acc loop seq
                  do idime = 1,ndime
                     K2mom(ipoin,idime) = K2mom(ipoin,idime)*dt
                  end do

                  K1mass(ipoin) = K1mass(ipoin) + b_i(istep)*K2mass(ipoin)
                  K1ener(ipoin) = K1ener(ipoin) + b_i(istep)*K2ener(ipoin)
                  K1eta(ipoin)  = K1eta(ipoin)  + b_i(istep)*K2eta(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     K1mom(ipoin,idime) = K1mom(ipoin,idime) + b_i(istep)*K2mom(ipoin,idime)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange
            end do
            call nvtxEndRange
            !
            ! RK update to variables
            !
            call nvtxStartRange("RK_UPDATE")
            !$acc parallel loop
            do ipoin = 1,npoin
               rho(ipoin,pos) = K1mass(ipoin)
               E(ipoin,pos) = K1ener(ipoin)
               !$acc loop seq
               do idime = 1,ndime
                  q(ipoin,idime,pos) = K1mom(ipoin,idime)
               end do
            end do
            !$acc end parallel loop
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
               eta(lpoin_w(ipoin),pos) = (rho(lpoin_w(ipoin),pos)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)**gamma_gas))
            end do
            !$acc end parallel loop
            call nvtxEndRange

            call nvtxStartRange("Entropy residual")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               K2eta(lpoin_w(ipoin)) =  -b_i(flag_rk_ls_stages)*K2eta(lpoin_w(ipoin))-(eta(lpoin_w(ipoin),2)-K1eta(lpoin_w(ipoin)))
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,K2eta)
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

            call nvtxStartRange("Entropy viscosity evaluation")
            !
            ! Compute entropy viscosity
            !
            call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,K2eta,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,pos),u(:,:,pos),csound,Tem(:,pos),eta(:,pos),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
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

      end module time_integ_ls
