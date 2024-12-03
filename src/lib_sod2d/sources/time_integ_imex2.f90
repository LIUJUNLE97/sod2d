module time_integ_imex2

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
   use mod_operators
   use mod_solver
   use mod_solver_imex
   use time_integ, only :  updateBuffer   
   use elem_stab, only : comp_tau

   implicit none

   real(rp), allocatable, dimension(:,:,:) :: Rmom
   real(rp), allocatable, dimension(:,:)   :: Rsource,Rwmles,Rener,Rmass
   real(rp), allocatable, dimension(:,:)   :: f_eta,Reta,Rflux
   real(rp), allocatable, dimension(:)     :: auxReta,aux_h
   real(rp), allocatable, dimension(:)     :: beta, alpha
   real(rp) :: gamma0

   contains
   subroutine init_imex2_solver(npoin)

      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      call nvtxStartRange("Init Incomp solver")

      allocate(Rmom(npoin,ndime,3))
      !$acc enter data create(Rmom(:,:,:))

      allocate(aux_h(npoin),Rsource(npoin,ndime+2),Rwmles(npoin,ndime),Rmass(npoin,3),Rener(npoin,3))
      !$acc enter data create(aux_h(:),Rsource(:,:),Rwmles(:,:),Rmass(:,:),Rener(:,:))

      allocate(auxReta(npoin),f_eta(npoin,ndime),Reta(npoin,3))
      !$acc enter data create(auxReta(:),f_eta(:,:),Reta(:,:))

      allocate(alpha(3),beta(3))
      !$acc enter data create(alpha(:),beta(:))
      call nvtxEndRange

   end subroutine init_imex2_solver

   subroutine end_imex2_solver()
      implicit none

   end subroutine end_imex2_solver
 
         subroutine imex2_main(igtime,iltime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,&
            ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
            ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
            listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)  ! Optional args

            implicit none

            logical,              intent(in)    :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,save_logFile_next,iltime
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin),lelpn(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow,maskMapped(npoin)
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: mue_l(nelem,nnode)
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
            real(rp),             intent(in)    :: coord(npoin,ndime)
            real(rp),             intent(in)  ::  wgp(ngaus)
            integer(4),            intent(in)    :: numBoundsWM
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(in)      :: source_term(npoin,ndime+2)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: istep,ipoin,idime,icode,iPer,ipoin_w
            real(rp)                            :: umag

#if 0
            call nvtxStartRange("AB2 init")
            if(iltime .eq. 1) then  
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*pr(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2)
               end do
               !$acc end parallel loop
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),q(:,:,1),rho(:,1),pr(:,1),&
               aux_h(:),Rmass(:,1),Rmom(:,:,1),Rener(:,1))
               
               if(mpi_size.ge.2) then
                  call nvtxStartRange("AB2 halo update")
                  call mpi_halo_atomic_update_real_mass_ener_momentum(Rmass(:,1),Rener(:,1),Rmom(:,:,1))
                  call nvtxEndRange
               end if               

               !$acc parallel loop
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(ipoin_w,idime) = u(ipoin_w,idime,1)*eta(ipoin_w,1)
                  end do
               end do
               !$acc end parallel loop

               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta(:,1))

               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:,1))
               end if

               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,1))

               gamma0 = 1.0_rp
               beta(1) = 1.0_rp
               beta(2) = 0.0_rp
               beta(3) = 0.0_rp
               alpha(1) = 1.0_rp
               alpha(2) = 0.0_rp
               alpha(3) = 0.0_rp
               !$acc update device(alpha(:))
               !$acc update device(beta(:))

               call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,1),u(:,:,1),csound,Tem(:,1),eta(:,1),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)

            else 
               if(iltime .eq. 2) then
                  gamma0 = 3.0_rp/2.0_rp
                  alpha(1) = 2.0_rp
                  alpha(2) = -0.5_rp
                  alpha(3) = 0.0_rp
                  beta(1) = 2.0_rp
                  beta(2) = -1.0_rp
                  beta(3) = 0.0_rp
                  !$acc update device(alpha(:))
                  !$acc update device(beta(:))      
               else
                  gamma0 = 11.0_rp/6.0_rp
                  alpha(1) = 3.0_rp
                  alpha(2) = -3.0_rp/2.0_rp
                  alpha(3) = 1.0_rp/3.0_rp
                  beta(1) = 3.0_rp
                  beta(2) = -3.0_rp
                  beta(3) = 1.0_rp
                  !$acc update device(alpha(:))
                  !$acc update device(beta(:))              
               end if
            end if
            call nvtxEndRange

            if(present(source_term) .or.flag_bouyancy_effect) then
               call nvtxStartRange("AB2 source")
               !$acc kernels
               Rsource(1:npoin,1:(ndime+2)) = 0.0_rp
               !$acc end kernels
               if(flag_bouyancy_effect) then
                  call mom_source_bouyancy_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,rho(:,2),Rsource(:,3:ndime+2))
                  call ener_source_bouyancy(nelem,npoin,connec,Ngp,dNgp,He,gpvol,q(:,:,2),Rsource(:,2))
               else if(present(source_term)) then
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,1:ndime,1),source_term(:,3:ndime+2),Rsource(:,3:ndime+2))
                  call ener_source_const(nelem,npoin,connec,Ngp,dNgp,He,gpvol,source_term(:,2),Rsource(:,2))
               end if
               call nvtxEndRange

               if(mpi_size.ge.2) then
                  call nvtxStartRange("AB2 halo update")
                  call mpi_halo_atomic_update_real_arrays(ndime+2,Rsource(:,:))
                  call nvtxEndRange
               end if
            end if
                  
            if((isWallModelOn) ) then
                  call nvtxStartRange("AB2 wall model")
                  if((numBoundsWM .ne. 0)) then
                     !$acc kernels
                     Rwmles(1:npoin,1:ndime) = 0.0_rp
                     !$acc end kernels
                     if(flag_type_wmles == wmles_type_reichardt) then
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),tauw,Rwmles)
                     else if (flag_type_wmles == wmles_type_abl) then
                        call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),zo,tauw,Rwmles)
                     end if
                  end if
                  call nvtxEndRange
          
                  if(mpi_size.ge.2) then
                     call nvtxStartRange("AB2 halo update")
                     call mpi_halo_atomic_update_real_arrays(ndime,Rwmles(:,:))
                     call nvtxEndRange
                  end if                  
            end if

            call nvtxStartRange("AB2 momentum")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*pr(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2)
            end do
            !$acc end parallel loop
            call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),q(:,:,2),rho(:,2),pr(:,2),&
            aux_h(:),Rmass(:,2),Rmom(:,:,2),Rener(:,2))            
            call nvtxEndRange

            if(mpi_size.ge.2) then
               call nvtxStartRange("AB2 halo update")
               call mpi_halo_atomic_update_real_mass_ener_momentum(Rmass(:,2),Rener(:,2),Rmom(:,:,2))
               call nvtxEndRange
            end if

            call nvtxStartRange("AB2 update u(2) & Rmom(1)")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoin_w = lpoin_w(ipoin)

               rho(ipoin_w,2) = -(beta(1)*Rmass(ipoin_w,2)+beta(2)*Rmass(ipoin_w,1)+beta(3)*Rmass(ipoin_w,3)) &
                              - Rsource(ipoin_w,1)
               rho(ipoin_w,2) =  Ml(ipoin_w)*(dt*rho(ipoin_w,2)/Ml(ipoin_w) + alpha(1)*rho(ipoin_w,1) + alpha(2)*rho(ipoin_w,3) + alpha(3)*rho(ipoin_w,4))/gamma0
               Rmass(ipoin_w,3) = Rmass(ipoin_w,1)
               Rmass(ipoin_w,1) = Rmass(ipoin_w,2)

               E(ipoin_w,2) = -(beta(1)*Rener(ipoin_w,2)+beta(2)*Rener(ipoin_w,1)+beta(3)*Rener(ipoin_w,3)) &
                              - Rsource(ipoin_w,2)
               E(ipoin_w,2) =  Ml(ipoin_w)*(dt*E(ipoin_w,2)/Ml(ipoin_w) + alpha(1)*E(ipoin_w,1) + alpha(2)*E(ipoin_w,3) + alpha(3)*E(ipoin_w,4))/gamma0
               Rener(ipoin_w,3) = Rener(ipoin_w,1)
               Rener(ipoin_w,1) = Rener(ipoin_w,2)

               !$acc loop seq   
               do idime = 1,ndime
                  q(ipoin_w,idime,2) = -(beta(1)*Rmom(ipoin_w,idime,2)+beta(2)*Rmom(ipoin_w,idime,1)+beta(3)*Rmom(ipoin_w,idime,3)) &
                                       -Rsource(ipoin_w,idime+2)-Rwmles(ipoin_w,idime)
                  q(ipoin_w,idime,2) =  Ml(ipoin_w)*(dt*q(ipoin_w,idime,2)/Ml(ipoin_w) + alpha(1)*q(ipoin_w,idime,1) + alpha(2)*q(ipoin_w,idime,3) + alpha(3)*q(ipoin_w,idime,4))/gamma0
                  Rmom(ipoin_w,idime,3) = Rmom(ipoin_w,idime,1)
                  Rmom(ipoin_w,idime,1) = Rmom(ipoin_w,idime,2)
              end do
            end do
            call nvtxEndRange

            call conjGrad_imex(1.0_rp/gamma0,igtime,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,numBoundsWM,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                              dlxigp_ip,He,gpvol,Ngp,Ml,helem,gamma_gas,Rgas,Cp,Prt,csound,mu_fluid,mu_e,mu_sgs,rho(:,1),rho(:,2),E(:,1),E(:,2),q(:,:,1),q(:,:,2), &
                              ndof,nbnodes,ldof,lbnodes,lnbn_nodes,bound,bou_codes,bou_codes_nodes,listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer)
                        
            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,maskMapped,rho(:,2),q(:,:,2),E(:,2),u_buffer)

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet(q(:,:,2),u(:,:,2),rho(:,2),E(:,2),pr(:,2))
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,2),q(:,:,2),u(:,:,2),pr(:,2),E(:,2),u_buffer)
               call nvtxEndRange
            end if

            if(flag_force_2D) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  q(ipoin,3,2) =  0.0_rp
               end do
               !$acc end parallel loop
            end if

            !$acc parallel loop
            do ipoin = 1,npoin_w
               umag = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  u(lpoin_w(ipoin),idime,2) = q(lpoin_w(ipoin),idime,2)/rho(lpoin_w(ipoin),2)
                  umag = umag + u(lpoin_w(ipoin),idime,2)**2
               end do
               e_int(lpoin_w(ipoin),2) = (E(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2))- &
                  0.5_rp*umag
               pr(lpoin_w(ipoin),2) = rho(lpoin_w(ipoin),2)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),2)
               csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2))
               umag = sqrt(umag)
               machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
               Tem(lpoin_w(ipoin),2) = pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)*Rgas)

               eta(lpoin_w(ipoin),1) = eta(lpoin_w(ipoin),2)
               eta(lpoin_w(ipoin),2) = (rho(lpoin_w(ipoin),2)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)**gamma_gas))
               !$acc loop seq
               do idime = 1,ndime                  
                  f_eta(lpoin_w(ipoin),idime)  = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
               end do
            end do
            !$acc end parallel loop        

            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
            gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta(:,2))

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(Reta(:,2))
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,2))

            call nvtxStartRange("Entropy residual")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoin_w = lpoin_w(ipoin)
               auxReta(ipoin_w) =  (beta(1)*Reta(ipoin_w,2)+beta(2)*Reta(ipoin_w,1)+beta(3)*Reta(ipoin_w,3)) + &
                                   (gamma0*eta(ipoin_w,2)-alpha(1)*eta(ipoin_w,1)-alpha(2)*eta(ipoin_w,3)-alpha(3)*eta(ipoin_w,4))/dt
               Reta(ipoin_w,3) = Reta(ipoin_w,1)
               Reta(ipoin_w,1) = Reta(ipoin_w,2)          
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,auxReta)
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
            ! If using Sutherland viscosity model:
            !
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call nvtxStartRange("Sutherland viscosity")
               call sutherland_viscosity(npoin,Tem(:,2),mu_factor,mu_fluid)
               call nvtxEndRange
            end if
            !
            ! Compute subgrid viscosity if active
            !
            if(flag_les == 1) then
               call nvtxStartRange("MU_SGS")
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3,mue_l) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),Ml,mu_sgs,mue_l)
               end if
               call nvtxEndRange
            end if
#endif
         end subroutine imex2_main

      end module time_integ_imex2
