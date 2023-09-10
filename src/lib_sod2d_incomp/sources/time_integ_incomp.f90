module time_integ_incomp

   use mod_nvtx
   use elem_convec
   use elem_convec_incomp
   use elem_diffu_incomp
   use elem_source
   use mod_solver
   use mod_entropy_viscosity
   use mod_numerical_params
   use mod_fluid_viscosity
   use mod_sgs_viscosity
   use mod_sgs_ilsa_viscosity
   use mod_bc_routines_incomp
   use mod_wall_model
   use mod_operators
   use mod_solver_incomp

   implicit none

   real(rp), allocatable, dimension(:,:,:) :: Rmom
   real(rp), allocatable, dimension(:,:) :: aux_q,Rsource
   real(rp), allocatable, dimension(:,:) :: Rmom_sum,Rdiff_mom
   real(rp), allocatable, dimension(:,:) ::GradP

   contains
   subroutine init_rk4_solver_incomp(npoin)

      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      allocate(Rmom(npoin,ndime,2))
      !$acc enter data create(Rmom(:,:,:))

      allocate(aux_q(npoin,ndime),Rsource(npoin,ndime))
      !$acc enter data create(aux_q(:,:),Rsource(:,:))

      allocate(Rmom_sum(npoin,ndime),Rdiff_mom(npoin,ndime))
      !$acc enter data create(Rmom_sum(:,:))
      !$acc enter data create(Rdiff_mom(:,:))

      allocate(gradP(npoin,ndime))
      !$acc enter data create(gradP(:,:))
   
      !$acc kernels
      Rmom(1:npoin,1:ndime,1:2) = 0.0_rp
      !$acc end kernels

      call nvtxEndRange

   end subroutine init_rk4_solver_incomp

   subroutine end_rk4_solver_incomp()
      implicit none


      !$acc exit data delete(Rmom(:,:,:))
      deallocate(Rmom)

      !$acc exit data delete(aux_q(:,:))
      !$acc exit data delete(Rsource(:,:))
      deallocate(aux_q,Rsource)

      !$acc exit data delete(Rmom_sum(:,:))
      !$acc exit data delete(Rdiff_mom(:,:))
      deallocate(Rmom_sum,Rdiff_mom)

      !$acc exit data delete(gradP(:,:))
      deallocate(gradP)

   end subroutine end_rk4_solver_incomp
 
         subroutine ab_main_incomp(igtime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,save_logFile_next
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn(nboun,npbou),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,3)
            real(rp),             intent(inout) :: u(npoin,ndime,3)
            real(rp),             intent(inout) :: q(npoin,ndime,3)
            real(rp),             intent(inout) :: pr(npoin,3)
            real(rp),             intent(inout) :: E(npoin,3)
            real(rp),             intent(inout) :: Tem(npoin,2)
            real(rp),             intent(inout) :: e_int(npoin,2)
            real(rp),             intent(inout) :: eta(npoin,3)
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
            real(rp), optional, intent(in)      :: source_term(npoin,ndime)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime)
            integer(4)                          :: istep, ipoin, idime,icode

            !
            ! Loop over all AB+CN steps
            ! 

            if(present(source_term)) then
               !$acc kernels
               Rsource(1:npoin,1:ndime) = 0.0_rp
               !$acc end kernels
               call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(ipoin,idime,1),source_term,Rsource)
            end if


            call full_diffusion_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mom)
               
            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq
               do idime = 1,ndime
                  aux_q(ipoin,idime) = u(ipoin,idime,1)*rho(ipoin,1)
               end do
            end do
            !$acc end parallel loop
            call full_convec_ijk_incomp(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),aux_q,rho(:,1),Rmom(:,:,2))

            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq   
               do idime = 1,ndime
                  u(ipoin,idime,2) = -dt*(1.5_rp*Rmom(ipoin,idime,2)-0.5_rp*(Rmom(ipoin,idime,1))+0.5_rp*Rdiff_mom(ipoin,idime))+dt*Rsource(ipoin,idime)
                  Rmom(ipoin,idime,1) = Rmom(ipoin,idime,2)
              end do
            end do

            !$acc end parallel loop     
            if(mpi_size.ge.2) then
               do idime = 1,ndime
                  call mpi_halo_atomic_update_real(u(:,idime,2))
                  end do              
            end if
            
            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq   
               do idime = 1,ndime
                  u(ipoin,idime,2) =  u(ipoin,idime,2) + u(ipoin,idime,1)*Ml(ipoin)
              end do
            end do
            !$acc end parallel loop   

            call conjGrad_veloc_incomp(igtime,save_logFile_next,dt,nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,mu_fluid,mu_e,mu_sgs,u(:,:,1),u(:,:,2))

            call eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),pr(:,2))
            !$acc kernels
            pr(:,2) = -pr(:,2)/dt
            !$acc end kernels

            call conjGrad_pressure_incomp(igtime,save_logFile_next,nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,pr(:,1),pr(:,2))
            call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,pr(:,2),gradP,.true.)
                        
            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq
               do idime = 1,ndime
                  u(ipoin,idime,2) = u(ipoin,idime,2)-dt*gradP(ipoin,idime)
               end do
            end do
            !$acc end parallel loop

            if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,u(:,:,2),u_buffer)

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,u(:,:,2),u_buffer)
               call nvtxEndRange
            end if
            !
            ! Compute subgrid viscosity if active
            !
            if(flag_les == 1) then
               call nvtxStartRange("MU_SGS")
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),Ml,mu_sgs)
               end if
               call nvtxEndRange
            end if
         end subroutine ab_main_incomp

         subroutine updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,u,u_buffer)
            implicit none
            integer(4),           intent(in)    :: npoin
            integer(4),           intent(in)    :: npoin_w
            real(rp),             intent(in)    :: coord(npoin,ndime)
            integer(4),           intent(in)    :: lpoin_w(npoin_w)
            real(rp),             intent(inout) :: u(npoin,ndime)
            real(rp),             intent(in) :: u_buffer(npoin,ndime)
            integer(4) :: ipoin
            real(rp)   :: xs,xb,xi,c1,c2

            c1 = 0.01_rp
            c2 = 10.0_rp

            call nvtxStartRange("Update buffer")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               xi = 1.0_rp
               !east
               if(flag_buffer_on_east .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),1)
                  if(xs>flag_buffer_e_min) then
                     xb = (xs-flag_buffer_e_min)/flag_buffer_e_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !west
               if(flag_buffer_on_west .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),1)
                  if(xs<flag_buffer_w_min) then
                     xb = (flag_buffer_w_min-xs)/flag_buffer_w_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !north
               if(flag_buffer_on_north .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),2)
                  if(xs>flag_buffer_n_min) then
                     xb = (xs-flag_buffer_n_min)/flag_buffer_n_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !south
               if(flag_buffer_on_south .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),2)
                  if(xs<flag_buffer_s_min) then
                     xb = (flag_buffer_s_min-xs)/flag_buffer_s_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !north
               if(flag_buffer_on_top .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),3)
                  if(xs>flag_buffer_t_min) then
                     xb = (xs-flag_buffer_t_min)/flag_buffer_t_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if
               !bottom
               if(flag_buffer_on_bottom .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),3)
                  if(xs<flag_buffer_b_min) then
                     xb = (flag_buffer_b_min-xs)/flag_buffer_b_size
                     xi = min((1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))),xi)
                  end if
               end if

               u(lpoin_w(ipoin),1) = u_buffer(lpoin_w(ipoin),1) + xi*(u(lpoin_w(ipoin),1)-u_buffer(lpoin_w(ipoin),1))
               u(lpoin_w(ipoin),2) = u_buffer(lpoin_w(ipoin),2) + xi*(u(lpoin_w(ipoin),2)-u_buffer(lpoin_w(ipoin),2))
               u(lpoin_w(ipoin),3) = u_buffer(lpoin_w(ipoin),3) + xi*(u(lpoin_w(ipoin),3)-u_buffer(lpoin_w(ipoin),3))

            end do
            !$acc end parallel loop
            call nvtxEndRange
         end subroutine updateBuffer_incomp

         subroutine filter_presure(dt,nelem,npoin,npoin_w,connec,lpoin_w,p,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,helem)

            implicit none

            integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),   intent(inout)  :: p(npoin)
            real(rp),   intent(in)    :: dt,Ml(npoin),dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),helem(nelem,nnode)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4)               :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode),ipoin2
            integer(4)              :: convertIJK(0:porder+2),ii,jj,kk,mm,nn,ll,iter
            real(rp)                :: p_l(npoin),al(-1:1),am(-1:1),an(-1:1),aux1
            real(rp)                :: gradIsoP(ndime)
            real(rp)                :: gradP(ndime),divDp,nu_e
            real(rp)                :: ResP(npoin),pl(nnode),gradPl(nnode,ndime)

            !$acc kernels
             ResP(:) = 0.0_rp
            !$acc end kernels

             !$acc parallel loop gang  private(ipoin,pl,gradPl)
             do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                   ipoin(inode) = connec(ielem,inode)
                end do
                !$acc loop vector
                do inode = 1,nnode
                   pl(inode) = p(ipoin(inode))
                end do
                gradPl(:,:) = 0.0_rp

                !$acc loop vector private(gradIsoP,gradP)
                do igaus = 1,ngaus

                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   gradIsoP(:) = 0.0_rp
                   !$acc loop seq
                   do ii=1,porder+1
                      gradIsoP(1) = gradIsoP(1) + dlxigp_ip(igaus,1,ii)*pl(invAtoIJK(ii,isoJ,isoK))
                      gradIsoP(2) = gradIsoP(2) + dlxigp_ip(igaus,2,ii)*pl(invAtoIJK(isoI,ii,isoK))
                      gradIsoP(3) = gradIsoP(3) + dlxigp_ip(igaus,3,ii)*pl(invAtoIJK(isoI,isoJ,ii))
                   end do

                   gradP(:) = 0.0_rp
                   !$acc loop seq
                   do idime=1, ndime
                      !$acc loop seq
                      do jdime=1, ndime
                         gradP(idime) = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                      end do
                   end do

                   !$acc loop seq
                   do idime = 1,ndime
                      gradPl(igaus,idime) =  gradP(idime)
                   end do
                end do

                !$acc loop vector private(divDp) 
                do igaus = 1,ngaus
                   !nu_e = ((2.0_rp*helem(ielem,igaus))**2)/24.0_rp
                   nu_e = dt/24.0_rp

                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   divDp = 0.0_rp
                   
                   !$acc loop seq
                   do ii=1,porder+1
                      !$acc loop seq
                      do idime=1,ndime
                         divDp = divDp + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*gradPl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDp = divDp + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*gradPl(invAtoIJK(isoI,ii,isoK),idime)
                         divDp = divDp + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*gradPl(invAtoIJK(isoI,isoJ,ii),idime)
                      end do
                   end do

                   !$acc atomic update
                   ResP(ipoin(igaus)) = ResP(ipoin(igaus))+nu_e*divDp
                   !$acc end atomic
                end do
             end do
             !$acc end parallel loop

               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(ResP(:))
                  call nvtxEndRange
               end if

               call nvtxStartRange("Lumped mass solver on generic")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,ResP(:))
               call nvtxEndRange
               
               !$acc parallel loop
               do ipoin2 = 1,npoin
                  p(ipoin2) = p(ipoin2) - ResP(ipoin2)
               end do
               !$acc end parallel loop
         end subroutine filter_presure
      end module time_integ_incomp
