module time_integ_incomp

   use mod_nvtx
   use elem_convec
   use elem_convec_incomp
   use elem_diffu_incomp
   use elem_source
   use mod_solver
   use mod_entropy_viscosity_incomp
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
   real(rp), allocatable, dimension(:,:) :: aux_q,Rsource,Rwmles
   real(rp), allocatable, dimension(:,:) :: Rmom_sum,Rdiff_mom
   real(rp), allocatable, dimension(:,:) ::GradP,f_eta,Reta
   real(rp), allocatable, dimension(:) :: auxReta


   contains
   subroutine init_rk4_solver_incomp(npoin)

      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      allocate(Rmom(npoin,ndime,2))
      !$acc enter data create(Rmom(:,:,:))

      allocate(aux_q(npoin,ndime),Rsource(npoin,ndime),Rwmles(npoin,ndime))
      !$acc enter data create(aux_q(:,:),Rsource(:,:),Rwmles(:,:))

      allocate(Rmom_sum(npoin,ndime),Rdiff_mom(npoin,ndime))
      !$acc enter data create(Rmom_sum(:,:))
      !$acc enter data create(Rdiff_mom(:,:))

      allocate(gradP(npoin,ndime))
      !$acc enter data create(gradP(:,:))

      allocate(auxReta(npoin),f_eta(npoin,ndime),Reta(npoin,ndime))
      !$acc enter data create(auxReta(:),f_eta(:,:),Reta(:,:))
   
      !$acc kernels
      Rmom(1:npoin,1:ndime,1:2) = 0.0_rp
      Rsource(1:npoin,1:ndime) = 0.0_rp
      Rwmles(1:npoin,1:ndime) = 0.0_rp
      !$acc end kernels

      call nvtxEndRange

   end subroutine init_rk4_solver_incomp

   subroutine end_rk4_solver_incomp()
      implicit none


      !$acc exit data delete(Rmom(:,:,:))
      deallocate(Rmom)

      !$acc exit data delete(aux_q(:,:))
      !$acc exit data delete(Rsource(:,:))
      !$acc exit data delete(Rwmles(:,:))
      deallocate(aux_q,Rsource,Rwmles)

      !$acc exit data delete(Rmom_sum(:,:))
      !$acc exit data delete(Rdiff_mom(:,:))
      deallocate(Rmom_sum,Rdiff_mom)

      !$acc exit data delete(gradP(:,:))
      deallocate(gradP)

      !$acc exit data delete(f_eta(:,:),auxReta(:),Reta(:,:))
      deallocate(f_eta,auxReta,Reta)

   end subroutine end_rk4_solver_incomp
 
         subroutine ab_main_incomp(igtime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         mue_l,convertIJK,al_weights,am_weights,an_weights, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,save_logFile_next
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn(nboun,npbou),lnbn_nodes(npoin),lelpn(npoin),convertIJK(0:porder+2)
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
            real(rp),             intent(inout) :: mue_l(nelem,nnode)
            real(rp),             intent(in)    :: al_weights(-1:1),am_weights(-1:1),an_weights(-1:1)
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

            if(igtime .eq. 1) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_q(ipoin,idime) = u(ipoin,idime,1)*rho(ipoin,1)
                  end do
               end do
               !$acc end parallel loop
               call full_convec_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),aux_q,rho(:,1),Rmom(:,:,1))          
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
            
            end if

            if(present(source_term)) then
               !$acc kernels
               Rsource(1:npoin,1:ndime) = 0.0_rp
               !$acc end kernels
               call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,1:ndime,1),source_term,Rsource)
            end if

            call full_diffusion_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),&
                                          mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mom)
                  
            if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  !$acc kernels
                  Rwmles(1:npoin,1:ndime) = 0.0_rp
                  !$acc end kernels
                  call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                     rho(:,1),walave_u(:,:),tauw,Rwmles)
            end if

            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq
               do idime = 1,ndime
                  aux_q(ipoin,idime) = u(ipoin,idime,1)*rho(ipoin,1)
               end do
            end do
            !$acc end parallel loop
            call full_convec_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),aux_q,rho(:,1),Rmom(:,:,2))          

            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq   
               do idime = 1,ndime
                  u(ipoin,idime,2) = -dt*(1.5_rp*Rmom(ipoin,idime,2)-0.5_rp*(Rmom(ipoin,idime,1))+0.5_rp*Rdiff_mom(ipoin,idime))&
                                     -dt*Rsource(ipoin,idime)-dt*Rwmles(ipoin,idime)
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

            call conjGrad_veloc_incomp(igtime,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,numBoundsWM,connec,lpoin_w,invAtoIJK,gmshAtoI,&
                                       gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,mu_fluid,mu_e,mu_sgs,u(:,:,1),u(:,:,2), &
                                       ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&                      
                                       listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u)

            call eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),pr(:,2))
            !$acc kernels
            pr(:,2) = -pr(:,2)/dt
            !$acc end kernels

            call conjGrad_pressure_incomp(igtime,save_logFile_next,noBoundaries,nelem,npoin,npoin_w,connec,lpoin_w,lelpn,invAtoIJK,gmshAtoI,gmshAtoJ,&
                                          gmshAtoK,dlxigp_ip,He,gpvol,Ngp,dNgp,Ml,pr(:,1),pr(:,2), &
                                          nboun,bou_codes_nodes,normalsAtNodes)
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_pressure_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,pr(:,2))
            end if
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
               call temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,u(:,:,2),u_buffer)
            end if
            !
            ! Compute subgrid viscosity if active
            !

            !$acc parallel loop
            do ipoin = 1,npoin_w
               eta(lpoin_w(ipoin),1) = eta(lpoin_w(ipoin),2)
               eta(lpoin_w(ipoin),2) = 0.5*(u(lpoin_w(ipoin),1,2)**2 + u(lpoin_w(ipoin),2,2)**2 + u(lpoin_w(ipoin),3,2)**2)
               !$acc loop seq
               do idime = 1,ndime
                  f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,2)*eta(lpoin_w(ipoin),2)
               end do
            end do
            !$acc end parallel loop

            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,2),u(:,:,2),Reta(:,2))

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(Reta(:,2))
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,2))

            !$acc parallel loop
            do ipoin = 1,npoin_w
               auxReta(lpoin_w(ipoin)) = (1.5_rp*Reta(lpoin_w(ipoin),2)-0.5_rp*Reta(lpoin_w(ipoin),1)) + &
                                          (eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
               Reta(lpoin_w(ipoin),1) = Reta(lpoin_w(ipoin),2)
            end do
            !$acc end parallel loop

            if (noBoundaries .eqv. .false.) then
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,auxReta)
            end if

            call smart_visc_spectral_incomp(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
                                            rho(:,2),u(:,:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                                            mue_l,convertIJK,al_weights,am_weights,an_weights)


            if(flag_les == 1) then
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,&
                              rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3,mue_l,convertIJK,al_weights,am_weights,an_weights) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),&
                                 Ml,mu_sgs,mue_l,convertIJK,al_weights,am_weights,an_weights)
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
      end module time_integ_incomp
