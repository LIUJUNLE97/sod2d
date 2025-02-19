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

   real(rp), allocatable, dimension(:,:,:) :: Rmom,aux_omega
   real(rp), allocatable, dimension(:,:) :: aux_q,Rsource,Rwmles,u_flux_buffer
   real(rp), allocatable, dimension(:,:) ::GradP,f_eta,Reta
   real(rp), allocatable, dimension(:) :: auxReta, p_buffer, aux_p
   real(rp), allocatable, dimension(:)   :: beta, alpha
   real(rp) :: gamma0
   logical :: firstTimeStep = .true.


   contains
   subroutine init_rk4_solver_incomp(npoin)

      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      call nvtxStartRange("Init Incomp solver")

      allocate(Rmom(npoin,ndime,3),aux_omega(npoin,ndime,3))
      !$acc enter data create(Rmom(:,:,:),aux_omega(:,:,:))

      allocate(aux_q(npoin,ndime),Rsource(npoin,ndime),Rwmles(npoin,ndime),u_flux_buffer(npoin,ndime))
      !$acc enter data create(aux_q(:,:),Rsource(:,:),Rwmles(:,:),u_flux_buffer(:,:))

      allocate(gradP(npoin,ndime))
      !$acc enter data create(gradP(:,:))

      allocate(auxReta(npoin),f_eta(npoin,ndime),Reta(npoin,3),p_buffer(npoin),aux_p(npoin))
      !$acc enter data create(auxReta(:),f_eta(:,:),Reta(:,:),p_buffer(:),aux_p(:))
   
      !$acc kernels
      Rmom(1:npoin,1:ndime,1:3) = 0.0_rp
      Rsource(1:npoin,1:ndime) = 0.0_rp
      Rwmles(1:npoin,1:ndime) = 0.0_rp
      p_buffer(1:npoin) = nscbc_p_inf
      u_flux_buffer(1:npoin,1:ndime) = 0.0_rp
      !$acc end kernels

      allocate(alpha(3),beta(3))
      !$acc enter data create(alpha(:),beta(:))
      call nvtxEndRange

   end subroutine init_rk4_solver_incomp

   subroutine end_rk4_solver_incomp()
      implicit none


      !$acc exit data delete(Rmom(:,:,:))
      !$acc exit data delete(aux_omega(:,:,:))
      deallocate(Rmom,aux_omega)

      !$acc exit data delete(aux_q(:,:))
      !$acc exit data delete(Rsource(:,:))
      !$acc exit data delete(Rwmles(:,:))
      !$acc exit data delete(u_flux_buffer(:,:))
      deallocate(aux_q,Rsource,Rwmles,u_flux_buffer)


      !$acc exit data delete(gradP(:,:))
      deallocate(gradP)

      !$acc exit data delete(f_eta(:,:),auxReta(:),Reta(:,:),p_buffer(:))
      deallocate(f_eta,auxReta,Reta,p_buffer)

   end subroutine end_rk4_solver_incomp
 
         subroutine ab_main_incomp(igtime,iltime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,leviCivi,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,numBouCodes,bouCodes2BCType,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u,walave_pr,wmles_thinBL_fit_d,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,iltime,save_logFile_next
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin),lelpn(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow,maskMapped(npoin)
            real(rp),              intent(in)   :: leviCivi(ndime,ndime,ndime)
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
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*),numBouCodes
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*),bouCodes2BCType(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(in)      :: source_term(npoin,ndime)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_pr(npoin)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4), optional, intent(in)    :: wmles_thinBL_fit_d(npoin)
            integer(4)                          :: istep,ipoin,idime,icode,iPer,ipoin_w

            call nvtxStartRange("AB2 init")
            if(iltime .eq. 1) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_q(ipoin,idime) = u(ipoin,idime,1)*rho(ipoin,1)
                  end do
               end do
               !$acc end parallel loop
               call full_convec_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),aux_q,rho(:,1),Rmom(:,:,1))          

               if(mpi_size.ge.2) then
                  call nvtxStartRange("AB2 halo update")
                  call mpi_halo_atomic_update_real_arrays(ndime,Rmom(:,:,1))
                  call nvtxEndRange
               end if               

               gamma0 = 1.0_rp
               beta(1) = 1.0_rp
               beta(2) = 0.0_rp
               beta(3) = 0.0_rp
               alpha(1) = 1.0_rp
               alpha(2) = 0.0_rp
               alpha(3) = 0.0_rp
               !$acc update device(alpha(:))
               !$acc update device(beta(:))
               
               !$acc kernels
               aux_omega(:,:,3) = 0.0_rp
               aux_omega(:,:,2) = 0.0_rp
               aux_omega(:,:,1) = 0.0_rp
               !$acc end kernels

               
               
               if(flag_type_wmles == wmles_type_reichardt) then
                  call evalEXAtFace(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)
               else if (flag_type_wmles == wmles_type_abl) then
                  call evalEXAtFace(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)
               else if (flag_type_wmles == wmles_type_reichardt_hwm) then
                  call evalEXAtHWM(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)
               else if (flag_type_wmles == wmles_type_thinBL_fit) then
                  call evalEXAt1OffNode(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)
               else if (flag_type_wmles == wmles_type_thinBL_fit_hwm) then
                  call evalEXAtHWM(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)
               end if   

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

            if (noBoundaries .eqv. .false.) then

               call evalPAtOutlet(nelem,npoin,npoin_w,nboun,connec,bound,point2elem,bou_codes,bou_codes_nodes,lpoin_w, &
                  bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,mu_e,mu_sgs,rho,u(:,:,1),p_buffer,u_flux_buffer)

            end if

            if(present(source_term)) then
               call nvtxStartRange("AB2 source")
               !$acc kernels
               Rsource(1:npoin,1:ndime) = 0.0_rp
               !$acc end kernels
               call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,1:ndime,1),source_term(:,1:ndime),Rsource)
               call nvtxEndRange

               if(mpi_size.ge.2) then
                  call nvtxStartRange("AB2 halo update")
                  call mpi_halo_atomic_update_real_arrays(ndime,Rsource(:,:))
                  call nvtxEndRange
               end if
            end if
                  
            if((isWallModelOn) ) then
                  call nvtxStartRange("AB2 wall model")
                  if ((flag_type_wmles == wmles_type_thinBL_fit) .or. (flag_type_wmles == wmles_type_thinBL_fit_hwm)) then
                     call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,walave_pr(:),gradP,.true.)
                  end if
                  if((numBoundsWM .ne. 0)) then
                     !$acc kernels
                     Rwmles(1:npoin,1:ndime) = 0.0_rp
                     !$acc end kernels
                     if((flag_type_wmles == wmles_type_reichardt) .or. (flag_type_wmles == wmles_type_reichardt_hwm)) then
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),tauw,Rwmles)
                     else if (flag_type_wmles == wmles_type_abl) then
                        call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),zo,tauw,Rwmles)
                     else if ((flag_type_wmles == wmles_type_thinBL_fit) .or. (flag_type_wmles == wmles_type_thinBL_fit_hwm)) then
                        call evalWallModelThinBLFit(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),gradP,tauw,wmles_thinBL_fit_d,Rwmles)
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
            do ipoin = 1,npoin
               !$acc loop seq
               do idime = 1,ndime
                  aux_q(ipoin,idime) = u(ipoin,idime,1)*rho(ipoin,1)
               end do
            end do
            !$acc end parallel loop
            call full_convec_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),aux_q,rho(:,1),Rmom(:,:,2))
            call nvtxEndRange

            if(mpi_size.ge.2) then
               call nvtxStartRange("AB2 halo update")
               call mpi_halo_atomic_update_real_arrays(ndime,Rmom(:,:,2))
               call nvtxEndRange
            end if

            call nvtxStartRange("AB2 update u(2) & Rmom(1)")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoin_w = lpoin_w(ipoin)
               !$acc loop seq   
               do idime = 1,ndime
                  Rmom(ipoin_w,idime,2) = Rmom(ipoin_w,idime,2) - Rsource(ipoin_w,idime)
                  u(ipoin_w,idime,2) = -(beta(1)*Rmom(ipoin_w,idime,2)+beta(2)*Rmom(ipoin_w,idime,1)+beta(3)*Rmom(ipoin_w,idime,3)) 
                  Rmom(ipoin_w,idime,3) = Rmom(ipoin_w,idime,1)
                  Rmom(ipoin_w,idime,1) = Rmom(ipoin_w,idime,2)
              end do
            end do
            call nvtxEndRange
            
            call nvtxStartRange("AB2 update u(2)")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoin_w = lpoin_w(ipoin)
               !$acc loop seq   
               do idime = 1,ndime
                  u(ipoin_w,idime,2) =  (dt*u(ipoin_w,idime,2)/Ml(ipoin_w) + alpha(1)*u(ipoin_w,idime,1) + alpha(2)*u(ipoin_w,idime,3) + alpha(3)*u(ipoin_w,idime,4))/gamma0
              end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (noBoundaries .eqv. .false.) then
               if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet_incomp(u(:,:,2))
               call temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes_nodes,lnbn_nodes,normalsAtNodes,u(:,:,2),u_buffer,wmles_thinBL_fit_d) 
            end if
            if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u(:,:,2),u_buffer)

            call eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),pr(:,2))
            !$acc kernels
            pr(:,2) =-gamma0*pr(:,2)/dt
            !$acc end kernels

            if (noBoundaries .eqv. .false.) then
               !$acc parallel loop 
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  !$acc loop seq   
                  do idime = 1,ndime
                     aux_q(ipoin_w,idime) = -(beta(1)*aux_omega(ipoin_w,idime,2)+beta(2)*aux_omega(ipoin_w,idime,1)+beta(3)*aux_omega(ipoin_w,idime,3))
                   end do
               end do 
               !$acc end parallel loop     
               call bc_routine_pressure_flux(nelem,npoin,nboun,connec,bound,point2elem,bou_codes,bou_codes_nodes,numBoundCodes,bouCodes2BCType, &
                                              bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,mu_sgs,rho,aux_q,aux_p)
               !$acc parallel loop 
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  pr(ipoin_w,2) = pr(ipoin_w,2) - aux_p(ipoin_w)
               end do 
               !$acc end parallel loop 
            end if
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_pressure_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,pr(:,1),p_buffer)              
            end if
            call conjGrad_pressure_incomp(igtime,save_logFile_next,noBoundaries,nelem,npoin,npoin_w,connec,lpoin_w,lelpn,invAtoIJK,gmshAtoI,gmshAtoJ,&
                                          gmshAtoK,dlxigp_ip,He,gpvol,Ngp,dNgp,Ml,pr(:,1),pr(:,2), &
                                          nboun,bou_codes_nodes,normalsAtNodes)
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_pressure_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,pr(:,2),p_buffer)              
            end if
 
            call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,pr(:,2),gradP,.true.)
                        
            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoin_w = lpoin_w(ipoin)
               !$acc loop seq
               do idime = 1,ndime
                  u(ipoin_w,idime,2) = (u(ipoin_w,idime,2)-dt*gradP(ipoin_w,idime)/gamma0)*Ml(ipoin_w) & 
                                     -dt*Rwmles(ipoin_w,idime)/gamma0
               end do
            end do
            !$acc end parallel loop
                        
            if (noBoundaries .eqv. .false.) then
               if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet_incomp(u(:,:,2))
               call temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes_nodes,lnbn_nodes,normalsAtNodes,u(:,:,2),u_buffer,wmles_thinBL_fit_d) 
            end if

            call conjGrad_veloc_incomp(igtime,1.0_rp/gamma0,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,connec,lpoin_w,invAtoIJK,&
                                       gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem,mu_fluid,mu_e,mu_sgs,u(:,:,1),u(:,:,2),&
                                       bou_codes_nodes,normalsAtNodes,u_buffer,wmles_thinBL_fit_d)

            if (noBoundaries .eqv. .false.) then
               if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet_incomp(u(:,:,2))
               call temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes_nodes,lnbn_nodes,normalsAtNodes,u(:,:,2),u_buffer,wmles_thinBL_fit_d)
               !$acc kernels
               aux_omega(:,:,3) = aux_omega(:,:,1)
               aux_omega(:,:,1) = aux_omega(:,:,2)
               !$acc end kernels
               call compute_vorticity(nelem,npoin,npoin_w,lpoin_w,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),aux_q,.true.)
               call compute_vorticity(nelem,npoin,npoin_w,lpoin_w,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_q,aux_omega(:,:,2),.true.)
            end if
            if (flag_buffer_on .eqv. .true.) call updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u(:,:,2),u_buffer)
            
            !
            ! Compute subgrid viscosity if active
            !

            if(flag_les == 1) then
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,&
                              rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3,mue_l) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),&
                                 Ml,mu_sgs,mue_l)
               end if
               call nvtxEndRange
            end if
         end subroutine ab_main_incomp

         subroutine updateBuffer_incomp(npoin,npoin_w,coord,lpoin_w,maskMapped,u,u_buffer)
            implicit none
            integer(4),           intent(in)    :: npoin
            integer(4),           intent(in)    :: npoin_w
            real(rp),             intent(in)    :: coord(npoin,ndime)
            integer(4),           intent(in)    :: lpoin_w(npoin_w),maskMapped(npoin)
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
               if (maskMapped(lpoin_w(ipoin)) == 1) then
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
                  !top
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
               end if

               u(lpoin_w(ipoin),1) = u_buffer(lpoin_w(ipoin),1) + xi*(u(lpoin_w(ipoin),1)-u_buffer(lpoin_w(ipoin),1))
               u(lpoin_w(ipoin),2) = u_buffer(lpoin_w(ipoin),2) + xi*(u(lpoin_w(ipoin),2)-u_buffer(lpoin_w(ipoin),2))
               u(lpoin_w(ipoin),3) = u_buffer(lpoin_w(ipoin),3) + xi*(u(lpoin_w(ipoin),3)-u_buffer(lpoin_w(ipoin),3))
            end do
            !$acc end parallel loop
            call nvtxEndRange
         end subroutine updateBuffer_incomp
      end module time_integ_incomp
