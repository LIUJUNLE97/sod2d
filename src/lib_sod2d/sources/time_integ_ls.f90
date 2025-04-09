module time_integ_ls

   use mod_nvtx
   use elem_convec
   use elem_diffu
   use elem_stab
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
   use time_integ, only : updateBuffer, limit_rho
   use mod_operators , only:eval_divergence,compute_vorticity
   use mod_rough_elem

   implicit none

   real(rp), allocatable, dimension(:)   :: aux_h,divU_ls
   real(rp), allocatable, dimension(:,:) :: f_eta,Rwmles,f_eta2,vorti, phiD
   real(rp), allocatable, dimension(:)   :: Rmass,Rener
   real(rp), allocatable, dimension(:,:) :: Rmom,Reta
   real(rp), allocatable, dimension(:)   :: auxReta
   real(rp), allocatable, dimension(:)   :: tau_stab_ls
   real(rp), allocatable, dimension(:,:) :: ProjMass_ls,ProjEner_ls,ProjMX_ls,ProjMY_ls,ProjMZ_ls

   real(rp), allocatable, dimension(:,:) :: lambda_ij, gamma_ij
   logical :: firstTimeStep = .true.

   contains

   subroutine init_rk4_ls_solver(nelem,npoin)
      implicit none
      integer(4),intent(in) :: nelem,npoin

      call nvtxStartRange("Init RK_LS solver")

      allocate(aux_h(npoin), divU_ls(npoin))
      !$acc enter data create(aux_h(:),divU_ls(:))

      allocate(auxReta(npoin),Rmass(npoin),Rener(npoin))
      !$acc enter data create(Rmass(:))
      !$acc enter data create(Rener(:))
      !$acc enter data create(auxReta(:))

      allocate(f_eta(npoin,ndime),f_eta2(npoin,ndime),Reta(npoin,2),Rmom(npoin,ndime),Rwmles(npoin,ndime),vorti(npoin,ndime),phiD(nelem,nnode))
      !$acc enter data create(Rmom(:,:))
      !$acc enter data create(f_eta(:,:))
      !$acc enter data create(f_eta2(:,:))
      !$acc enter data create(Reta(:,:))
      !$acc enter data create(Rwmles(:,:))
      !$acc enter data create(vorti(:,:))
      !$acc enter data create(phiD(:,:))

      allocate(lambda_ij(flag_rk_ls_stages+1,flag_rk_ls_stages+1),gamma_ij(flag_rk_ls_stages+1,flag_rk_ls_stages+1))
      !$acc enter data create(lambda_ij(:,:))
      !$acc enter data create(gamma_ij(:,:))

      allocate(ProjMass_ls(npoin,ndime),ProjEner_ls(npoin,ndime),ProjMX_ls(npoin,ndime),ProjMY_ls(npoin,ndime),ProjMZ_ls(npoin,ndime),tau_stab_ls(nelem))
      !$acc enter data create(ProjMass_ls(:,:),ProjEner_ls(:,:),ProjMX_ls(:,:),ProjMY_ls(:,:),ProjMZ_ls(:,:),tau_stab_ls(:))

      !$acc kernels
      aux_h(:) = 0.0_rp
      Rmass(:) = 0.0_rp
      Rener(:) = 0.0_rp
      auxReta(:) = 0.0_rp
      Rmom(:,:) = 0.0_rp
      f_eta(:,:) = 0.0_rp
      f_eta2(:,:) = 0.0_rp
      Reta(:,:) = 0.0_rp
      Rwmles(:,:) = 0.0_rp
      ProjMass_ls(:,:) = 0.0_rp
      ProjEner_ls(:,:) = 0.0_rp
      ProjMX_ls(:,:) = 0.0_rp
      ProjMY_ls(:,:) = 0.0_rp
      ProjMZ_ls(:,:) = 0.0_rp
      tau_stab_ls(:) = 0.0_rp
      lambda_ij(:,:) = 0.0_rp
      gamma_ij(:,:)  = 0.0_rp
      phiD(:,:) = 0.0_rp
      !$acc end kernels

      if (flag_rk_ls_stages == 3) then         
         lambda_ij(:,:) = 0.0_rp         
         lambda_ij(2,1) = 1.0_rp
         lambda_ij(3,1) = 3.0_rp/4.0_rp
         lambda_ij(3,2) = 1.0_rp/4.0_rp
         lambda_ij(4,1) = 1.0_rp/3.0_rp
         lambda_ij(4,3) = 2.0_rp/3.0_rp
      
         gamma_ij(:,:) = 0.0_rp
         gamma_ij(2,1) = 1.0_rp 
         gamma_ij(3,2) = 0.25_rp 
         gamma_ij(4,3) = 2.0_rp/3.0_rp
      else if (flag_rk_ls_stages == 4) then         
         lambda_ij(:,:) = 0.0_rp         
         lambda_ij(2,1) = 1.0_rp
         lambda_ij(3,2) = 1.0_rp
         lambda_ij(4,1) = 2.0_rp/3.0_rp
         lambda_ij(4,3) = 1.0_rp/3.0_rp
         lambda_ij(5,4) = 1.0_rp
         
         gamma_ij(:,:) = 0.0_rp
         gamma_ij(2,1) = 0.5_rp 
         gamma_ij(3,2) = 0.5_rp 
         gamma_ij(4,3) = 1.0_rp/6.0_rp
         gamma_ij(5,4) = 0.5_rp
      else if (flag_rk_ls_stages == 5) then
         lambda_ij(:,:) = 0.0_rp         
         lambda_ij(2,1) = 1.0_rp
         lambda_ij(3,2) = 1.0_rp
         lambda_ij(4,3) = 1.0_rp
         lambda_ij(5,4) = 1.0_rp
         lambda_ij(6,5) = 1.0_rp
         gamma_ij(:,:) = 0.0_rp
         select case (flag_rk_ls_n)
            case (1)
               !SSP53_2N_1
               lambda_ij(5,1) = real(0.571403511494104d0,rp)
               lambda_ij(5,4) = real(0.428596488505896d0,rp)

               gamma_ij(2,1) = real(0.443568244942995d0,rp)
               gamma_ij(3,2) = real(0.291111420073766d0,rp)
               gamma_ij(4,3) = real(0.270612601278217d0,rp)
               gamma_ij(5,4) = real(0.110577759392786d0,rp)
               gamma_ij(6,5) = real(0.458557505351052d0,rp)
            case (2)
               !SSP53_2N_2
               lambda_ij(4,1) = real(0.682342861037239d0,rp) ! lambda_41
               lambda_ij(4,3) = real(0.317657138962761d0,rp) ! lambda_43
               lambda_ij(6,1) = real(0.045230974482400d0,rp) ! lambda_61
               lambda_ij(6,5) = real(0.954769025517600d0,rp) ! lambda_65

               gamma_ij(2,1) = real(0.465388589249323,rp) ! gamma_21
               gamma_ij(3,2) = real(0.465388589249323,rp) ! gamma_32
               gamma_ij(4,3) = real(0.124745797313998,rp) ! gamma_43
               gamma_ij(5,4) = real(0.465388589249323,rp) ! gamma_54
               gamma_ij(6,5) = real(0.154263303748666,rp) ! gamma_65
            case (3)
               !SSP53_2N_3
               lambda_ij(5,1) = real(0.592032910942121,rp)
               lambda_ij(5,4) = 1.0_rp - lambda_ij(5,1)
               
               gamma_ij(2,1) = real(0.266541020678955,rp)
               gamma_ij(3,2) = real(0.548560709048532,rp)
               gamma_ij(4,3) = real(0.289517014154401,rp)
               gamma_ij(5,4) = real(0.086408328057923,rp)
               gamma_ij(6,5) = real(0.462943578481813,rp)
            case (4)
               !SSP53_2N_4
               lambda_ij(4,1) = real(0.707858560931430,rp)
               lambda_ij(4,3) = 1.0_rp - lambda_ij(4,1)
               lambda_ij(6,1) = real(0.222853615080669,rp)
               lambda_ij(6,5) = 1.0_rp - lambda_ij(6,1)
               
               gamma_ij(2,1) = real(0.292845746913355,rp)
               gamma_ij(3,2) = real(0.339532793976408,rp)
               gamma_ij(4,3) = real(0.200532330324672,rp)
               gamma_ij(5,4) = real(0.701676169006879,rp)
               gamma_ij(6,5) = real(0.155278812461877,rp)
            case default
               print *, 'wrong option bro!!!'
         end select
      else
         write(1,*) "--| NOT CODED FOR RK_LS stages > 5 YET!"
         stop 1
      end if
      !$acc update device(lambda_ij(:,:))
      !$acc update device(gamma_ij(:,:))

      call nvtxEndRange

   end subroutine init_rk4_ls_solver

   subroutine end_rk4_ls_solver()
      implicit none

      !$acc exit data delete(aux_h(:))
      deallocate(aux_h)

      !$acc exit data delete(Rmass(:))
      !$acc exit data delete(Rener(:))
      !$acc exit data delete(auxReta(:))
      deallocate(auxReta,Rmass,Rener)

      !$acc exit data delete(Rmom(:,:))
      !$acc exit data delete(f_eta(:,:))
      !$acc exit data delete(Reta(:,:))
      !$acc exit data delete(Rwmles(:,:))
      deallocate(f_eta,Reta,Rmom,Rwmles)

      !$acc exit data delete(lambda_ij(:,:))
      !$acc exit data delete(gamma_ij(:,:))
      deallocate(lambda_ij,gamma_ij)

   end subroutine end_rk4_ls_solver

         subroutine rk_4_ls_main(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,leviCivi,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,invMl,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,u_mapped,tauw,source_term,walave_u,walave_pr,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin),lelpn(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow,maskMapped(npoin)
            real(rp),              intent(in)   :: leviCivi(ndime,ndime,ndime)
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin), invMl(npoin)
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
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime), u_mapped(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(inout)      :: source_term(npoin,ndime+2)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_pr(npoin)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: pos, ipoin_w, ielem, inode
            integer(4)                          :: istep, ipoin, idime,icode
            real(rp)                            :: umag, rho_min, rho_avg, divUl

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value

            if(flag_lps_stab) then
               call comp_tau(nelem,npoin,connec,csound,u(:,:,pos),helem,dt,tau_stab_ls)
            end if

            if(firstTimeStep .eqv. .true.) then
               firstTimeStep = .false.

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
         
               call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Reta(:,1))      
                  
               call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,Reta(:,1),Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,1),u(:,:,1),csound,Tem(:,1),eta(:,1),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)                  

               if((isWallModelOn) ) then
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
               end if
               
               call updateF(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,&
                        ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,invMl,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                        rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                        ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&            
                        listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,u_mapped,tauw,source_term,walave_u,walave_pr,zo)
                              
            end if           
            

            !$acc parallel loop
            do ipoin = 1,npoin_w
               rho(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos) + gamma_ij(2,1)*dt*Rmass(lpoin_w(ipoin))
               E(lpoin_w(ipoin),pos) = E(lpoin_w(ipoin),pos) + gamma_ij(2,1)*dt*Rener(lpoin_w(ipoin))
               !$acc loop seq
               do idime = 1,ndime
                  q(lpoin_w(ipoin),idime,pos) = q(lpoin_w(ipoin),idime,pos) + gamma_ij(2,1)*dt*Rmom(lpoin_w(ipoin),idime)
               end do
            end do
            !$acc end parallel loop

            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")

            do istep = 2,flag_rk_ls_stages

               call updateF(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,&
                        ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,invMl,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                        rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                        ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&            
                        listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,u_mapped,tauw,source_term,walave_u,walave_pr,zo)
               !
               ! Accumulate the residuals
               !
               call nvtxStartRange("Accumulate residuals")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  rho(lpoin_w(ipoin),pos) = lambda_ij(istep+1,1)*rho(lpoin_w(ipoin),1) + lambda_ij(istep+1,istep)*rho(lpoin_w(ipoin),pos) + gamma_ij(istep+1,istep)*dt*Rmass(lpoin_w(ipoin))
                  E(lpoin_w(ipoin),pos) = lambda_ij(istep+1,1)*E(lpoin_w(ipoin),1) + lambda_ij(istep+1,istep)*E(lpoin_w(ipoin),pos) + gamma_ij(istep+1,istep)*dt*Rener(lpoin_w(ipoin))
                  !$acc loop seq
                  do idime = 1,ndime
                     q(lpoin_w(ipoin),idime,pos) = lambda_ij(istep+1,1)*q(lpoin_w(ipoin),idime,1) + lambda_ij(istep+1,istep)*q(lpoin_w(ipoin),idime,pos) + gamma_ij(istep+1,istep)*dt*Rmom(lpoin_w(ipoin),idime)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange                              
            end do
            call nvtxEndRange



            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,maskMapped,rho(:,pos),q(:,:,pos),E(:,pos),u_buffer)

            !
            ! Apply bcs after update
            !
            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet(q(:,:,2),u(:,:,2),rho(:,2),E(:,2),pr(:,2))
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer,u_mapped)
               call nvtxEndRange
            end if

            if(flag_force_2D) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  q(ipoin,3,pos) =  0.0_rp
               end do
               !$acc end parallel loop
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
                  f_eta(lpoin_w(ipoin),idime)  = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
                  f_eta2(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,2)*eta(lpoin_w(ipoin),2)
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if(flag_use_ducros .eqv. .true.) then            
   
               call eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),divU_ls)
               call compute_vorticity(nelem,npoin,npoin_w,lpoin_w,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),vorti,.true.) 

               call ducros_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,coord,dNgp,gpvol,wgp, &               
                  rho(:,2),u(:,:,2),divU_ls,vorti,csound,helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)

            else

               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta(:,2))

               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:,2))
               end if

               call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Reta(:,2))

               call nvtxStartRange("Entropy residual")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  auxReta(lpoin_w(ipoin)) = (1.5_rp*Reta(lpoin_w(ipoin),2)-0.5_rp*Reta(lpoin_w(ipoin),1)) + &
                                                factor_comp*(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
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
               ! Compute entropy viscosity
               !
               call nvtxStartRange("Entropy viscosity evaluation")
               call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
                  gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
               call nvtxEndRange   

               call eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),divU_ls)
               call compute_vorticity(nelem,npoin,npoin_w,lpoin_w,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),vorti,.true.) 

               call ducros_visc_limit(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,coord,dNgp,gpvol,wgp, &               
                  rho(:,2),u(:,:,2),divU_ls,vorti,csound,helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)

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

         subroutine updateF(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,&
            ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,invMl,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
            rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
            ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
            listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,u_mapped,tauw,source_term,walave_u,walave_pr,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow,maskMapped(npoin)
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin), invMl(npoin)
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
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime),u_mapped(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(inout)      :: source_term(npoin,ndime+2)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_pr(npoin)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode
            real(rp)                            :: umag, rho_min, rho_avg

            call nvtxStartRange("UpdateF")

            pos = 2

            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,maskMapped,rho(:,pos),q(:,:,pos),E(:,pos),u_buffer)

            !
            ! Apply bcs after update
            !
            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet(q(:,:,2),u(:,:,2),rho(:,2),E(:,2),pr(:,2))
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer,u_mapped)
               call nvtxEndRange
            end if

            if(flag_force_2D) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  q(ipoin,3,pos) =  0.0_rp
               end do
               !$acc end parallel loop
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
            call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),q(:,:,pos),rho(:,pos),&
                                    pr(:,pos),aux_h,Rmass,Rmom,Rener,.true.,-1.0_rp)    
            if(flag_lps_stab) then
               call nvtxStartRange("Stab + Diffu")
               call full_proj_ijk(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),u(:,:,pos),&
                                 Tem(:,pos),Ml,invMl,ProjMass_ls,ProjEner_ls,ProjMX_ls,ProjMY_ls,ProjMZ_ls) 
               call full_diff_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),rho(:,pos),u(:,:,pos),&
                                 Tem(:,pos),mu_fluid,mu_e,mu_sgs,Ml,ProjMass_ls,ProjEner_ls,ProjMX_ls,ProjMY_ls,ProjMZ_ls,tau_stab_ls,Rmass,Rmom,Rener,.false.,-1.0_rp)
               call nvtxEndRange
            else 
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),rho(:,pos),u(:,:,pos),&
                                    Tem(:,pos),mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener,.false.,-1.0_rp)                                 
            end if
            call nvtxEndRange
            !
            ! Call source term if applicable
            !
            
            if(present(source_term) .or.flag_bouyancy_effect) then
               call nvtxStartRange("Extra terms")

               if(flag_trip_element) then
                  call mom_source_rough_elem(npoin,coord,rho(:,pos),u(:,:,pos),source_term(:,3:ndime+2))
                  call ener_source_rough_elem(npoin,coord,rho(:,pos),q(:,:,pos),source_term(:,:))
               end if

               if(flag_bouyancy_effect) then
                  call mom_source_bouyancy_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,rho(:,pos),Rmom,1.0_rp)
                  call ener_source_bouyancy(nelem,npoin,connec,Ngp,dNgp,He,gpvol,q(:,:,pos),Rener,1.0_rp)
               else if(present(source_term)) then
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,1:ndime,pos),source_term(:,3:ndime+2),Rmom,1.0_rp) 
                  call ener_source_const(nelem,npoin,connec,Ngp,dNgp,He,gpvol,source_term(:,2),Rener,1.0_rp)
               end if               
               call nvtxEndRange
            end if

            !
            ! Evaluate wall models

            if((isWallModelOn) ) then
               call nvtxStartRange("WALL MODEL")
               if ((flag_type_wmles == wmles_type_thinBL_fit) .or. (flag_type_wmles == wmles_type_thinBL_fit_hwm)) then
                  call nvtxStartRange("thinBL_WM_init")
                  call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,walave_pr(:),f_eta2,.true.)
                  call nvtxEndRange
               end if
               if((numBoundsWM .ne. 0)) then
                  if((flag_type_wmles == wmles_type_reichardt) .or. (flag_type_wmles == wmles_type_reichardt_hwm)) then
                     call nvtxStartRange("Reichardt_WM")
                     call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,rho(:,pos),walave_u(:,:),tauw,Rmom,-1.0_rp)
                     call nvtxEndRange
                  else if (flag_type_wmles == wmles_type_abl) then
                     call nvtxStartRange("ABL_WM")
                     call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                                          bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                                          rho(:,pos),walave_u(:,:),zo,tauw,Rmom,-1.0_rp)
                     call nvtxEndRange
                  else if ((flag_type_wmles == wmles_type_thinBL_fit) .or. (flag_type_wmles == wmles_type_thinBL_fit_hwm)) then
                     call nvtxStartRange("thinBL_WM")
                     call evalWallModelThinBLFit(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                        rho(:,1),walave_u(:,:),f_eta2,tauw,Rmom,-1.0_rp)
                     call nvtxEndRange
                  end if
               end if
               call nvtxEndRange
            end if

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real_massEnerMom(Rmass(:),Rener(:),Rmom(:,:))
               call nvtxEndRange
            end if

            !
            ! Call lumped mass matrix solver
            !
            call nvtxStartRange("Call solver")
            call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Rmass(:))
            call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Rener(:))
            call lumped_solver_vect_opt(npoin,npoin_w,lpoin_w,invMl,Rmom(:,:))
            call nvtxEndRange

            call nvtxEndRange
            
         end subroutine updateF                 


      end module time_integ_ls
