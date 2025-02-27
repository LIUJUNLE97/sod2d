module time_integ_imex

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
   use mod_operators
   use time_integ, only :  updateBuffer   
   use elem_stab, only : comp_tau



   implicit none

   real(rp), allocatable, dimension(:,:,:) :: Rmom_imex
   real(rp), allocatable, dimension(:,:) :: Rmass_imex, Rener_imex,Reta_imex,Rmom_stab
   real(rp), allocatable, dimension(:,:) :: Rsource_imex,Rwmles_imex
   real(rp), allocatable, dimension(:,:,:) :: Rdiff_mom_imex
   real(rp), allocatable, dimension(:,:) :: f_eta_imex, f_eta_imex2
   real(rp), allocatable, dimension(:,:)   :: Rdiff_mass_imex,Rdiff_ener_imex
   real(rp), allocatable, dimension(:) :: auxReta_imex,aux_h,Rmass_stab,Rener_stab
   real(rp)  , allocatable, dimension(:) 	:: tau_stab_imex
   real(rp)  , allocatable, dimension(:,:) :: ProjMass_imex,ProjEner_imex,ProjMX_imex,ProjMY_imex,ProjMZ_imex
   real(rp), dimension(4,4) :: aij_e, aij_i
   real(rp), dimension(4) :: bij_e, bij_i
   logical :: firstTimeStep = .true.

  contains

   subroutine init_imex_solver(npoin,nelem)

      implicit none
      integer(4),intent(in) :: npoin, nelem

      call nvtxStartRange("Init IMEX solver")

      allocate(Rmom_imex(npoin,ndime,flag_imex_stages),Rmass_stab(npoin),Rener_stab(npoin),Rmom_stab(npoin,ndime))
      !$acc enter data create(Rmom_imex(:,:,:),Rmass_stab(:),Rener_stab(:),Rmom_stab(:,:))

      allocate(Rmass_imex(npoin,flag_imex_stages),Rener_imex(npoin,flag_imex_stages),Reta_imex(npoin,2))
      !$acc enter data create(Rmass_imex(:,:),Rener_imex(:,:),Reta_imex(:,:))

      allocate(Rsource_imex(npoin,ndime+2),Rwmles_imex(npoin,ndime))
      !$acc enter data create(Rsource_imex(:,:),Rwmles_imex(:,:))

      allocate(Rdiff_mom_imex(npoin,ndime,flag_imex_stages))
      !$acc enter data create(Rdiff_mom_imex(:,:,:))

      allocate(auxReta_imex(npoin),f_eta_imex(npoin,ndime),f_eta_imex2(npoin,ndime),aux_h(npoin))
      !$acc enter data create(auxReta_imex(:),f_eta_imex(:,:),f_eta_imex2(:,:),aux_h(:))

      allocate(Rdiff_mass_imex(npoin,flag_imex_stages),Rdiff_ener_imex(npoin,flag_imex_stages))
      !$acc enter data create(Rdiff_mass_imex(:,:),Rdiff_ener_imex(:,:))

      allocate(ProjMass_imex(npoin,ndime),ProjEner_imex(npoin,ndime),ProjMX_imex(npoin,ndime),ProjMY_imex(npoin,ndime),ProjMZ_imex(npoin,ndime),tau_stab_imex(nelem))
      !$acc enter data create(ProjMass_imex(:,:),ProjEner_imex(:,:),ProjMX_imex(:,:),ProjMY_imex(:,:),ProjMZ_imex(:,:),tau_stab_imex(:))
      if(flag_imex_stages .eq. 4) then
         bij_i(1) = 4.0_rp/15.0_rp 
         bij_i(2) = 1.0_rp/3.0_rp 
         bij_i(3) = 7.0_rp/30.0_rp 
         bij_i(4) = 1.0_rp/6.0_rp 

         bij_e(1) = 1.0_rp/4.0_rp 
         bij_e(2) = 0.0_rp 
         bij_e(3) = 3.0_rp/4.0_rp 
         bij_e(4) = 0.0_rp 

         aij_i(1,1) = 0.0_rp
         aij_i(1,2) = 0.0_rp
         aij_i(1,3) = 0.0_rp
         aij_i(1,4) = 0.0_rp
         aij_i(2,1) = bij_i(1)
         aij_i(2,2) = bij_i(1)
         aij_i(2,3) = 0.0_rp
         aij_i(2,4) = 0.0_rp
         aij_i(3,1) = bij_i(1)
         aij_i(3,2) = bij_i(2)
         aij_i(3,3) = 1.0_rp/15.0_rp
         aij_i(3,4) = 0.0_rp
         aij_i(4,1) = bij_i(1)
         aij_i(4,2) = bij_i(2)
         aij_i(4,3) = bij_i(3)
         aij_i(4,4) = bij_i(4)

         aij_e(1,1) = 0.0_rp
         aij_e(1,2) = 0.0_rp
         aij_e(1,3) = 0.0_rp
         aij_e(1,4) = 0.0_rp
         aij_e(2,1) = 8.0_rp/15.0_rp
         aij_e(2,2) = 0.0_rp
         aij_e(2,3) = 0.0_rp
         aij_e(2,4) = 0.0_rp
         aij_e(3,1) = bij_e(1)
         aij_e(3,2) = 5.0_rp/12.0_rp
         aij_e(3,3) = 0.0_rp
         aij_e(3,4) = 0.0_rp
         aij_e(4,1) = bij_e(1)
         aij_e(4,2) = bij_e(2)
         aij_e(4,3) = bij_e(3)
         aij_e(4,4) = bij_e(4)
      else
         bij_i(1) = real(0.398930808264688d0,rp) 
         bij_i(2) = real(0.345755244189623d0,rp)  
         bij_i(3) = real(0.255313947545689d0,rp) 

         bij_e(1) = real(0.398930808264688d0,rp) 
         bij_e(2) = real(0.345755244189623d0,rp)  
         bij_e(3) = real(0.255313947545689d0,rp) 

         aij_i(1,1) = 0.0_rp
         aij_i(1,2) = 0.0_rp
         aij_i(1,3) = 0.0_rp
         aij_i(2,1) = real(0.353842865099275d0,rp)
         aij_i(2,2) = real(0.353842865099275d0,rp)
         aij_i(2,3) = 0.0_rp
         aij_i(3,1) = real(0.398930808264689d0,rp)
         aij_i(3,2) = real(0.345755244189622d0,rp)
         aij_i(3,3) = real(0.255313947545689d0,rp)

         aij_e(1,1) = 0.0_rp
         aij_e(1,2) = 0.0_rp
         aij_e(1,3) = 0.0_rp
         aij_e(2,1) = real(0.711664700366941d0,rp)
         aij_e(2,2) = 0.0_rp
         aij_e(2,3) = 0.0_rp
         aij_e(3,1) = real(0.077338168947683d0,rp)
         aij_e(3,2) = real(0.917273367886007d0,rp)
         aij_e(3,3) = 0.0_rp
      endif

      !$acc enter data copyin(bij_i(:))
      !$acc enter data copyin(bij_e(:))
      !$acc enter data copyin(aij_i(:,:))
      !$acc enter data copyin(aij_e(:,:))

      !$acc kernels
      Rsource_imex(1:npoin,1:ndime+2) = 0.0_rp
      Rwmles_imex(1:npoin,1:ndime) = 0.0_rp
      Rmass_stab(1:npoin) = 0.0_rp
      Rener_stab(1:npoin) = 0.0_rp
      Rmom_stab(1:npoin,1:ndime) = 0.0_rp
      !$acc end kernels

      call nvtxEndRange

   end subroutine init_imex_solver

   subroutine end_imex_solver()
      implicit none


      !$acc exit data delete(Rmom_imex(:,:,:))
      deallocate(Rmom_imex)

      !$acc exit data delete(Rmass_imex(:,:),Rener_imex(:,:),Reta_imex(:,:))
      deallocate(Rmass_imex,Rener_imex,Reta_imex)

      !$acc exit data delete(Rsource_imex(:,:),Rwmles_imex(:,:))
      deallocate(Rsource_imex,Rwmles_imex)

      !$acc exit data delete(auxReta_imex(:),f_eta_imex(:,:))
      deallocate(auxReta_imex,f_eta_imex)

      !$acc exit data delete(Rdiff_mass_imex(:,:),Rdiff_ener_imex(:,:))
      deallocate(Rdiff_mass_imex,Rdiff_ener_imex)

   end subroutine end_imex_solver
        subroutine imex_main(igtime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,u_mapped,tauw,source_term,walave_u,walave_pr,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,save_logFile_next
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
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime), u_mapped(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(in)      :: source_term(npoin,ndime+2)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_pr(npoin)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: istep, ipoin, idime,icode,jstep,ipoin_w
            real(rp)                            :: umag

            if(firstTimeStep .eqv. .true.) then
               firstTimeStep = .false.
               
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta_imex(ipoin_w,idime) = u(ipoin_w,idime,1)*eta(ipoin_w,1)
                  end do
               end do
               !$acc end parallel loop
               
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta_imex,eta(:,1),u(:,:,1),Reta_imex(:,1))
               
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta_imex(:,1))
               end if
               
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta_imex(:,1))     
               
               call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,Reta_imex(:,1),Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,1),u(:,:,1),csound,Tem(:,1),eta(:,1),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)   

               !$acc parallel loop
               do ipoin = 1,npoin_w
                  aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*pr(lpoin_w(ipoin),1)/rho(lpoin_w(ipoin),1)
               end do
               !$acc end parallel loop
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),q(:,:,1),rho(:,1),pr(:,1),&
                                    aux_h(:),Rmass_imex(:,1),Rmom_imex(:,:,1),Rener_imex(:,1))
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,1),rho(:,1),u(:,:,1),&
                                       Tem(:,1),mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass_imex(:,1),Rdiff_mom_imex(:,:,1),Rdiff_ener_imex(:,1))

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
               !$acc parallel loop
               do ipoin = 1,npoin
                  Rmass_imex(ipoin,1)      =  Rmass_imex(ipoin,flag_imex_stages)
                  Rdiff_mass_imex(ipoin,1) =  Rdiff_mass_imex(ipoin,flag_imex_stages)
                  Rener_imex(ipoin,1)      =  Rener_imex(ipoin,flag_imex_stages)
                  Rdiff_ener_imex(ipoin,1) =  Rdiff_ener_imex(ipoin,flag_imex_stages)                  
                  !$acc loop seq   
                  do idime = 1,ndime
                     Rmom_imex(ipoin,idime,1)      =  Rmom_imex(ipoin,idime,flag_imex_stages)
                     Rdiff_mom_imex(ipoin,idime,1) =  Rdiff_mom_imex(ipoin,idime,flag_imex_stages)  
                  end do
               end do
               !$acc end parallel loop
            end if

            if(flag_lps_stab) then
               call comp_tau(nelem,npoin,connec,csound,u(:,:,2),helem,dt,tau_stab_imex)
               call full_proj_ijk(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,1),u(:,:,1),&
                              Tem(:,1),Ml,ProjMass_imex,ProjEner_imex,ProjMX_imex,ProjMY_imex,ProjMZ_imex)
               call full_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,1),rho(:,1),u(:,:,1),&
                              Tem(:,1),Ml,ProjMass_imex,ProjEner_imex,ProjMX_imex,ProjMY_imex,ProjMZ_imex,tau_stab_imex,Rmass_stab,Rmom_stab,Rener_stab,.true.,-1.0_rp) 
            end if

            if(isWallModelOn) then
               !$acc kernels
               Rwmles_imex(1:npoin,1:ndime) = 0.0_rp
               !$acc end kernels
               if ((flag_type_wmles == wmles_type_thinBL_fit) .or. (flag_type_wmles == wmles_type_thinBL_fit_hwm)) then
                  call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,walave_pr(:),f_eta_imex,.true.)
               end if
               if(numBoundsWM .ne. 0) then
                  if((flag_type_wmles == wmles_type_reichardt) .or. (flag_type_wmles == wmles_type_reichardt_hwm)) then
                     call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                        rho(:,2),walave_u(:,:),tauw,Rwmles_imex)
                  else if (flag_type_wmles == wmles_type_abl) then
                     call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                        rho(:,2),walave_u(:,:),zo,tauw,Rwmles_imex)
                  else if ((flag_type_wmles == wmles_type_thinBL_fit) .or. (flag_type_wmles == wmles_type_thinBL_fit_hwm)) then
                     call evalWallModelThinBLFit(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                        rho(:,1),walave_u(:,:),f_eta_imex,tauw,Rwmles_imex)
                  end if  
               end if
            end if

            do istep = 2,flag_imex_stages 

               if(present(source_term) .or.flag_bouyancy_effect) then
                  !$acc kernels
                  Rsource_imex(1:npoin,1:ndime+2) = 0.0_rp
                  !$acc end kernels
                  if(flag_bouyancy_effect) then
                     call mom_source_bouyancy_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,rho(:,2),Rsource_imex(:,3:ndime+2),-1.0_rp)
                     call ener_source_bouyancy(nelem,npoin,connec,Ngp,dNgp,He,gpvol,q(:,:,2),Rsource_imex(:,2),-1.0_rp)
                  else if(present(source_term)) then
                     call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,1:ndime,1),source_term(:,3:ndime+2),Rsource_imex(:,3:ndime+2),-1.0_rp)
                     call ener_source_const(nelem,npoin,connec,Ngp,dNgp,He,gpvol,source_term(:,2),Rsource_imex(:,2),-1.0_rp)
                  end if
               end if

               !$acc parallel loop
               do ipoin = 1,npoin
                  rho(ipoin,2) =  -dt*Rsource_imex(ipoin,1)-dt*Rmass_stab(ipoin)
                  E(ipoin,2)   =  -dt*Rsource_imex(ipoin,2)-dt*Rener_stab(ipoin)
                  !$acc loop seq   
                  do idime = 1,ndime
                     q(ipoin,idime,2) =  -dt*(Rsource_imex(ipoin,idime+2)+Rwmles_imex(ipoin,idime)+Rmom_stab(ipoin,idime))
                  end do
                  do jstep = 1, istep-1
                     rho(ipoin,2) = rho(ipoin,2) -dt*aij_e(istep,jstep)*(Rmass_imex(ipoin,jstep)+Rsource_imex(ipoin,1))-dt*aij_i(istep,jstep)*Rdiff_mass_imex(ipoin,jstep)
                     E(ipoin,2)   = E(ipoin,2)   -dt*aij_e(istep,jstep)*(Rener_imex(ipoin,jstep)+Rsource_imex(ipoin,2))-dt*aij_i(istep,jstep)*Rdiff_ener_imex(ipoin,jstep)
                     do idime = 1,ndime
                        q(ipoin,idime,2) = q(ipoin,idime,2) -dt*aij_e(istep,jstep)*(Rmom_imex(ipoin,idime,jstep)+Rsource_imex(ipoin,idime+2)+Rwmles_imex(ipoin,idime)) &
                                                            -dt*aij_i(istep,jstep)*Rdiff_mom_imex(ipoin,idime,jstep)                    
                     end do
                  end do
               end do
               !$acc end parallel loop
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real_massEnerMom(rho(:,2),E(:,2),q(:,:,2))
               end if
               
               !$acc parallel loop
               do ipoin = 1,npoin
                  rho(ipoin,2) =  rho(ipoin,2) + rho(ipoin,1)*Ml(ipoin)
                  E(ipoin,2) =  E(ipoin,2) + E(ipoin,1)*Ml(ipoin)
                  !$acc loop seq   
                  do idime = 1,ndime
                     q(ipoin,idime,2) =  q(ipoin,idime,2) + q(ipoin,idime,1)*Ml(ipoin)
                  end do
               end do
               !$acc end parallel loop
               !if(mpi_rank.eq.0) write(111,*)   " before cg"
               call nvtxStartRange("IMEX_CONJ_GRAD")
               call conjGrad_imex(aij_i(istep,istep),igtime,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,numBoundsWM,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                                 dlxigp_ip,He,gpvol,Ngp,Ml,helem,gamma_gas,Rgas,Cp,Prt,csound,mu_fluid,mu_e,mu_sgs,rho(:,1),rho(:,2),E(:,1),E(:,2),q(:,:,1),q(:,:,2), &
                                 ndof,nbnodes,ldof,lbnodes,lnbn_nodes,bound,bou_codes,bou_codes_nodes,listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer)
               call nvtxEndRange
               !if(mpi_rank.eq.0) write(111,*)   " after cg"
               
               !call limit_rho(nelem,npoin,connec,rho(:,2),epsilon(umag))

               if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,maskMapped,rho(:,2),q(:,:,2),E(:,2),u_buffer)

               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  if(isMappedFaces.and.isMeshPeriodic) call copy_periodicNodes_for_mappedInlet(q(:,:,2),u(:,:,2),rho(:,2),E(:,2),pr(:,2))
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,2),q(:,:,2),u(:,:,2),pr(:,2),E(:,2),u_buffer,u_mapped)
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
                  aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*pr(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2)
                  csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2))
                  umag = sqrt(umag)
                  machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
                  Tem(lpoin_w(ipoin),2) = pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)*Rgas)
               end do
               !$acc end parallel loop
               
               !call comp_tau(nelem,npoin,connec,csound,u(:,:,2),helem,dt,tau_stab_imex)

               if(flag_total_enthalpy .eqv. .true.) then
                  call full_convec_ijk_H(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),q(:,:,2),rho(:,2),pr(:,2),&
                                    E(:,2),Rmass_imex(:,istep),Rmom_imex(:,:,istep),Rener_imex(:,istep))
               else
                  call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),q(:,:,2),rho(:,2),pr(:,2),&
                                    aux_h(:),Rmass_imex(:,istep),Rmom_imex(:,:,istep),Rener_imex(:,istep))
               end if
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,2),rho(:,2),u(:,:,2),&
                                    Tem(:,2),mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass_imex(:,istep),Rdiff_mom_imex(:,:,istep),Rdiff_ener_imex(:,istep))
            end do
            !if(mpi_rank.eq.0) write(111,*)   " after in"
      
            !$acc parallel loop
            do ipoin = 1,npoin_w
               eta(lpoin_w(ipoin),1) = eta(lpoin_w(ipoin),2)
               eta(lpoin_w(ipoin),2) = (rho(lpoin_w(ipoin),2)/(gamma_gas-1.0_rp))* &
                  log(max(pr(lpoin_w(ipoin),2),0.0_rp)/(rho(lpoin_w(ipoin),2)**gamma_gas))
               !$acc loop seq
               do idime = 1,ndime
                  f_eta_imex(lpoin_w(ipoin),idime)  = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
                  f_eta_imex2(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,2)*eta(lpoin_w(ipoin),2)
               end do
            end do
            !$acc end parallel loop
#if 1
            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta_imex,eta(:,1),u(:,:,1),Reta_imex(:,2))

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(Reta_imex(:,2))
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta_imex(:,2))

            call nvtxStartRange("Entropy residual")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               auxReta_imex(lpoin_w(ipoin)) = (1.5_rp*Reta_imex(lpoin_w(ipoin),2)-0.5_rp*Reta_imex(lpoin_w(ipoin),1)) + &
                                             factor_comp*(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
               Reta_imex(lpoin_w(ipoin),1) = Reta_imex(lpoin_w(ipoin),2)            
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,auxReta_imex)
               call nvtxEndRange
            end if

            !
            ! Compute entropy viscosity
            !
            call nvtxStartRange("Entropy viscosity evaluation")
            call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,auxReta_imex,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
            call nvtxEndRange  
#else
            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
            gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta_imex2,eta(:,2),u(:,:,2),Reta_imex(:,2))

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(Reta_imex(:,2))
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta_imex(:,2))

            call generic_scalar_convec_projection_residual_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta_imex,eta(:,1),u(:,:,1),Reta_imex(:,2),auxReta_imex)

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(auxReta_imex)
            end if
            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,auxReta_imex)    

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,auxReta_imex)
               call nvtxEndRange
            end if
            !
            ! Compute entropy viscosity
            !
            call nvtxStartRange("Entropy viscosity evaluation")
            call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,auxReta_imex,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
            call nvtxEndRange   

#endif
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
         end subroutine imex_main
      end module time_integ_imex
