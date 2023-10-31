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
   use time_integ, only :  updateBuffer


   implicit none

   real(rp), allocatable, dimension(:,:,:) :: Rmom_imex
   real(rp), allocatable, dimension(:,:) :: Rmass_imex, Rener_imex,Reta_imex
   real(rp), allocatable, dimension(:,:) :: Rsource_imex,Rwmles_imex
   real(rp), allocatable, dimension(:,:,:) :: Rdiff_mom_imex
   real(rp), allocatable, dimension(:,:) :: f_eta_imex
   real(rp), allocatable, dimension(:,:)   :: Rdiff_mass_imex,Rdiff_ener_imex
   real(rp), allocatable, dimension(:) :: auxReta_imex
   integer(4), parameter :: numSteps = 5
   !integer(4), parameter :: numSteps = 4
   real(rp), dimension(numSteps,numSteps) :: aij_e, aij_i
   real(rp), dimension(numSteps) :: bij_e, bij_i
   logical :: firstTimeStep = .true.

  contains

   subroutine init_imex_solver(npoin)

      implicit none
      integer(4),intent(in) :: npoin

      allocate(Rmom_imex(npoin,ndime,numSteps))
      !$acc enter data create(Rmom_imex(:,:,:))

      allocate(Rmass_imex(npoin,numSteps),Rener_imex(npoin,numSteps),Reta_imex(npoin,2))
      !$acc enter data create(Rmass_imex(:,:),Rener_imex(:,:),Reta_imex(:,:))

      allocate(Rsource_imex(npoin,ndime),Rwmles_imex(npoin,ndime))
      !$acc enter data create(Rsource_imex(:,:),Rwmles_imex(:,:))

      allocate(Rdiff_mom_imex(npoin,ndime,numSteps))
      !$acc enter data create(Rdiff_mom_imex(:,:,:))

      allocate(auxReta_imex(npoin),f_eta_imex(npoin,ndime))
      !$acc enter data create(auxReta_imex(:),f_eta_imex(:,:))

      allocate(Rdiff_mass_imex(npoin,numSteps),Rdiff_ener_imex(npoin,numSteps))
      !$acc enter data create(Rdiff_mass_imex(:,:),Rdiff_ener_imex(:,:))
#if 1
      bij_i(1) = 0.0_rp
      bij_i(2) = 1.5_rp
      bij_i(3) = -1.5_rp 
      bij_i(4) = 0.5_rp
      bij_i(5) = 0.5_rp

      bij_e(1) = 0.25_rp
      bij_e(2) = 1.75_rp
      bij_e(3) = 0.75_rp
      bij_e(4) = -1.75_rp
      bij_e(5) = 0.0_rp


      aij_i(1,1) = 0.0_rp
      aij_i(1,2) = 0.0_rp
      aij_i(1,3) = 0.0_rp
      aij_i(1,4) = 0.0_rp
      aij_i(1,5) = 0.0_rp
      aij_i(2,1) = 0.0_rp
      aij_i(2,2) = 0.5_rp
      aij_i(2,3) = 0.0_rp
      aij_i(2,4) = 0.0_rp
      aij_i(2,5) = 0.0_rp
      aij_i(3,1) = 0.0_rp
      aij_i(3,2) = 0.1666666666_rp
      aij_i(3,3) = 0.5_rp
      aij_i(3,4) = 0.0_rp
      aij_i(3,5) = 0.0_rp
      aij_i(4,1) = 0.0_rp
      aij_i(4,2) = -0.5_rp
      aij_i(4,3) = 0.5_rp
      aij_i(4,4) = 0.5_rp
      aij_i(4,5) = 0.0_rp
      aij_i(5,1) = 0.0_rp
      aij_i(5,2) = 1.5_rp
      aij_i(5,3) = -1.5_rp
      aij_i(5,4) = 0.5_rp
      aij_i(5,5) = 0.5_rp

      aij_e(1,1) = 0.0_rp
      aij_e(1,2) = 0.0_rp
      aij_e(1,3) = 0.0_rp
      aij_e(1,4) = 0.0_rp
      aij_e(1,5) = 0.0_rp
      aij_e(2,1) = 0.5_rp
      aij_e(2,2) = 0.0_rp
      aij_e(2,3) = 0.0_rp
      aij_e(2,4) = 0.0_rp
      aij_e(2,5) = 0.0_rp
      aij_e(3,1) = 0.6111111111_rp
      aij_e(3,2) = 0.0555555555_rp
      aij_e(3,3) = 0.0_rp
      aij_e(3,4) = 0.0_rp
      aij_e(3,5) = 0.0_rp
      aij_e(4,1) = 0.8333333333_rp
      aij_e(4,2) = -0.8333333333_rp
      aij_e(4,3) = 0.5_rp
      aij_e(4,4) = 0.0_rp
      aij_e(4,5) = 0.0_rp
      aij_e(5,1) = 0.25_rp
      aij_e(5,2) = 1.75_rp
      aij_e(5,3) = 0.75_rp
      aij_e(5,4) = -1.75_rp
      aij_e(5,5) = 0.0_rp
#else
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
#endif

      !$acc update device(bij_i(:))
      !$acc update device(bij_e(:))
      !$acc update device(aij_i(:,:))
      !$acc update device(aij_e(:,:))

      !$acc kernels
      Rsource_imex(1:npoin,1:ndime) = 0.0_rp
      Rwmles_imex(1:npoin,1:ndime) = 0.0_rp
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
        subroutine imex_main(igtime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,lelpn,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u,zo)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,save_logFile_next
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin),lelpn(npoin)
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
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: istep, ipoin, idime,icode,jstep
            real(rp)                            :: umag

            !if(mpi_rank.eq.0) write(111,*)   " all in"

            if(firstTimeStep .eqv. .true.) then
               firstTimeStep = .false.

               !$acc parallel loop
               do ipoin = 1,npoin_w
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta_imex(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
                  end do
               end do
               !$acc end parallel loop
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta_imex,eta(:,1),u(:,:,1),Reta_imex(:,1))

               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta_imex(:,1))
               end if
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta_imex(:,1))   
             
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,1),q(:,:,1),rho(:,1),pr(:,1),E(:,1),Rmass_imex(:,1),Rmom_imex(:,:,1),Rener_imex(:,1))     
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,1),rho(:,1),u(:,:,1),&
                                       Tem(:,1),mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass_imex(:,1),Rdiff_mom_imex(:,:,1),Rdiff_ener_imex(:,1))
            else 
               !$acc parallel loop
               do ipoin = 1,npoin
                  Rmass_imex(ipoin,1)      =  Rmass_imex(ipoin,numSteps)
                  Rdiff_mass_imex(ipoin,1) =  Rdiff_mass_imex(ipoin,numSteps)
                  Rener_imex(ipoin,1)      =  Rener_imex(ipoin,numSteps)
                  Rdiff_ener_imex(ipoin,1) =  Rdiff_ener_imex(ipoin,numSteps)                  
                  !$acc loop seq   
                  do idime = 1,ndime
                     Rmom_imex(ipoin,idime,1)      =  Rmom_imex(ipoin,idime,numSteps)
                     Rdiff_mom_imex(ipoin,idime,1) =  Rdiff_mom_imex(ipoin,idime,numSteps)  
                  end do
               end do
               !$acc end parallel loop
            end if

            if(present(source_term)) then
               !$acc kernels
               Rsource_imex(1:npoin,1:ndime) = 0.0_rp
               !$acc end kernels
               call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,1:ndime,1),source_term,Rsource_imex)
            end if

            if(isWallModelOn) then
                  !$acc kernels
                  Rwmles_imex(1:npoin,1:ndime) = 0.0_rp
                  !$acc end kernels
                  if(numBoundsWM .ne. 0) then
                     if(flag_type_wmles == wmles_type_reichardt) then
                        call evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),tauw,Rwmles_imex)
                     else if (flag_type_wmles == wmles_type_abl) then
                        call evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,&
                           rho(:,1),walave_u(:,:),zo,tauw,Rwmles_imex)
                     end if
                  end if
            end if
            !if(mpi_rank.eq.0) write(111,*)   " before in"
            do istep = 2,numSteps 
               !$acc parallel loop
               do ipoin = 1,npoin
                  rho(ipoin,2) =  0.0_rp
                  E(ipoin,2) =  0.0_rp
                  !$acc loop seq   
                  do idime = 1,ndime
                     q(ipoin,idime,2) =  -dt*(Rsource_imex(ipoin,idime)+Rwmles_imex(ipoin,idime))
                  end do
                  do jstep = 1, istep-1
                     rho(ipoin,2) = rho(ipoin,2) -dt*aij_e(istep,jstep)*Rmass_imex(ipoin,jstep)-dt*aij_i(istep,jstep)*Rdiff_mass_imex(ipoin,jstep)
                     E(ipoin,2)   = E(ipoin,2)   -dt*aij_e(istep,jstep)*Rener_imex(ipoin,jstep)-dt*aij_i(istep,jstep)*Rdiff_ener_imex(ipoin,jstep)
                     do idime = 1,ndime
                        q(ipoin,idime,2) = q(ipoin,idime,2) -dt*aij_e(istep,jstep)*Rmom_imex(ipoin,idime,jstep) &
                                                            -dt*aij_i(istep,jstep)*Rdiff_mom_imex(ipoin,idime,jstep)
                                           
                     end do
                  end do
               end do
               !$acc end parallel loop
               !$acc end parallel loop     
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(rho(:,2))
                  call mpi_halo_atomic_update_real(E(:,2))
                  do idime = 1,ndime
                     call mpi_halo_atomic_update_real(q(:,idime,2))
                  end do
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
               call conjGrad_imex(aij_i(istep,istep),igtime,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,numBoundsWM,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                                 dlxigp_ip,He,gpvol,Ngp,Ml,gamma_gas,Rgas,Cp,Prt,mu_fluid,mu_e,mu_sgs,rho(:,1),rho(:,2),E(:,1),E(:,2),q(:,:,1),q(:,:,2), &
                                 ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer) 
               !if(mpi_rank.eq.0) write(111,*)   " after cg"

               if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,2),q(:,:,2),u_buffer)

               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,2),q(:,:,2),u(:,:,2),pr(:,2),E(:,2),u_buffer)
                  call nvtxEndRange
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
               end do
               !$acc end parallel loop

               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,2),q(:,:,2),rho(:,2),pr(:,2),&
                                    E(:,2),Rmass_imex(:,istep),Rmom_imex(:,:,istep),Rener_imex(:,istep))
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,2),rho(:,2),u(:,:,2),&
                                    Tem(:,2),mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass_imex(:,istep),Rdiff_mom_imex(:,:,istep),Rdiff_ener_imex(:,istep))
         
            end do
            !if(mpi_rank.eq.0) write(111,*)   " after in"
          
            !$acc parallel loop
            do ipoin = 1,npoin_w
               eta(lpoin_w(ipoin),1) = eta(lpoin_w(ipoin),2)
               eta(lpoin_w(ipoin),2) = (rho(lpoin_w(ipoin),2)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)**gamma_gas))
               !$acc loop seq
               do idime = 1,ndime
                  f_eta_imex(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,2)*eta(lpoin_w(ipoin),2)
               end do
            end do
            !$acc end parallel loop

            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta_imex,eta(:,2),u(:,:,2),Reta_imex(:,2))

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(Reta_imex(:,2))
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta_imex(:,2))

            call nvtxStartRange("Entropy residual")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               auxReta_imex(lpoin_w(ipoin)) = (1.5_rp*Reta_imex(lpoin_w(ipoin),2)-0.5_rp*Reta_imex(lpoin_w(ipoin),1)) !+ &
                                              !(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
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
            ! If using Sutherland viscosity model:
            !
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call nvtxStartRange("Sutherland viscosity")
               call sutherland_viscosity(npoin,Tem(:,2),mu_factor,mu_fluid)
               call nvtxEndRange
            end if

            call nvtxStartRange("Entropy viscosity evaluation")
            !
            ! Compute entropy viscosity
            !
            call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,auxReta_imex,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
            call nvtxEndRange
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