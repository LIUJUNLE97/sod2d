module time_integ_species_imex

   use mod_nvtx
   use elem_convec_species
   use elem_convec , only : generic_scalar_convec_ijk, generic_scalar_convec_projection_residual_ijk
   use elem_diffu_species
   use elem_stab_species
   use mod_solver_species
   use mod_entropy_viscosity_species
   use mod_numerical_params
   use mod_bc_routines_species
   use mod_operators
   use time_integ, only :  limit_rho

   implicit none

   real(rp), allocatable, dimension(:,:,:) :: RYk,f_eta
   real(rp), allocatable, dimension(:,:) ::Reta,Rflux,gradYk
   real(rp), allocatable, dimension(:) :: auxReta,Rstab,tau
   real(rp), allocatable, dimension(:)   :: beta, alpha
   real(rp) :: gamma0

   contains
   subroutine init_imex_species_solver(npoin,nelem)

      implicit none
      integer(4),intent(in) :: npoin,nelem
      integer(4) :: numSteps

      call nvtxStartRange("Init species solver")

      allocate(RYk(npoin,nspecies,3),f_eta(npoin,ndime,2))
      !$acc enter data create(RYk(:,:,:),f_eta(:,:,:))

      allocate(auxReta(npoin),Reta(npoin,3),gradYk(npoin,ndime),Rstab(npoin),tau(nelem))
      !$acc enter data create(auxReta(:),Reta(:,:),gradYk(:,:),Rstab(:),tau)
   
      !$acc kernels
      RYk(1:npoin,1:nspecies,1:3) = 0.0_rp
      Rstab(1:npoin) = 0.0_rp
      !$acc end kernels

      allocate(alpha(3),beta(3))
      !$acc enter data create(alpha(:),beta(:))
      call nvtxEndRange

   end subroutine init_imex_species_solver

   subroutine end_imex_species_solver()
      implicit none


      !TODO

   end subroutine end_imex_species_solver
 
         subroutine imex_species_main(ispc,igtime,iltime,save_logFile_next,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
                         rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer)  ! Optional args

            implicit none

            integer(4),           intent(in)    :: ispc
            logical,              intent(in)    :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime,iltime,save_logFile_next
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
            real(rp),             intent(in)    :: Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,4)
            real(rp),             intent(inout) :: u(npoin,ndime,4)
            real(rp),             intent(inout) :: Yk(npoin,nspecies,4)
            real(rp),             intent(inout) :: eta_Yk(npoin,nspecies,4)
            real(rp),             intent(inout) :: mu_fluid(npoin)
            real(rp),             intent(inout) :: mu_e_Yk(nelem,ngaus,nspecies)
            real(rp),             intent(inout) :: mu_sgs(nelem,ngaus)
            real(rp),             intent(inout) :: mue_l(nelem,nnode)
            real(rp),             intent(in)    :: coord(npoin,ndime)
            real(rp),             intent(in)  ::  wgp(ngaus)
            integer(4),            intent(in)    :: numBoundsWM
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: Yk_buffer(npoin,nspecies)
            integer(4)                          :: istep,ipoin,idime,icode,iPer,ipoin_w

            call nvtxStartRange("AB2 init")
            if(iltime .eq. 1) then
               call species_convec_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,1),Yk(:,ispc,1),u(:,:,1),RYk(:,ispc,1))               

               if(mpi_size.ge.2) then
                  call nvtxStartRange("AB2 halo update")
                  call mpi_halo_atomic_update_real(RYk(:,ispc,1))
                  call nvtxEndRange
               end if    
               if(flag_entropy_stab_in_species .eqv. .true.) then 
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     ipoin_w = lpoin_w(ipoin)
                     !$acc loop seq
                     do idime = 1,ndime
                        f_eta(ipoin_w,idime,1) = u(ipoin_w,idime,1)*eta_Yk(ipoin_w,ispc,1)
                     end do
                  end do
                  !$acc end parallel loop

                  call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                     gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:,1),eta_Yk(:,ispc,1),u(:,:,1),Reta(:,1))

                  if(mpi_size.ge.2) then
                     call mpi_halo_atomic_update_real(Reta(:,1))
                  end if

                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,1))
                  
                  call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta(:,1),Ngp,coord,dNgp,gpvol,wgp, &
                     rho(:,1),u(:,:,1),eta_Yk(:,ispc,1),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
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

            call species_tau(nelem,npoin,connec,u(:,:,1),helem,dt,tau)

            call nvtxStartRange("AB2 species")

            call species_convec_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,1),Yk(:,ispc,1),u(:,:,1),RYk(:,ispc,2))               

            call nvtxEndRange

            if(mpi_size.ge.2) then
               call nvtxStartRange("AB2 halo update")
               call mpi_halo_atomic_update_real(RYk(:,ispc,2))
               call nvtxEndRange
            end if

            call nvtxStartRange("AB2 update Yk(2) & RYk(1)")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoin_w = lpoin_w(ipoin)
               Yk(ipoin_w,ispc,2) = -(beta(1)*RYk(ipoin_w,ispc,2)+beta(2)*RYk(ipoin_w,ispc,1)+beta(3)*RYk(ipoin_w,ispc,3))
               Yk(ipoin_w,ispc,2) = Cp*rho(ipoin_w,2)*Ml(ipoin_w)*(dt*Yk(ipoin_w,ispc,2)/Ml(ipoin_w) + alpha(1)*Yk(ipoin_w,ispc,1) + alpha(2)*Yk(ipoin_w,ispc,3) + alpha(3)*Yk(ipoin_w,ispc,4))/gamma0
               RYk(ipoin_w,ispc,3) = RYk(ipoin_w,ispc,1)
               RYk(ipoin_w,ispc,1) = RYk(ipoin_w,ispc,2)
            end do
            !$acc end parallel loop

            call nvtxEndRange
            
            call conjGrad_species(ispc,igtime,1.0_rp/gamma0,dt,save_logFile_next,noBoundaries,nelem,npoin,npoin_w,nboun,connec,lpoin_w,invAtoIJK,&
                             gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,mu_fluid,mu_e_Yk(:,:,ispc),mu_sgs,tau,Cp,Prt,rho(:,2),Yk(:,ispc,1),Yk(:,ispc,2),&
                             bou_codes,bound,nbnodes,lbnodes,lnbn_nodes,bou_codes_nodes,normalsAtNodes,Yk_buffer)

            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Yk(:,ispc,2),Yk_buffer(:,ispc))
            end if

            !
            ! Compute subgrid viscosity if active
            !

            if(flag_entropy_stab_in_species .eqv. .true.) then 
#if 1            
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  eta_Yk(ipoin_w,ispc,2) = 0.5_rp*Yk(ipoin_w,ispc,2)*Yk(ipoin_w,ispc,2)
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(ipoin_w,idime,1) = u(ipoin_w,idime,1)*eta_Yk(ipoin_w,ispc,1)
                  end do
               end do
               !$acc end parallel loop

               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:,1),eta_Yk(:,ispc,1),u(:,:,1),Reta(:,2))
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:,2))
               end if

               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,2))

               !$acc parallel loop
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  auxReta(ipoin_w) = (beta(1)*Reta(ipoin_w,2)+beta(2)*Reta(ipoin_w,1)+beta(3)*Reta(ipoin_w,3)) !+ &
                                     !(gamma0*eta_Yk(lpoin_w(ipoin),ispc,2)-alpha(1)*eta_Yk(lpoin_w(ipoin),ispc,1)-alpha(2)*eta_Yk(lpoin_w(ipoin),ispc,3)-alpha(3)*eta_Yk(lpoin_w(ipoin),ispc,4))/dt
                  Reta(ipoin_w,3) = Reta(ipoin_w,1)
                  Reta(ipoin_w,1) = Reta(ipoin_w,2)
               end do
               !$acc end parallel loop

               if (noBoundaries .eqv. .false.) then
                  call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,auxReta)
               end if

               call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
                                             rho(:,2),u(:,:,2),eta_Yk(:,ispc,2),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
#else
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  ipoin_w = lpoin_w(ipoin)
                  eta_Yk(ipoin_w,ispc,2) = 0.5_rp*Yk(ipoin_w,ispc,2)*Yk(ipoin_w,ispc,2)
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(ipoin_w,idime,2) = u(ipoin_w,idime,2)*eta_Yk(ipoin_w,ispc,2)
                     f_eta(ipoin_w,idime,1) = u(ipoin_w,idime,1)*eta_Yk(ipoin_w,ispc,1)
                  end do
               end do
               !$acc end parallel loop

               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:,2),eta_Yk(:,ispc,2),u(:,:,2),Reta(:,2))
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:,2))
               end if
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,2))

               call generic_scalar_convec_projection_residual_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:,1),eta_Yk(:,ispc,1),u(:,:,1),Reta(:,2),auxReta)
               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(auxReta)
               end if
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,auxReta)               

               if (noBoundaries .eqv. .false.) then
                  call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,auxReta)
               end if

               call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
                                             rho(:,2),u(:,:,2),eta_Yk(:,ispc,2),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)  
#endif
            end if
         end subroutine imex_species_main

       end module time_integ_species_imex
