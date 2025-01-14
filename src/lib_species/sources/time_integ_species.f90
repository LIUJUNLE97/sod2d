
module time_integ_species

   use mod_nvtx
   use elem_convec_species
   use elem_convec , only : generic_scalar_convec_ijk, generic_scalar_convec_projection_residual_ijk
   use elem_diffu_species
   use mod_solver
   use mod_entropy_viscosity_species
   use mod_numerical_params
   use mod_bc_routines_species
   use mod_operators
   use elem_stab_species
   use mod_wall_model_species

   implicit none

   real(rp), allocatable, dimension(:,:) :: Reta
   real(rp), allocatable, dimension(:,:) :: f_eta,gradYk
   real(rp), allocatable, dimension(:) :: Rspc
   real(rp), allocatable, dimension(:) :: auxReta,tau

   real(rp), allocatable, dimension(:,:) :: lambda_ij, gamma_ij
   logical :: firstTimeStep = .true.

   contains

   subroutine init_rk4_ls_species_solver(npoin,nelem)
      implicit none
      integer(4),intent(in) :: npoin,nelem

      write(*,*) "--| Allocating species"

      allocate(auxReta(npoin))
      !$acc enter data create(auxReta(:))

      allocate(f_eta(npoin,ndime),Reta(npoin,2),Rspc(npoin),tau(nelem),gradYk(npoin,ndime))
      !$acc enter data create(Rspc(:))
      !$acc enter data create(tau(:))
      !$acc enter data create(f_eta(:,:))
      !$acc enter data create(gradYk(:,:))
      !$acc enter data create(Reta(:,:))

      allocate(lambda_ij(6,6),gamma_ij(6,6))
      !$acc enter data create(lambda_ij(:,:))
      !$acc enter data create(gamma_ij(:,:))

      lambda_ij(:,:) = 0.0_rp         
      lambda_ij(2,1) = 1.0_rp
      lambda_ij(3,2) = 1.0_rp
      lambda_ij(4,3) = 1.0_rp
      lambda_ij(5,1) = real(0.571403511494104d0,rp)
      lambda_ij(5,4) = real(0.428596488505896d0,rp)
      lambda_ij(6,5) = 1.0_rp

      gamma_ij(:,:) = 0.0_rp
      gamma_ij(2,1) = real(0.443568244942995d0,rp)
      gamma_ij(3,2) = real(0.291111420073766d0,rp)
      gamma_ij(4,3) = real(0.270612601278217d0,rp)
      gamma_ij(5,4) = real(0.110577759392786d0,rp)
      gamma_ij(6,5) = real(0.458557505351052d0,rp)
     
      !$acc update device(lambda_ij(:,:))
      !$acc update device(gamma_ij(:,:))

      call nvtxEndRange

   end subroutine init_rk4_ls_species_solver

   subroutine end_rk4_ls_species_solver()
      implicit none

      !$acc exit data delete(auxReta(:))
      deallocate(auxReta)

      !$acc exit data delete(Rspc(:))
      !$acc exit data delete(f_eta(:,:))
      !$acc exit data delete(Reta(:,:))
      deallocate(f_eta,Reta,Rspc)

      !$acc exit data delete(lambda_ij(:,:))
      !$acc exit data delete(gamma_ij(:,:))
      deallocate(lambda_ij,gamma_ij)

   end subroutine end_rk4_ls_species_solver

   subroutine rk_4_ls_species_main(ispc,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                   ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
                   rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
                   ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                   listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer,walave_u,walave_t,zo)  ! Optional args

      implicit none

      integer(4),           intent(in)    :: ispc
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
      real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_t(npoin)
      real(rp), optional, intent(in)      :: zo(npoin)
      integer(4)                          :: pos
      integer(4)                          :: istep, ipoin, idime,icode
      real(rp)                            :: umag, rho_min, rho_avg

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! New version of RK4 using loops                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      pos = 2 ! Set correction as default value

      call species_tau(nelem,npoin,connec,u(:,:,pos),helem,dt,tau)
      
      if(firstTimeStep .eqv. .true.) then
         firstTimeStep = .false.
         call updateFspecies(ispc,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                              ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
                              rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
                              ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                              listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer,walave_u,walave_t,zo)
      end if           

      !$acc parallel loop
      do ipoin = 1,npoin_w
         Yk(lpoin_w(ipoin),ispc,pos) = Yk(lpoin_w(ipoin),ispc,pos) + gamma_ij(2,1)*dt*Rspc(lpoin_w(ipoin))              
      end do
      !$acc end parallel loop

      call nvtxStartRange("Loop over RK steps")
      do istep = 2,5
         call updateFspecies(ispc,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                              ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
                              rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
                              ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                              listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer,walave_u,walave_t,zo)

         !
         ! Accumulate the residuals
         !
         call nvtxStartRange("Accumulate residuals")
         !$acc parallel loop
         do ipoin = 1,npoin_w
            Yk(lpoin_w(ipoin),ispc,pos) = lambda_ij(istep+1,1)*Yk(lpoin_w(ipoin),ispc,1) + lambda_ij(istep+1,istep)*Yk(lpoin_w(ipoin),ispc,pos) + gamma_ij(istep+1,istep)*dt*Rspc(lpoin_w(ipoin))                  
         end do
         !$acc end parallel loop
         call nvtxEndRange
      end do
      call nvtxEndRange
      !
      ! Apply bcs after update
      !
      if (noBoundaries .eqv. .false.) then
         call nvtxStartRange("BCS_AFTER_UPDATE")
         call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Yk(:,ispc,pos),Yk_buffer(:,ispc))
         call nvtxEndRange
      end if

      !
      ! Update velocity and equations of state
      !
      call nvtxStartRange("Update u and EOS")

      !$acc parallel loop
      do ipoin = 1,npoin_w
         eta_Yk(lpoin_w(ipoin),ispc,1) = eta_Yk(lpoin_w(ipoin),ispc,2)
         eta_Yk(lpoin_w(ipoin),ispc,2) = 0.5_rp*Yk(lpoin_w(ipoin),ispc,pos)*Yk(lpoin_w(ipoin),ispc,pos)
         !$acc loop seq
         do idime = 1,ndime
            f_eta(lpoin_w(ipoin),idime)  = u(lpoin_w(ipoin),idime,1)*eta_Yk(lpoin_w(ipoin),ispc,1)
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

      call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
         gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:),eta_Yk(:,ispc,1),u(:,:,1),Reta(:,2))
      if(mpi_size.ge.2) then
         call mpi_halo_atomic_update_real(auxReta)
      end if
      call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,2))

      !$acc parallel loop
      do ipoin = 1,npoin_w
         auxReta(lpoin_w(ipoin)) = (1.5_rp*Reta(lpoin_w(ipoin),2)-0.5_rp*Reta(lpoin_w(ipoin),1)) + &
                                 (eta_Yk(lpoin_w(ipoin),ispc,2)-eta_Yk(lpoin_w(ipoin),ispc,1))/dt
         Reta(lpoin_w(ipoin),1) = Reta(lpoin_w(ipoin),2)            
      end do
      !$acc end parallel loop

      call nvtxStartRange("Entropy residual")               

      !
      ! Compute entropy viscosity
      !
      call nvtxStartRange("Entropy viscosity evaluation")
      call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
         rho(:,2),u(:,:,2),eta_Yk(:,ispc,2),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
      call nvtxEndRange

   end subroutine rk_4_ls_species_main

   subroutine updateFspecies(ispc,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                              ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
                              rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
                              ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&                ! Optional args
                              listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer,walave_u,walave_t,zo)                       ! Optional args

         implicit none

         integer(4),           intent(in)    :: ispc
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
         real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_t(npoin)
         real(rp), optional, intent(in)      :: zo(npoin)
         integer(4)                          :: pos
         integer(4)                          :: istep, ipoin, idime,icode
         real(rp)                            :: umag, rho_min, rho_avg

         
         pos = 2 ! Set correction as default value
        
         !
         ! Apply bcs after update
         !
         if (noBoundaries .eqv. .false.) then
           call nvtxStartRange("BCS_AFTER_UPDATE")
           call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Yk(:,ispc,pos),Yk_buffer(:,ispc))
           call nvtxEndRange
         end if

               

         ! Compute diffusion terms with values at current substep
         !
         call nvtxStartRange("CONVDIFFS")

         call species_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),Yk(:,ispc,pos),mu_fluid,mu_e_Yk(:,:,ispc),mu_sgs,Ml,Rspc,.true.,-1.0_rp)
         call species_convec_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),Yk(:,ispc,pos),u(:,:,pos),Rspc,.false.,-1.0_rp)               
         call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,Yk(:,ispc,pos),gradYk,.true.)
         call species_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Yk(:,ispc,pos),gradYk,Cp,Prt,rho,tau,Ml,Rspc,.false.,1.0_rp)

         if((isWallModelOn) ) then
            call nvtxStartRange("AB2 wall model")
            if((numBoundsWM .ne. 0)) then
               call evalWallModelABLtemp(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,Cp,mu_fluid,&
                     rho(:,1),walave_u(:,:),walave_t(:),Yk_buffer(:,ispc),zo,Rspc,-1.0_rp)
            end if              
         end if  

         call nvtxEndRange


         !TESTING NEW LOCATION FOR MPICOMMS
         if(mpi_size.ge.2) then
            call nvtxStartRange("MPI_comms_tI")
            call mpi_halo_atomic_update_real(Rspc(:))
            call nvtxEndRange
         end if

         !
         ! Call lumped mass matrix solver
         !
         call nvtxStartRange("Call solver")
         call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rspc(:))
         call nvtxEndRange

   end subroutine updateFspecies
end module time_integ_species
