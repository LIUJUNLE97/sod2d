
module time_integ_euler

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
   use time_integ, only : updateBuffer

   implicit none


   real(rp)  , allocatable, dimension(:,:)   :: sigMass,sigEner
   real(rp)  , allocatable, dimension(:,:,:) :: sigMom
   real(rp), allocatable,   dimension(:,:)   :: aijKjMass,aijKjEner,pt
   real(rp), allocatable,   dimension(:,:,:) :: aijKjMom
   real(rp), allocatable, dimension(:)   :: tau_stab_ls
   real(rp), allocatable, dimension(:,:) :: ProjMass_ls,ProjEner_ls,ProjMX_ls,ProjMY_ls,ProjMZ_ls
   real(rp), allocatable, dimension(:)   :: b_i, b_i2
   real(rp), allocatable, dimension(:,:) ::a_ij


   real(rp), allocatable, dimension(:)   :: Reta, Rrho
   real(rp), allocatable, dimension(:,:) :: aux_u, aux_q
   real(rp), allocatable, dimension(:)   :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta,aux_h
   real(rp), allocatable, dimension(:)   :: Rmass, Rener, Rmass_sum, Rener_sum, alpha,dt_min
   real(rp), allocatable, dimension(:,:) :: Rmom, Rmom_sum, f_eta

   contains

   subroutine init_rk_pseudo_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps,i

      call nvtxStartRange("Init RK4 solver")

      allocate(sigMass(npoin,2), sigEner(npoin,2), sigMom(npoin,ndime,2))
      allocate(aijKjMass(npoin,11),aijKjEner(npoin,11),pt(npoin,11))
      allocate(aijKjMom(npoin,ndime,11))

      !$acc enter data create(sigMass(:,:),sigEner(:,:),sigMom(:,:,:))
      !$acc enter data create(aijKjMass(:,:),aijKjEner(:,:),pt(:,:))
      !$acc enter data create(aijKjMom(:,:,:))

      allocate(ProjMass_ls(npoin,ndime),ProjEner_ls(npoin,ndime),ProjMX_ls(npoin,ndime),ProjMY_ls(npoin,ndime),ProjMZ_ls(npoin,ndime),tau_stab_ls(npoin))
      !$acc enter data create(ProjMass_ls(:,:),ProjEner_ls(:,:),ProjMX_ls(:,:),ProjMY_ls(:,:),ProjMZ_ls(:,:),tau_stab_ls(:))

      allocate(Reta(npoin), Rrho(npoin), aux_u(npoin,ndime), aux_q(npoin,ndime), aux_rho(npoin), aux_pr(npoin), aux_E(npoin), aux_Tem(npoin), aux_e_int(npoin),aux_eta(npoin),aux_h(npoin))
      !$acc enter data create(Reta(:), Rrho(:), aux_u(:,:), aux_q(:,:), aux_rho(:), aux_pr(:), aux_E(:), aux_Tem(:), aux_e_int(:),aux_eta(:),aux_h(:))      

      allocate(Rmass(npoin), Rener(npoin), Rmass_sum(npoin), Rener_sum(npoin), alpha(npoin),dt_min(npoin),Rmom(npoin,ndime), Rmom_sum(npoin,ndime), f_eta(npoin,ndime))
      !$acc enter data create(Rmass(:), Rener(:), Rmass_sum(:), Rener_sum(:), alpha(:),dt_min(:),Rmom(:,:), Rmom_sum(:,:), f_eta(:,:))


      !$acc kernels
      sigMass(1:npoin,1:2) = 0.0_rp
      sigEner(1:npoin,1:2) = 0.0_rp
      sigMom(1:npoin,1:ndime,1:2) = 0.0_rp
      ProjMass_ls(:,:) = 0.0_rp
      ProjEner_ls(:,:) = 0.0_rp
      ProjMX_ls(:,:) = 0.0_rp
      ProjMY_ls(:,:) = 0.0_rp
      ProjMZ_ls(:,:) = 0.0_rp
      tau_stab_ls(:) = 0.0_rp      
      !$acc end kernels

      allocate(a_ij(11,11),b_i(11),b_i2(11))
      !$acc enter data create(a_ij(:,:))
      !$acc enter data create(b_i(:))
      !$acc enter data create(b_i2(:))

      call nvtxEndRange

   end subroutine init_rk_pseudo_solver

   subroutine end_rk_pseudo_solver()
      implicit none

   end subroutine end_rk_pseudo_solver

   subroutine rk_pseudo_main(save_logFile_step,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,maskMapped,&
      ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,invMl,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
      rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor,mue_l, &
      ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
      listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,u_mapped,tauw,source_term,walave_u,walave_pr,zo)

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: save_logFile_step 
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
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime), u_mapped(npoin,ndime)
            real(rp), optional,   intent(inout) :: tauw(npoin,ndime)
            real(rp), optional, intent(in)      :: source_term(npoin,ndime+2)
            real(rp), optional, intent(in)      :: walave_u(npoin,ndime),walave_pr(npoin)
            real(rp), optional, intent(in)      :: zo(npoin)
            integer(4)                          :: pos,maxIterL, save_logFile_next
            integer(4)                          :: istep, ipoin, idime,icode,itime,jstep,inode,ielem,npoin_w_g
            real(rp)                            :: alfa_pt(5)
            real(rp)                            :: umag,aux,vol_rank,kappa=1e-6,phi=0.4,xi=0.7,f_save=1.0,f_max=2.0_rp,f_min=0.5_rp,errMax
            real(8)                             :: auxN(5),auxN2(5),vol_tot_d, res(2),aux2,res_ini

            kappa = sqrt(epsilon(kappa))


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            ! Choose between updating prediction or correction
            !
            pos = 2 ! Set correction as default value

            call nvtxStartRange("Updating local dt")
            call adapt_local_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,pseudo_cfl,dt_min,pseudo_cfl,mu_fluid,mu_sgs,rho(:,2))
            call nvtxEndRange()

               call nvtxStartRange("Initialize pt")
               !$acc kernels
               pt(:,1) = dt_min(:)
               pt(:,2) = dt_min(:)
               pt(:,3) = dt_min(:)
               pt(:,4) = dt_min(:)
               pt(:,5) = dt_min(:)
               !$acc end kernels
               call nvtxEndRange()
            !
            ! Butcher tableau
            !
            call nvtxStartRange("Create tableau")
            if (pseudo_steps == 4) then               
               a_ij(:,:) = 0.0_rp

               a_ij(2,1) = 0.32416573882874605_rp  ! real(11847461282814,rp)/real(36547543011857,rp)

               a_ij(3,1) = 0.10407986927510238_rp  ! real(1017324711453,rp)/real(9774461848756,rp)
               a_ij(3,2) = 0.5570978645055429_rp   ! real(3943225443063,rp)/real(7078155732230,rp)

               a_ij(4,1) = 0.10407986927510238_rp  ! real(1017324711453,rp)/real(9774461848756,rp)
               a_ij(4,2) = 0.6019391368822611_rp   ! real(8237718856693,rp)/real(13685301971492,rp)
               a_ij(4,3) = -0.08605491431272755_rp !-real(346793006927,rp)/real(4029903576067,rp)

               b_i(:) = 0.0_rp

               b_i(1) = 0.10407986927510238_rp ! real(1017324711453,rp)/real(9774461848756,rp)
               b_i(2) = 0.6019391368822611_rp  ! real(8237718856693,rp)/real(13685301971492,rp)
               b_i(3) = 2.9750900268840206_rp  ! real(57731312506979,rp)/real(19404895981398,rp)
               b_i(4) = -2.681109033041384_rp  !-real(101169746363290,rp)/real(37734290219643,rp)

               b_i2(:) = 0.0_rp

               b_i2(1) = 0.3406814840808433_rp  ! real(15763415370699,rp)/real(46270243929542,rp)
               b_i2(2) = 0.09091523008632837_rp ! real(514528521746,rp)/real(5659431552419,rp)
               b_i2(3) = 2.866496742725443_rp   ! real(27030193851939,rp)/real(9429696342944,rp)
               b_i2(4) = -2.298093456892615_rp  !-real(69544964788955,rp)/real(30262026368149,rp)
             else if (pseudo_steps == 10) then  
               a_ij(:,:) = 0.0_rp

               a_ij(2,1) = 0.1111111111_rp

               a_ij(3,1) = 0.1900199097_rp
               a_ij(3,2) = 0.0322023124_rp

               a_ij(4,1) = 0.2810938259_rp
               a_ij(4,3) = 0.0522395073_rp

               a_ij(5,1) = 0.3683599872_rp
               a_ij(5,4) = 0.0760844571_rp

               a_ij(6,1) = 0.4503724121_rp
               a_ij(6,5) = 0.1051831433_rp

               a_ij(7,1) = 0.5247721825_rp
               a_ij(7,6) = 0.1418944841_rp
               
               a_ij(8,1) = 0.5874505094_rp
               a_ij(8,7) = 0.1903272683_rp

               a_ij(9,1) = 0.6304783975_rp
               a_ij(9,8) = 0.2584104913_rp

               a_ij(10,1) = 0.6358199324_rp
               a_ij(10,9) = 0.3641800675_rp

               b_i(:) = 0.0_rp

               b_i(1)  = 0.4988192238_rp
               b_i(10)  = 0.5011807761_rp

               b_i2(:) = 0.0_rp

               b_i2(1)  = 0.3906281980_rp
               b_i2(2)  = 0.1179848341_rp
               b_i2(3)  = 1.7353065354_rp
               b_i2(4)  = -7.9567462555_rp
               b_i2(5)  = 17.3753753701_rp
               b_i2(6)  = -23.4057667136_rp
               b_i2(7)  = 20.5007152462_rp
               b_i2(8)  = -11.4042315893_rp
               b_i2(9)  = 3.6467343745_rp

             else
               write(1,*) "--| NOT CODED YET!"
               stop 1
            end if
            call nvtxEndRange

            !$acc update device(a_ij(:,:))
            !$acc update device(b_i(:))
            !$acc update device(b_i2(:))

            !
            ! Initialize variables to zero
            !
            call nvtxStartRange("Initialize variables")
            !$acc kernels
            aux_rho(1:npoin) = 0.0_rp
            aux_u(1:npoin,1:ndime) = 0.0_rp
            aux_q(1:npoin,1:ndime) = 0.0_rp
            aux_pr(1:npoin) = 0.0_rp
            aux_E(1:npoin) = 0.0_rp
            aux_Tem(1:npoin) = 0.0_rp
            aux_e_int(1:npoin) = 0.0_rp
            aux_eta(1:npoin) = 0.0_rp
            Rmass(1:npoin) = 0.0_rp
            Rmom(1:npoin,1:ndime) = 0.0_rp
            Rener(1:npoin) = 0.0_rp
            Reta(1:npoin) = 0.0_rp
            Rmass_sum(1:npoin) = 0.0_rp
            Rener_sum(1:npoin) = 0.0_rp
            Rmom_sum(1:npoin,1:ndime) = 0.0_rp
            aijKjMass(1:npoin,1:11) = 0.0_rp
            aijKjEner(1:npoin,1:11) = 0.0_rp
            aijKjMom(1:npoin,1:ndime,1:11) = 0.0_rp
            !$acc end kernels
            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")
            maxIterL = maxIterNonLineal
            save_logFile_next = 1
            do itime =1, maxIterL
               !$acc kernels
               Rmass_sum(1:npoin) = 0.0_rp
               Rener_sum(1:npoin) = 0.0_rp
               Rmom_sum(1:npoin,1:ndime) = 0.0_rp
               Rmass(1:npoin) = 0.0_rp
               Rmom(1:npoin,1:ndime) = 0.0_rp
               Rener(1:npoin) = 0.0_rp
               aijKjMass(1:npoin,1:11) = 0.0_rp
               aijKjEner(1:npoin,1:11) = 0.0_rp
               aijKjMom(1:npoin,1:ndime,1:11) = 0.0_rp
               sigMass(1:npoin,1) = sigMass(1:npoin,2)
               sigEner(1:npoin,1) = sigEner(1:npoin,2)
               sigMom(1:npoin,1:ndime,1) = sigMom(1:npoin,1:ndime,2)
               sigMass(1:npoin,2) = 0.0_rp
               sigEner(1:npoin,2) = 0.0_rp
               sigMom(1:npoin,1:ndime,2) = 0.0_rp
               !$acc end kernels
               call nvtxStartRange("Loop over pseudo-steps")
               if(flag_lps_stab) call comp_tau(nelem,npoin,connec,csound,u(:,:,pos),helem,dt,tau_stab_ls)
               do istep = 1,pseudo_steps
                  !
                  ! Compute variable at substep (y_i = y_n+dt*A_ij*R_j)
                  !
                  call nvtxStartRange("Update aux_*")
                  !$acc parallel loop
                  do ipoin = 1,npoin
                     aux_rho(ipoin) = rho(ipoin,pos) + pt(ipoin,1)*aijKjMass(ipoin,istep)
                     aux_E(ipoin)   = E(ipoin,pos)   + pt(ipoin,2)*aijKjEner(ipoin,istep)
                     !$acc loop seq
                     do idime = 1,ndime
                        aux_q(ipoin,idime) = q(ipoin,idime,pos) + pt(ipoin,idime+2)*aijKjMom(ipoin,idime,istep)
                     end do
                  end do
                  !$acc end parallel loop
                  call nvtxEndRange

                  if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,maskMapped,aux_rho,aux_q,aux_E,u_buffer)
                  
                  call nvtxStartRange("Update u and EOS")
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     !$acc loop seq
                     do idime = 1,ndime
                        aux_u(lpoin_w(ipoin),idime) = aux_q(lpoin_w(ipoin),idime)/aux_rho(lpoin_w(ipoin))
                     end do
                     aux_e_int(lpoin_w(ipoin)) = (aux_E(lpoin_w(ipoin))/aux_rho(lpoin_w(ipoin)))- &
                        0.5_rp*dot_product(aux_u(lpoin_w(ipoin),:),aux_u(lpoin_w(ipoin),:))
                     aux_pr(lpoin_w(ipoin)) = aux_rho(lpoin_w(ipoin))*(gamma_gas-1.0_rp)*aux_e_int(lpoin_w(ipoin))
                     aux_Tem(lpoin_w(ipoin)) = aux_pr(lpoin_w(ipoin))/(aux_rho(lpoin_w(ipoin))*Rgas)
                     aux_h(lpoin_w(ipoin)) = (gamma_gas/(gamma_gas-1.0_rp))*aux_pr(lpoin_w(ipoin))/aux_rho(lpoin_w(ipoin))
                  end do
                  !$acc end parallel loop
                  call nvtxEndRange
                 
                  call nvtxStartRange("CONVECTIONS")
                  call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_h,Rmass,Rmom,Rener)    
                  call nvtxEndRange
                  if(flag_lps_stab) then           
                     call full_proj_ijk(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_u,&
                                       aux_Tem,Ml,invMl,ProjMass_ls,ProjEner_ls,ProjMX_ls,ProjMY_ls,ProjMZ_ls)
                     call full_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_rho,aux_u,&
                                       aux_Tem,Ml,ProjMass_ls,ProjEner_ls,ProjMX_ls,ProjMY_ls,ProjMZ_ls,tau_stab_ls,Rmass,Rmom,Rener,.false.,-1.0_rp) 
                  end if
                  call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_rho,aux_u,&
                                    aux_Tem,mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener,.false.)                     

                  !TESTING NEW LOCATION FOR MPICOMMS
                  if(mpi_size.ge.2) then
                     call nvtxStartRange("MPI_comms_tI")
                     call mpi_halo_atomic_update_real(Rmass)
                     call mpi_halo_atomic_update_real(Rener)
                     do idime = 1,ndime
                        call mpi_halo_atomic_update_real(Rmom(:,idime))
                     end do
                     call nvtxEndRange
                  end if

                  !
                  ! Call lumped mass matrix solver
                  !
                  call nvtxStartRange("Call solver")
                  call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Rmass)
                  call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Rener)
                  call lumped_solver_vect_opt(npoin,npoin_w,lpoin_w,invMl,Rmom)
                  call nvtxEndRange
                  !
                  ! Accumulate the residuals
                  !
                  call nvtxStartRange("Accumulate residuals")

                  !$acc parallel loop 
                  do ipoin = 1,npoin
                     Rmass(ipoin) = -Rmass(ipoin)
                     Rmass_sum(ipoin) = Rmass_sum(ipoin) + b_i(istep)*Rmass(ipoin)
                     sigMass(ipoin,2) = sigMass(ipoin,2) + abs(pt(ipoin,1)*(b_i(istep)-b_i2(istep))*Rmass(ipoin))/kappa
                     Rener(ipoin) = -Rener(ipoin)
                     Rener_sum(ipoin) = Rener_sum(ipoin) + b_i(istep)*Rener(ipoin)
                     sigEner(ipoin,2) = sigEner(ipoin,2) + abs(pt(ipoin,2)*(b_i(istep)-b_i2(istep))*Rener(ipoin))/kappa
                     !$acc loop seq
                     do idime = 1,ndime
                        Rmom(ipoin,idime) = -Rmom(ipoin,idime)
                        Rmom_sum(ipoin,idime) = Rmom_sum(ipoin,idime) + b_i(istep)*Rmom(ipoin,idime)
                        sigMom(ipoin,idime,2) = sigMom(ipoin,idime,2) + abs(pt(ipoin,idime+2)*(b_i(istep)-b_i2(istep))*Rmom(ipoin,idime))/kappa
                     end do
                     !$acc loop seq
                     do jstep=istep+1,pseudo_steps
                        aijKjMass(ipoin,jstep) = aijKjMass(ipoin,jstep) + a_ij(jstep,istep)*Rmass(ipoin)
                        aijKjEner(ipoin,jstep) = aijKjEner(ipoin,jstep) + a_ij(jstep,istep)*Rener(ipoin)
                        !$acc loop seq
                        do idime = 1,ndime
                           aijKjMom(ipoin,idime,jstep) = aijKjMom(ipoin,idime,jstep) + a_ij(jstep,istep)*RMom(ipoin,idime)
                        end do
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
               aux2 = 0.0_rp
               !$acc parallel loop reduction(+:aux2)
               do ipoin = 1,npoin
                  rho(ipoin,pos) = rho(ipoin,pos)+pt(ipoin,1)*Rmass_sum(ipoin)
                  aux2 = aux2 + real(Rmass_sum(ipoin)**2,8)
                  E(ipoin,pos) = (E(ipoin,pos)+pt(ipoin,2)*Rener_sum(ipoin))
                  aux2 = aux2 + real(Rener_sum(ipoin)**2,8)
                  !$acc loop seq
                  do idime = 1,ndime
                     q(ipoin,idime,pos) = q(ipoin,idime,pos)+pt(ipoin,idime+2)*Rmom_sum(ipoin,idime)
                     aux2 = aux2 + real(Rmom_sum(ipoin,idime)**2,8)
                  end do
                  
                 ! pseudo stepping
                  aux = ((sigMass(ipoin,2))**(-phi/3.0_rp))*((sigMass(ipoin,1))**(-xi/3.0_rp))
                  aux = min(f_max,max(f_min,f_save*aux))
                  pt(ipoin,1) = max(dt_min(ipoin),min(dt_min(ipoin)*pseudo_ftau,aux*pt(ipoin,1)))
                  aux = ((sigEner(ipoin,2))**(-phi/3.0_rp))*((sigEner(ipoin,2))**(-xi/3.0_rp))
                  aux = min(f_max,max(f_min,f_save*aux))
                  pt(ipoin,2) = max(dt_min(ipoin),min(dt_min(ipoin)*pseudo_ftau,aux*pt(ipoin,2)))
                  !$acc loop seq
                  do idime = 1,ndime
                     aux = ((sigMom(ipoin,idime,2))**(-phi/3.0_rp))*((sigMom(ipoin,idime,2))**(-xi/3.0_rp))
                     aux = min(f_max,max(f_min,f_save*aux))
                     pt(ipoin,idime+2) = max(dt_min(ipoin),min(dt_min(ipoin)*pseudo_ftau,aux*pt(ipoin,idime+2)))
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               if (flag_buffer_on .eqv. .true.) then
                  call nvtxStartRange("Apply buffer")
                  call updateBuffer(npoin,npoin_w,coord,lpoin_w,maskMapped,rho(:,pos),q(:,:,pos),E(:,pos),u_buffer)
                  call nvtxEndRange
               end if

               !
               ! Apply bcs after update
               !
               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer,u_mapped)
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
                  umag = sqrt(umag)
                  e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))- &
                     0.5_rp*dot_product(u(lpoin_w(ipoin),:,pos),u(lpoin_w(ipoin),:,pos))
                  pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),pos)
                  csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))
                  machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
                  Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)
                  eta(lpoin_w(ipoin),pos) = (rho(lpoin_w(ipoin),pos)/(gamma_gas-1.0_rp))* &
                     log(pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)**gamma_gas))
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),pos)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               call nvtxStartRange("Update generic convection")
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta)
               call nvtxEndRange


               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(Reta)
                  call nvtxEndRange
               end if

               call nvtxStartRange("Lumped mass solver on generic")
               call lumped_solver_scal_opt(npoin,npoin_w,lpoin_w,invMl,Reta)
               call nvtxEndRange

               call nvtxStartRange("Update sign Reta")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin))-(3.0_rp*eta(lpoin_w(ipoin),2)-4.0_rp*eta(lpoin_w(ipoin),1)+eta(lpoin_w(ipoin),3))/(2.0_rp*dt)
               end do
               !$acc end parallel loop
               call nvtxEndRange

               !
               ! Compute entropy viscosity
               !
               call nvtxStartRange("Entropy viscosity evaluation")
               call smart_visc_spectral_imex(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Ngp,coord,dNgp,gpvol,wgp, &
                  gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
               call nvtxEndRange

               call nvtxStartRange("Accumullatee auxN in auxN2")
               call MPI_Allreduce(aux2,res(1),1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange

               res(1) = sqrt(res(1))

               if(itime .lt. 3) then
                  res_ini = res(1)
               endif
               !errMax = abs(res(1))/abs(res_ini)    
               errMax = abs(res(1)-res(2))/abs(res_ini)
               res(2) = res(1)
           

               if(errMax .lt. tol) exit
               if(itime==save_logFile_next.and.mpi_rank.eq.0) then
                  write(111,*) " it ",itime," err ",errMax
                  call flush(111)
                  save_logFile_next = save_logFile_next + save_logFile_step
               end if
            end do
            call nvtxEndRange

            call nvtxStartRange("Last update")
            !$acc kernels
            rho(:,3) = rho(:,1)
            E(:,3) = E(:,1)
            q(:,:,3) = q(:,:,1)
            eta(:,3) = eta(:,1)
            !$acc end kernels
            call nvtxEndRange

         end subroutine rk_pseudo_main


end module time_integ_euler
