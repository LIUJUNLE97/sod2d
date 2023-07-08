module time_integ

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

   implicit none

   real(rp), allocatable, dimension(:)   :: Rmass,Rener,Reta
   real(rp), allocatable, dimension(:,:) :: Rmom
   real(rp), allocatable, dimension(:,:)   :: sigMass,sigEner
   real(rp), allocatable, dimension(:,:,:) :: sigMom
   real(rp), allocatable, dimension(:,:)   :: aijKjMass,aijKjEner,pt
   real(rp), allocatable, dimension(:,:,:) :: aijKjMom
   real(rp), allocatable, dimension(:)     :: dt_min,alfa_pt

   real(rp), allocatable, dimension(:)   :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta
   real(rp), allocatable, dimension(:,:) :: aux_u,aux_q
   real(rp), allocatable, dimension(:)   :: Rmass_sum,Rener_sum,Reta_sum,alpha,Rdiff_mass,Rdiff_ener
   real(rp), allocatable, dimension(:,:) :: Rmom_sum,Rdiff_mom,f_eta

   real(rp), allocatable, dimension(:)   :: a_i, b_i, c_i,b_i2
   real(rp), allocatable, dimension(:,:)   :: a_ij

   contains

   subroutine init_rk4_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      if(flag_implicit == 1) then
         if(implicit_solver == implicit_solver_bdf2_rk10) then
            allocate(sigMass(npoin,2), sigEner(npoin,2), sigMom(npoin,ndime,2))
            !$acc enter data create(sigMass(:,:))
            !$acc enter data create(sigEner(:,:))
            !$acc enter data create(sigMom(:,:,:))
            allocate(aijKjMass(npoin,11),aijKjEner(npoin,11),pt(npoin,11))
            !$acc enter data create(aijKjMass(:,:))
            !$acc enter data create(aijKjEner(:,:))
            !$acc enter data create(pt(:,:))
            allocate(aijKjMom(npoin,ndime,11))
            !$acc enter data create(aijKjMom(:,:,:))
            allocate(dt_min(npoin))
            !$acc enter data create(dt_min(:))
            allocate(alfa_pt(5))
            !$acc enter data create(alfa_pt(:))
         endif
      end if

      allocate(Rmass(npoin),Rener(npoin),Reta(npoin),Rmom(npoin,ndime))
      !$acc enter data create(Rmass(:))
      !$acc enter data create(Rener(:))
      !$acc enter data create(Reta(:))
      !$acc enter data create(Rmom(:,:))

      allocate(aux_rho(npoin),aux_pr(npoin),aux_E(npoin),aux_Tem(npoin),aux_e_int(npoin),aux_eta(npoin))
      !$acc enter data create(aux_rho(:))
      !$acc enter data create(aux_pr(:))
      !$acc enter data create(aux_E(:))
      !$acc enter data create(aux_Tem(:))
      !$acc enter data create(aux_e_int(:))
      !$acc enter data create(aux_eta(:))
      allocate(aux_u(npoin,ndime),aux_q(npoin,ndime))
      !$acc enter data create(aux_u(:,:))
      !$acc enter data create(aux_q(:,:))

      allocate(Rmass_sum(npoin),Rener_sum(npoin),Reta_sum(npoin),alpha(npoin),Rdiff_mass(npoin),Rdiff_ener(npoin))
      !$acc enter data create(Rmass_sum(:))
      !$acc enter data create(Rener_sum(:))
      !$acc enter data create(Reta_sum(:))
      !$acc enter data create(alpha(:))
      !$acc enter data create(Rdiff_mass(:))
      !$acc enter data create(Rdiff_ener(:))
      allocate(Rmom_sum(npoin,ndime),Rdiff_mom(npoin,ndime),f_eta(npoin,ndime))
      !$acc enter data create(Rmom_sum(:,:))
      !$acc enter data create(Rdiff_mom(:,:))
      !$acc enter data create(f_eta(:,:))

      if(flag_implicit == 1) then
         allocate(a_ij(11,11),b_i(11),b_i2(11))
         !$acc enter data create(a_ij(:,:))
         !$acc enter data create(b_i(:))
         !$acc enter data create(b_i2(:))
         if (pseudo_steps == 10) then  
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
            write(1,*) "--| ONLY CODED RK10 !"
            stop 1
         end if
         !$acc update device(a_ij(:,:))
         !$acc update device(b_i(:))
         !$acc update device(b_i2(:))
      else
         allocate(a_i(4),b_i(4),c_i(4))
         !$acc enter data create(a_i(:))
         !$acc enter data create(b_i(:))
         !$acc enter data create(c_i(:))

         if (flag_rk_order == 1) then
            a_i = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
            c_i = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
            b_i = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
         else if (flag_rk_order == 2) then
            a_i = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
            c_i = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
            b_i = [0.5_rp, 0.5_rp, 0.0_rp, 0.0_rp]
         else if (flag_rk_order == 3) then
            write(1,*) "--| NOT CODED FOR RK3 YET!"
            stop 1
         else if (flag_rk_order == 4) then
            a_i = [0.0_rp, 0.5_rp, 0.5_rp, 1.0_rp]
            c_i = [0.0_rp, 0.5_rp, 0.5_rp, 1.0_rp]
            b_i = [1.0_rp/6.0_rp, 1.0_rp/3.0_rp, 1.0_rp/3.0_rp, 1.0_rp/6.0_rp]
         else
            write(1,*) "--| NOT CODED FOR RK > 4 YET!"
            stop 1
         end if
         !$acc update device(a_i(:))
         !$acc update device(b_i(:))
         !$acc update device(c_i(:))
      end if

      call nvtxEndRange

   end subroutine init_rk4_solver

   subroutine end_rk4_solver()
      implicit none

      !$acc exit data delete(Rmass(:))
      !$acc exit data delete(Rener(:))
      !$acc exit data delete(Reta(:))
      !$acc exit data delete(Rmom(:,:))
      deallocate(Rmass,Rener,Reta,Rmom)

      !$acc exit data delete(aux_rho(:))
      !$acc exit data delete(aux_pr(:))
      !$acc exit data delete(aux_E(:))
      !$acc exit data delete(aux_Tem(:))
      !$acc exit data delete(aux_e_int(:))
      !$acc exit data delete(aux_eta(:))
      deallocate(aux_rho,aux_pr,aux_E,aux_Tem,aux_e_int,aux_eta)
      !$acc exit data delete(aux_u(:,:))
      !$acc exit data delete(aux_q(:,:))
      deallocate(aux_u,aux_q)

      !$acc exit data delete(Rmass_sum(:))
      !$acc exit data delete(Rener_sum(:))
      !$acc exit data delete(Reta_sum(:))
      !$acc exit data delete(alpha(:))
      !$acc exit data delete(Rdiff_mass(:))
      !$acc exit data delete(Rdiff_ener(:))
      deallocate(Rmass_sum,Rener_sum,Reta_sum,alpha,Rdiff_mass,Rdiff_ener)

      !$acc exit data delete(Rmom_sum(:,:))
      !$acc exit data delete(Rdiff_mom(:,:))
      !$acc exit data delete(f_eta(:,:))
      deallocate(Rmom_sum,Rdiff_mom,f_eta)

      if(flag_implicit == 1) then

         if(implicit_solver == implicit_solver_bdf2_rk10) then
            !$acc exit data delete(sigMass(:,:))
            !$acc exit data delete(sigEner(:,:))
            !$acc exit data delete(sigMom(:,:,:))
            !$acc exit data delete(aijKjMass(:,:))
            !$acc exit data delete(aijKjEner(:,:))
            !$acc exit data delete(pt(:,:))
            !$acc exit data delete(aijKjMom(:,:,:))
            !$acc exit data delete(dt_min(:))
            !$acc exit data delete(alfa_pt(:))
            deallocate(sigMass,sigEner,sigMom)
            deallocate(aijKjMass,aijKjEner,pt)
            deallocate(aijKjMom)
            deallocate(dt_min)
            deallocate(alfa_pt)
         endif

         !$acc exit data delete(a_ij(:,:))
         !$acc exit data delete(b_i(:))
         !$acc exit data delete(b_i2(:))
         deallocate(a_ij,b_i,b_i2)
      else
         !$acc exit data delete(a_i(:))
         !$acc exit data delete(b_i(:))
         !$acc exit data delete(c_i(:))
         deallocate(a_i,b_i,c_i)
      end if

   end subroutine end_rk4_solver

         subroutine rk_implicit_bdf2_rk10_main(igtime,save_logFile_next,currIter,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u)  ! Optional args

            implicit none

            logical,              intent(in)    :: noBoundaries,isWallModelOn
            integer(4),           intent(out)   :: currIter
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
            real(rp),             intent(inout) :: u(npoin,ndime,2)
            real(rp),             intent(inout) :: q(npoin,ndime,3)
            real(rp),             intent(inout) :: pr(npoin,2)
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
            integer(4)                          :: pos,maxIterL
            integer(4)                          :: istep, ipoin, idime,icode,itime,jstep,inode,ielem,npoin_w_g
            real(rp),    dimension(npoin)       :: Rrho
            real(rp)                            :: umag,aux,kappa=1e-6,phi=0.4_rp,xi=0.7_rp,f_save=1.0_rp,f_max=1.01_rp,f_min=0.98_rp
            real(8)                             :: aux2,res_ini,res(2),errMax

            !
            ! Choose between updating prediction or correction
            !
            pos = 2 ! Set correction as default value
            kappa = sqrt(epsilon(kappa))

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
            call nvtxEndRange

            !
            ! Initialize variables to zero
            !
            call nvtxStartRange("Initialize variables")
            !$acc parallel loop
            do ipoin = 1,npoin
               aux_rho(ipoin) = 0.0_rp
               aux_pr(ipoin) = 0.0_rp
               aux_E(ipoin) = 0.0_rp
               aux_Tem(ipoin) = 0.0_rp
               aux_e_int(ipoin) = 0.0_rp
               aux_eta(ipoin) = 0.0_rp
               Rdiff_mass(ipoin) = 0.0_rp
               Rdiff_ener(ipoin) = 0.0_rp
               Rmass(ipoin) = 0.0_rp
               Rener(ipoin) = 0.0_rp
               Reta(ipoin) = 0.0_rp
               Rmass_sum(ipoin) = 0.0_rp
               Rener_sum(ipoin) = 0.0_rp
               aijKjMass(ipoin,1:11) = 0.0_rp
               aijKjEner(ipoin,1:11) = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  aux_u(ipoin,idime) = 0.0_rp
                  aux_q(ipoin,idime) = 0.0_rp
                  Rdiff_mom(ipoin,idime) = 0.0_rp
               
                  Rmom(ipoin,idime) = 0.0_rp
                  Rmom_sum(ipoin,idime) = 0.0_rp
                  aijKjMom(ipoin,idime,1:11) = 0.0_rp
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")
            maxIterL = maxIterNonLineal
            res(:) = 0.0d0
            do itime =1, maxIterL
               !$acc parallel loop
               do ipoin = 1,npoin
                  Rmass_sum(ipoin) = 0.0_rp
                  Rener_sum(ipoin) = 0.0_rp
                  sigMass(ipoin,1) = sigMass(ipoin,2)
                  sigEner(ipoin,1) = sigEner(ipoin,2)
                  sigMass(ipoin,2) = 0.0_rp
                  sigEner(ipoin,2) = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     Rmom_sum(ipoin,idime) = 0.0_rp
                     sigMom(ipoin,idime,1) = sigMom(ipoin,idime,2)
                     sigMom(ipoin,idime,2) = 0.0_rp
                     aijKjMom(ipoin,idime,1:11) = 0.0_rp
                  end do
                  aijKjMass(ipoin,1:11) = 0.0_rp
                  aijKjEner(ipoin,1:11) = 0.0_rp
               end do
               !$acc end parallel loop
               call nvtxStartRange("Loop over pseudo-steps")
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

                  if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,aux_rho,aux_q,u_buffer)

                  !
                  ! Apply bcs after update
                  !
                  if (noBoundaries .eqv. .false.) then
                     call nvtxStartRange("BCS_AFTER_UPDATE")
                     call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,aux_rho(:),aux_q(:,:),aux_u(:,:),aux_pr(:),aux_E(:),u_buffer)
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
                        aux_u(lpoin_w(ipoin),idime) = aux_q(lpoin_w(ipoin),idime)/aux_rho(lpoin_w(ipoin))
                        umag = umag + (aux_u(lpoin_w(ipoin),idime)*aux_u(lpoin_w(ipoin),idime))
                     end do
                     aux_e_int(lpoin_w(ipoin)) = (aux_E(lpoin_w(ipoin))/aux_rho(lpoin_w(ipoin)))-0.5_rp*umag
                     aux_pr(lpoin_w(ipoin)) = aux_rho(lpoin_w(ipoin))*(gamma_gas-1.0_rp)*aux_e_int(lpoin_w(ipoin))
                     aux_Tem(lpoin_w(ipoin)) = aux_pr(lpoin_w(ipoin))/(aux_rho(lpoin_w(ipoin))*Rgas)
                  end do
                  !$acc end parallel loop
                  call nvtxEndRange

                  !
                  ! Compute viscosities and diffusion
                  !
                  !
                  ! Update viscosity if Sutherland's law is active
                  !
                  if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
                     call nvtxStartRange("MU_SUT")
                     call sutherland_viscosity(npoin,aux_Tem,mu_factor,mu_fluid)
                     call nvtxEndRange
                  end if
                  !
                  ! Compute diffusion terms with values at current substep
                  !
                  call nvtxStartRange("DIFFUSIONS")
                  call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_u,aux_Tem,mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass,Rdiff_mom,Rdiff_ener)
                  call nvtxEndRange
                  !
                  ! Call source term if applicable
                  !
                  if(present(source_term)) then
                     call nvtxStartRange("SOURCE TERM")
                     call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,Rdiff_mom)
                     call nvtxEndRange
                  end if
                  !
                  ! Evaluate wall models
      
                  if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                     call nvtxStartRange("WALL MODEL")
                     if(flag_walave == 0) then
                        call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),aux_u(:,:),tauw,Rdiff_mom)
                     else
                        call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                           bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),walave_u(:,:),tauw,Rdiff_mom)
                     end if
                     call nvtxEndRange
                  end if
                  !
                  !
                  ! Compute convective terms
                  !
                  call nvtxStartRange("CONVECTIONS")
                  call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_E,Rmass(:),Rmom(:,:),Rener(:))
                  call nvtxEndRange

                  ! entropy advection
                  !
                  ! Add convection and diffusion terms (Rdiff_* is zero during prediction)
                  !
                  call nvtxStartRange("Add convection and diffusion")
                  !$acc kernels
                  Rmass(:) = Rmass(:) + Rdiff_mass(:)
                  Rener(:) = Rener(:) + Rdiff_ener(:)
                  Rmom(:,:) = Rmom(:,:) + Rdiff_mom(:,:)
                  !$acc end kernels
                  call nvtxEndRange

                  !TESTING NEW LOCATION FOR MPICOMMS
                  if(mpi_size.ge.2) then
                     call nvtxStartRange("MPI_comms_tI")
                     call mpi_halo_atomic_update_real(Rmass(:))
                     call mpi_halo_atomic_update_real(Rener(:))
                     do idime = 1,ndime
                        call mpi_halo_atomic_update_real(Rmom(:,idime))
                     end do
                     call nvtxEndRange
                  end if

                  !
                  ! Call lumped mass matrix solver
                  !
                  call nvtxStartRange("Call solver")
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass(:))
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener(:))
                  call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom(:,:))
                  call nvtxEndRange
                  !
                  ! Accumulate the residuals
                  !
                  call nvtxStartRange("Accumulate residuals")

                  !$acc parallel loop 
                  do ipoin = 1,npoin
                     alfa_pt(1) = 1.0_rp+b_i(istep)*pt(ipoin,1)/(dt*2.0_rp/3.0_rp)
                     alfa_pt(2) = 1.0_rp+b_i(istep)*pt(ipoin,2)/(dt*2.0_rp/3.0_rp)
                     alfa_pt(3) = 1.0_rp+b_i(istep)*pt(ipoin,3)/(dt*2.0_rp/3.0_rp)
                     alfa_pt(4) = 1.0_rp+b_i(istep)*pt(ipoin,4)/(dt*2.0_rp/3.0_rp)
                     alfa_pt(5) = 1.0_rp+b_i(istep)*pt(ipoin,5)/(dt*2.0_rp/3.0_rp)

                     Rmass(ipoin) = -Rmass(ipoin)-(rho(ipoin,2)-4.0_rp*rho(ipoin,1)/3.0_rp + rho(ipoin,3)/3.0_rp)/(dt*2.0_rp/3.0_rp)
                     Rmass_sum(ipoin) = Rmass_sum(ipoin) + b_i(istep)*Rmass(ipoin)/alfa_pt(1)
                     sigMass(ipoin,2) = sigMass(ipoin,2) + abs(pt(ipoin,1)*(b_i(istep)-b_i2(istep))*Rmass(ipoin))/kappa
                     Rener(ipoin) = -Rener(ipoin)-(E(ipoin,2)-4.0_rp*E(ipoin,1)/3.0_rp + E(ipoin,3)/3.0_rp)/(dt*2.0_rp/3.0_rp)
                     Rener_sum(ipoin) = Rener_sum(ipoin) + b_i(istep)*Rener(ipoin)/alfa_pt(2)
                     sigEner(ipoin,2) = sigEner(ipoin,2) + abs(pt(ipoin,2)*(b_i(istep)-b_i2(istep))*Rener(ipoin))/kappa
                     !$acc loop seq
                     do idime = 1,ndime
                        Rmom(ipoin,idime) = -Rmom(ipoin,idime)-(q(ipoin,idime,2)-4.0_rp*q(ipoin,idime,1)/3.0_rp + q(ipoin,idime,3)/3.0_rp)/(dt*2.0_rp/3.0_rp)
                        Rmom_sum(ipoin,idime) = Rmom_sum(ipoin,idime) + b_i(istep)*Rmom(ipoin,idime)/alfa_pt(idime+2)
                        sigMom(ipoin,idime,2) = sigMom(ipoin,idime,2) + abs(pt(ipoin,idime+2)*(b_i(istep)-b_i2(istep))*Rmom(ipoin,idime))/kappa
                     end do
                     if(istep .eq. 1) then
                        !$acc loop seq
                        do jstep=istep+1,pseudo_steps
                           aijKjMass(ipoin,jstep) = aijKjMass(ipoin,jstep) + a_ij(jstep,istep)*Rmass(ipoin)
                           aijKjEner(ipoin,jstep) = aijKjEner(ipoin,jstep) + a_ij(jstep,istep)*Rener(ipoin)
                           !$acc loop seq
                           do idime = 1,ndime
                              aijKjMom(ipoin,idime,jstep) = aijKjMom(ipoin,idime,jstep) + a_ij(jstep,istep)*RMom(ipoin,idime)
                           end do
                        end do
                     else 
                        aijKjMass(ipoin,istep+1) = aijKjMass(ipoin,istep+1) + a_ij(istep+1,istep)*Rmass(ipoin)
                        aijKjEner(ipoin,istep+1) = aijKjEner(ipoin,istep+1) + a_ij(istep+1,istep)*Rener(ipoin)
                        !$acc loop seq
                        do idime = 1,ndime
                           aijKjMom(ipoin,idime,istep+1) = aijKjMom(ipoin,idime,istep+1) + a_ij(istep+1,istep)*RMom(ipoin,idime)
                        end do
                     end if
                  end do
                  !$acc end parallel loop
                  call nvtxEndRange
               end do
               call nvtxEndRange
               !
               ! RK update to variables
               !
               call nvtxStartRange("RK_UPDATE AND PT UPDATE")
               aux2 = 0.d0
               !$acc parallel loop reduction(+:aux2)
               do ipoin = 1,npoin
                  rho(ipoin,pos) = rho(ipoin,pos)+pt(ipoin,1)*Rmass_sum(ipoin)
                  aux2 = aux2 + real(Rmass_sum(ipoin),8)**2
                  E(ipoin,pos) = (E(ipoin,pos)+pt(ipoin,2)*Rener_sum(ipoin))
                  aux2 = aux2 + real(Rener_sum(ipoin),8)**2
                  !$acc loop seq
                  do idime = 1,ndime
                     q(ipoin,idime,pos) = q(ipoin,idime,pos)+pt(ipoin,idime+2)*Rmom_sum(ipoin,idime)
                     aux2 = aux2 + real(Rmom_sum(ipoin,idime),8)**2
                  end do
                  
                 ! pseudo stepping
                  aux = ((sigMass(ipoin,2))**(-phi))*((sigMass(ipoin,1))**(-xi))
                  aux = min(f_max,max(f_min,f_save*aux))
                  pt(ipoin,1) = max(dt_min(ipoin),min(dt_min(ipoin)*pseudo_ftau,aux*pt(ipoin,1)))
                  aux = ((sigEner(ipoin,2))**(-phi))*((sigEner(ipoin,2))**(-xi))
                  aux = min(f_max,max(f_min,f_save*aux))
                  pt(ipoin,2) = max(dt_min(ipoin),min(dt_min(ipoin)*pseudo_ftau,aux*pt(ipoin,2)))
                  !$acc loop seq
                  do idime = 1,ndime
                     aux = ((sigMom(ipoin,idime,2))**(-phi))*((sigMom(ipoin,idime,2))**(-xi))
                     aux = min(f_max,max(f_min,f_save*aux))
                     pt(ipoin,idime+2) = max(dt_min(ipoin),min(dt_min(ipoin)*pseudo_ftau,aux*pt(ipoin,idime+2)))
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               if (flag_buffer_on .eqv. .true.) then
                  call nvtxStartRange("Apply buffer")
                  call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,pos),q(:,:,pos),u_buffer)
                  call nvtxEndRange
               end if

               !
               ! Apply bcs after update
               !
               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer)
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
                  e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))- &
                     0.5_rp*umag!dot_product(u(lpoin_w(ipoin),:,pos),u(lpoin_w(ipoin),:,pos))
                  pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),pos)
                  csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))
                  umag = sqrt(umag)
                  machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
                  Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)
                  aux_eta(lpoin_w(ipoin)) = eta(lpoin_w(ipoin),pos) 
                  eta(lpoin_w(ipoin),pos) = (rho(lpoin_w(ipoin),pos)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)**gamma_gas))
                  !$acc loop seq
                     do idime = 1,ndime
                        f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,pos)*eta(lpoin_w(ipoin),pos)
                     end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               call nvtxStartRange("Update generic convection")
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,pos),u(:,:,pos),Reta(:),alpha)
               call nvtxEndRange

               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(Reta(:))
                  call nvtxEndRange
               end if

               call nvtxStartRange("Lumped mass solver on generic")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:))
               call nvtxEndRange

               call nvtxStartRange("Update sign Reta")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin)) &
                                          -(3.0_rp*eta(lpoin_w(ipoin),2)-4.0_rp*eta(lpoin_w(ipoin),1)+eta(lpoin_w(ipoin),3))/(2.0_rp*dt) &
                                          -(eta(lpoin_w(ipoin),2)-aux_eta(lpoin_w(ipoin)))/pt(lpoin_w(ipoin),1)
               end do
               !$acc end parallel loop
               call nvtxEndRange
               
               if (noBoundaries .eqv. .false.) then
                  call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,Reta)
               end if
               !
               ! Compute entropy viscosity
               !
               call nvtxStartRange("Entropy viscosity evaluation")
               call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Rrho,Ngp,coord,dNgp,gpvol,wgp, &
                  gamma_gas,rho(:,pos),u(:,:,pos),csound,Tem(:,pos),eta(:,pos),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)
               call nvtxEndRange

               call nvtxStartRange("Accumullate aux2 in res(1)")
               call MPI_Allreduce(aux2,res(1),1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange

               res(1) = sqrt(res(1))

               if(itime .lt. 2) then
                  res_ini = res(1)
               endif
               errMax = abs(res(1)-res(2))/abs(res_ini)

               res(2) = res(1)            

               if(errMax .lt. tol) exit
               !if(mpi_rank.eq.0) write(111,*)"   non lineal residual ",errMax," non lineal iterations ",itime
            end do
            call nvtxEndRange

            if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*)"   non lineal residual ",errMax," non lineal iterations ",itime

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
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,pos),u(:,:,pos),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),u(:,:,pos),Ml,mu_sgs)
               end if
               call nvtxEndRange
            end if

            currIter = itime

         end subroutine rk_implicit_bdf2_rk10_main

         subroutine rk_4_main(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term,walave_u)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
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
            real(rp),             intent(inout) :: u(npoin,ndime,2)
            real(rp),             intent(inout) :: q(npoin,ndime,3)
            real(rp),             intent(inout) :: pr(npoin,2)
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
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode
            real(rp),    dimension(npoin)       :: Rrho
            real(rp)                            :: umag


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value

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
            Rdiff_mass(1:npoin) = 0.0_rp
            Rdiff_mom(1:npoin,1:ndime) = 0.0_rp
            Rdiff_ener(1:npoin) = 0.0_rp
            Rmass(1:npoin) = 0.0_rp
            Rmom(1:npoin,1:ndime) = 0.0_rp
            Rener(1:npoin) = 0.0_rp
            Reta(1:npoin) = 0.0_rp
            Rmass_sum(1:npoin) = 0.0_rp
            Rener_sum(1:npoin) = 0.0_rp
            Reta_sum(1:npoin) = 0.0_rp
            Rmom_sum(1:npoin,1:ndime) = 0.0_rp
            !$acc end kernels

            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")

            do istep = 1,flag_rk_order
               !
               ! Compute variable at substep (y_i = y_n+dt*A_ij*R_j)
               !
               call nvtxStartRange("Update aux_*")
               !$acc parallel loop
               do ipoin = 1,npoin
                  aux_rho(ipoin) = rho(ipoin,pos) - dt*a_i(istep)*Rmass(ipoin)
                  aux_E(ipoin)   = E(ipoin,pos)   - dt*a_i(istep)*Rener(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_q(ipoin,idime) = q(ipoin,idime,pos) - dt*a_i(istep)*Rmom(ipoin,idime)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,aux_rho,aux_q,u_buffer)

               if (noBoundaries .eqv. .false.) then
                     call nvtxStartRange("BCS_AFTER_UPDATE")
                     call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,aux_rho(:),aux_q(:,:),aux_u(:,:),aux_pr(:),aux_E(:),u_buffer)
                     call nvtxEndRange
               end if

               !
               ! Apply bcs after update
               !
               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,aux_rho(:),aux_q(:,:),aux_u(:,:),aux_pr(:),aux_E(:),u_buffer)
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
                     aux_u(lpoin_w(ipoin),idime) = aux_q(lpoin_w(ipoin),idime)/aux_rho(lpoin_w(ipoin))
                     umag = umag + (aux_u(lpoin_w(ipoin),idime)*aux_u(lpoin_w(ipoin),idime))
                  end do
                  aux_e_int(lpoin_w(ipoin)) = (aux_E(lpoin_w(ipoin))/aux_rho(lpoin_w(ipoin)))-0.5_rp*umag
                  aux_pr(lpoin_w(ipoin)) = aux_rho(lpoin_w(ipoin))*(gamma_gas-1.0_rp)*aux_e_int(lpoin_w(ipoin))
                  aux_Tem(lpoin_w(ipoin)) = aux_pr(lpoin_w(ipoin))/(aux_rho(lpoin_w(ipoin))*Rgas)
               end do
               !$acc end parallel loop
               call nvtxEndRange

               !
               ! Compute viscosities and diffusion
               !
               !
               ! Update viscosity if Sutherland's law is active
               !
               if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
                  call nvtxStartRange("MU_SUT")
                  call sutherland_viscosity(npoin,aux_Tem,mu_factor,mu_fluid)
                  call nvtxEndRange
               end if
               !
               ! Compute diffusion terms with values at current substep
               !
               call nvtxStartRange("DIFFUSIONS")
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_u,aux_Tem,mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass,Rdiff_mom,Rdiff_ener)
               call nvtxEndRange
               !
               ! Call source term if applicable
               !
               if(present(source_term)) then
                  call nvtxStartRange("SOURCE TERM")
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,Rdiff_mom)
                  call nvtxEndRange
               end if
               !
               ! Evaluate wall models
     
               if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  call nvtxStartRange("WALL MODEL")
                  if(flag_walave == 0) then
                     call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),aux_u(:,:),tauw,Rdiff_mom)
                  else
                     call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),walave_u(:,:),tauw,Rdiff_mom)
                  end if
                  call nvtxEndRange
               end if
               !
               !
               ! Compute convective terms
               !
               call nvtxStartRange("CONVECTIONS")
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_E,Rmass(:),Rmom(:,:),Rener(:))
               call nvtxEndRange

               call nvtxStartRange("Add convection and diffusion")
               !$acc kernels
               Rmass(:) = Rmass(:) + Rdiff_mass(:)
               Rener(:) = Rener(:) + Rdiff_ener(:)
               Rmom(:,:) = Rmom(:,:) + Rdiff_mom(:,:)
               !$acc end kernels
               call nvtxEndRange

               !TESTING NEW LOCATION FOR MPICOMMS
               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(Rmass(:))
                  call mpi_halo_atomic_update_real(Rener(:))
                  do idime = 1,ndime
                     call mpi_halo_atomic_update_real(Rmom(:,idime))
                  end do
                  call nvtxEndRange
               end if

               !
               ! Call lumped mass matrix solver
               !
               call nvtxStartRange("Call solver")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass(:))
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener(:))
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom(:,:))
               call nvtxEndRange
               !
               ! Accumulate the residuals
               !
               call nvtxStartRange("Accumulate residuals")
               !$acc parallel loop
               do ipoin = 1,npoin
                  Rmass_sum(ipoin) = Rmass_sum(ipoin) + b_i(istep)*Rmass(ipoin)
                  Rener_sum(ipoin) = Rener_sum(ipoin) + b_i(istep)*Rener(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     Rmom_sum(ipoin,idime) = Rmom_sum(ipoin,idime) + b_i(istep)*Rmom(ipoin,idime)
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
            !$acc parallel loop
            do ipoin = 1,npoin
               rho(ipoin,pos) = rho(ipoin,pos)-dt*Rmass_sum(ipoin)
               E(ipoin,pos) = E(ipoin,pos)-dt*Rener_sum(ipoin)
               !$acc loop seq
               do idime = 1,ndime
                  q(ipoin,idime,pos) = q(ipoin,idime,pos)-dt*Rmom_sum(ipoin,idime)
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,pos),q(:,:,pos),u_buffer)

            !
            ! Apply bcs after update
            !
            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,rho(:,pos),q(:,:,pos),u(:,:,pos),pr(:,pos),E(:,pos),u_buffer)
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
               e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))- &
                  0.5_rp*umag!dot_product(u(lpoin_w(ipoin),:,pos),u(lpoin_w(ipoin),:,pos))
               pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),pos)
               csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))
               umag = sqrt(umag)
               machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
               Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)
               eta(lpoin_w(ipoin),pos) = (rho(lpoin_w(ipoin),pos)/(gamma_gas-1.0_rp))* &
                  log(pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)**gamma_gas))
               !$acc loop seq
               do idime = 1,ndime
                  f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,pos)*eta(lpoin_w(ipoin),pos)
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            call nvtxStartRange("Entropy convection")
            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,pos),u(:,:,pos),Reta(:),alpha)
            call nvtxEndRange


            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(Reta(:))
               call nvtxEndRange
            end if

            call nvtxStartRange("Entropy Residual")
            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)

            !$acc parallel loop
            do ipoin = 1,npoin_w
               Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin))!-(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
            end do
            !$acc end parallel loop
            call nvtxEndRange
            if (noBoundaries .eqv. .false.) then
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,Reta)
            end if
            !
            ! If using Sutherland viscosity model:
            !
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call nvtxStartRange("Sutherland viscosity")
               call sutherland_viscosity(npoin,Tem(:,pos),mu_factor,mu_fluid)
               call nvtxEndRange
            end if

            call nvtxStartRange("Entropy viscosity evaluation")
            !
            ! Compute entropy viscosity
            !
            call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Rrho,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,pos),u(:,:,pos),csound,Tem(:,pos),eta(:,pos),helem_l,helem,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)
            call nvtxEndRange
            !
            ! Compute subgrid viscosity if active
            !
            if(flag_les == 1) then
               call nvtxStartRange("MU_SGS")
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,pos),u(:,:,pos),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),u(:,:,pos),Ml,mu_sgs)
               end if
               call nvtxEndRange
            end if

         end subroutine rk_4_main

         subroutine updateBuffer(npoin,npoin_w,coord,lpoin_w,rho,q,u_buffer)
            integer(4),           intent(in)    :: npoin
            integer(4),           intent(in)    :: npoin_w
            real(rp),             intent(in)    :: coord(npoin,ndime)
            integer(4),           intent(in)    :: lpoin_w(npoin_w)
            real(rp),             intent(inout) :: rho(npoin)
            real(rp),             intent(inout) :: q(npoin,ndime)
            real(rp),             intent(in) :: u_buffer(npoin,ndime)
            integer(4) :: ipoin
            real(rp)   :: xs,xb,xi,c1,c2

            c1 = 0.01_rp
            c2 = 10.0_rp

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

               q(lpoin_w(ipoin),1) = u_buffer(lpoin_w(ipoin),1)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),1)-u_buffer(lpoin_w(ipoin),1)*rho(lpoin_w(ipoin)))
               q(lpoin_w(ipoin),2) = u_buffer(lpoin_w(ipoin),2)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),2)-u_buffer(lpoin_w(ipoin),2)*rho(lpoin_w(ipoin)))
               q(lpoin_w(ipoin),3) = u_buffer(lpoin_w(ipoin),3)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),3)-u_buffer(lpoin_w(ipoin),3)*rho(lpoin_w(ipoin)))

            end do
            !$acc end parallel loop
         end subroutine updateBuffer

      end module time_integ
