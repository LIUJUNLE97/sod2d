module time_integ

   use mod_nvtx
   use elem_convec
   use elem_diffu
   use elem_source
   use mod_solver
   use mod_entropy_viscosity
   use mod_constants
   use mod_fluid_viscosity
   use mod_sgs_viscosity
   use mod_sgs_ilsa_viscosity
   use mod_bc_routines
   use mod_wall_model

      implicit none

      real(rp)  , allocatable, dimension(:,:)   :: sigMass,sigEner
      real(rp)  , allocatable, dimension(:,:,:) :: sigMom
      real(rp), allocatable,   dimension(:,:)     :: aijKjMass,aijKjEner,pt
      real(rp), allocatable,   dimension(:,:,:) :: aijKjMom
      real(rp)                                  :: volT
      logical                                   :: flag_mem_alloc=.true.
      logical                                   :: initialize_pt=.true.

      contains

         subroutine implicit_rosenbrock_main(noBoundaries,isWallModelOn,flag_predic,flag_emac,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: flag_predic, flag_emac
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn(nboun,npbou),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem) !helem_l(npoin) TO REVIEW I THINK IS BUG!
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,2)
            real(rp),             intent(inout) :: u(npoin,ndime,2)
            real(rp),             intent(inout) :: q(npoin,ndime,2)
            real(rp),             intent(inout) :: pr(npoin,2)
            real(rp),             intent(inout) :: E(npoin,2)
            real(rp),             intent(inout) :: Tem(npoin,2)
            real(rp),             intent(inout) :: e_int(npoin,2)
            real(rp),             intent(inout) :: eta(npoin,2)
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
            integer(4)                          :: nstep
            integer(4)                          :: istep, ipoin,idime,icode,jstep
            real(rp),    dimension(npoin)       :: Reta, Reta2, Rrho
            real(rp),    dimension(8)           :: m_i
            real(rp),    dimension(8,8)         :: a_ij, c_ij
            real(rp),    dimension(npoin,ndime) :: aux_u, aux_q,aux_u_wall,u_eta
            real(rp),    dimension(npoin)       :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta
            real(rp),    dimension(npoin)       :: Rmass, Rener,alpha,cMass,cEner
            real(rp),    dimension(npoin,ndime) :: Rmom, f_eta,cMom
            real(rp),    dimension(npoin,8)       :: Yrho,YE
            real(rp),    dimension(npoin,ndime,8) :: Yq
            real(rp)                              :: Rdiff_mass(npoin), Rdiff_mom(npoin,ndime), Rdiff_ener(npoin)
            real(rp)                              :: Aemac(npoin,ndime), Femac(npoin), umag,gamma_RK

            !
            ! Butcher tableau
            !
            call nvtxStartRange("Create tableau")

            if (flag_rk_order == 2) then
               nstep = 2
               gamma_RK = 1.0_rp - sqrt(2.0_rp)*0.5_rp

               a_ij(:,:) = 0.0_rp

               a_ij(2,1) = 4.0_rp -8.0_rp*gamma_RK

               c_ij(:,:) = 0.0_rp

               m_i(:) = 0.0_rp

               m_i(1) = (1.0_rp/gamma_RK)*(1.0_rp-1.0_rp/(8.0_rp*gamma_RK))
               m_i(2) = 1.0_rp/(8.0_rp*gamma_RK*gamma_RK)
            else if (flag_rk_order == 3) then
               nstep = 3
               gamma_RK = real(7.886751345948129e-01,rp)

               a_ij(:,:) = 0.0_rp

               a_ij(2,1) = real(1.267949192431123e+00,rp)
               a_ij(3,1) = real(1.267949192431123e+00,rp)

               c_ij(:,:) = 0.0_rp

               c_ij(2,1) = real(-1.607695154586736e+00,rp) 
               c_ij(3,1) = real(-3.464101615137755e+00,rp) 
               c_ij(3,2) = real(-1.732050807568877e+00,rp) 

               m_i(:) = 0.0_rp

               m_i(1) = real(2.000000000000000e+00,rp) 
               m_i(2) = real(5.773502691896258e-01,rp)
               m_i(3) = real(4.226497308103742e-01,rp)
            else if (flag_rk_order == 4) then
               nstep = 4
               gamma_RK = 0.5_rp

               a_ij(:,:) = 0.0_rp

               a_ij(2,1) = 2.0_rp
               a_ij(3,1) = 48.0_rp/25.0_rp
               a_ij(3,2) = 6.0_rp/25.0_rp
               a_ij(4,1) = 48.0_rp/25.0_rp
               a_ij(4,2) = 6.0_rp/25.0_rp
               a_ij(4,3) = 0.0_rp

               c_ij(:,:) = 0.0_rp

               c_ij(2,1) = -8.0_rp
               c_ij(3,1) = 372.0_rp/25.0_rp
               c_ij(3,2) = 12.0_rp/5.0_rp
               c_ij(4,1) = -112.0_rp/125.0_rp
               c_ij(4,2) = -54.0_rp/125.0_rp
               c_ij(4,3) = -2.0_rp/5.0_rp
               
               m_i(:) = 0.0_rp

               m_i(1) = 19.0_rp/9.0_rp
               m_i(2) = 0.5_rp
               m_i(3) = 25.0_rp/108.0_rp
               m_i(4) = 125.0_rp/108.0_rp
            else if (flag_rk_order == 5) then
               nstep = 8

               gamma_RK = 1.9e-1

               a_ij(:,:) = 0.0_rp

               a_ij(2,1) =real(2.0,rp)
               a_ij(3,1) =real(3.040894194418781,rp)
               a_ij(3,2) =real(1.041747909077569,rp)
               a_ij(4,1) =real(2.576417536461461,rp)
               a_ij(4,2) =real(1.622083060776640,rp)
               a_ij(4,3) =real(-9.089668560264532e-1,rp)
               a_ij(5,1) =real(2.760842080225597,rp)
               a_ij(5,2) =real(1.446624659844071,rp)
               a_ij(5,3) =real(-3.036980084553738e-1,rp)
               a_ij(5,4) =real(2.877498600325443e-1,rp)
               a_ij(6,1) =real(-1.409640773051259e1,rp)
               a_ij(6,2) =real(6.925207756232704,rp)
               a_ij(6,3) =real(-4.147510893210728e1,rp)
               a_ij(6,4) =real(2.343771018586405,rp)
               a_ij(6,5) =real(2.413215229196062e1,rp)
               a_ij(7,1) =real(-1.409640773051259e1,rp)
               a_ij(7,2) =real(6.925207756232704,rp)
               a_ij(7,3) =real(-4.147510893210728e1,rp)
               a_ij(7,4) =real(2.343771018586405,rp)
               a_ij(7,5) =real(2.413215229196062e1,rp)
               a_ij(7,6) =real(1.0,rp)
               a_ij(8,1) =real(-1.409640773051259e1,rp)
               a_ij(8,2) =real(6.925207756232704,rp)
               a_ij(8,3) =real(-4.147510893210728e1,rp)
               a_ij(8,4) =real(2.343771018586405,rp)
               a_ij(8,5) =real(2.413215229196062e1,rp)
               a_ij(8,6) =real(1.0,rp)
               a_ij(8,7) =real(1.0,rp)

               c_ij(:,:) = 0.0_rp

               c_ij(2,1) =real(-1.031323885133993e1,rp)
               c_ij(3,1) =real(-2.104823117650003e1,rp)
               c_ij(3,2) =real(-7.234992135176716,rp)
               c_ij(4,1) =real(3.222751541853323e1,rp)
               c_ij(4,2) =real(-4.943732386540191,rp)
               c_ij(4,3) =real(1.944922031041879e1,rp)
               c_ij(5,1) =real(-2.069865579590063e1,rp)
               c_ij(5,2) =real(-8.816374604402768,rp)
               c_ij(5,3) =real(1.260436877740897,rp)
               c_ij(5,4) =real(-7.495647613787146e-1,rp)
               c_ij(6,1) =real(-4.622004352711257e1,rp)
               c_ij(6,2) =real(-1.749534862857472e1,rp)
               c_ij(6,3) =real(-2.896389582892057e2,rp)
               c_ij(6,4) =real(9.360855400400906e1,rp)
               c_ij(6,5) =real(3.183822534212147e2,rp)
               c_ij(7,1) =real(3.420013733472935e1,rp)
               c_ij(7,2) =real(-1.415535402717690e1,rp)
               c_ij(7,3) =real(5.782335640988400e1,rp)
               c_ij(7,4) =real(2.583362985412365e1,rp)
               c_ij(7,5) =real(1.408950972071624,rp)
               c_ij(7,6) =real(-6.551835421242162,rp)
               c_ij(8,1) =real(4.257076742291101e1,rp)
               c_ij(8,2) =real(-1.380770672017997e1,rp)
               c_ij(8,3) =real(9.398938432427124e1,rp)
               c_ij(8,4) =real(1.877919633714503e1,rp)
               c_ij(8,5) =real(-3.158359187223370,rp)
               c_ij(8,6) =real(-6.685968952921985,rp)
               c_ij(8,7) =real(-5.810979938412932,rp)
               
               m_i(:) = 0.0_rp

               m_i(1) =real(-1.409640773051259e1,rp)
               m_i(2) =real(6.925207756232704,rp)
               m_i(3) =real(-4.147510893210728e1,rp)
               m_i(4) =real(2.343771018586405,rp)
               m_i(5) =real(2.413215229196062e1,rp)
               m_i(6) =real(1.0,rp)
               m_i(7) =real(1.0,rp)
               m_i(8) =real(1.0,rp)

            else
               write(1,*) "--| NOT CODED FOR RK > 5 YET!"
               stop 1
            end if
            call nvtxEndRange
            !
            ! Initialize variables to zero
            !
            call nvtxStartRange("Initialize variables")
            !$acc kernels
            aux_rho(1:npoin) = 0.0_rp
            aux_u(1:npoin,1:ndime) = u(:,:,1)
            u_eta(1:npoin,1:ndime) = u(:,:,1)
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
            Yq(1:npoin,1:ndime,1:8) = 0.0_rp
            YE(1:npoin,1:8) = 0.0_rp
            Yrho(1:npoin,1:8) = 0.0_rp
            cMass(1:npoin) = 0.0_rp
            cMom(1:npoin,1:ndime) = 0.0_rp
            cEner(1:npoin) = 0.0_rp
            !$acc end kernels

            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")

            do istep = 1,nstep
               !
               ! Compute variable at substep (y_i = y_n+dt*A_ij*R_j)
               !
               call nvtxStartRange("Update aux_*")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  aux_rho(lpoin_w(ipoin)) = rho(lpoin_w(ipoin),1)
                  aux_E(lpoin_w(ipoin))   = E(lpoin_w(ipoin),1)
                  cMass(lpoin_w(ipoin)) = 0.0_rp
                  cEner(lpoin_w(ipoin))   = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_q(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,1)
                     cMom(lpoin_w(ipoin),idime) = 0.0_rp
                  end do
                  !$acc loop seq
                  do jstep=1, istep-1
                     aux_rho(lpoin_w(ipoin)) = aux_rho(lpoin_w(ipoin)) + a_ij(istep,jstep)*Yrho(lpoin_w(ipoin),jstep)
                     aux_E(lpoin_w(ipoin))   = aux_E(lpoin_w(ipoin))   + a_ij(istep,jstep)*YE(lpoin_w(ipoin),jstep)
                     cMass(lpoin_w(ipoin)) = cMass(lpoin_w(ipoin)) + c_ij(istep,jstep)*Yrho(lpoin_w(ipoin),jstep)/dt
                     cEner(lpoin_w(ipoin))   = cEner(lpoin_w(ipoin))   + c_ij(istep,jstep)*YE(lpoin_w(ipoin),jstep)/dt
                     !$acc loop seq
                     do idime = 1,ndime
                        aux_q(lpoin_w(ipoin),idime) = aux_q(lpoin_w(ipoin),idime) + a_ij(istep,jstep)*Yq(lpoin_w(ipoin),idime,jstep)
                        cMom(lpoin_w(ipoin),idime) = cMom(lpoin_w(ipoin),idime) + c_ij(istep,jstep)*Yq(lpoin_w(ipoin),idime,jstep)/dt
                     end do
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,aux_rho,aux_q,u_buffer)

               !
               ! Update velocity and equations of state
               !
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
               end do
               !$acc end parallel loop

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
               
               ! Compute diffusion terms with values at current substep
               !
               call nvtxStartRange("DIFFUSIONS")
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,aux_rho,aux_u,aux_Tem,mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass,Rdiff_mom,Rdiff_ener)
               call nvtxEndRange
               !
               ! Call source term if applicable
               !
               if(present(source_term)) then
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,Rdiff_mom)
               end if
               !
               ! Evaluate wall models
     
               if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),aux_u(:,:),tauw,Rdiff_mom)
               end if
               !
               !
               ! Compute convective terms
               !
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_E,Rmass,Rmom,Rener)

               ! entropy advection
               !
               ! Add convection and diffusion terms (Rdiff_* is zero during prediction)
               !
               call nvtxStartRange("Add convection and diffusion")
               !$acc kernels
               Rmass(:) =  Rmass(:) + Rdiff_mass(:) 
               Rener(:) =  Rener(:) + Rdiff_ener(:) 
               Rmom(:,:) =  Rmom(:,:) + Rdiff_mom(:,:) 
               !$acc end kernels
               call nvtxEndRange

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
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass)
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener)
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom)
               call nvtxEndRange

               !$acc kernels
               Rmass(:) =  -Rmass(:)   + cMass(:) 
               Rener(:) =  -Rener(:)   + cEner(:)
               Rmom(:,:) =  -Rmom(:,:) + cMom(:,:)
               !$acc end kernels
               call nvtxEndRange

#if 0
             call jacobi_full(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                    atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                    !aux_rho(:),aux_u(:,:),aux_q(:,:),aux_pr(:),aux_E(:),aux_Tem(:),Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                    rho(:,1),u(:,:,1),q(:,:,1),pr(:,1),E(:,1),Tem(:,1),Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                    gamma_RK,dt,Rmass,Rmom,Rener,Yrho(:,istep),Yq(:,:,istep),YE(:,istep),istep)
#else
              call gmres_full(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                    atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                    aux_rho(:),aux_u(:,:),aux_q(:,:),aux_pr(:),aux_E(:),aux_Tem(:),Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                    !rho(:,1),u(:,:,1),q(:,:,1),pr(:,1),E(:,1),Tem(:,1),Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                    gamma_RK,dt,Rmass,Rmom,Rener,Yrho(:,istep),Yq(:,:,istep),YE(:,istep),istep)
#endif
               !
               ! RK update to variables
               !
               call nvtxStartRange("RK_UPDATE")
            end do
            call nvtxEndRange
            !$acc parallel loop
            do ipoin = 1,npoin_w
               rho(lpoin_w(ipoin),2) = rho(lpoin_w(ipoin),1)
               E(lpoin_w(ipoin),2) = E(lpoin_w(ipoin),1)
               !$acc loop seq
               do idime = 1,ndime
                  q(lpoin_w(ipoin),idime,2) = q(lpoin_w(ipoin),idime,1)
               end do
               do istep=1,nstep
                  rho(lpoin_w(ipoin),2) = rho(lpoin_w(ipoin),2)+m_i(istep)*Yrho(lpoin_w(ipoin),istep)
                  E(lpoin_w(ipoin),2) = E(lpoin_w(ipoin),2)+m_i(istep)*YE(lpoin_w(ipoin),istep)
                  !$acc loop seq
                  do idime = 1,ndime
                     q(lpoin_w(ipoin),idime,2) = q(lpoin_w(ipoin),idime,2)+m_i(istep)*Yq(lpoin_w(ipoin),idime,istep)
                  end do
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (flag_buffer_on .eqv. .true.) call updateBuffer(npoin,npoin_w,coord,lpoin_w,rho(:,2),q(:,:,2),u_buffer)

            !
            ! Apply bcs after update
            !
            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,rho(:,2),q(:,:,2),u(:,:,2),pr(:,2),E(:,2),u_buffer)
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
                  u(lpoin_w(ipoin),idime,2) = q(lpoin_w(ipoin),idime,2)/rho(lpoin_w(ipoin),2)
                  umag = umag + u(lpoin_w(ipoin),idime,2)**2
               end do
               umag = sqrt(umag)
               e_int(lpoin_w(ipoin),2) = (E(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2))- &
                  0.5_rp*dot_product(u(lpoin_w(ipoin),:,2),u(lpoin_w(ipoin),:,2))
               pr(lpoin_w(ipoin),2) = rho(lpoin_w(ipoin),2)*(gamma_gas-1.0_rp)*e_int(lpoin_w(ipoin),2)
               csound(lpoin_w(ipoin)) = sqrt(gamma_gas*pr(lpoin_w(ipoin),2)/rho(lpoin_w(ipoin),2))
               machno(lpoin_w(ipoin)) = umag/csound(lpoin_w(ipoin))
               Tem(lpoin_w(ipoin),2) = pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)*Rgas)
               eta(lpoin_w(ipoin),2) = (rho(lpoin_w(ipoin),2)/(gamma_gas-1.0_rp))* &
                  log(abs(pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)**gamma_gas)))
               !$acc loop seq
               do idime = 1,ndime
                  f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
               end do
            end do
            !$acc end parallel loop

            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta,alpha)


            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(Reta)
               call nvtxEndRange
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)

            !$acc parallel loop
            do ipoin = 1,npoin_w
               Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin))!-(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
            end do
            !$acc end parallel loop

            call nvtxStartRange("Entropy viscosity evaluation")
            !
            ! Compute entropy viscosity
            !
            call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Rrho,Ngp,coord,dNgp,gpvol,wgp, &
               gamma_gas,rho(:,2),u(:,:,2),csound,Tem(:,2),eta(:,2),helem_l,helem,Ml,mu_e)
            call nvtxEndRange

            !
            ! If using Sutherland viscosity model:
            !
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call nvtxStartRange("Sutherland viscosity")
               call sutherland_viscosity(npoin,Tem(:,2),mu_factor,mu_fluid)
               call nvtxEndRange
            end if

            if(flag_les == 1) then
               call nvtxStartRange("MU_SGS")
               if(flag_les_ilsa == 1) then
                  call sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt,rho(:,2),u(:,:,2),mu_sgs,mu_fluid,mu_e,kres,etot,au,ax1,ax2,ax3) 
               else
                  call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,2),u(:,:,2),Ml,mu_sgs)
               end if
               call nvtxEndRange
            end if

         end subroutine implicit_rosenbrock_main

         subroutine rk_pseudo_main(noBoundaries,isWallModelOn,flag_predic,flag_emac,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: flag_predic, flag_emac
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn(nboun,npbou),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem) !helem_l(npoin) TO REVIEW I THINK IS BUG!
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,2)
            real(rp),             intent(inout) :: u(npoin,ndime,2)
            real(rp),             intent(inout) :: q(npoin,ndime,2)
            real(rp),             intent(inout) :: pr(npoin,2)
            real(rp),             intent(inout) :: E(npoin,2)
            real(rp),             intent(inout) :: Tem(npoin,2)
            real(rp),             intent(inout) :: e_int(npoin,2)
            real(rp),             intent(inout) :: eta(npoin,2)
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
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode,itime,jstep,inode,ielem,npoin_w_g
            real(rp),    dimension(npoin)       :: Reta, Rrho
            real(rp),    dimension(11)           :: b_i, b_i2
            real(rp),    dimension(11,11)         :: a_ij
            real(rp),    dimension(npoin,ndime) :: aux_u, aux_q,aux_u_wall
            real(rp),    dimension(npoin)       :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta
            real(rp),    dimension(npoin)       :: Rmass, Rener, Rmass_sum, Rener_sum, alpha,dt_min
            real(rp),    dimension(npoin,ndime) :: Rmom, Rmom_sum, f_eta
            real(rp)                            :: Rdiff_mass(npoin), Rdiff_mom(npoin,ndime), Rdiff_ener(npoin)
            real(rp)                            :: umag,aux,vol_rank,kappa=1e-6,phi=0.4,xi=0.7,f_save=1.0,f_max=1.02_rp,f_min=0.98_rp,errMax
            real(8)                             :: auxN(5),auxN2(5),vol_tot_d

            !kappa = sqrt(epsilon(kappa))

            if (flag_mem_alloc .eqv. .true.) then
               allocate(sigMass(npoin,2), sigEner(npoin,2), sigMom(npoin,ndime,2))
               allocate(aijKjMass(npoin,11),aijKjEner(npoin,11),pt(npoin,11))
               allocate(aijKjMom(npoin,ndime,11))
               flag_mem_alloc = .false.
               !$acc kernels
               sigMass(1:npoin,1:2) = 0.0_rp
               sigEner(1:npoin,1:2) = 0.0_rp
               sigMom(1:npoin,1:ndime,1:2) = 0.0_rp
               !$acc end kernels

               vol_rank  = 0.0
               vol_tot_d = 0.0
               call nvtxStartRange("Computingg ttotal volume per rank")
               !$acc parallel loop collapse(2) reduction(+:vol_rank)
               do ielem = 1,numElemsRankPar
                  do inode = 1,nnode
                     vol_rank = vol_rank+gpvol(1,inode,ielem)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange()

               call nvtxStartRange("Computing total volume")
               call MPI_Allreduce(vol_rank,vol_tot_d,1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange()
               call nvtxStartRange("Computing total number of points")
               call MPI_Allreduce(npoin_w,npoin_w_g,1,mpi_datatype_int,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange()

               volT = real(vol_tot_d,rp)/npoin_w_g

            end if


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

            !if(initialize_pt .eqv. .true.) then
               call nvtxStartRange("Initialize pt")
               !$acc kernels
               pt(:,1) = dt_min(:)
               pt(:,2) = dt_min(:)
               pt(:,3) = dt_min(:)
               pt(:,4) = dt_min(:)
               pt(:,5) = dt_min(:)
               !$acc end kernels
               call nvtxEndRange()
               !initialize_pt = .true.
            !end if
            !
            ! Butcher tableau
            !
            call nvtxStartRange("Create tableau")
            if (flag_pseudo_steps == 4) then               
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
               b_i(4) = -2.681109033041384_rp  ! real(101169746363290,rp)/real(37734290219643,rp)

               b_i2(:) = 0.0_rp

               b_i2(1) = 0.3406814840808433_rp  ! real(15763415370699,rp)/real(46270243929542,rp)
               b_i2(2) = 0.09091523008632837_rp ! real(514528521746,rp)/real(5659431552419,rp)
               b_i2(3) = 2.866496742725443_rp   ! real(27030193851939,rp)/real(9429696342944,rp)
               b_i2(4) = -2.298093456892615_rp  ! real(69544964788955,rp)/real(30262026368149,rp)
             else if (flag_pseudo_steps == 10) then  
               a_ij(:,:) = 0.0_rp

               a_ij(2,1) = 0.1111111111

               a_ij(3,1) = 0.1900199097
               a_ij(3,2) = 0.0322023124

               a_ij(4,1) = 0.2810938259
               a_ij(4,3) = 0.0522395073

               a_ij(5,1) = 0.3683599872
               a_ij(5,4) = 0.0760844571

               a_ij(6,1) = 0.4503724121
               a_ij(6,5) = 0.1051831433

               a_ij(7,1) = 0.5247721825
               a_ij(7,6) = 0.1418944841
               
               a_ij(8,1) = 0.5874505094
               a_ij(8,7) = 0.1903272683

               a_ij(9,1) = 0.6304783975
               a_ij(9,8) = 0.2584104913

               a_ij(10,1) = 0.6358199324
               a_ij(10,9) = 0.3641800675

               b_i(:) = 0.0_rp

               b_i(1)  = 0.4988192238
               b_i(10)  = 0.5011807761

               b_i2(:) = 0.0_rp

               b_i2(1)  = 0.3906281980
               b_i2(2)  = 0.1179848341
               b_i2(3)  = 1.7353065354
               b_i2(4)  = -7.9567462555
               b_i2(5)  = 17.3753753701
               b_i2(6)  = -23.4057667136
               b_i2(7)  = 20.5007152462
               b_i2(8)  = -11.4042315893
               b_i2(9)  = 3.6467343745

             else
               write(1,*) "--| NOT CODED YET!"
               stop 1
            end if
            call nvtxEndRange
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
            Rmom_sum(1:npoin,1:ndime) = 0.0_rp
            aijKjMass(1:npoin,1:5) = 0.0_rp
            aijKjEner(1:npoin,1:5) = 0.0_rp
            aijKjMom(1:npoin,1:ndime,1:5) = 0.0_rp
            !$acc end kernels
            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")
            do itime =1, maxIter
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
               do istep = 1,flag_pseudo_steps
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
                  ! Update velocity and equations of state
                  !
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
                     call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                        bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),aux_u(:,:),tauw,Rdiff_mom)
                     call nvtxEndRange
                  end if
                  !
                  !
                  ! Compute convective terms
                  !
                  call nvtxStartRange("CONVECTIONS")
                  call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_E,Rmass,Rmom,Rener)
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
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass)
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener)
                  call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom)
                  call nvtxEndRange
                  !
                  ! Accumulate the residuals
                  !
                  call nvtxStartRange("Accumulate residuals")
                  !$acc parallel loop
                  do ipoin = 1,npoin
                     Rmass(ipoin) = -Rmass(ipoin)-(rho(ipoin,2)-4.0_rp*rho(ipoin,1)/3.0_rp + rho(ipoin,3)/3.0_rp)/(dt*2.0_rp/3.0_rp)
                     Rmass_sum(ipoin) = Rmass_sum(ipoin) + b_i(istep)*Rmass(ipoin)
                     sigMass(ipoin,2) = sigMass(ipoin,2) + abs(pt(ipoin,1)*(b_i(istep)-b_i2(istep))*Rmass(ipoin))/kappa
                     Rener(ipoin) = -Rener(ipoin)-(E(ipoin,2)-4.0_rp*E(ipoin,1)/3.0_rp + E(ipoin,3)/3.0_rp)/(dt*2.0_rp/3.0_rp)
                     Rener_sum(ipoin) = Rener_sum(ipoin) + b_i(istep)*Rener(ipoin)
                     sigEner(ipoin,2) = sigEner(ipoin,2) + abs(pt(ipoin,2)*(b_i(istep)-b_i2(istep))*Rener(ipoin))/kappa
                     !$acc loop seq
                     do idime = 1,ndime
                        Rmom(ipoin,idime) = -Rmom(ipoin,idime)-(q(ipoin,idime,2)-4.0_rp*q(ipoin,idime,1)/3.0_rp + q(ipoin,idime,3)/3.0_rp)/(dt*2.0_rp/3.0_rp)
                        Rmom_sum(ipoin,idime) = Rmom_sum(ipoin,idime) + b_i(istep)*Rmom(ipoin,idime)
                        sigMom(ipoin,idime,2) = sigMom(ipoin,idime,2) + abs(pt(ipoin,idime+2)*(b_i(istep)-b_i2(istep))*Rmom(ipoin,idime))/kappa
                     end do
                     !$acc loop seq
                     do jstep=istep+1,flag_pseudo_steps
                        aijKjMass(ipoin,jstep) = aijKjMass(ipoin,jstep) + a_ij(jstep,istep)*Rmass(ipoin)
                        aijKjEner(ipoin,jstep) = aijKjEner(ipoin,jstep) + a_ij(jstep,istep)*Rener(ipoin)
                        !$acc loop seq
                        do idime = 1,ndime
                           aijKjMom(ipoin,idime,jstep) = aijKjMom(ipoin,idime,jstep) + a_ij(jstep,istep)*RMom(ipoin,idime)
                        end do
                     end do
                     !if(istep==4) then
                     !   if(mpi_rank.eq.0)print*, " sig ",sigMass(ipoin,2)," kappa ",kappa," dt ",pt(ipoin,1)
                     !end if
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
                  aux = rho(ipoin,pos) 
                  rho(ipoin,pos) = rho(ipoin,pos)+pt(ipoin,1)*Rmass_sum(ipoin)
                  Rmass(ipoin) = (rho(ipoin,pos)-aux)*Ml(ipoin)/volT
                  aux = E(ipoin,pos)
                  E(ipoin,pos) = (E(ipoin,pos)+pt(ipoin,2)*Rener_sum(ipoin))
                  Rener(ipoin) = (E(ipoin,pos)-aux)*Ml(ipoin)/volT
                  !$acc loop seq
                  do idime = 1,ndime
                     aux = q(ipoin,idime,pos)
                     q(ipoin,idime,pos) = q(ipoin,idime,pos)+pt(ipoin,idime+2)*Rmom_sum(ipoin,idime)
                     Rmom(ipoin,idime) = (q(ipoin,idime,pos)-aux)*Ml(ipoin)/volT
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
                  gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,pos),u(:,:,pos),Reta,alpha)
               call nvtxEndRange


               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(Reta)
                  call nvtxEndRange
               end if

               call nvtxStartRange("Lumped mass solver on generic")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)
               call nvtxEndRange

               call nvtxStartRange("Update sign Reta")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin))!-(3.0_rp*eta(lpoin_w(ipoin),2)-4.0_rp*eta(lpoin_w(ipoin),1)+eta(lpoin_w(ipoin),3))/(2.0_rp*dt)
               end do
               !$acc end parallel loop
               call nvtxEndRange

               !
               ! Compute entropy viscosity
               !
               call nvtxStartRange("Entropy viscosity evaluation")
               call smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Rrho,Ngp,coord,dNgp,gpvol,wgp, &
                  gamma_gas,rho(:,pos),u(:,:,pos),csound,Tem(:,pos),eta(:,pos),helem_l,helem,Ml,mu_e)
               call nvtxEndRange

               !!! Dot products forr errNorm
               call nvtxStartRange("Dot products for errNorm")
               aux = 0.0_rp
               !$acc parallel loop reduction(+:aux)
               do ipoin = 1,npoin_w
                  aux = aux + real(rho(lpoin_w(ipoin),pos)**2,8)
               end do
               !$acc end parallel loop
               auxN(1) = aux

               aux = 0.0_rp
               !$acc parallel loop reduction(+:aux)
               do ipoin = 1,npoin_w
                  aux = aux + real(E(lpoin_w(ipoin),pos)**2,8)
               end do
               !$acc end parallel loop
               auxN(2) = aux

               do idime = 1,ndime
                  aux = 0.0_rp
                  !$acc parallel loop reduction(+:aux)
                  do ipoin = 1,npoin_w
                     aux = aux + real(q(lpoin_w(ipoin),idime,pos)**2,8)
                  end do
                  !$acc end parallel loop
                  auxN(idime+2) = aux
               end do
               call nvtxEndRange

               call nvtxStartRange("Accumullatee auxN in auxN2")
               call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange

               if(auxN2(1)<epsilon(aux)) auxN2(1) = 1.0_rp
               if(auxN2(2)<epsilon(aux)) auxN2(2) = 1.0_rp
               if(auxN2(3)<epsilon(aux)) auxN2(3) = 1.0_rp
               if(auxN2(4)<epsilon(aux)) auxN2(4) = 1.0_rp
               if(auxN2(5)<epsilon(aux)) auxN2(5) = 1.0_rp

               epsR = 1.0_rp/(sqrt(real(auxN2(1),rp)))
               epsE = 1.0_rp/(sqrt(real(auxN2(2),rp)))
               epsQ(1) = 1.0_rp/(sqrt(real(auxN2(3),rp)))
               epsQ(2) = 1.0_rp/(sqrt(real(auxN2(4),rp)))
               epsQ(3) = 1.0_rp/(sqrt(real(auxN2(5),rp)))

               call nvtxStartRange("Dot products of the residuals")
               aux = 0.0_rp
               !$acc parallel loop reduction(+:aux)
               do ipoin = 1,npoin_w
                  aux = aux + real(Rmass(lpoin_w(ipoin)),8)
               end do
               !$acc end parallel loop
               auxN(1) = aux

               aux = 0.0_rp
               !$acc parallel loop reduction(+:aux)
               do ipoin = 1,npoin_w
                  aux = aux + real(Rener(lpoin_w(ipoin)),8)
               end do
               !$acc end parallel loop
               auxN(2) = aux

               do idime = 1,ndime
                  aux = 0.0_rp
                  !$acc parallel loop reduction(+:aux)
                  do ipoin = 1,npoin_w
                     aux = aux + real(Rmom(lpoin_w(ipoin),idime),8)
                  end do
                  !$acc end parallel loop
                  auxN(idime+2) = aux
               end do
               call nvtxEndRange

               call nvtxStartRange("Accumullatee auxN in auxN2")
               call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange

               epsR = (sqrt(real(auxN2(1),rp)))*epsR
               epsE = (sqrt(real(auxN2(2),rp)))*epsE
               epsQ(1) = (sqrt(real(auxN2(3),rp)))*epsQ(1)
               epsQ(2) = (sqrt(real(auxN2(4),rp)))*epsQ(2)
               epsQ(3) = (sqrt(real(auxN2(5),rp)))*epsQ(3)

               errMax = 0.0_rp
               errMax = max(errMax, epsR)
               errMax = max(errMax, epsE)
               errMax = max(errMax, epsQ(1))
               errMax = max(errMax, epsQ(2))
               errMax = max(errMax, epsQ(3))

               if(errMax .lt. tol) exit
               !if(mpi_rank.eq.0)print*, " err ",errMax," it ",itime,' emass ',epsR," eener ",epsE," emom ", epsQ(1)," ", epsQ(2)," ",epsQ(3)
            end do
            call nvtxEndRange

            if(mpi_rank.eq.0)print*, " err ",errMax," it ",itime,' emass ',epsR," eener ",epsE," emom ", epsQ(1)," ", epsQ(2)," ",epsQ(3)

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

            call nvtxStartRange("Last update")
            !$acc kernels
            rho(:,3) = rho(:,1)
            E(:,3) = E(:,1)
            q(:,:,3) = q(:,:,1)
            eta(:,3) = eta(:,1)
            !$acc end kernels
            call nvtxEndRange

         end subroutine rk_pseudo_main

         subroutine rk_4_main(noBoundaries,isWallModelOn,flag_predic,flag_emac,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: flag_predic, flag_emac
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn(nboun,npbou),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem) !helem_l(npoin) TO REVIEW I THINK IS BUG!
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: mu_factor(npoin)
            real(rp),             intent(in)    :: Rgas, gamma_gas, Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,2)
            real(rp),             intent(inout) :: u(npoin,ndime,2)
            real(rp),             intent(inout) :: q(npoin,ndime,2)
            real(rp),             intent(inout) :: pr(npoin,2)
            real(rp),             intent(inout) :: E(npoin,2)
            real(rp),             intent(inout) :: Tem(npoin,2)
            real(rp),             intent(inout) :: e_int(npoin,2)
            real(rp),             intent(inout) :: eta(npoin,2)
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
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode
            real(rp),    dimension(npoin)       :: Reta, Rrho
            real(rp),    dimension(4)           :: a_i, b_i, c_i
            real(rp),    dimension(npoin,ndime) :: aux_u, aux_q,aux_u_wall
            real(rp),    dimension(npoin)       :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta
            real(rp),    dimension(npoin)       :: Rmass, Rener, Rmass_sum, Rener_sum, alpha,Reta_sum
            real(rp),    dimension(npoin,ndime) :: Rmom, Rmom_sum, f_eta
            real(rp)                            :: Rdiff_mass(npoin), Rdiff_mom(npoin,ndime), Rdiff_ener(npoin)
            real(rp)                            :: Aemac(npoin,ndime), Femac(npoin), umag

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            ! Choose between updating prediction or correction
            !
            pos = 2 ! Set correction as default value
            if (flag_predic == 1) then
               pos = 1 ! Change to prediction update
            end if
            !
            ! Butcher tableau
            !
            call nvtxStartRange("Create tableau")
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
            call nvtxEndRange
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

               !
               ! Update velocity and equations of state
               !
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
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,Rdiff_mom)
               end if
               !
               ! Evaluate wall models
     
               if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,aux_rho(:),aux_u(:,:),tauw,Rdiff_mom)
               end if
               !
               !
               ! Compute convective terms
               !
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,aux_u,aux_q,aux_rho,aux_pr,aux_E,Rmass,Rmom,Rener)

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
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass)
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener)
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom)
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
                  f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*eta(lpoin_w(ipoin),1)
               end do
            end do
            !$acc end parallel loop

            call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
               gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta(:,1),u(:,:,1),Reta,alpha)


            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(Reta)
               call nvtxEndRange
            end if

            call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)

            !$acc parallel loop
            do ipoin = 1,npoin_w
               Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin))!-(eta(lpoin_w(ipoin),2)-eta(lpoin_w(ipoin),1))/dt
            end do
            !$acc end parallel loop

            call nvtxEndRange
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
               gamma_gas,rho(:,pos),u(:,:,pos),csound,Tem(:,pos),eta(:,pos),helem_l,helem,Ml,mu_e)
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
                     xi = (1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))) 
                  end if
               end if
               !west 
               if(flag_buffer_on_west .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),1)
                  if(xs<flag_buffer_w_min) then
                     xb = (flag_buffer_w_min-xs)/flag_buffer_w_size
                     xi = (1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))) 
                  end if
               end if
               !north 
               if(flag_buffer_on_north .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),2)
                  if(xs>flag_buffer_n_min) then
                     xb = (xs-flag_buffer_n_min)/flag_buffer_n_size
                     xi = (1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))) 
                  end if
               end if
               !south
               if(flag_buffer_on_south .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),2)
                  if(xs<flag_buffer_s_min) then
                     xb = (flag_buffer_s_min-xs)/flag_buffer_s_size
                     xi = (1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))) 
                  end if
               end if
               !north 
               if(flag_buffer_on_top .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),3)
                  if(xs>flag_buffer_t_min) then
                     xb = (xs-flag_buffer_t_min)/flag_buffer_t_size
                     xi = (1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))) 
                  end if
               end if
               !bottom
               if(flag_buffer_on_bottom .eqv. .true.) then
                  xs = coord(lpoin_w(ipoin),3)
                  if(xs<flag_buffer_b_min) then
                     xb = (flag_buffer_b_min-xs)/flag_buffer_b_size
                     xi = (1.0_rp-c1*xb*xb)*(1.0_rp-(1.0_rp-exp(c2*xb*xb))/(1.0_rp-exp(c2))) 
                  end if
               end if

               q(lpoin_w(ipoin),1) = u_buffer(lpoin_w(ipoin),1)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),1)-u_buffer(lpoin_w(ipoin),1)*rho(lpoin_w(ipoin)))
               q(lpoin_w(ipoin),2) = u_buffer(lpoin_w(ipoin),2)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),2)-u_buffer(lpoin_w(ipoin),2)*rho(lpoin_w(ipoin)))
               q(lpoin_w(ipoin),3) = u_buffer(lpoin_w(ipoin),3)*rho(lpoin_w(ipoin)) + xi*(q(lpoin_w(ipoin),3)-u_buffer(lpoin_w(ipoin),3)*rho(lpoin_w(ipoin)))

            end do
            !$acc end parallel loop
         end subroutine updateBuffer

      end module time_integ
