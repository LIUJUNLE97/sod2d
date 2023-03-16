module time_integ

#define PSEUDO_RK 0

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

      real(rp)                                  :: volT
      logical                                   :: flag_mem_alloc=.true.
      logical                                   :: initialize_pt=.true.

      contains

         subroutine rk_pseudo_main(igtime,noBoundaries,isWallModelOn,flag_predic,flag_emac,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas,Cp,Prt, &
                         rho,u,q,pr,E,Tem,csound,machno,e_int,eta,mu_e,mu_sgs,kres,etot,au,ax1,ax2,ax3,lpoin_w,mu_fluid,mu_factor, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer,tauw,source_term)  ! Optional args

            implicit none

            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: igtime, flag_predic, flag_emac
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
            integer(4)                          :: pos,maxIterL
            integer(4)                          :: istep, ipoin, idime,icode,itime,jstep,inode,ielem,npoin_w_g
            real(rp),    dimension(npoin)       :: Reta, Rrho
            real(rp),    dimension(npoin,ndime) :: aux_u, aux_q,aux_u_wall
            real(rp),    dimension(npoin)       :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int,aux_eta
            real(rp),    dimension(npoin)       :: Rmass, Rener, Rmass_sum, Rener_sum, alpha,dt_min
            real(rp),    dimension(npoin,ndime) :: Rmom, Rmom_sum, f_eta
            real(rp)                            :: Rdiff_mass(npoin), Rdiff_mom(npoin,ndime), Rdiff_ener(npoin),alfa_pt(5), f_alpha,pt_g,dt_min_g
            real(rp)                            :: umag,aux,vol_rank,errMax
            real(8)                             :: auxN(5),auxN2(5),vol_tot_d, res(2),aux2,res_ini


            if (flag_mem_alloc .eqv. .true.) then
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
               flag_mem_alloc = .false.
            end if


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            ! Choose between updating prediction or correction
            !
            pos = 2 ! Set correction as default value
            
            call expl_adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,0.1_rp,pt_g,0.1_rp,mu_fluid,mu_sgs,rho(:,2))
            dt_min_g = pt_g

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
            !$acc end kernels
            call nvtxEndRange
            !
            ! Loop over all RK steps
            !

            res(:) = 0.0
            call nvtxStartRange("Loop over RK steps")
           ! if(igtime .lt. 5) then 
            !   maxIterL = 5
            !else              
               maxIterL = maxIterNonLineal
            !end if
            do itime =1, maxIterL
               !
               ! Compute diffusion terms with values at current substep
               !
               call nvtxStartRange("DIFFUSIONS")
               call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),u(:,:,pos),Tem(:,pos),mu_fluid,mu_e,mu_sgs,Ml,Rdiff_mass,Rdiff_mom,Rdiff_ener)
               call nvtxEndRange
               !
               ! Call source term if applicable
               !
               if(present(source_term)) then
                  call nvtxStartRange("SOURCE TERM")
                  call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u(:,:,pos),source_term,Rdiff_mom)
                  call nvtxEndRange
               end if
               !
               ! Evaluate wall models
      
               if((isWallModelOn) .and. (numBoundsWM .ne. 0)) then
                  call nvtxStartRange("WALL MODEL")
                  call evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK,bou_codes,&
                     bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol, mu_fluid,rho(:,pos),u(:,:,pos),tauw,Rdiff_mom)
                  call nvtxEndRange
               end if
               !
               !
               ! Compute convective terms
               !
               call nvtxStartRange("CONVECTIONS")
               call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u(:,:,pos),q(:,:,pos),rho(:,pos),pr(:,pos),E(:,pos),Rmass,Rmom,Rener)
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

               aux2 = 0.0_rp
               !$acc parallel loop reduction(+:aux2)
               do ipoin = 1,npoin
                  Rmass(ipoin) = -(Rmass(ipoin)+(rho(ipoin,2)-4.0_rp*rho(ipoin,1)/3.0_rp + rho(ipoin,3)/3.0_rp)/(dt*2.0_rp/3.0_rp))
                  aux2 = aux2 + real(Rmass(ipoin)**2,8)
                  Rener(ipoin) = -(Rener(ipoin)+(E(ipoin,2)-4.0_rp*E(ipoin,1)/3.0_rp + E(ipoin,3)/3.0_rp)/(dt*2.0_rp/3.0_rp))
                  aux2 = aux2 + real(Rener(ipoin)**2,8)
                  !$acc loop seq
                  do idime = 1,ndime
                     Rmom(ipoin,idime) = -(Rmom(ipoin,idime)+(q(ipoin,idime,2)-4.0_rp*q(ipoin,idime,1)/3.0_rp + q(ipoin,idime,3)/3.0_rp)/(dt*2.0_rp/3.0_rp))
                     aux2 = aux2 + real(Rmom(ipoin,idime)**2,8)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               !$acc parallel loop
               do ipoin = 1,npoin
                  aux_rho(ipoin) = 0.0_rp
                  aux_E(ipoin)   = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_q(ipoin,idime) = 0.0_rp
                  end do
               end do
               !$acc end parallel loop
               call nvtxStartRange("GMRES")
               call gmres_full(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                              atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                              rho(:,pos),u(:,:,pos),q(:,:,pos),pr(:,pos),E(:,pos),Tem(:,pos),Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                              2.0_rp/3.0_rp,dt,pt_g,Rmass,Rmom,Rener,aux_rho,aux_q,aux_E,1,0.001_rp)
               call nvtxEndRange

               !$acc parallel loop
               do ipoin = 1,npoin
                  Rmass(ipoin) = aux_rho(ipoin)*Ml(ipoin)/volT
                  rho(ipoin,pos) = rho(ipoin,pos)+aux_rho(ipoin)
                  Rener(ipoin) = aux_E(ipoin)*Ml(ipoin)/volT
                  E(ipoin,pos)   = E(ipoin,pos)+aux_E(ipoin)
                  !$acc loop seq
                  do idime = 1,ndime
                     Rmom(ipoin,idime) = aux_q(ipoin,idime)*Ml(ipoin)/volT
                     q(ipoin,idime,pos) = q(ipoin,idime,pos)+aux_q(ipoin,idime)
                  end do
               end do
               !$acc end parallel loop

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
                     f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,pos)*eta(lpoin_w(ipoin),pos)
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
                  Reta(lpoin_w(ipoin)) = -Reta(lpoin_w(ipoin))-(3.0_rp*eta(lpoin_w(ipoin),2)-4.0_rp*eta(lpoin_w(ipoin),1)+eta(lpoin_w(ipoin),3))/(2.0_rp*dt)
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
                  aux = aux + real(Rmass(lpoin_w(ipoin))**2,8)
               end do
               !$acc end parallel loop
               auxN(1) = aux

               aux = 0.0_rp
               !$acc parallel loop reduction(+:aux)
               do ipoin = 1,npoin_w
                  aux = aux + real(Rener(lpoin_w(ipoin))**2,8)
               end do
               !$acc end parallel loop
               auxN(2) = aux

               do idime = 1,ndime
                  aux = 0.0_rp
                  !$acc parallel loop reduction(+:aux)
                  do ipoin = 1,npoin_w
                     aux = aux + real(Rmom(lpoin_w(ipoin),idime)**2,8)
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


               call nvtxStartRange("Accumullatee auxN in auxN2")
               call MPI_Allreduce(aux2,res(1),1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
               call nvtxEndRange

               res(1) = sqrt(res(1))

               if(itime .gt. 1) then
                  pt_g = min(pt_g*(res(2)/res(1)),pseudo_max_dt)
               else
                  res_ini = res(1)
               endif

               errMax = abs(res(1)-res(2))/abs(res_ini)

               res(2) = res(1)


               if(errMax .lt. tol) exit  
               !if(mpi_rank.eq.0)print*, " #### err ",errMax," it ",itime,' emass ',epsR," eener ",epsE," emom ", epsQ(1)," ", epsQ(2)," ",epsQ(3)," pt ",pt_g
        
            end do
            call nvtxEndRange

            if(mpi_rank.eq.0)print*,"time ",igtime, "err ",errMax," it ",itime,' emass ',epsR," eener ",epsE," emom ", epsQ(1)," ", epsQ(2)," ",epsQ(3)," pt ",pt_g

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
