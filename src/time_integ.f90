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
   use mod_bc_routines

      contains

         subroutine rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                         ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                         rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional arg

            implicit none

            integer(4),           intent(in)    :: flag_predic, flag_emac
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w)
            integer(4),           intent(in)    :: ppow
            real(8),              intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(8),              intent(in)    :: He(ndime,ndime,ngaus,nelem)
            real(8),              intent(in)    :: gpvol(1,ngaus,nelem)
            real(8),              intent(in)    :: dt, helem(nelem), helem_l(npoin)
            real(8),              intent(in)    :: Ml(npoin)
            real(8),              intent(in)    :: Rgas, gamma_gas
            real(8),              intent(inout) :: rho(npoin,2)
            real(8),              intent(inout) :: u(npoin,ndime,2)
            real(8),              intent(inout) :: q(npoin,ndime,2)
            real(8),              intent(inout) :: pr(npoin,2)
            real(8),              intent(inout) :: E(npoin,2)
            real(8),              intent(inout) :: Tem(npoin,2)
            real(8),              intent(inout) :: e_int(npoin,2)
            real(8),              intent(inout) :: mu_fluid(npoin)
            real(8),              intent(out)   :: mu_e(nelem,ngaus)
            real(8),              intent(out)   :: mu_sgs(nelem,ngaus)
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(ndof), lbnodes(nbnodes)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun,2)
            real(8),    optional, intent(in)    :: source_term(ndime)
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime
            real(8),    dimension(npoin)        :: Reta, Rrho
            real(8),    dimension(4)            :: a_i, b_i, c_i
            real(8),    dimension(npoin,ndime)  :: aux_u, aux_q
            real(8),    dimension(npoin)        :: aux_rho, aux_pr, aux_E, aux_Tem, aux_e_int
            real(8),    dimension(npoin)        :: Rmass, Rener, Rmass_sum, Rener_sum
            real(8),    dimension(npoin,ndime)  :: Rmom, Rmom_sum
            real(8)                             :: Rdiff_mass(npoin), Rdiff_mom(npoin,ndime), Rdiff_ener(npoin)
            real(8)                             :: Aemac(npoin,ndime), Femac(npoin)
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
               a_i = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
               c_i = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
               b_i = [1.0d0, 0.0d0, 0.0d0, 0.0d0]
            else if (flag_rk_order == 2) then
               a_i = [0.0d0, 1.0d0, 0.0d0, 0.0d0]
               c_i = [0.0d0, 1.0d0, 0.0d0, 0.0d0]
               b_i = [0.5d0, 0.5d0, 0.0d0, 0.0d0]
            else if (flag_rk_order == 3) then
               write(1,*) "--| NOT CODED FOR RK3 YET!"
               STOP(1)
            else if (flag_rk_order == 4) then
               a_i = [0.0d0, 0.5d0, 0.5d0, 1.0d0]
               c_i = [0.0d0, 0.5d0, 0.5d0, 1.0d0]
               b_i = [1.0d0/6.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/6.0d0]
            else
               write(1,*) "--| NOT CODED FOR RK > 4 YET!"
               STOP(1)
            end if
            call nvtxEndRange
            !
            ! Initialize variables to zero
            !
            call nvtxStartRange("Initialize variables")
            !$acc kernels
            aux_rho(1:npoin) = 0.0d0
            aux_u(1:npoin,1:ndime) = 0.0d0
            aux_q(1:npoin,1:ndime) = 0.0d0
            aux_pr(1:npoin) = 0.0d0
            aux_E(1:npoin) = 0.0d0
            aux_Tem(1:npoin) = 0.0d0
            aux_e_int(1:npoin) = 0.0d0
            Rdiff_mass(1:npoin) = 0.0d0
            Rdiff_mom(1:npoin,1:ndime) = 0.0d0
            Rdiff_ener(1:npoin) = 0.0d0
            Rmass(1:npoin) = 0.0d0
            Rmom(1:npoin,1:ndime) = 0.0d0
            Rener(1:npoin) = 0.0d0
            Rmass_sum(1:npoin) = 0.0d0
            Rener_sum(1:npoin) = 0.0d0
            Rmom_sum(1:npoin,1:ndime) = 0.0d0
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
               !
               ! Impose boundary conditions
               !
               if (nboun .ne. 0) then
                  call nvtxStartRange("Boundary conditions")
                  call temporary_bc_routine(npoin,nboun,bou_codes,bound,nbnodes,lbnodes,aux_rho,aux_q)
                  call nvtxEndRange
               end if
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
                     0.5d0*dot_product(aux_u(lpoin_w(ipoin),:),aux_u(lpoin_w(ipoin),:))
                  aux_pr(lpoin_w(ipoin)) = aux_rho(lpoin_w(ipoin))*(gamma_gas-1.0d0)*aux_e_int(lpoin_w(ipoin))
                  aux_Tem(lpoin_w(ipoin)) = aux_pr(lpoin_w(ipoin))/(aux_rho(lpoin_w(ipoin))*Rgas)
               end do
               !$acc end parallel loop
               call nvtxEndRange
               !
               ! If updating the correction, compute viscosities and diffusion
               !
               if (flag_predic == 0) then
                  !
                  ! Update residuals for envit
                  !
                  call nvtxStartRange("ENVIT")
                  call residuals(nelem,npoin,npoin_w,lpoin_w, &
                                 ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                 dt, aux_rho, aux_u, aux_pr, aux_q, &
                                 rho, u, pr, q, gamma_gas, &
                                 Reta, Rrho)
                  !
                  ! Compute entropy viscosity
                  !
                  if (flag_SpectralElem == 1) then
                     call smart_visc_spectral(nelem,npoin,connec,Reta,Rrho,Ngp, &
                        gamma_gas,aux_rho,aux_u,aux_Tem,helem_l,helem,mu_e)
                  else
                     call smart_visc(nelem,npoin,connec,Reta,Rrho,Ngp, &
                        gamma_gas,aux_rho,aux_u,aux_Tem,helem,mu_e)
                  end if
                  call nvtxEndRange
                  !
                  ! Compute subgrid viscosity if active
                  !
                  if(flag_les == 1) then
                     call nvtxStartRange("MU_SGS")
                     call sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_rho,aux_u,mu_sgs)
                     call nvtxEndRange
                  end if
                  !
                  ! Update viscosity if Sutherland's law is active
                  !
                  if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
                     call nvtxStartRange("MU_SUT")
                     call sutherland_viscosity(npoin,aux_Tem,mu_fluid)
                     call nvtxEndRange
                  end if
                  !
                  ! Compute diffusion terms with values at current substep
                  !
                  call nvtxStartRange("DIFFUSIONS")
                  call mass_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_rho,mu_e,Rdiff_mass)
                  call mom_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,mu_fluid,mu_e,mu_sgs,Rdiff_mom)
                  call ener_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,aux_Tem,mu_fluid,mu_e,mu_sgs,Rdiff_ener)
                  call nvtxEndRange
                  !
                  ! Call source term if applicable
                  !
                  if(present(source_term)) then
                     call mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,source_term,Rdiff_mom)
                  end if
               end if
               !
               ! Compute convective terms
               !
               call mass_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_q,aux_rho,aux_u,Rmass)
               call ener_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,aux_pr,aux_E,aux_rho,aux_q,Rener)
               if (flag_emac .eq. 0) then
                  !
                  ! Conservation momentum convection div(qi*uj)+grad(p)
                  !
                  call mom_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,aux_u,aux_q,aux_pr,Rmom)
               else if (flag_emac .eq. 1) then
                  !
                  ! EMAC term
                  !
                  !$acc parallel loop
                  do ipoin = 1,npoin
                     !$acc loop seq
                     do idime = 1,ndime
                        Aemac(ipoin,idime) = aux_u(ipoin,idime)*sqrt(aux_rho(ipoin))
                     end do
                     Femac(ipoin) = dot_product(Aemac(ipoin,:),Aemac(ipoin,:))
                  end do
                  !$acc end parallel loop
                  call mom_convec_emac(nelem,npoin,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,aux_pr,Rmom)
               end if
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
               !
               ! Call lumped mass matrix solver
               !
               call nvtxStartRange("Call solver")
               if (flag_solver_type == 1) then ! Lumped mass
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass)
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener)
                  call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom)
               else if (flag_solver_type == 2) then ! Lumped+apinv
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rmass)
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rener)
                  call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,Rmom)
                  call approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,Rmass)
                  call approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,Rener)
                  call approx_inverse_vect(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,Rmom)
               else if (flag_solver_type == 3) then ! ConjGrad
                  call conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,Rmass)
                  call conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,Rener)
                  call conjGrad_vector(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,Rmom)
               else 
                  write(1,*) "--| SOLVER NOT CODED YET!"
                  STOP(1)
               end if
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
            !
            ! Apply bcs after update
            !
            if (nboun .ne. 0) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call temporary_bc_routine(npoin,nboun,bou_codes,bound,nbnodes,lbnodes,rho(:,pos),q(:,:,pos))
               call nvtxEndRange
            end if
            !
            ! Update velocity and equations of state
            !
            call nvtxStartRange("Update u and EOS")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               !$acc loop seq
               do idime = 1,ndime
                  u(lpoin_w(ipoin),idime,pos) = q(lpoin_w(ipoin),idime,pos)/rho(lpoin_w(ipoin),pos)
               end do
               e_int(lpoin_w(ipoin),pos) = (E(lpoin_w(ipoin),pos)/rho(lpoin_w(ipoin),pos))- &
                  0.5d0*dot_product(u(lpoin_w(ipoin),:,pos),u(lpoin_w(ipoin),:,pos))
               pr(lpoin_w(ipoin),pos) = rho(lpoin_w(ipoin),pos)*(gamma_gas-1.0d0)*e_int(lpoin_w(ipoin),pos)
               Tem(lpoin_w(ipoin),pos) = pr(lpoin_w(ipoin),pos)/(rho(lpoin_w(ipoin),pos)*Rgas)
            end do
            !$acc end parallel loop
            call nvtxEndRange
            !
            ! If using Sutherland viscosity model:
            !
            if (flag_real_diff == 1 .and. flag_diff_suth == 1) then
               call nvtxStartRange("Sutherland viscosity")
               call sutherland_viscosity(npoin,Tem(:,pos),mu_fluid)
               call nvtxEndRange
            end if
         end subroutine rk_4_main
end module time_integ
