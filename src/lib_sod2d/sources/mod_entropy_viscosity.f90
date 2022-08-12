module mod_entropy_viscosity

   use mod_constants
   use mod_nvtx
   use mod_veclen

      ! TODO: Finish module and create unit tests

      contains

              subroutine residuals(nelem,npoin,npoin_w,lpoin_w, &
                                   ppow, connec, Ngp, dNgp, He, gpvol, Ml, &
                                   dt, rhok, uk, prk, qk, &
                                   rho, u, pr, q, gamma_gas, &
                                   Reta, Rrho)

                      use mass_matrix
                      use mod_solver
                      use elem_convec

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin, ppow, npoin_w
                      integer(4), intent(in)  :: connec(nelem,nnode), lpoin_w(npoin_w)
                      real(rp),    intent(in)  :: Ngp(nnode,ngaus), dNgp(ndime,nnode,ngaus)
                      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(rp),    intent(in)  :: gpvol(1,ngaus,nelem), Ml(npoin)
                      real(rp),    intent(in)  :: dt, gamma_gas
                      real(rp),    intent(in)  :: rhok(npoin), uk(npoin,ndime), prk(npoin), qk(npoin,ndime)     ! From substep
                      real(rp),    intent(in)  :: rho(npoin,2), u(npoin,ndime,2), pr(npoin,2), q(npoin,ndime,2) ! From prediction
                      real(rp),    intent(out) :: Reta(npoin), Rrho(npoin)
                      integer(4)              :: ipoin, idime
                      real(rp)                 :: eta(npoin), eta_p(npoin), alpha(npoin), alpha_p(npoin)
                      real(rp)                 :: f_eta(npoin,ndime), f_rho(npoin,ndime), R1(npoin), R2(npoin)
                      real(rp)                 :: aux1(npoin), maxEta, maxRho,vol

                       !
                       ! Entropy function and temporal terms
                       !
                       call nvtxStartRange("Entropy transport")
                       !$acc parallel loop
                       do ipoin = 1,npoin_w

                          eta(lpoin_w(ipoin)) = (rhok(lpoin_w(ipoin))/(gamma_gas-1.0_rp))* &
                             log(prk(lpoin_w(ipoin))/(rhok(lpoin_w(ipoin))**gamma_gas))
                          eta_p(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),2)/(gamma_gas-1.0_rp))* &
                             log(pr(lpoin_w(ipoin),2)/(rho(lpoin_w(ipoin),2)**gamma_gas))
                          !eta_p(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),1)/(gamma_gas-1.0_rp))* &
                          !   log(pr(lpoin_w(ipoin),1)/(rho(lpoin_w(ipoin),1)**gamma_gas))

                          !$acc loop vector
                          do idime = 1,ndime
                             f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,2)*eta_p(lpoin_w(ipoin))
!                             f_rho(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime,2)
                          end do

                          R1(lpoin_w(ipoin)) = (eta(lpoin_w(ipoin))-eta_p(lpoin_w(ipoin)))/dt  ! Temporal entropy
                          !R1(lpoin_w(ipoin)) = (eta_p(lpoin_w(ipoin))-eta(lpoin_w(ipoin)))/dt  ! Temporal entropy
                          !Reta(lpoin_w(ipoin)) = 0.0_rp
                          alpha(lpoin_w(ipoin)) = 1.0_rp

                       end do
                       !$acc end parallel loop

                       !
                       ! Entropy residual
                       !
                       call generic_scalar_convec(nelem,npoin,connec,Ngp,dNgp,He, &
                                                  gpvol,f_eta,eta_p,u(:,:,2),Reta,alpha)
                       !
                       ! Alter Reta with inv(Mc)
                       !
                       if (flag_solver_type == 1) then
                          call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)
                       else if (flag_solver_type == 2) then
                          call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)
                          call approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,Reta)
                       else if (flag_solver_type == 3) then
                          call conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,Reta)
                       else
                          write(1,*) "--| SOLVER NOT CODED!"
                          STOP(1)
                       end if
                       !
                       ! Update Reta
                       !
                      !!$acc parallel loop
                      !do ipoin = 1,npoin_w
                      !   Reta(lpoin_w(ipoin)) = Reta(lpoin_w(ipoin))+1.0_rp*R1(lpoin_w(ipoin))

                      !   R2(lpoin_w(ipoin)) = (rhok(lpoin_w(ipoin))-rho(lpoin_w(ipoin),2))/dt
                      !   !R2(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),1)-rhok(lpoin_w(ipoin)))/dt
                      !   alpha(lpoin_w(ipoin)) = eta(lpoin_w(ipoin))/rhok(lpoin_w(ipoin))
                      !end do
                      !!$acc end parallel loop
                      !!
                      !! Alter R2 with Mcw
                      !!
                      !if (flag_solver_type == 2 .or. flag_solver_type == 3) then
                      !   call wcmass_times_vector(nelem,npoin,connec,gpvol,Ngp,R2,aux1,alpha)
                      !else if (flag_solver_type == 1) then
                      !   !call wlmass_times_vector()
                      !   !write(1,*) "--| NOT CODED YET!"
                      !   !STOP(1)
                      !end if
                      !!
                      !! Compute weighted mass convec
                      !!
                      !! oriol: is no needed generic_scalar_conv already does
                      !! this
                      !call generic_scalar_convec(nelem,npoin,connec,Ngp, &
                      !                           dNgp,He,gpvol,f_rho,rho,u,Rrho,alpha)
                      !!
                      !! Update Rrho with both terms
                      !!
                      !!$acc parallel loop
                      !do ipoin = 1,npoin_w
                      !   Rrho(lpoin_w(ipoin)) = Rrho(lpoin_w(ipoin))+1.0_rp*aux1(lpoin_w(ipoin))
                      !end do
                      !!$acc end parallel loop
                      !!
                      !! Apply solver
                      !!
                      !if (flag_solver_type == 1) then
                      !   call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rrho)
                      !else if (flag_solver_type == 2) then
                      !   call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rrho)
                      !   call approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,Rrho)
                      !else if (flag_solver_type == 3) then
                      !   call conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,Rrho)
                      !else
                      !   write(1,*) "--| SOLVER NOT CODED!"
                      !   STOP(1)
                      !end if
                      !call nvtxEndRange

                  !!
                  !! Normalize
                  !!
               !   maxEta = 0.0_rp !maxval(abs(eta(lpoin_w(:))))
               !   vol = 0.0_rp

               !   !$acc parallel loop reduction(+:maxEta,vol)
               !   do ipoin = 1,npoin_w
               !      maxEta = maxEta + eta(lpoin_w(ipoin))*Ml(lpoin_w(ipoin))
               !      vol    = vol + Ml(lpoin_w(ipoin))
               !   end do
               !   !$acc end parallel loop
               !   maxEta = maxEta/vol

                !  !$acc parallel loop 
                !  do ipoin = 1,npoin_w
                !     Reta(lpoin_w(ipoin)) = abs(Reta(lpoin_w(ipoin)))/abs(maxEta)
                !     !Reta(lpoin_w(ipoin)) = abs(Reta(lpoin_w(ipoin)))/abs(maxEta-eta(lpoin_w(ipoin)))
                !  end do

              end subroutine residuals

              subroutine smart_visc(nelem,npoin,connec,Reta,Rrho,Ngp, &
                                    gamma_gas,rho,u,Tem,helem,mu_e)
              
                      ! TODO: Compute element size h
              
                      implicit none

                      integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
                      real(rp),    intent(in)  :: Reta(npoin), Rrho(npoin), Ngp(ngaus,nnode),helem(nelem), gamma_gas
                      real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime), Tem(npoin)
                      real(rp),    intent(out) :: mu_e(nelem,ngaus)
                      integer(4)              :: ielem, inode, igaus
                      real(rp)                 :: R1, R2, Ve
                      real(rp)                 :: betae
                      real(rp)                 :: L3, aux1, aux2, aux3

                      !$acc parallel loop gang vector_length(vecLength)
                      do ielem = 1,nelem
                         !$acc loop worker
                         do igaus = 1,ngaus
                            R1 = 0.0_rp
                            R2 = 0.0_rp
                            !$acc loop vector reduction(+:R1,R2)
                            do inode = 1,nnode
                               R1 = R1+Ngp(igaus,inode)*abs(Reta(connec(ielem,inode))) ! Linf norm of Reta on element
                               R2 = R2+Ngp(igaus,inode)*abs(Rrho(connec(ielem,inode))) ! Linf norm of Rrho on element
                            end do
                            Ve = ce*max(R1,R2)*((helem(ielem)/real(porder,rp))**2) ! Normalized residual for element
                            !
                            ! Max. Wavespeed at element
                            !
                            aux1 = 0.0_rp
                            !$acc loop vector reduction(+:aux1)
                            do inode = 1,nnode
                               aux2 = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
                               aux3 = sqrt(gamma_gas*Tem(connec(ielem,inode)))     ! Speed of sound at node
                               aux1 = aux1+Ngp(igaus,inode)*(aux2+aux3)
                            end do
                            !
                            ! Select against Upwind viscosity
                            !
                            betae = cmax*(helem(ielem)/real(porder,rp))*aux1
                            mu_e(ielem,igaus) = cglob*min(Ve,betae) ! Dynamic viscosity
                         end do
                      end do
                      !$acc end parallel loop

              end subroutine smart_visc

              subroutine smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Rrho,Ngp, &
                                    gamma_gas,rho,u,Tem,eta,helem,helem_k,Ml,mu_e)
              
                      ! TODO: Compute element size h
              
                      implicit none

                      integer(4), intent(in)  :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
                      real(rp),    intent(in)  :: Reta(npoin), Rrho(npoin), Ngp(ngaus,nnode),gamma_gas
                      real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime), Tem(npoin), eta(npoin),helem(nelem,nnode),helem_k(nelem),Ml(npoin)
                      real(rp),    intent(out) :: mu_e(nelem,ngaus)
                      integer(4)              :: ielem, inode,igaus,ipoin
                      real(rp)                 :: R1, R2, Ve
                      real(rp)                 :: betae,mu,vol,vol2
                      real(rp)                 :: L3, aux1, aux2, aux3
                      real(rp)                 ::  maxEta, maxRho, norm, Rgas

                      !Rgas = 1.0_rp*1.0_rp/(1.4_rp*1.0_rp*1.25_rp*1.25_rp)
                      Rgas = nscbc_Rgas_inf

                      norm = 1.0_rp

                      if(flag_normalise_entropy .eq. 1) then
                         maxEta = 0.0_rp
                         !$acc parallel loop reduction(+:maxEta)
                         do ipoin = 1,npoin_w
                            maxEta = maxEta + eta(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                         maxEta = maxEta/dble(npoin_w)

                         norm = 0.0_rp
                         !$acc parallel loop reduction(max:norm)
                         do ipoin = 1,npoin_w
                            norm = max(norm, abs(eta(lpoin_w(ipoin))-maxEta))
                         end do
                         !$acc end parallel loop
                      end if

                      !$acc parallel loop gang vector_length(vecLength)
                      do ielem = 1,nelem
                         mu = 0.0_rp
                         betae = 0.0_rp
                         !$acc loop vector reduction(max:betae)
                         do inode = 1,nnode
                            aux2 = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
                            aux3 = sqrt(Rgas*gamma_gas*Tem(connec(ielem,inode)))     ! Speed of sound at node
                            aux1 = aux2+aux3
                            betae = max(betae,(rho(connec(ielem,inode))*helem_k(ielem))*(cmax/real(porder,rp))*aux1)
                         end do
                         !$acc loop vector reduction(+:mu)
                         do inode = 1,nnode
                            R1 = rho(connec(ielem,inode))*abs(Reta(connec(ielem,inode)))/norm
                            Ve = ce*R1*(helem(ielem,inode))**2
                            mu = mu + cglob*min(Ve,betae)
                         end do
                         mu = mu/real(nnode,rp)
                         !$acc loop vector
                         do inode = 1,nnode
                            mu_e(ielem,inode) = mu
                         end do
                      end do
                      !$acc end parallel loop

              end subroutine smart_visc_spectral
end module mod_entropy_viscosity
