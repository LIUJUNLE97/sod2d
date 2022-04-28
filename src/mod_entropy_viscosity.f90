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
                      real(8),    intent(in)  :: Ngp(nnode,ngaus), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem), Ml(npoin)
                      real(8),    intent(in)  :: dt, gamma_gas
                      real(8),    intent(in)  :: rhok(npoin), uk(npoin,ndime), prk(npoin), qk(npoin,ndime)     ! From substep
                      real(8),    intent(in)  :: rho(npoin,2), u(npoin,ndime,2), pr(npoin,2), q(npoin,ndime,2) ! From prediction
                      real(8),    intent(out) :: Reta(npoin), Rrho(npoin)
                      integer(4)              :: ipoin, idime
                      real(8)                 :: eta(npoin), eta_p(npoin), alpha(npoin), alpha_p(npoin)
                      real(8)                 :: f_eta(npoin,ndime), f_rho(npoin,ndime), R1(npoin), R2(npoin)
                      real(8)                 :: aux1(npoin), maxEta, maxRho

                       !
                       ! Entropy function and temporal terms
                       !
                       call nvtxStartRange("Entropy transport")
                       !$acc parallel loop
                       do ipoin = 1,npoin_w

                          eta(lpoin_w(ipoin)) = (rhok(lpoin_w(ipoin))/(gamma_gas-1.0d0))* &
                             log(prk(lpoin_w(ipoin))/(rhok(lpoin_w(ipoin))**gamma_gas))
                          eta_p(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),1)/(gamma_gas-1.0d0))* &
                             log(pr(lpoin_w(ipoin),1)/(rho(lpoin_w(ipoin),1)**gamma_gas))

                          !$acc loop vector
                          do idime = 1,ndime
                             f_eta(lpoin_w(ipoin),idime) = uk(lpoin_w(ipoin),idime)*eta(lpoin_w(ipoin))
                             f_rho(lpoin_w(ipoin),idime) = qk(lpoin_w(ipoin),idime)
                          end do

                          R1(lpoin_w(ipoin)) = (eta_p(lpoin_w(ipoin))-eta(lpoin_w(ipoin)))/dt  ! Temporal entropy
                          !Reta(lpoin_w(ipoin)) = 0.0d0
                          alpha(lpoin_w(ipoin)) = 1.0d0

                       end do
                       !$acc end parallel loop

                       !
                       ! Entropy residual
                       !
                       call generic_scalar_convec(nelem,npoin,connec,Ngp,dNgp,He, &
                                                  gpvol,f_eta,Reta,alpha)
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
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          Reta(lpoin_w(ipoin)) = Reta(lpoin_w(ipoin))+1.0d0*R1(lpoin_w(ipoin))

                          R2(lpoin_w(ipoin)) = (rho(lpoin_w(ipoin),1)-rhok(lpoin_w(ipoin)))/dt
                          alpha(lpoin_w(ipoin)) = eta(lpoin_w(ipoin))/rhok(lpoin_w(ipoin))
                       end do
                       !$acc end parallel loop
                       !
                       ! Alter R2 with Mcw
                       !
                       if (flag_solver_type == 2 .or. flag_solver_type == 3) then
                          call wcmass_times_vector(nelem,npoin,connec,gpvol,Ngp,R2,aux1,alpha)
                       else if (flag_solver_type == 1) then
                          !call wlmass_times_vector()
                          write(1,*) "--| NOT CODED YET!"
                          STOP(1)
                       end if
                       !
                       ! Compute weighted mass convec
                       !
                       ! oriol: is no needed generic_scalar_conv already does
                       ! this
                       call generic_scalar_convec(nelem,npoin,connec,Ngp, &
                                                  dNgp,He,gpvol,f_rho,Rrho,alpha)
                       !
                       ! Update Rrho with both terms
                       !
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          Rrho(lpoin_w(ipoin)) = Rrho(lpoin_w(ipoin))+1.0d0*aux1(lpoin_w(ipoin))
                       end do
                       !$acc end parallel loop
                       !
                       ! Apply solver
                       !
                       if (flag_solver_type == 1) then
                          call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rrho)
                       else if (flag_solver_type == 2) then
                          call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rrho)
                          call approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,Rrho)
                       else if (flag_solver_type == 3) then
                          call conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,Rrho)
                       else
                          write(1,*) "--| SOLVER NOT CODED!"
                          STOP(1)
                       end if
                       call nvtxEndRange

                       !!
                       !! Normalize
                       !!
                       !maxEta = maxval(abs(eta(lpoin_w(:))))
                       !maxRho = maxval(abs(rhok(lpoin_w(:))))

                       !!$acc parallel loop
                       !do ipoin = 1,npoin_w
                       !   Reta(lpoin_w(ipoin)) = Reta(lpoin_w(ipoin))/maxEta
                       !   Rrho(lpoin_w(ipoin)) = Rrho(lpoin_w(ipoin))/maxRho
                       !end do
                       !!$acc end parallel loop

              end subroutine residuals

              subroutine smart_visc(nelem,npoin,connec,Reta,Rrho,Ngp, &
                                    gamma_gas,rho,u,Tem,helem,mu_e)
              
                      ! TODO: Compute element size h
              
                      implicit none

                      integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
                      real(8),    intent(in)  :: Reta(npoin), Rrho(npoin), Ngp(ngaus,nnode),helem(nelem), gamma_gas
                      real(8),    intent(in)  :: rho(npoin), u(npoin,ndime), Tem(npoin)
                      real(8),    intent(out) :: mu_e(nelem,ngaus)
                      integer(4)              :: ielem, inode, igaus
                      real(8)                 :: R1, R2, Ve
                      real(8)                 :: betae
                      real(8)                 :: L3, aux1, aux2, aux3

                      !$acc parallel loop gang vector_length(vecLength)
                      do ielem = 1,nelem
                         !$acc loop worker
                         do igaus = 1,ngaus
                            R1 = 0.0d0
                            R2 = 0.0d0
                            !$acc loop vector reduction(+:R1,R2)
                            do inode = 1,nnode
                               R1 = R1+Ngp(igaus,inode)*abs(Reta(connec(ielem,inode))) ! Linf norm of Reta on element
                               R2 = R2+Ngp(igaus,inode)*abs(Rrho(connec(ielem,inode))) ! Linf norm of Rrho on element
                            end do
                            Ve = ce*max(R1,R2)*((helem(ielem)/dble(porder))**2) ! Normalized residual for element
                            !
                            ! Max. Wavespeed at element
                            !
                            aux1 = 0.0d0
                            !$acc loop vector reduction(+:aux1)
                            do inode = 1,nnode
                               aux2 = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
                               aux3 = sqrt(gamma_gas*Tem(connec(ielem,inode)))     ! Speed of sound at node
                               aux1 = aux1+Ngp(igaus,inode)*(aux2+aux3)
                            end do
                            !
                            ! Select against Upwind viscosity
                            !
                            betae = cmax*(helem(ielem)/dble(porder))*aux1
                            mu_e(ielem,igaus) = cglob*min(Ve,betae) ! Dynamic viscosity
                         end do
                      end do
                      !$acc end parallel loop

              end subroutine smart_visc
end module mod_entropy_viscosity
