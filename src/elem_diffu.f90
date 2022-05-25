module elem_diffu

      use mod_nvtx
      use mod_constants
      use mod_veclen

      ! TODO: Create unit tests for all subroutines

      contains

              subroutine mass_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,rho,mu_e,Rmass)

                      ! TODO: Add stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: rho(npoin), mu_e(nelem,ngaus)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode), nu_e
                      real(8)                 :: tmp1, gpcar(ndime,nnode), tmp2

                      call nvtxStartRange("Mass diffusion")
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang  private(gpcar,Re) vector_length(vecLength)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                           !$acc loop seq
                           do idime = 1,ndime
                              tmp1 = 0.0d0
                              tmp2 = 0.0d0
                              !$acc loop vector reduction(+:tmp1,tmp2)
                              do inode = 1,nnode
                                 tmp1 = tmp1+gpcar(idime,inode)*rho(connec(ielem,inode))
                                 tmp2 = tmp2+Ngp(igaus,inode)*rho(connec(ielem,inode)) 
                              end do
                              nu_e = c_rho*mu_e(ielem,igaus)/tmp2
                              !$acc loop vector
                              do inode = 1,nnode
                                 Re(inode) = Re(inode)+nu_e*gpvol(1,igaus,ielem)* &
                                             gpcar(idime,inode)*tmp1
                              end do
                           end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mass_diffusion

              subroutine mom_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,mu_fluid,mu_e,mu_sgs,Rmom)

                      ! TODO: Add. stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime),mu_e(nelem,ngaus),mu_sgs(nelem,ngaus)
                      real(8),    intent(in)  :: mu_fluid(npoin)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      integer(4)              :: ielem, igaus, idime, jdime, inode, jnode
                      real(8)                 :: Re(nnode,ndime), twoThirds, gpcar(ndime,nnode), tau(ndime,ndime)
                      real(8)                 :: divU, gradU(ndime,ndime), mu_fgp, aux

                      twoThirds = 2.0d0/3.0d0
                      call nvtxStartRange("Momentum diffusion")
                      !$acc kernels
                      Rmom(:,:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang  private(Re,gpcar,tau,gradU) vector_length(vecLength)
                      do ielem = 1,nelem
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               Re(inode,idime) = 0.0d0
                            end do
                         end do
                         !
                         ! Gradient structure:
                         !
                         !         | u1,1 u1,2 u1,3 |
                         ! u_i,j = | u2,1 u2,2 u2,3 |
                         !         | u3,1 u3,2 u3,3 |
                         !
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !
                            ! Compute gpcar
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !
                            ! Compute combined viscosity at Gauss point
                            !
                            mu_fgp = 0.0d0
                            !$acc loop vector reduction(+:mu_fgp)
                            do inode = 1,nnode
                               mu_fgp = mu_fgp+(Ngp(igaus,inode)*mu_fluid(connec(ielem,inode)))
                            end do
                            mu_fgp = mu_fgp+mu_e(ielem,igaus)+mu_sgs(ielem,igaus)
                            !
                            ! Compute grad(u)
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  aux = 0.0d0
                                  !$acc loop vector reduction(+:aux)
                                  do inode = 1,nnode
                                     aux = aux+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                                  end do
                                  gradU(idime,jdime) = aux
                               end do
                            end do
                            !
                            ! Compute div(u)
                            !
                            divU = 0.0d0
                            ! divU = gradU(1,1)+gradU(2,2)+gradU(3,3)
                            !$acc loop seq
                            do idime = 1,ndime
                               divU = divU+gradU(idime,idime)
                            end do

                            ! TODO: Remove this
                            !divU = 0.0d0
                            !!$acc loop vector collapse(2) reduction(+:divU)
                            !do idime = 1,ndime
                            !   do inode = 1,nnode
                            !      divU = divU+gpcar(idime,inode)*u(connec(ielem,inode),idime)
                            !   end do
                            !end do

                            !
                            ! Finish computing tau_ij = u_i,j + u_j,i - (2/3)*div(u)*d_ij
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  tau(idime,jdime) = gradU(idime,jdime)+gradU(jdime,idime)
                               end do
                            end do
                            !$acc loop seq
                            do idime = 1,ndime
                               tau(idime,idime) = tau(idime,idime)-twoThirds*divU
                            end do
                            !
                            ! Compute div(tau) at the Gauss point
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     Re(inode,idime) = Re(inode,idime)+gpvol(1,igaus,ielem)* &
                                        gpcar(jdime,inode)*mu_fgp*tau(idime,jdime)
                                  end do
                               end do
                            end do
                         end do
                         !
                         ! Assembly
                         !
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               !$acc atomic update
                               Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+1.0d0*Re(inode,idime)
                               !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mom_diffusion

              subroutine ener_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,Tem,mu_fluid,mu_e,mu_sgs,Rener)

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
                      real(8),    intent(in)  :: mu_fluid(npoin)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode), kappa_e, mu_fgp, aux, divU, tauU(ndime), twoThirds
                      real(8)                 :: gpcar(ndime,nnode), gradU(ndime,ndime), gradT(ndime)

                      call nvtxStartRange("Energy diffusion")
                      twoThirds = 2.0d0/3.0d0
                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang  private(Re,gpcar,gradU,gradT,tauU)  vector_length(vecLength)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !
                            ! Compute viscosity and conductivity at Gauss point
                            !
                            mu_fgp = 0.0d0
                            !$acc loop vector reduction(+:mu_fgp)
                            do inode = 1,nnode
                               mu_fgp = mu_fgp+(Ngp(igaus,inode)*mu_fluid(connec(ielem,inode)))
                            end do
                            kappa_e =mu_fgp*1004.0d0/0.71d0+c_ener*mu_e(ielem,igaus)/0.4d0 + mu_sgs(ielem,igaus)/0.9d0
                            mu_fgp = mu_fgp+mu_e(ielem,igaus)
                            !
                            ! Compute grad(U)
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  aux = 0.0d0
                                  !$acc loop vector reduction(+:aux)
                                  do inode = 1,nnode
                                     aux = aux+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                                  end do
                                  gradU(idime,jdime) = aux
                               end do
                            end do
                            !
                            ! Compute div(u)
                            divU = 0.0d0
                            !$acc loop seq
                            do idime = 1,ndime
                               divU = divU+gradU(idime,idime)
                            end do

                            ! TODO: Remove this lines
                            !divU = 0.0d0
                            !!$acc loop vector collapse(2) reduction(+:divU)
                            !do idime = 1,ndime
                            !   do inode = 1,nnode
                            !      divU = divU+gpcar(idime,inode)*u(connec(ielem,inode),idime)
                            !   end do
                            !end do

                            !
                            ! Compute tau*u
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               tauU(idime) = 0.0d0
                               !$acc loop seq
                               do jdime = 1,ndime
                                  aux = 0.0d0
                                  !$acc loop vector reduction(+:aux)
                                  do inode = 1,nnode
                                     aux = aux+Ngp(igaus,inode)*u(connec(ielem,inode),jdime)
                                  end do
                                  tauU(idime) = tauU(idime) + &
                                     gradU(idime,jdime)*aux + gradU(jdime,idime)*aux
                               end do
                               aux = 0.0d0
                               !$acc loop vector reduction(+:aux)
                               do inode = 1,nnode
                                  aux = aux+Ngp(igaus,inode)*u(connec(ielem,inode),idime)
                               end do
                               tauU(idime) = tauU(idime)-twoThirds*divU*aux
                            end do
                            !
                            ! Compute grad(T)
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               aux = 0.0d0
                               !$acc loop vector reduction(+:aux)
                               do inode = 1,nnode
                                  aux = aux+gpcar(idime,inode)*Tem(connec(ielem,inode))
                               end do
                               gradT(idime) = aux
                            end do
                            !
                            ! Gaussian integ
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*gpcar(idime,inode)* &
                                     (mu_fgp*tauU(idime)+kappa_e*gradT(idime))
                               end do
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+1.0d0*Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine ener_diffusion

              subroutine full_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Ml,Rmass,Rmom,Rener)
                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: rho(npoin),u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus),Ml(npoin)
                      real(8),    intent(in)  :: mu_fluid(npoin)
                      real(8),    intent(out) :: Rmass(npoin)
                      real(8),    intent(out) :: Rmom(npoin,ndime)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re_mass(nnode),Re_mom(nnode,ndime),Re_ener(nnode)
                      real(8)                 :: kappa_e, mu_fgp, mu_egp,aux, divU, tauU(ndime), twoThirds,nu_e,tau(ndime,ndime)
                      real(8)                 :: gpcar(ndime,nnode), gradU(ndime,ndime), gradT(ndime),tmp1,ugp(ndime),vol,arho

                      call nvtxStartRange("Full diffusion")
                      twoThirds = 2.0d0/3.0d0
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      Rmom(:,:) = 0.0d0
                      Rener(:) = 0.0d0
                      !$acc end kernels

                      !$acc parallel loop gang  private(gpcar,Re_mass,Re_mom,tau,gradU,Re_ener,gradT,tauU,ugp) vector_length(vecLength)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re_mass(inode) = 0.0d0
                            Re_ener(inode) = 0.0d0
                         end do
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               Re_mom(inode,idime) = 0.0d0
                            end do
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                               ugp(idime) = u(connec(ielem,igaus),idime)
                            end do
                            nu_e = c_rho*mu_e(ielem,igaus)/rho(connec(ielem,igaus))
                            mu_fgp = mu_fluid(connec(ielem,igaus))+rho(connec(ielem,igaus))*mu_sgs(ielem,igaus)
                            mu_egp = mu_e(ielem,igaus)
                            kappa_e =mu_fluid(connec(ielem,igaus))*1.6d0/0.71d0+c_ener*mu_e(ielem,igaus)/0.4d0 + rho(connec(ielem,igaus))*mu_sgs(ielem,igaus)/0.9d0
                            !kappa_e =mu_fluid(connec(ielem,igaus))*1004.0d0/0.71d0+c_ener*mu_e(ielem,igaus)/0.4d0 + rho(connec(ielem,igaus))*mu_sgs(ielem,igaus)/0.9d0
                            !
                            ! Compute grad(u)
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  aux = 0.0d0
                                  !$acc loop vector reduction(+:aux)
                                  do inode = 1,nnode
                                     aux = aux+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                                  end do
                                  gradU(idime,jdime) = aux
                               end do
                            end do
                            !
                            ! Compute div(u)
                            !
                            divU = gradU(1,1)+gradU(2,2)+gradU(3,3)
                            !
                            ! Finish computing tau_ij = u_i,j + u_j,i - (2/3)*div(u)*d_ij
                            ! Compute tau*u
                            !
                            !$acc loop seq
                            do idime = 1,ndime
                               tauU(idime) = 0.0d0
                               !$acc loop seq
                               do jdime = 1,ndime
                                  tauU(idime) = tauU(idime) + &
                                     (mu_fgp+mu_egp)*(gradU(idime,jdime)+ gradU(jdime,idime))*ugp(jdime)
                                  tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
                               end do
                               tauU(idime) = tauU(idime)-(mu_fgp)*twoThirds*divU*ugp(idime)
                               tau(idime,idime) = tau(idime,idime)-(mu_fgp)*twoThirds*divU
                            end do

                            ! Dif rho
                            ! Dif T
                            ! Compute div(tau) at the Gauss point
                            !$acc loop seq
                            do idime = 1,ndime
                               tmp1 = 0.0d0
                               aux = 0.0d0
                               !$acc loop vector reduction(+:tmp1,aux)
                               do inode = 1,nnode
                                  tmp1 = tmp1+gpcar(idime,inode)*rho(connec(ielem,inode))
                                  aux = aux+gpcar(idime,inode)*Tem(connec(ielem,inode))
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re_mass(inode) = Re_mass(inode)+nu_e*gpvol(1,igaus,ielem)* &
                                     gpcar(idime,inode)*tmp1
                                  Re_ener(inode) = Re_ener(inode)+gpvol(1,igaus,ielem)*gpcar(idime,inode)* &
                                                   (tauU(idime)+kappa_e*aux)
                               end do
                               !$acc loop seq
                               do jdime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     Re_mom(inode,idime) = Re_mom(inode,idime)+gpvol(1,igaus,ielem)* &
                                        gpcar(jdime,inode)*tau(idime,jdime)
                                  end do
                               end do
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re_mass(inode)
                            !$acc end atomic
                            !$acc atomic update
                            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re_ener(inode)
                            !$acc end atomic
                         end do
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               !$acc atomic update
                               Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                               !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine full_diffusion

end module elem_diffu
