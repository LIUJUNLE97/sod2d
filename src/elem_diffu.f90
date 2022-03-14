module elem_diffu

      use mod_nvtx
      use mod_constants

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
                      real(8),    intent(in)  :: rho(npoin), mu_e(nelem)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode), nu_e
                      real(8)                 :: tmp1, gpcar(ndime,nnode)

                      call nvtxStartRange("Mass diffusion")
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(gpcar,Re) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         nu_e = mu_e(ielem)/maxval(abs(rho(connec(ielem,:))))
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
                               !$acc loop vector reduction(+:tmp1)
                               do inode = 1,nnode
                                  tmp1 = tmp1+gpcar(idime,inode)*rho(connec(ielem,inode))
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                              gpcar(idime,inode)*tmp1
                               end do
                            end do
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            !$acc atomic update
                            Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+nu_e*1.0d0*Re(inode)
                            !$acc end atomic
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine mass_diffusion

              subroutine mom_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,mu_fluid,mu_e,Rmom)

                      ! TODO: Add. stab. viscosity

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), mu_e(nelem)
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
                      !$acc parallel loop gang private(Re,gpcar,tau,gradU) vector_length(32)
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
                            mu_fgp = mu_fgp+mu_e(ielem)
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
                            !$acc loop vector collapse(2) reduction(+:divU)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  divU = divU+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                               end do
                            end do
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
                            do idime = 1,ndime
                               tau(idime,idime) = tau(idime,idime)-twoThirds*divU
                            end do
                            !
                            ! Compute div(tau) at the Gauss point
                            !
                            !$acc loop vector collapse(3)
                            do idime = 1,ndime
                               do jdime = 1,ndime
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

              subroutine ener_diffusion(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,Tem,mu_fluid,mu_e,Rener)

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), Tem(npoin), mu_e(nelem)
                      real(8),    intent(in)  :: mu_fluid(npoin)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode), kappa_e, mu_fgp
                      real(8)                 :: el_Ke(nnode)
                      real(8)                 :: gradT, gradKe, gpcar(ndime,nnode)

                      call nvtxStartRange("Energy diffusion")
                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      !$acc parallel loop gang private(el_Ke,Re,gpcar) vector_length(32)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            Re(inode) = 0.0d0
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            el_Ke(inode) = dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))/2.0d0
                         end do
                         !$acc loop seq
                         do igaus = 1,ngaus
                            mu_fgp = dot_product(Ngp(igaus,:),mu_fluid(connec(ielem,:)))
                            kappa_e = (mu_fgp+mu_e(ielem))*1004.0d0/0.71d0
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop vector
                               do inode = 1,nnode
                                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                               end do
                            end do
                            !$acc loop seq
                            do idime = 1,ndime
                               gradT = kappa_e*dot_product(gpcar(idime,:),Tem(connec(ielem,:)))
                               gradKe = (mu_fgp+mu_e(ielem))*dot_product(gpcar(idime,:),el_Ke(:))
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                     gpcar(idime,inode)*(gradT+gradKe)
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

end module elem_diffu
