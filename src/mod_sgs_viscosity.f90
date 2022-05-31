module mod_sgs_viscosity

   use mod_constants
   use mod_nvtx
   use mod_veclen

      ! TODO: Finish module and create unit tests

      contains
               
               ! it implents the Vreman SGS

              subroutine sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol, &
                                    rho,u,Ml,mu_sgs)
              
                      implicit none

                      integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: rho(npoin), u(npoin,ndime),Ml(npoin)
                      real(8),    intent(out) :: mu_sgs(nelem,ngaus)
                      integer(4)              :: ielem, inode, igaus, kdime, idime, jdime
                      real(8)                 :: hLES, evol,gpcar(ndime,nnode), aux, mue(npoin), ave(npoin)
                      real(8)                 :: gradU(ndime,ndime),gradV2(ndime,ndime),Bbeta, alpha, tau


                      if (flag_SpectralElem == 0) then

                         !$acc parallel loop gang  private(gpcar,gradU,gradV2) vector_length(vecLength)
                         do ielem = 1,nelem
                            evol = 0.0d0
                            !$acc loop vector reduction(+:evol)
                            do igaus = 1,ngaus
                               evol = evol + gpvol(1,igaus,ielem)
                            end do
                            hLES = (evol**(1.0d0/3.0d0))/dble(porder)
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
                               ! Compute SGS
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
                               alpha = 0.0d0
                               !$acc loop vector collapse(2)
                               do idime = 1, ndime
                                  do jdime = 1, ndime
                                     gradV2(idime,jdime) = 0.0d0
                                  end do
                               end do

                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop seq
                                  do jdime = 1,ndime
                                     !$acc loop seq
                                     do kdime = 1,ndime
                                        gradV2(idime,jdime) = gradV2(idime,jdime) + gradU(kdime,idime)*gradU(kdime,jdime)*hLES*hLES
                                     end do
                                     alpha = alpha + gradU(idime,jdime)*gradU(idime,jdime)
                                  end do
                               end do
                               Bbeta = gradV2(1,1)*gradV2(2,2) + gradV2(1,1)*gradV2(3,3) + gradV2(2,2)*gradV2(3,3) &
                                  - gradV2(1,2)*gradV2(1,2) - gradV2(1,3)*gradV2(1,3) - gradV2(2,3)*gradV2(2,3)
                               if(alpha > 1.0e-10) then
                                  mu_sgs(ielem,igaus) = c_sgs*sqrt(max(Bbeta,1.0e-10)/alpha)
                               else
                                  mu_sgs(ielem,igaus) = 0.0d0
                               end if
                            end do
                         end do
                         !$acc end parallel loop

                      else
                         !!$acc kernels
                         !mue(:) = 0.0d0
                         !ave(:) = 0.0d0
                         !!$acc end kernels
                         !$acc parallel loop gang  private(gpcar,gradU,gradV2) vector_length(vecLength)
                         do ielem = 1,nelem
                            evol = 0.0d0
                            !$acc loop vector reduction(+:evol)
                            do igaus = 1,ngaus
                               evol = evol + gpvol(1,igaus,ielem)
                            end do
                            hLES = (evol**(1.0d0/3.0d0))/dble(porder+1)
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
                               ! Compute SGS
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
                               alpha = 0.0d0
                               !$acc loop vector collapse(2)
                               do idime = 1, ndime
                                  do jdime = 1, ndime
                                     gradV2(idime,jdime) = 0.0d0
                                  end do
                               end do

                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop seq
                                  do jdime = 1,ndime
                                     !$acc loop seq
                                     do kdime = 1,ndime
                                        gradV2(idime,jdime) = gradV2(idime,jdime) + gradU(kdime,idime)*gradU(kdime,jdime)*hLES*hLES
                                     end do
                                     alpha = alpha + gradU(idime,jdime)*gradU(idime,jdime)
                                  end do
                               end do
                               Bbeta = gradV2(1,1)*gradV2(2,2) + gradV2(1,1)*gradV2(3,3) + gradV2(2,2)*gradV2(3,3) &
                                  - gradV2(1,2)*gradV2(1,2) - gradV2(1,3)*gradV2(1,3) - gradV2(2,3)*gradV2(2,3)
                               aux = 0.0d0
                               if(alpha > 1.0e-10) then
                                  aux = c_sgs*sqrt(max(Bbeta,1.0e-10)/alpha)
                               end if
                               mu_sgs(ielem,igaus) = aux
                              !!$acc atomic update
                              !mue(connec(ielem,igaus)) = mue(connec(ielem,igaus))+aux
                              !!$acc end atomic
                              !!$acc atomic update
                              !ave(connec(ielem,igaus)) = ave(connec(ielem,igaus))+1.0d0
                              !!$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop

                       ! !$acc parallel loop gang vector_length(vecLength)
                       ! do ielem = 1,nelem
                       !    !$acc loop vector
                       !    do inode = 1,nnode
                       !       mu_sgs(ielem,inode) = mue(connec(ielem,inode))/ave(connec(ielem,inode))
                       !    end do
                       ! end do
                       ! !$acc end parallel loop

                      end if

                   end subroutine sgs_visc
                end module mod_sgs_viscosity
