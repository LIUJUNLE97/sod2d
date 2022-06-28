module elem_convec

      use mod_nvtx
      use mod_constants
      use mod_veclen
      use mod_maths

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Computes convective term for Euler/NS equation system, as well             !
      ! as for any generic scalar transport that might occur. Based                !
      ! on Ljunkvist matrix-free implementation (assembles only rhs vector).       !
      ! This module can be passed to CUDA in order to do fine-grained parallelism. !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

              subroutine mass_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,q,rho,u,Rmass)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(rho*u) using a standard Continuous  !
                      ! Galerkin formulation. In the above, rho is the scalar sensity and !
                      ! u is the vector of velocities in all dimensions.                  !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(in)  :: rho(npoin)
                      real(8),    intent(in)  :: u(npoin,ndime)
                      real(8),    intent(out) :: Rmass(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: Re(nnode)
                      real(8)                 :: tmp1,gpcar(ndime,nnode),tmp2,tmp3,tmp4

                      call nvtxStartRange("Mass Convection")
                      !$acc kernels
                      Rmass(:) = 0.0d0
                      !$acc end kernels

                      if (flag_SpectralElem == 0) then
                         !$acc parallel loop gang  private(Re,gpcar) vector_length(vecLength)
                         do ielem = 1,nelem
                            !$acc loop vector
                            do inode = 1,nnode
                               Re(inode) = 0.0d0
                            end do
                            !
                            ! Quadrature
                            !
                            !$acc loop seq
                            do igaus = 1,ngaus
                               tmp1 = 0.0d0
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                                  end do
                               end do
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop vector reduction(+:tmp1)
                                  do inode = 1,nnode
                                     tmp1 = tmp1+(gpcar(idime,inode)*q(connec(ielem,inode),idime))
                                  end do
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                     Ngp(igaus,inode)*tmp1
                               end do
                            end do
                            !
                            ! Assembly
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               !$acc atomic update
                               Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re(inode)
                               !$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop
                      else
                         !$acc parallel loop gang  private(Re,gpcar) vector_length(vecLength)
                         do ielem = 1,nelem
                            !
                            ! Quadrature
                            !
                            !$acc loop seq
                            do igaus = 1,ngaus
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                                  end do
                               end do
                               tmp1 = 0.0d0
                               tmp2 = 0.0d0
                               tmp4 = 0.0d0
                               !$acc loop seq
                               do idime = 1,ndime
                                  tmp3 = 0.0d0
                                  !$acc loop vector reduction(+:tmp1,tmp2,tmp3)
                                  do inode = 1,nnode
                                     tmp1 = tmp1+(gpcar(idime,inode)*q(connec(ielem,inode),idime))
                                     tmp2 = tmp2+(gpcar(idime,inode)*u(connec(ielem,inode),idime))
                                     tmp3 = tmp3+(gpcar(idime,inode)*rho(connec(ielem,inode)))
                                  end do
                                  tmp4 = tmp4 + u(connec(ielem,igaus),idime)*tmp3
                               end do
                               !Re(igaus) = gpvol(1,igaus,ielem)*tmp1
                               Re(igaus) = gpvol(1,igaus,ielem)*0.5d0*(tmp1+tmp2*rho(connec(ielem,igaus))+tmp4)
                            end do
                            !
                            ! Assembly
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               !$acc atomic update
                               Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re(inode)
                               !$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop
                      end if
                      call nvtxEndRange

              end subroutine mass_convec

              subroutine mom_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,q,rho,pr,Rmom)

                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Subroutine to compute R = div(q*u) using a standard Continuous   !
                 ! Galerkin formulation. In the above, q is the momentum vector and !
                 ! u is the vector of velocities in all dimensions. The product     !
                 ! inside div() is an outer tensor product between the 2 vectors.   !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                 implicit none

                 integer(4), intent(in)  :: nelem, npoin
                 integer(4), intent(in)  :: connec(nelem,nnode)
                 real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                 real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                 real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                 real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin)
                 real(8),    intent(out) :: Rmom(npoin,ndime)
                 integer(4)              :: ielem, igaus, idime, jdime, inode
                 real(8)                 :: Re(nnode,ndime), aux, divU,gradQ(ndime, ndime)
                 real(8)                 :: tmp1(ndime), tmp2(ndime), gpcar(ndime,nnode)
                 real(8)                 :: aux2,aux3,aux4,aux5

                 call nvtxStartRange("Momentum convection")
                 !$acc kernels
                 Rmom(:,:) = 0.0d0
                 !$acc end kernels

                 if (flag_SpectralElem == 0) then
                    !$acc parallel loop gang  private(Re,gpcar,tmp1,tmp2,gradQ,divU) vector_length(vecLength)
                    do ielem = 1,nelem
                       !$acc loop vector collapse(2)
                       do inode = 1,nnode
                          do idime = 1,ndime
                             Re(inode,idime) = 0.0d0
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
                          end do
                          !$acc loop seq
                          do idime = 1,ndime
                             aux = 0.0d0
                             !$acc loop seq
                             do jdime = 1,ndime
                                !$acc loop vector reduction(+:aux)
                                do inode = 1,nnode
                                   aux = aux+gpcar(jdime,inode)*(q(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime))
                                end do
                             end do
                             tmp1(idime) = aux
                          end do
                          !$acc loop seq
                          do idime = 1,ndime
                             aux = 0.0d0
                             !$acc loop vector reduction(+:aux)
                             do inode = 1,nnode
                                aux = aux+gpcar(idime,inode)*pr(connec(ielem,inode))
                             end do
                             tmp2(idime) = aux
                          end do
                          !$acc loop vector
                          do inode = 1,nnode
                             Re(inode,1) = Re(inode,1)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*(tmp1(1)+tmp2(1))
                             Re(inode,2) = Re(inode,2)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*(tmp1(2)+tmp2(2))
                             Re(inode,3) = Re(inode,3)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*(tmp1(3)+tmp2(3))
                          end do
                       end do
                       !
                       ! Final assembly
                       !
                       !$acc loop vector collapse(2)
                       do idime = 1,ndime
                          do inode = 1,nnode
                             !$acc atomic update
                             Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re(inode,idime)
                             !$acc end atomic
                          end do
                       end do
                    end do
                    !$acc end parallel loop
                 else
                    !$acc parallel loop gang  private(Re,gpcar,tmp1,tmp2,gradQ,divU) vector_length(vecLength)
                    do ielem = 1,nelem
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
                             aux  = 0.0d0
                             aux3  = 0.0d0
                             aux4  = 0.0d0
                             aux5  = 0.0d0
                             !$acc loop seq
                             do jdime = 1,ndime
                                aux2  = 0.0d0
                                !$acc loop vector reduction(+:aux,aux2,aux3)
                                do inode = 1,nnode
                                   aux = aux  +gpcar(jdime,inode)*(q(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime))
                                   aux2 = aux2+gpcar(jdime,inode)*u(connec(ielem,inode),idime)
                                   aux3 = aux3+gpcar(jdime,inode)*u(connec(ielem,inode),jdime)
                                   aux5 = aux5+u(connec(ielem,igaus),jdime)*gpcar(jdime,inode)*rho(connec(ielem,inode))
                                end do
                                aux4 = aux4+q(connec(ielem,igaus),jdime)*aux2
                             end do
                             tmp1(idime) = 0.5d0*(aux+q(connec(ielem,igaus),idime)*aux3+aux4+u(connec(ielem,igaus),idime)*aux5)
                          end do
                          !$acc loop seq
                          do idime = 1,ndime
                             aux = 0.0d0
                             !$acc loop vector reduction(+:aux)
                             do inode = 1,nnode
                                aux = aux+gpcar(idime,inode)*pr(connec(ielem,inode))
                             end do
                             tmp2(idime) = aux
                          end do
                          Re(igaus,1) = gpvol(1,igaus,ielem)*(tmp1(1)+tmp2(1))
                          Re(igaus,2) = gpvol(1,igaus,ielem)*(tmp1(2)+tmp2(2))
                          Re(igaus,3) = gpvol(1,igaus,ielem)*(tmp1(3)+tmp2(3))
                       end do
                       !
                       ! Final assembly
                       !
                       !$acc loop vector collapse(2)
                       do idime = 1,ndime
                          do inode = 1,nnode
                             !$acc atomic update
                             Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re(inode,idime)
                             !$acc end atomic
                          end do
                       end do
                    end do
                    !$acc end parallel loop
                 end if

                 call nvtxEndRange

              end subroutine mom_convec

              subroutine mom_convec_emac(nelem,npoin,connec,Ngp,dNgp,He,gpvol,Aemac,Femac,pr,Rmom)

                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! Subroutine to compute R = EMAC(A) using a standard Continuous    !
                 ! Galerkin model. A_i is q_i/sqrt(rho), and the EMAC term is       !
                 ! defined as:                                                      !
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                 implicit none

                 integer(4), intent(in)  :: nelem,npoin
                 integer(4), intent(in)  :: connec(nelem,nnode)
                 real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                 real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                 real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                 real(8),    intent(in)  :: Aemac(npoin,ndime), Femac(npoin), pr(npoin)
                 real(8),    intent(out) :: Rmom(npoin,ndime)
                 integer(4)              :: ielem, igaus, idime, jdime, inode,jnode, kdime
                 real(8)                 :: Re(nnode,ndime), tmp(ndime),gpA(ndime),aux2
                 real(8)                 :: gpcar(ndime,nnode), gradA(ndime,ndime), divA, aux

                 call nvtxStartRange("EMAC Momentum convection")

                 !
                 ! Start global RHS
                 !
                 !$acc kernels
                 Rmom(:,:) = 0.0d0
                 !$acc end kernels

                 if (flag_SpectralElem == 0) then
                    ! Start elemental ops
                    !
                    !$acc parallel loop gang  private(Re,gpcar,gradA,tmp,gpA) vector_length(vecLength)
                    do ielem = 1,nelem
                       !
                       ! Initialize element vector to 0
                       !
                       !$acc loop vector collapse(2)
                       do idime = 1,ndime
                          do inode = 1,nnode
                             Re(inode,idime) = 0.0d0
                          end do
                       end do
                       !
                       ! Loop over Gauss points
                       !
                       !$acc loop seq
                       do igaus = 1,ngaus
                          !
                          ! Create GPCAR(ndime,nnode) for each element at each Gauss point
                          !
                          !$acc loop seq
                          do idime = 1,ndime
                             !$acc loop vector
                             do inode = 1,nnode
                                gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                             end do
                          end do
                          !
                          ! grad(A) gpA
                          !
                          !$acc loop seq
                          do idime = 1,ndime
                             !$acc loop seq
                             do jdime = 1,ndime
                                aux = 0.0d0
                                !$acc loop vector reduction(+:aux)
                                do inode = 1,nnode
                                   aux = aux+gpcar(jdime,inode)*Aemac(connec(ielem,inode),idime)
                                end do
                                gradA(idime,jdime) = aux
                             end do
                             aux = 0.0d0
                             !$acc loop vector reduction(+:aux)
                             do inode = 1,nnode
                                aux = aux+Ngp(igaus,inode)*Aemac(connec(ielem,inode),idime)
                             end do
                             gpA(idime) = aux
                          end do
                          !
                          ! div(A)
                          !

                          divA = gradA(1,1)+gradA(2,2)+gradA(3,3)

                          !
                          ! grad(A)*A + gradT(A)*A + div(A)*A
                          !
                               !$acc loop seq
                               do idime = 1,ndime
                                  tmp(idime) = 0.0d0
                                  !$acc loop seq
                                  do jdime = 1,ndime
                                     tmp(idime) = tmp(idime) + (gradA(idime,jdime)+gradA(jdime,idime))*gpA(jdime)
                                  end do
                                  tmp(idime) = tmp(idime)+divA*gpA(idime)
                               end do
                               !
                               ! Subtract -0.5*grad(F), where F = (A.A)
                               !
                                !$acc loop seq
                                do idime = 1,ndime
                                   aux = 0.0d0
                                   !$acc loop vector reduction(+:aux)
                                   do inode = 1,nnode
                                      aux = aux+gpcar(idime,inode)*Femac(connec(ielem,inode))
                                   end do
                                   tmp(idime) = tmp(idime)-0.5d0*aux
                                end do

                               !
                               ! Addd pressure
                               !
                               !$acc loop seq
                               do idime = 1,ndime
                                  aux = 0.0d0
                                  !$acc loop vector reduction(+:aux)
                                  do inode = 1,nnode
                                     aux = aux+gpcar(idime,inode)*pr(connec(ielem,inode))
                                  end do
                                  tmp(idime) = tmp(idime) + aux
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode,1) = Re(inode,1)+gpvol(1,igaus,ielem)*(Ngp(igaus,inode)*tmp(1))
                                  Re(inode,2) = Re(inode,2)+gpvol(1,igaus,ielem)*(Ngp(igaus,inode)*tmp(2))
                                  Re(inode,3) = Re(inode,3)+gpvol(1,igaus,ielem)*(Ngp(igaus,inode)*tmp(3))
                               end do
                            end do
                            !
                            ! Final assembly
                            !
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  !$acc atomic update
                                  Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re(inode,idime)
                                  !$acc end atomic
                               end do
                            end do
                         end do
                         !$acc end parallel loop
                      else !spectral
                         ! Start elemental ops
                         !
                         !$acc parallel loop gang  private(Re,gpcar,gradA,tmp,gpA) vector_length(vecLength)
                         do ielem = 1,nelem
                            !
                            ! Loop over Gauss points
                            !
                            !$acc loop seq
                            do igaus = 1,ngaus
                               !
                               ! Create GPCAR(ndime,nnode) for each element at each Gauss point
                               !
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                                  end do
                               end do
                               !
                               ! grad(A) gpA
                               !
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop seq
                                  do jdime = 1,ndime
                                     aux = 0.0d0
                                     !$acc loop vector reduction(+:aux)
                                     do inode = 1,nnode
                                        aux = aux+gpcar(jdime,inode)*Aemac(connec(ielem,inode),idime)
                                     end do
                                     gradA(idime,jdime) = aux
                                  end do
                                  gpA(idime) = Aemac(connec(ielem,igaus),idime)
                               end do
                               !
                               ! div(A)
                               !

                               divA = gradA(1,1)+gradA(2,2)+gradA(3,3)

                               !
                               ! grad(A)*A + gradT(A)*A + div(A)*A
                               !
                               !$acc loop seq
                               do idime = 1,ndime
                                  tmp(idime) = 0.0d0
                                  !$acc loop seq
                                  do jdime = 1,ndime
                                     tmp(idime) = tmp(idime) + (gradA(idime,jdime)+gradA(jdime,idime))*gpA(jdime)
                                  end do
                                  tmp(idime) = tmp(idime)+divA*gpA(idime)
                               end do

                               !
                               ! Addd pressure
                               !
                               !$acc loop seq
                               do idime = 1,ndime
                                  aux = 0.0d0
                                  !$acc loop vector reduction(+:aux)
                                  do inode = 1,nnode
                                     aux = aux+gpcar(idime,inode)*(pr(connec(ielem,inode))-0.5d0*Femac(connec(ielem,inode)))
                                  end do
                                  tmp(idime) = tmp(idime) + aux
                               end do
                               Re(igaus,1) = gpvol(1,igaus,ielem)*tmp(1)
                               Re(igaus,2) = gpvol(1,igaus,ielem)*tmp(2)
                               Re(igaus,3) = gpvol(1,igaus,ielem)*tmp(3)
                            end do
                            !
                            ! Final assembly
                            !
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  !$acc atomic update
                                  Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re(inode,idime)
                                  !$acc end atomic
                               end do
                            end do
                         end do
                         !$acc end parallel loop

                      end if

                      call nvtxEndRange

                   end subroutine mom_convec_emac

                   subroutine ener_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,pr,E,rho,q,Rener)

                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ! Subroutine to compute R = div(u*(E+p)) using a standard Continuous !
                      ! Galerkin formulation. In the above, E is the scalar total energy,  !
                      ! p is the scalar pressure and u is the vector of velocities in all  !
                      ! dimensions.                                                        !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: u(npoin,ndime), pr(npoin), E(npoin), rho(npoin),q(npoin,ndime)
                      real(8),    intent(out) :: Rener(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, ipoin, jdime
                      real(8)                 :: Re(nnode), aux
                      real(8)                 :: tmp1,gpcar(ndime,nnode),tmp2,tmp3,tmp4,tmp5

                      call nvtxStartRange("Energy Convection")
                      !$acc kernels
                      Rener(:) = 0.0d0
                      !$acc end kernels
                      if (flag_SpectralElem == 0) then
                         !$acc parallel loop gang  private(Re,gpcar) vector_length(vecLength)
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
                               tmp1 = 0.0d0
                               !$acc loop vector collapse(2) reduction(+:tmp1)
                               do idime = 1,ndime
                                  do inode = 1,nnode
                                     tmp1 = tmp1+(gpcar(idime,inode)* &
                                        (u(connec(ielem,inode),idime)*(E(connec(ielem,inode))+pr(connec(ielem,inode)))))
                                  end do
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)*Ngp(igaus,inode)*tmp1
                               end do
                            end do
                            !
                            ! Assembly
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               !$acc atomic update
                               Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re(inode)
                               !$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop
                      else

                         !$acc parallel loop gang  private(Re,gpcar) vector_length(vecLength)
                         do ielem = 1,nelem
                            !$acc loop seq
                            do igaus = 1,ngaus
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop vector
                                  do inode = 1,nnode
                                     gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                                  end do
                               end do
                               tmp1 = 0.0d0
                               tmp2 = 0.0d0
                               tmp4 = 0.0d0
                               tmp5 = 0.0d0
                               !$acc loop seq
                               do idime = 1,ndime
                                  tmp3 = 0.0d0
                                  !$acc loop vector reduction(+:tmp1,tmp5,tmp2,tmp3)
                                  do inode = 1,nnode
                                     !tmp1 = tmp1+gpcar(idime,inode)*u(connec(ielem,inode),idime)*(E(connec(ielem,inode))+pr(connec(ielem,inode)))
                                     tmp1 = tmp1+gpcar(idime,inode)*u(connec(ielem,inode),idime)*E(connec(ielem,inode))
                                     tmp5 = tmp5+gpcar(idime,inode)*u(connec(ielem,inode),idime)*pr(connec(ielem,inode))
                                     tmp2 = tmp2+gpcar(idime,inode)*u(connec(ielem,inode),idime)
                                     tmp3 = tmp3+gpcar(idime,inode)*E(connec(ielem,inode))
                                  end do
                                  tmp4 = tmp4 + u(connec(ielem,igaus),idime)*tmp3
                               end do
                               Re(igaus) = gpvol(1,igaus,ielem)*(0.5d0*(tmp1+E(connec(ielem,igaus))*tmp2+tmp4)+tmp5)
                               !Re(igaus) = gpvol(1,igaus,ielem)*tmp1
                            end do
                            !
                            ! Assembly
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               !$acc atomic update
                               Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re(inode)
                               !$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop

                      end if
                      call nvtxEndRange

                   end subroutine ener_convec

                   subroutine generic_scalar_convec(nelem,npoin,connec,Ngp, &
                         dNgp,He,gpvol,q,eta,u,Rconvec,alpha)

                      implicit none

                      integer(4), intent(in)  :: nelem, npoin
                      integer(4), intent(in)  :: connec(nelem,nnode)
                      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                      real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                      real(8),    intent(in)  :: q(npoin,ndime)
                      real(8),    intent(in)  :: eta(npoin)
                      real(8),    intent(in)  :: u(npoin,ndime)
                      real(8),    intent(in)  :: alpha(npoin)
                      real(8),    intent(out) :: Rconvec(npoin)
                      integer(4)              :: ielem, igaus, inode, idime, jdime
                      real(8)                 :: tmp1, tmp2, Re(nnode), gpcar(ndime,nnode),tmp3,tmp4,tmp5
                      real(8)                 :: ul(nnode,ndime), ql(nnode,ndime), etal(nnode)

                      call nvtxStartRange("Generic Convection")
                      !$acc kernels
                      Rconvec(:) = 0.0d0
                      !$acc end kernels
                      if (flag_SpectralElem == 0) then
                         !$acc parallel loop gang  private(Re,gpcar) vector_length(vecLength)
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
                               tmp1 = 0.0d0
                               !$acc loop vector reduction(+:tmp1)
                               do inode = 1,nnode
                                  tmp1 = tmp1+(Ngp(igaus,inode)*alpha(connec(ielem,inode)))
                               end do
                               tmp2 = 0.0d0
                               !$acc loop vector collapse(2) reduction(+:tmp2)
                               do idime = 1,ndime
                                  do inode = 1,nnode
                                     tmp2 = tmp2+(gpcar(idime,inode)*q(connec(ielem,inode),idime))
                                  end do
                               end do
                               !$acc loop vector
                               do inode = 1,nnode
                                  Re(inode) = Re(inode)+gpvol(1,igaus,ielem)* &
                                     Ngp(igaus,inode)*tmp1*tmp2
                               end do
                            end do
                            !
                            ! Assembly
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               !$acc atomic update
                               Rconvec(connec(ielem,inode)) = Rconvec(connec(ielem,inode))+Re(inode)
                               !$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop
                      else
                         !$acc parallel loop gang  private(Re,ul,ql,etal) !!vector_length(vecLength)
                         do ielem = 1,nelem
                            !$acc loop vector collapse(2)
                            do idime = 1,ndime
                               do inode = 1,nnode
                                  ul(inode,idime) = u(connec(ielem,inode),idime)
                                  ql(inode,idime) = q(connec(ielem,inode),idime)
                               end do
                            end do
                            !$acc loop vector
                            do inode = 1,nnode
                               etal(inode) = eta(connec(ielem,inode))
                            end do
                            !$acc loop worker private(gpcar)
                            do igaus = 1,ngaus
                               !$acc loop vector collapse(2)
                               do idime = 1,ndime
                                  do inode = 1,nnode
                                     tmp1 = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                                     gpcar(idime,inode) = tmp1
                                  end do
                               end do
                               tmp1 = 0.0d0
                               tmp2 = 0.0d0
                               tmp4 = 0.0d0
                               !$acc loop seq
                               do idime = 1,ndime
                                  !$acc loop vector reduction(+:tmp1)
                                  do inode = 1,nnode
                                     tmp1 = tmp1+gpcar(idime,inode)*ql(inode,idime)&
                                                 +etal(inode)*gpcar(idime,inode)*ul(inode,idime)
                                  end do
                               end do
                               !$acc loop seq
                               do idime = 1,ndime
                                  tmp3 = 0.0d0
                                  !$acc loop vector reduction(+:tmp3)
                                  do inode = 1,nnode
                                     tmp3 = tmp3+gpcar(idime,inode)*etal(inode)
                                  end do
                                  tmp4 = tmp4 + ul(igaus,idime)*tmp3
                               end do
                               Re(igaus) = gpvol(1,igaus,ielem)*0.5d0*(tmp1+tmp4)
                            end do
                            !
                            ! Assembly
                            !
                            !$acc loop vector
                            do inode = 1,nnode
                               !$acc atomic update
                               Rconvec(connec(ielem,inode)) = Rconvec(connec(ielem,inode))+Re(inode)
                               !$acc end atomic
                            end do
                         end do
                         !$acc end parallel loop
                      end if
                      call nvtxEndRange

              end subroutine generic_scalar_convec

              subroutine full_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,q,rho,pr,E,Rmass,Rmom,Rener)

                 implicit none

                 integer(4), intent(in)  :: nelem, npoin
                 integer(4), intent(in)  :: connec(nelem,nnode)
                 real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                 real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem)
                 real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                 real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
                 real(8),    intent(out) :: Rmass(npoin)
                 real(8),    intent(out) :: Rmom(npoin,ndime)
                 real(8),    intent(out) :: Rener(npoin)
                 integer(4)              :: ielem, igaus, idime, jdime, inode
                 real(8)                 :: Re_mom(nnode,ndime), aux_mom
                 real(8)                 :: Re_mass(nnode), Re_ener(nnode)
                 real(8)                 :: tmp1_mom(ndime), tmp2_mom(ndime), gpcar(ndime,nnode)
                 real(8)                 :: aux2_mom,aux3_mom,aux4_mom,aux5_mom
                 real(8)                 :: tmp1_mass,tmp2_mass,tmp3_mass,tmp4_mass
                 real(8)                 :: tmp1_ener,tmp2_ener,tmp3_ener,tmp4_ener, tmp5_ener
                 real(8)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode)


                 call nvtxStartRange("Full convection")
                 !$acc kernels
                 Rmom(:,:) = 0.0d0
                 Rmass(:) = 0.0d0
                 Rener(:) = 0.0d0
                 !$acc end kernels

                 !$acc parallel loop gang private(Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El) !!num_workers(4) vector_length(vecLength)
                 do ielem = 1,nelem
                    !$acc loop vector collapse(2)
                    do idime = 1,ndime
                       do inode = 1,nnode
                          ul(inode,idime) = u(connec(ielem,inode),idime)
                          ql(inode,idime) = q(connec(ielem,inode),idime)
                       end do
                    end do
                    !$acc loop vector
                    do inode = 1,nnode
                       rhol(inode) = rho(connec(ielem,inode))
                       El(inode) = E(connec(ielem,inode))
                       prl(inode) = pr(connec(ielem,inode))
                    end do
                    !$acc loop worker private(tmp1_mom,tmp2_mom,gpcar)
                    do igaus = 1,ngaus
                       !$acc loop vector collapse(2)
                       do idime = 1,ndime
                          do inode = 1,nnode
                             aux_mom =  dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                             gpcar(idime,inode) = aux_mom
                          end do
                       end do
                       tmp1_mass = 0.0d0
                       tmp4_mass = 0.0d0
                       tmp1_ener = 0.0d0
                       tmp4_ener = 0.0d0
                       tmp5_ener = 0.0d0
                       !$acc loop seq
                       do idime = 1,ndime
                          aux_mom  = 0.0d0
                          aux3_mom  = 0.0d0
                          aux4_mom  = 0.0d0
                          aux5_mom  = 0.0d0
                          !$acc loop vector collapse(2) reduction(+:aux_mom,aux3_mom)
                          do jdime = 1,ndime
                             do inode = 1,nnode
                                aux_mom = aux_mom+gpcar(jdime,inode)*(ql(inode,idime)*ul(inode,jdime)) &
                                                 +ul(igaus,idime)*ul(igaus,jdime)*gpcar(jdime,inode)*rhol(inode)
                                aux3_mom = aux3_mom+gpcar(jdime,inode)*ul(inode,jdime)
                             end do
                          end do
                          !$acc loop seq
                          do jdime = 1,ndime
                             aux2_mom  = 0.0d0
                             !$acc loop vector reduction(+:aux2_mom)
                             do inode = 1,nnode
                                aux2_mom = aux2_mom+gpcar(jdime,inode)*ul(inode,idime)
                             end do
                             aux4_mom = aux4_mom+ql(igaus,jdime)*aux2_mom
                          end do
                          tmp1_mom(idime) = 0.5d0*(aux_mom+ql(igaus,idime)*aux3_mom+aux4_mom)

                          aux_mom = 0.0d0
                          tmp3_mass = 0.0d0
                          tmp3_ener = 0.0d0
                          !$acc loop vector reduction(+:tmp1_mass,aux_mom,tmp3_mass,tmp1_ener,tmp5_ener,tmp3_ener)
                          do inode = 1,nnode
                             tmp1_mass = tmp1_mass+(gpcar(idime,inode)*ql(inode,idime))
                             tmp3_mass = tmp3_mass+(gpcar(idime,inode)*rhol(inode))
                             aux_mom = aux_mom+gpcar(idime,inode)*prl(inode)
                             tmp1_ener = tmp1_ener+gpcar(idime,inode)*ul(inode,idime)*El(inode)
                             tmp5_ener = tmp5_ener+gpcar(idime,inode)*ul(inode,idime)*prl(inode)
                             tmp3_ener = tmp3_ener+gpcar(idime,inode)*El(inode)
                          end do
                          tmp4_mass = tmp4_mass + ul(igaus,idime)*tmp3_mass
                          tmp4_ener = tmp4_ener + ul(igaus,idime)*tmp3_ener
                          tmp2_mom(idime) = aux_mom
                       end do
                       Re_mom(igaus,1) = gpvol(1,igaus,ielem)*(tmp1_mom(1)+tmp2_mom(1))
                       Re_mom(igaus,2) = gpvol(1,igaus,ielem)*(tmp1_mom(2)+tmp2_mom(2))
                       Re_mom(igaus,3) = gpvol(1,igaus,ielem)*(tmp1_mom(3)+tmp2_mom(3))

                       Re_mass(igaus) = gpvol(1,igaus,ielem)*0.5d0*(tmp1_mass+aux3_mom*rhol(igaus)+tmp4_mass)
                       Re_ener(igaus) = gpvol(1,igaus,ielem)*(0.5d0*(tmp1_ener+El(igaus)*aux3_mom+tmp4_ener)+tmp5_ener)
                    end do
                    !
                    ! Final assembly
                    !
                    !$acc loop vector collapse(2)
                    do idime = 1,ndime
                       do inode = 1,nnode
                          !$acc atomic update
                          Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                          !$acc end atomic
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
                 end do
                 !$acc end parallel loop

                 call nvtxEndRange

              end subroutine full_convec

              subroutine full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener)

                 implicit none

                 integer(4), intent(in)  :: nelem, npoin
                 integer(4), intent(in)  :: connec(nelem,nnode)
                 real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                 real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
                 real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                 integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                 real(8),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
                 real(8),    intent(out) :: Rmass(npoin)
                 real(8),    intent(out) :: Rmom(npoin,ndime)
                 real(8),    intent(out) :: Rener(npoin)
                 integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
                 real(8)                 :: Re_mom(nnode,ndime)
                 real(8)                 :: Re_mass(nnode), Re_ener(nnode)
                 real(8)                 :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime)
                 real(8)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ
                 real(8)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime)
                 real(8), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

                 call nvtxStartRange("Full convection")
                 !$acc kernels
                 Rmom(:,:) = 0.0d0
                 Rmass(:) = 0.0d0
                 Rener(:) = 0.0d0
                 !$acc end kernels

                 !$acc parallel loop gang private(Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,fl,fel) !!vector_length(vecLength)
                 do ielem = 1,nelem
                    !$acc loop vector collapse(2)
                    do idime = 1,ndime
                       do inode = 1,nnode
                          ul(inode,idime) = u(connec(ielem,inode),idime)
                          ql(inode,idime) = q(connec(ielem,inode),idime)
                          fel(inode,idime) = (E(connec(ielem,inode))+pr(connec(ielem,inode)))*u(connec(ielem,inode),idime)
                       end do
                    end do
                    !$acc loop vector collapse(3)
                    do idime = 1,ndime
                       do jdime = 1,ndime
                          do inode = 1,nnode
                             fl(inode,idime,jdime)  = q(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime)
                          end do
                       end do
                    end do
                    !$acc loop vector
                    do inode = 1,nnode
                       rhol(inode) = rho(connec(ielem,inode))
                       El(inode) = E(connec(ielem,inode))
                       prl(inode) = pr(connec(ielem,inode))
                    end do
                    !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoP, gradIsoE,gradIsoU, gradIsoF, gradIsoQ, gradIsoFe,gradRho,gradP,gradE,gradU,divF,divU,divQ,divFe)
                    do igaus = 1,ngaus
                       !$acc loop seq
                       do ii=1,porder+1
                          dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                          dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                          dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                       end do
                       isoI = gmshAtoI(igaus) 
                       isoJ = gmshAtoJ(igaus) 
                       isoK = gmshAtoK(igaus) 

                       gradIsoRho(:) = 0.0d0
                       gradIsoP(:) = 0.0d0
                       gradIsoE(:) = 0.0d0
                       gradIsoU(:,:) = 0.0d0
                       gradIsoF(:,:,:) = 0.d0
                       gradIsoQ(:,:) = 0.d0
                       gradIsoFe(:,:) = 0.d0
                       !$acc loop seq
                       do ii=1,porder+1
                          gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                          gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                          gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                          gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*prl(invAtoIJK(ii,isoJ,isoK))
                          gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*prl(invAtoIJK(isoI,ii,isoK))
                          gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*prl(invAtoIJK(isoI,isoJ,ii))

                          gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*(prl(invAtoIJK(ii,isoJ,isoK))  + El(invAtoIJK(ii,isoJ,isoK)))
                          gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*(prl(invAtoIJK(isoI,ii,isoK)) + El(invAtoIJK(isoI,ii,isoK)))
                          gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*(prl(invAtoIJK(isoI,isoJ,ii))+ El(invAtoIJK(isoI,isoJ,ii)))
                          
                          !$acc loop seq
                          do idime=1,ndime
                             gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                             gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                             gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                             gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                             gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                             gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

                             gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                             gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                             gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                             
                             !$acc loop seq
                             do jdime=1, ndime
                                 gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*fl(invAtoIJK(ii,isoJ,isoK),idime,jdime)
                                 gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*fl(invAtoIJK(isoI,ii,isoK),idime,jdime)
                                 gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*fl(invAtoIJK(isoI,isoJ,ii),idime,jdime)
                             end do
                          end do
                       end do

                       gradRho(:) = 0.0d0
                       gradP(:) = 0.0d0
                       gradE(:) = 0.0d0
                       gradU(:,:) = 0.0d0
                       divF(:) = 0.0d0
                       divQ = 0.0d0
                       divFe = 0.0d0
                       !$acc loop seq
                       do idime=1, ndime
                          !$acc loop seq
                          do jdime=1, ndime
                              gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                              gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                              gradE(idime)   = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                              divQ = divQ + He(idime,jdime,igaus,ielem) * gradIsoQ(idime,jdime)
                              divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                              !$acc loop seq
                              do kdime=1,ndime
                                 gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                                 divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                              end do
                          end do
                       end do
                       divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3) 
                       Re_mass(igaus) = 0.5d0*(divQ+rhol(igaus)*divU) 
                       Re_ener(igaus) = 0.5d0*(divFe+(El(igaus)+prl(igaus))*divU)
                       !$acc loop seq
                       do idime=1, ndime
                          Re_mom(igaus,idime) = 0.5d0*(divF(idime)+ql(igaus,idime)*divU)+ gradP(idime)
                          Re_mass(igaus) = Re_mass(igaus) + 0.5d0*(ul(igaus,idime)*gradRho(idime))
                          Re_ener(igaus) = Re_ener(igaus) + 0.5d0*(ul(igaus,idime)*gradE(idime))
                          !$acc loop seq
                          do jdime=1, ndime
                             Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.5d0*(ul(igaus,idime)*ul(igaus,jdime)*gradRho(jdime) &
                                                                       + ql(igaus,jdime)*gradU(idime,jdime))
                          end do
                          Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                       end do
                       Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
                       Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
                    end do
                    !
                    ! Final assembly
                    !
                    !$acc loop vector collapse(2)
                    do idime = 1,ndime
                       do inode = 1,nnode
                          !$acc atomic update
                          Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                          !$acc end atomic
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
                 end do
                 !$acc end parallel loop

                 call nvtxEndRange

              end subroutine full_convec_ijk

              subroutine generic_scalar_convec_ijk(nelem,npoin,connec,Ngp, &
                    dNgp,He,gpvol,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,q,eta,u,Rconvec,alpha)

                 implicit none

                 integer(4), intent(in)  :: nelem, npoin
                 integer(4), intent(in)  :: connec(nelem,nnode)
                 real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                 real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
                 real(8),    intent(in)  :: gpvol(1,ngaus,nelem)
                 integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                 real(8),    intent(in)  :: q(npoin,ndime)
                 real(8),    intent(in)  :: eta(npoin)
                 real(8),    intent(in)  :: u(npoin,ndime)
                 real(8),    intent(in)  :: alpha(npoin)
                 real(8),    intent(out) :: Rconvec(npoin)
                 integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ii
                 real(8)                 :: gradIsoE(ndime),gradIsoU(ndime,ndime),gradIsoFe(ndime,ndime)
                 real(8)                 :: gradE(ndime),gradU(ndime,ndime),divU,divFe
                 real(8)                 :: ul(nnode,ndime), fel(nnode,ndime), etal(nnode), Re(nnode)
                 real(8), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

                 call nvtxStartRange("Generic Convection")
                 !$acc kernels
                 Rconvec(:) = 0.0d0
                 !$acc end kernels
                 !$acc parallel loop gang  private(Re,ul,fel,etal) !!vector_length(vecLength)
                 do ielem = 1,nelem
                    !$acc loop vector collapse(2)
                    do idime = 1,ndime
                       do inode = 1,nnode
                          ul(inode,idime) = u(connec(ielem,inode),idime)
                          fel(inode,idime) = q(connec(ielem,inode),idime)
                       end do
                    end do
                    !$acc loop vector
                    do inode = 1,nnode
                       etal(inode) = eta(connec(ielem,inode))
                    end do
                    !$acc loop vector private(gradIsoE,gradIsoU,gradIsoFe,gradE,gradU,divU,divFe, dlxi_ip, dleta_ip, dlzeta_ip)
                    do igaus = 1,ngaus
                       !$acc loop seq
                       do ii=1,porder+1
                          dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                          dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                          dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                       end do
                       isoI = gmshAtoI(igaus) 
                       isoJ = gmshAtoJ(igaus) 
                       isoK = gmshAtoK(igaus) 

                       gradIsoE(:) = 0.0d0
                       gradIsoU(:,:) = 0.0d0
                       gradIsoFe(:,:) = 0.d0
                       !$acc loop seq
                       do ii=1,porder+1
                          gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*etal(invAtoIJK(ii,isoJ,isoK))
                          gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*etal(invAtoIJK(isoI,ii,isoK))
                          gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*etal(invAtoIJK(isoI,isoJ,ii))
                          !$acc loop seq
                          do idime=1,ndime
                             gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                             gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                             gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                             gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                             gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                             gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                          end do
                       end do

                       gradE(:) = 0.0d0
                       divU = 0.0d0
                       divFe = 0.0d0
                       !$acc loop seq
                       do idime=1, ndime
                          !$acc loop seq
                          do jdime=1, ndime
                              gradE(idime) = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                              divU = divU + He(idime,jdime,igaus,ielem) * gradIsoU(idime,jdime)
                              divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                          end do
                       end do

                       Re(igaus) = 0.5d0*(divFe+etal(igaus)*divU)
                       !$acc loop seq
                       do idime=1, ndime
                          Re(igaus) = Re(igaus) + 0.5d0*(ul(igaus,idime)*gradE(idime))
                       end do
                       Re(igaus) = gpvol(1,igaus,ielem)*Re(igaus)
                    end do
                    !
                    ! Assembly
                    !
                    !$acc loop vector
                    do inode = 1,nnode
                       !$acc atomic update
                       Rconvec(connec(ielem,inode)) = Rconvec(connec(ielem,inode))+Re(inode)
                       !$acc end atomic
                    end do
                 end do
                 !$acc end parallel loop
                 call nvtxEndRange

              end subroutine generic_scalar_convec_ijk
end module elem_convec
