module mod_solver

      use mod_constants
      use mod_nvtx

      contains

              subroutine lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,R)

                      implicit none

                      integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
                      real(8),    intent(in)    :: Ml(npoin)
                      real(8),    intent(inout) :: R(npoin)
                      integer(4)                :: ipoin

                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         R(lpoin_w(ipoin)) = R(lpoin_w(ipoin))/Ml(lpoin_w(ipoin))
                      end do
                      !$acc end parallel

              end subroutine lumped_solver_scal

              subroutine lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,R)

                      implicit none

                      integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
                      real(8),    intent(in)    :: Ml(npoin)
                      real(8),    intent(inout) :: R(npoin,ndime)
                      integer(4)                :: idime, ipoin

                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            R(lpoin_w(ipoin),idime) = R(lpoin_w(ipoin),idime)/Ml(lpoin_w(ipoin))
                         end do
                      end do
                      !$acc end  parallel loop

              end subroutine lumped_solver_vect

              subroutine approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,R)

                      use mass_matrix

                      implicit none

                      integer(4), intent(in)       :: nelem,npoin,ppow,npoin_w
                      integer(4), intent(in)       :: connec(nelem,nnode),lpoin_w(npoin_w)
                      real(8),    intent(in)       :: gpvol(1,ngaus,nelem),Ngp(ngaus,nnode),Ml(npoin)
                      real(8),    intent(inout)    :: R(npoin)
                      integer(4)                   :: ipoin, ipow
                      real(8),    dimension(npoin) :: b, v, x

                      call nvtxStartRange("Scalar APINV")

                      !
                      ! Initialize series at k=0
                      !
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         b(lpoin_w(ipoin)) = R(lpoin_w(ipoin))
                         v(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                         x(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !
                      ! Step over sucessive powers
                      !
                      do ipow = 1,ppow
                         call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,v,b)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            v(lpoin_w(ipoin)) = v(lpoin_w(ipoin))- &
                               (b(lpoin_w(ipoin))/Ml(lpoin_w(ipoin)))
                            x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+v(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end do
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         R(lpoin_w(ipoin)) = x(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine approx_inverse_scalar

              subroutine approx_inverse_vect(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,R)

                      use mass_matrix

                      implicit none

                      integer(4), intent(in)       :: nelem,npoin,ppow,npoin_w
                      integer(4), intent(in)       :: connec(nelem,nnode),lpoin_w(npoin_w)
                      real(8),    intent(in)       :: gpvol(1,ngaus,nelem),Ngp(ngaus,nnode),Ml(npoin)
                      real(8),    intent(inout)    :: R(npoin,ndime)
                      integer(4)                   :: ipoin, idime, ipow
                      real(8),    dimension(npoin) :: b, v, x

                      call nvtxStartRange("Vector APINV")

                      do idime = 1,ndime
                         !
                         ! Initialize series at k=0
                         !
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            b(lpoin_w(ipoin)) = R(lpoin_w(ipoin),idime)
                            v(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                            x(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop

                         !
                         ! Step over sucessive powers
                         !
                         do ipow = 1,ppow
                            call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,v,b)
                            !$acc parallel loop
                            do ipoin = 1,npoin_w
                               v(lpoin_w(ipoin)) = v(lpoin_w(ipoin))- &
                                  (b(lpoin_w(ipoin))/Ml(lpoin_w(ipoin)))
                               x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+v(lpoin_w(ipoin))
                            end do
                            !$acc end parallel loop
                         end do
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            R(lpoin_w(ipoin),idime) = x(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end do
                      call nvtxEndRange

              end subroutine approx_inverse_vect

              subroutine CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Mc,v,u)

                      implicit none

                      integer(4), intent(in)  :: npoin, nzdom
                      integer(4), intent(in)  :: rdom(npoin+1), cdom(nzdom)
                      real(8),    intent(in)  :: Mc(nzdom), v(npoin)
                      real(8),    intent(out) :: u(npoin)
                      integer(4)              :: ipoin, izdom, jpoin, rowb, rowe

                      call nvtxStartRange("SPMV")
                      !$acc kernels
                      u(:) = 0.0d0
                      !$acc end kernels
                      !
                      !$acc parallel loop
                      do ipoin = 1,npoin
                         !
                         ! Get CSR section for row ipoin
                         !
                         rowb = rdom(ipoin)+1
                         rowe = rdom(ipoin+1)
                         !
                         ! Loop inside CSR section
                         !
                         !!$acc loop seq
                         do izdom = rowb,rowe
                            jpoin = cdom(izdom) ! Col. index
                            u(ipoin) = u(ipoin)+Mc(izdom)*v(jpoin) ! Dot product
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine CSR_SpMV_scal

              subroutine conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,R)

                 use mass_matrix

                 implicit none

                 integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
                 real(8)   , intent(in)    :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                 real(8)   , intent(inout) :: R(npoin)
                 integer(4)                :: ipoin, iter
                 real(8), dimension(npoin) :: x, r0, p0, q
                 real(8)                   :: Q1, Q2, T1, alpha, beta
                 !
                 ! Initialize solver
                 !
                 !$acc kernels
                 x(:) = R(:)
                 !$acc end kernels
                 call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,q,x)
                 !$acc parallel loop
                 do ipoin = 1,npoin_w
                    r0(lpoin_w(ipoin)) = x(lpoin_w(ipoin))-q(lpoin_w(ipoin))
                    p0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))
                 end do
                 !$acc end parallel loop
                 !
                 ! Start iterations
                 !
                 do iter = 1,maxIter
                    call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,q,p0)
                    Q1 = dot_product(r0,r0)
                    Q2 = dot_product(p0,q)
                    alpha = Q1/Q2
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+alpha*p0(lpoin_w(ipoin))
                       r0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))-alpha*q(lpoin_w(ipoin))
                    end do
                    !$acc end parallel loop
                    T1 = dot_product(r0,r0)
                    !
                    ! Stop cond
                    !
                    if (sqrt(T1) .gt. tol) then
                       exit
                    end if
                    !
                    ! Update p
                    !
                    beta = T1/Q1
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       p0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))+beta*p0(lpoin_w(ipoin))
                    end do
                    !$acc end parallel loop
                 end do
                 !$acc kernels
                 R(:) = x(:)
                 !$acc end kernels

              end subroutine conjGrad_scalar

end module mod_solver
