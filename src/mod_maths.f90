!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for introducing mathematical operations that can be          !
! reutilized, such as interpolation. Included routines:               !
!                                                                     !
!   - Interpolation                                                   !
!   - Chebyshev roots                                                 !
!   - Polynomial evaluation                                           !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_maths

   use mod_constants

   contains

      pure subroutine chebyshev_roots(xi_chb)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Computes the roots of Tn(cos(y)) = cos(ny)       !
         ! in the interval [-1,1], and includes the         !
         ! endpoints. Roots are given in the following      !
         ! order:                                           !
         !  - xi(1) = -1                                    !
         !  - xi(2) =  1                                    !
         !    xi(3:porder+1) = -cos(pi*[1:porder-1]/porder) !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         real(rp), intent(out) :: xi_chb(porder+1)
         integer(4)           :: i

         !do i = 1,porder+1
         !   xi_chb(i) = -cos(v_pi*dble(i-1)/dble(porder))
         !end do

         xi_chb(1) = -1.0_rp
         xi_chb(2) = -0.447213595499958_rp
         xi_chb(3) =  0.447213595499958_rp
         xi_chb(4) =  1.0_rp


      end subroutine chebyshev_roots

      pure subroutine eval_chebyshevPoly1(xi,Tn)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Evaluates all type 1 Chebyshev polys from        !
         ! order 0 to 1 using the recursion expression.     !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         real(rp), intent(in)  :: xi
         real(rp), intent(out) :: Tn(porder+1)
         integer(4)           :: n

         Tn(1) = 1.0_rp
         Tn(2) = xi
         do n = 3,porder+1
            Tn(n) = 2.0_rp*xi*Tn(n-1)-Tn(n-2)
         end do

      end subroutine eval_chebyshevPoly1

      pure subroutine eval_chebyshevPoly2(xi,Un)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Evaluates all type 2 Chebyshev polys from        !
         ! order 0 to 1 using the recursion expression.     !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         real(rp), intent(in)  :: xi
         real(rp), intent(out) :: Un(porder+1)
         integer(4)           :: n

         Un(1) = 1.0_rp
         Un(2) = 2.0_rp*xi
         do n = 3,porder+1
            Un(n) = 2.0_rp*xi*Un(n-1)-Un(n-2)
         end do

      end subroutine eval_chebyshevPoly2

      pure subroutine lagrange_roots(xi_lag)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Computes the equispaced loci for the Lagrange    !
         ! poly in the interval [-1,1], and includes the    !
         ! endpoints. Roots are given in the following      !
         ! order:                                           !
         !  - xi(1) = -1                                    !
         !  - xi(2) =  1                                    !
         !    xi(3:porder+1) = xi(1)+[1:N-1]*(2/N)          !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         implicit none

         real(rp), intent(out) :: xi_lag(porder+1)
         integer(4)           :: i

         do i = 1,porder+1
            xi_lag(i) = -1.0_rp+(2.0_rp*real(i-1,rp)/real(porder,rp))
         end do

      end subroutine lagrange_roots

      pure subroutine eval_lagrangePoly(xi,xi_p,l_ip)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Evaluates the Lagrange poly of order N at a      !
         ! location xi(p), where xi E [-1,1]. Returns an    !
         ! array with all possible values l_ip can assume.  !
         ! Remark that l_ip = 1 if p==i and 0 if            !
         ! p==j~=i.                                         !
         ! Expects a grid of the form:                      !
         !                                                  !
         ! -1                        1 xi                   !
         !  o----v----v----v----v----o                      !
         !  0    2    3    4    5    1 i                    !
         !                                                  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none


         real(rp),    intent(in)  :: xi(porder+1), xi_p
         real(rp),    intent(out) :: l_ip(porder+1)
         integer(4)              :: i, j, lorder(porder+1)

         lorder(1) = 1
         lorder(2) = porder+1
         do i = 3,porder+1
            lorder(i) = i-1
         end do
         do i = 0,porder ! Nodal loop
            l_ip(i+1) = 1.0_rp
            do j = 0,porder ! Product series
               if (j .ne. (lorder(i+1)-1)) then
                  l_ip(i+1) = l_ip(i+1)*((xi_p-xi(j+1))/(xi(lorder(i+1))-xi(j+1)))
               end if
            end do
         end do

      end subroutine eval_lagrangePoly

      pure subroutine eval_lagrangePolyDeriv(xi,xi_p,dl_ip)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Evaluates the derivative of the Lagrange poly    !
         ! of order N at a location xi(p), where            !
         ! xi E [-1,1]. Returns an array with all possible  !
         ! values dl_ip can assume.                         !
         ! Expects a grid of the form:                      !
         !                                                  !
         ! -1                        1 xi                   !
         !  o----v----v----v----v----o                      !
         !  0    2    3    4    5    1 i                    !
         !                                                  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none


         real(rp),    intent(in)  :: xi(porder+1), xi_p
         real(rp),    intent(out) :: dl_ip(porder+1)
         integer(4)              :: i, j, m, lorder(porder+1)
         real(rp)                 :: aux

         lorder(1) = 1
         lorder(2) = porder+1
         do i = 3,porder+1
            lorder(i) = i-1
         end do
         do i = 0,porder ! Nodal loop
            dl_ip(i+1) = 0.0_rp
            do j = 0,porder ! Product series
               aux = 1.0_rp
               if (j .ne. (lorder(i+1)-1)) then
                  do m = 0,porder
                     if (m .ne. j .and. m .ne. lorder(i+1)-1) then
                        aux = aux * (xi_p-xi(m+1))/(xi(lorder(i+1))-xi(m+1))
                     end if
                  end do
                  dl_ip(i+1) = dl_ip(i+1) + (1.0_rp/(xi(lorder(i+1))-xi(j+1)))*aux
               end if
            end do
         end do

      end subroutine eval_lagrangePolyDeriv

      pure subroutine TripleTensorProduct(xi_grid,s,t,z,atoIJK,N,dN,dlxigp_ip)

         implicit none

         integer(4), intent(in)       :: atoIJK(nnode)
         real(rp), intent(in)          :: s, t, z, xi_grid(porder+1)
         real(rp), intent(out)         :: N(nnode), dN(ndime,nnode)
         real(rp),    optional, intent(out) :: dlxigp_ip(ndime,porder+1)
         integer(4)                   :: i, j, k, c
         real(rp), dimension(porder+1) :: lxi_ip, leta_ip, lzeta_ip
         real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

         call eval_lagrangePoly(xi_grid,s,lxi_ip)
         call eval_lagrangePoly(xi_grid,t,leta_ip)
         call eval_lagrangePoly(xi_grid,z,lzeta_ip)
         call eval_lagrangePolyDeriv(xi_grid,s,dlxi_ip)
         call eval_lagrangePolyDeriv(xi_grid,t,dleta_ip)
         call eval_lagrangePolyDeriv(xi_grid,z,dlzeta_ip)

         if(present(dlxigp_ip)) then
            dlxigp_ip(1,:) = dlxi_ip(:)
            dlxigp_ip(2,:) = dleta_ip(:)
            dlxigp_ip(3,:) = dlzeta_ip(:)
         end if

         c = 0
         do k = 1,porder+1
            do i = 1,porder+1
               do j = 1,porder+1
                  c = c+1
                  N(atoIJK(c)) = lxi_ip(i)*leta_ip(j)*lzeta_ip(k)
                  dN(1,atoIJK(c)) = dlxi_ip(i)*leta_ip(j)*lzeta_ip(k)
                  dN(2,atoIJK(c)) = lxi_ip(i)*dleta_ip(j)*lzeta_ip(k)
                  dN(3,atoIJK(c)) = lxi_ip(i)*leta_ip(j)*dlzeta_ip(k)
               end do
            end do
         end do

      end subroutine TripleTensorProduct

      pure subroutine var_interpolate(var,Neval,var_a)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Interpolates a variable using element shape      !
         ! functions N(xi,eta,zeta). Given a coordinate     !
         ! set X = (x,y,z), v(X) = N_a*v_a, where v_a are   !
         ! element nodal values.
         ! Expects a grid of the form:                      !
         !                                                  !
         ! -1                        1 xi                   !
         !  o----v----v----v----v----o                      !
         !  0    2    3    4    5    1 i                    !
         !                                                  !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         implicit none

         real(rp), intent(in)  :: var(nnode), Neval(nnode)
         real(rp), intent(out) :: var_a
         integer(4)           :: inode

         var_a = 0.0_rp
         do inode = 1,nnode
            var_a = var_a+Neval(inode)*var(inode)
         end do

      end subroutine var_interpolate

end module mod_maths
