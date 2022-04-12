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

         real(8), intent(out) :: xi_chb(porder+1)
         integer(4)           :: i

         do i = 1,porder+1
            xi_chb(i) = -cos(v_pi*dble(i-1)/dble(porder))
         end do

      end subroutine chebyshev_roots

      pure subroutine eval_chebyshevPoly1(xi,Tn)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Evaluates all type 1 Chebyshev polys from        !
         ! order 0 to 1 using the recursion expression.     !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         real(8), intent(in)  :: xi
         real(8), intent(out) :: Tn(porder+1)
         integer(4)           :: n

         Tn(1) = 1.0d0
         Tn(2) = xi
         do n = 3,porder+1
            Tn(n) = 2.0d0*xi*Tn(n-1)-Tn(n-2)
         end do

      end subroutine eval_chebyshevPoly1

      pure subroutine eval_chebyshevPoly2(xi,Un)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Evaluates all type 2 Chebyshev polys from        !
         ! order 0 to 1 using the recursion expression.     !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         real(8), intent(in)  :: xi
         real(8), intent(out) :: Un(porder+1)
         integer(4)           :: n

         Un(1) = 1.0d0
         Un(2) = 2.0d0*xi
         do n = 3,porder+1
            Un(n) = 2.0d0*xi*Un(n-1)-Un(n-2)
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

         real(8), intent(out) :: xi_lag(porder+1)
         integer(4)           :: i

         do i = 1,porder+1
            xi_lag(i) = -1.0d0+(2.0d0*dble(i-1)/dble(porder))
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

         real(8),    intent(in)  :: xi(porder+1), xi_p
         real(8),    intent(out) :: l_ip(porder+1)
         integer(4)              :: i, j

         do i = 0,porder ! Nodal loop
            l_ip(i+1) = 1.0d0
            do j = 0,porder ! Product series
               if (j .ne. i) then
                  l_ip(i+1) = l_ip(i+1)*((xi_p-xi(j+1))/(xi(i+1)-xi(j+1)))
               end if
            end do
         end do

      end subroutine eval_lagrangePoly

      !!pure subroutine var_interpolate()

      !!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!   ! Interpolates a variable using element shape      !
      !!   ! functions N(xi,eta,zeta). Given a coordinate     !
      !!   ! set X = (x,y,z), v(X) = N_a*v_a, where v_a are   !
      !!   ! element nodal values.
      !!   ! Expects a grid of the form:                      !
      !!   !                                                  !
      !!   ! -1                        1 xi                   !
      !!   !  o----v----v----v----v----o                      !
      !!   !  0    2    3    4    5    1 i                    !
      !!   !                                                  !
      !!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!   !
      !!   implicit none

      !!   var_a = 0.0d0
      !!   do inode = 1,nnode
      !!      var_a = var_a+N(inode)*var(inode)
      !!   end do

      !!end subroutine var_interpolate

end module mod_maths
