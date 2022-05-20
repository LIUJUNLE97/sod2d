module quadrature_rules

   use mod_constants
   use mod_maths

   contains

      subroutine chebyshev_hex(atoIJK,xgp,wgp)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Closed rule quadrature that uses the     !
         ! Chebyshev grid+endpoints as abcissas.   !
         ! Weights are obtained by evaluating       !
         ! w_j = int(l^n_i(xi_j),-1,1).             !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         integer(4), intent(in)  :: atoIJK(ngaus)
         real(8),    intent(out) :: xgp(ngaus,ndime), wgp(ngaus)
         integer(4)              :: inode, i, j, k, lorder(porder+1)
         real(8)                 :: xi(porder+1),w1d(porder+1),w0,w1,w2,w3

         call chebyshev_roots(xi)
         lorder(1) = 1
         lorder(2) = porder+1
         do i = 3,porder+1
            lorder(i) = i-1
         end do
         if (porder == 3) then
            w1d(1:4) = [1.0d0/9.0d0, 8.0d0/9.0d00, 8.0d0/9.0d0, 1.0d0/9.0d0]
         else
            STOP(1)
         end if
         inode = 0
         do k = 1,porder+1
            do i = 1,porder+1
               do j = 1,porder+1
                  inode = inode + 1
                  xgp(atoIJK(inode),1:3) = [xi(lorder(i)), xi(lorder(j)), xi(lorder(k))]
                  wgp(atoIJK(inode)) = w1d(lorder(i))*w1d(lorder(j))*w1d(lorder(k))
               end do
            end do
         end do

      end subroutine chebyshev_hex

      subroutine lagrange_hex(atoIJK,xgp,wgp)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Closed rule quadrature that uses the     !
         ! Lagrange grid+endpoints as abcissas.     s!
         ! Weights are obtained by evaluating       !
         ! w_j = int(l^n_i(xi_j),-1,1).             !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         implicit none

         integer(4), intent(in)  :: atoIJK(ngaus)
         real(8),    intent(out) :: xgp(ngaus,ndime), wgp(ngaus)
         integer(4)              :: inode, i, j, k, lorder(porder+1)
         real(8)                 :: xi(porder+1),w1d(porder+1),w0,w1,w2,w3

         call lagrange_roots(xi)
         lorder(1) = 1
         lorder(2) = porder+1
         do i = 3,porder+1
            lorder(i) = i-1
         end do
         if (porder == 3) then
            w1d(1:4) = [1.0d0/4.0d0, 3.0d0/4.0d00, 3.0d0/4.0d0, 1.0d0/4.0d0]
         else
            STOP(1)
         end if
         inode = 0
         do k = 1,porder+1
            do i = 1,porder+1
               do j = 1,porder+1
                  inode = inode + 1
                  xgp(atoIJK(inode),1:3) = [xi(lorder(i)), xi(lorder(j)), xi(lorder(k))]
                  wgp(atoIJK(inode)) = w1d(lorder(i))*w1d(lorder(j))*w1d(lorder(k))
               end do
            end do
         end do
         
      end subroutine lagrange_hex

end module quadrature_rules
