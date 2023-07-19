!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module which has all the functions to use witness points in sod2d                                                                                                       !
! Includes the following subroutines:                                                                                                                                     !
!  - Read an input .txt file with all the witness points                                                                                                                  !
!  - Compute the element to which the point belongs and its isoparametric coordinates (xi_1, xi_2, xi_3) ** Possible additional subroutine for domain splitting**         !
!  - Interpolation of the magnitude to the witness point considering the value at the nodes of the element                                                                !
!  - Output of the magnitude values at the witness points in HDF5 format.                                                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

module mod_witness_points
   
   use mod_constants
   use elem_hex

   implicit none
   contains
      subroutine read_points(fname, np, xyz)
         !
         ! Subroutine which reads an ASCII file for the input of the witness points.
         ! First row of the file contains the number of points and the rest the coordinates of the wintess points as of: X Y Z
         !
         implicit none
         character(512), intent(in)  :: fname           ! Input 1: path to the witness points file   
         integer(4),    intent(in)  :: np              ! Output 1: number of witness points
         real(rp),       intent(out) :: xyz(np,ndime) ! Output 2: coordinates of the witness points in a 1D array as xyz = [x1, y1, z1, ..., xn, yn, zn]
         integer(4)                 :: ii

         open(unit=99, file=fname, status='old', action='read')

         do ii = 1, np
            read(99,*) xyz(ii, :)
         end do
         close(99)
      end subroutine read_points

      subroutine isocoords(elpoints, wit, xi, isinside, Niwit)
         !
         ! Subroutine which computes the isoparametric coordinates of a point in an HEX64 element.
         ! If any of them is outside the bounds of -1 and 1 it means that the point is outside of the element.
         !
         implicit none
         real(rp), intent(in)   :: elpoints(nnode, ndime)   ! Input 1: coordinates of the element nodes following the gmesh ordering
         real(rp), intent(in)   :: wit(ndime)               ! Input 2: coordinates of the point we are looking the isoparametric coordinates from
         real(rp), intent(out)  :: xi(ndime)                ! Output 1: isoparametric coordinates of the point
         logical,  intent(out)  :: isinside
         real(rp), intent(out)  :: Niwit(nnode) 
         real(rp)               :: xi_0(ndime), xi_n(ndime)
         real(rp)               :: N(nnode), N_lagrange(nnode) 
         integer(4)             :: atoIJK(nnode)
         real(rp)               :: dlxigp_ip(ndime, porder+1)
         real(rp)               :: dN(ndime, nnode), dN_lagrange(ndime, nnode)
         real(rp)               :: f(ndime)
         real(rp)               :: a(ndime*ndime), b(ndime*ndime)
         real(rp)               :: j(ndime, ndime), k(ndime, ndime)
         real(rp)               :: detJ
         integer(4)            :: ii, ip
         real(rp), parameter    :: tol = 1e-10, alpha = 1, div = 100
         integer(4), parameter :: maxite = 50

         xi_0(:) = 0
         xi(:)   = xi_0(:)
         write(*,*) 'TODO mod_witness_points.f90 line 62'
         !call set_hex64_lists(atoIJK)
         isinside = .false.

         do ii = 1, maxite
            write(*,*) 'TODO mod_witness_points.f90 line 67'
            !call hex_highorder(xi(1), xi(2), xi(3), atoIJK, N, dN, N_lagrange, dN_lagrange, dlxigp_ip)
            f(:)   = wit(:)
            j(:,:) = 0
            do ip = 1, nnode
               f(:)   = f(:)   - N(ip)*elpoints(ip,:)
               j(1,1) = j(1,1) - dN(1,ip)*elpoints(ip,1)
               j(1,2) = j(1,2) - dN(2,ip)*elpoints(ip,1)
               j(1,3) = j(1,3) - dN(3,ip)*elpoints(ip,1)
               j(2,1) = j(2,1) - dN(1,ip)*elpoints(ip,2)
               j(2,2) = j(2,2) - dN(2,ip)*elpoints(ip,2)
               j(2,3) = j(2,3) - dN(3,ip)*elpoints(ip,2)
               j(3,1) = j(3,1) - dN(1,ip)*elpoints(ip,3)
               j(3,2) = j(3,2) - dN(2,ip)*elpoints(ip,3)
               j(3,3) = j(3,3) - dN(3,ip)*elpoints(ip,3)
            end do
            detJ = j(1,1)*j(2,2)*j(3,3)+j(1,2)*j(2,3)*j(3,1)+j(1,3)*j(2,1)*j(3,2)-j(3,1)*j(2,2)*j(1,3)&
            -j(3,2)*j(2,3)*j(1,1)-j(3,3)*j(2,1)*j(1,2)
            !
            ! Inverse computation
            !
            ! Minors for inverse
            !
            a(1) = j(2,2)*j(3,3)-j(3,2)*j(2,3)
            a(2) = j(2,1)*j(3,3)-j(3,1)*j(2,3)
            a(3) = j(2,1)*j(3,2)-j(3,1)*j(2,2)
            a(4) = j(1,2)*j(3,3)-j(3,2)*j(1,3)
            a(5) = j(1,1)*j(3,3)-j(3,1)*j(1,3)
            a(6) = j(1,1)*j(3,2)-j(3,1)*j(1,2)
            a(7) = j(1,2)*j(2,3)-j(2,2)*j(1,3)
            a(8) = j(1,1)*j(2,3)-j(2,1)*j(1,3)
            a(9) = j(1,1)*j(2,2)-j(2,1)*j(1,2)
            !
            ! Sign changes
            !
            a(2) = -a(2)
            a(4) = -a(4)
            a(6) = -a(6)
            a(8) = -a(8)
            !
            ! Transpose a into b
            !
            do ip = 1,9
               b(ip) = a(ip)
            end do
            b(2) = a(4)
            b(3) = a(7)
            b(4) = a(2)
            b(6) = a(8)
            b(7) = a(3)
            b(8) = a(6)
            !
            ! Divide by detj
            !
            do ip = 1,9
               b(ip) = 1.0_rp/detJ*b(ip)
            end do
            !
            ! Organize into inverse
            !
            k(1,1) = b(1)
            k(1,2) = b(2)
            k(1,3) = b(3)
            k(2,1) = b(4)
            k(2,2) = b(5)
            k(2,3) = b(6)
            k(3,1) = b(7)
            k(3,2) = b(8)
            k(3,3) = b(9)
            !
            ! Newton-Raphson method to find the new xi values
            !
            xi_n(:) = xi(:) - alpha*matmul(k, f)
            xi(:)   = xi_n(:)
            if (dot_product(f, f) < tol) then
               isinside = .true.
               Niwit    = N
               exit
            end if
            if (dot_product(f, f) > div) then
               exit
            end if
         end do
      end subroutine isocoords

      subroutine wit_interpolation(xiwit,elvalues,N,witval)
         !
         ! Subroutine which interpolates the values from the element nodes to any point inside the element
         !
         implicit none
         real(rp),intent(in)   :: xiwit(ndime)    ! Input 1: isoparametric coordinates of the point we want to interpolate to
         real(rp),intent(in)   :: elvalues(nnode) ! Input 2: values of the magnitude at the element nodes 
         real(rp),intent(in)   :: N(nnode)
         real(rp),intent(out)  :: witval          ! Output 1: value interpolated at the point
         
         witval = 0.0_rp

         !call set_hex64_lists(atoIJK)
         !call hex_highorder(xiwit(1), xiwit(2), xiwit(3), atoIJK, N, dN, N_lagrange, dN_lagrange, dlxigp_ip)
         call var_interpolate(nnode,elvalues,N,witval)
      end subroutine wit_interpolation

end module mod_witness_points
