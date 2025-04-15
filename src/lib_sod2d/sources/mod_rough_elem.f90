module mod_rough_elem
   use mod_constants
   use mod_numerical_params
   use mod_nvtx

   implicit none

   contains

   subroutine mom_source_rough_elem(npoin,coord,rho,u,source_term)


      implicit none

      integer(4), intent(in)          :: npoin
      real(rp),    intent(in)         :: coord(npoin,ndime),rho(npoin), u(npoin,ndime)
      real(rp),    intent(inout)      :: source_term(npoin,ndime)
      integer(4) ::  inode      
      real(rp)   :: cd = 1.0_rp, xmin , xmax, ymin, ymax


      xmin = x_trip_o
      xmax = xmin+l_trip_x
      ymin = y_trip_o
      ymax = ymin+l_trip_y

      !$acc parallel loop
      do inode = 1,npoin
         if((coord(inode,2)<ymax)  .and. (coord(inode,2)>ymin)) then
            if((coord(inode,1)<xmax)  .and. (coord(inode,1)>xmin)) then
               source_term(iNode,1) = -0.5_rp*rho(iNode)*cd*u(iNode,1)*abs(u(iNode,1))/l_trip_x
               source_term(iNode,2) = 0.00_rp
               source_term(iNode,3) = 0.00_rp
            end if
         end if
      end do
      !$acc end parallel loop

   end subroutine mom_source_rough_elem

   
   subroutine ener_source_rough_elem(npoin,coord,rho,q,source_term)


      implicit none

      integer(4), intent(in)          :: npoin
      real(rp),    intent(in)         :: coord(npoin,ndime),rho(npoin), q(npoin,ndime)
      real(rp),    intent(inout)      :: source_term(npoin,5)
      integer(4) ::  inode      
      real(rp)   :: cd = 1.0_rp, xmin , xmax, ymin, ymax


      xmin = x_trip_o
      xmax = xmin+l_trip_x
      ymin = y_trip_o
      ymax = ymin+l_trip_y

      !$acc parallel loop
      do inode = 1,npoin
         if((coord(inode,2)<ymax)  .and. (coord(inode,2)>ymin)) then
            if((coord(inode,1)<xmax)  .and. (coord(inode,1)>xmin)) then
               source_term(iNode,2) = q(inode,1)*(source_term(inode,3))
            end if
         end if
      end do
      !$acc end parallel loop

   end subroutine ener_source_rough_elem


    



end module mod_rough_elem