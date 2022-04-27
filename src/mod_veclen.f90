module mod_veclen

   use mod_constants

   implicit none

   integer(4) :: vecLength
   PUBLIC     :: vecLength

   contains

      subroutine define_veclen()

         implicit none

         !*********************************************************************!
         ! Define vector length to be used                                     !
         !*********************************************************************!

         if (nnode .gt. 32 .and. mod(nnode,32) .ne. 0) then
            vecLength = ((nnode/32)+1)*32 ! Ensures enough threads are present
         else if (nnode .gt. 32 .and. mod(nnode,32) == 0) then
            vecLength = nnode ! One thread per node
         else if (nnode .le. 32) then
            vecLength = 32 ! Default for small elements
         end if

      end subroutine define_veclen

end module mod_veclen
