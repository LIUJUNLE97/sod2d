module mod_constants

      implicit none

      !
      ! Dimensions
      !
      integer(4), parameter :: ndime=3
      
      !
      ! Element characteristics
      !
      integer(4), parameter :: nnode=27
      integer(4), parameter :: porder=2
      integer(4), parameter :: npbou=9
      integer(4), parameter :: ngaus=27

      !integer(4), parameter :: nnode=8
      !integer(4), parameter :: porder=1
      !integer(4), parameter :: npbou=4
      !integer(4), parameter :: ngaus=8

      !
      ! Flags
      !
      integer(4), parameter :: flag_real_diff=0

end module mod_constants
