module mod_constants

      implicit none

      !
      ! Time-step to read, matching the vtk file
      !
      integer(4), parameter :: istep=0

      !
      ! Dimensions
      !
      integer(4), parameter :: ndime=3
      
      !
      ! Linear mesh characteristics
      !
      integer(4), parameter :: lin_nnode=8
      integer(4), parameter :: lin_porder=1
      integer(4), parameter :: lin_npbou=4
      integer(4), parameter :: lin_ngaus=8

      !
      ! High order mesh characteristics
      !
      integer(4), parameter :: nnode=64
      integer(4), parameter :: porder=3
      integer(4), parameter :: npbou=16
      integer(4), parameter :: ngaus=64

      !
      ! Common characteristics
      !
      integer(4), parameter :: nncorner=8

      !
      ! Constants
      !
      real(8), parameter :: v_pi=3.1415926535897932385d0

end module mod_constants
