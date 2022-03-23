module mod_constants

      implicit none

      !
      ! Dimensions
      !
      integer(4), parameter :: ndime=3
      
      !
      ! Element characteristics
      !
      !integer(4), parameter :: nnode=27
      !integer(4), parameter :: porder=2
      !integer(4), parameter :: npbou=9
      !integer(4), parameter :: ngaus=27

      integer(4), parameter :: nnode=8
      integer(4), parameter :: porder=1
      integer(4), parameter :: npbou=4
      integer(4), parameter :: ngaus=8

      !
      ! Flags
      !
      integer(4), parameter :: flag_real_diff=1
      integer(4), parameter :: flag_diff_suth=1
      integer(4), parameter :: flag_rk_order=4

      real(8), parameter :: ce = 1.0d0
      real(8), parameter :: cmax = 0.5d0
      real(8), parameter :: cglob = 1.0d0
      real(8), parameter :: c_rho = 0.1d0
      real(8), parameter :: c_ener = 0.1d0

      real(8) :: flag_mu_factor=1.0d0

end module mod_constants
