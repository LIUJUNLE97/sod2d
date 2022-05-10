module mod_constants

      implicit none

      !
      ! Dimensions
      !
      integer(4), parameter :: ndime=3
      
      !
      ! Element characteristics
      !
      integer(4), parameter :: nnode=64
      integer(4), parameter :: porder=3
      integer(4), parameter :: npbou=16
      integer(4), parameter :: ngaus=64

      !integer(4), parameter :: nnode=8
      !integer(4), parameter :: porder=1
      !integer(4), parameter :: npbou=4
      !integer(4), parameter :: ngaus=8

      !
      ! Flags
      !
      integer(4), parameter :: flag_real_diff=1
      integer(4), parameter :: flag_diff_suth=1
      integer(4), parameter :: flag_rk_order=4
      integer(4), parameter :: flag_les=1
      integer(4), parameter :: flag_solver_type=1    ! 1 = Lumped, 2 = APINV, 3 = CG
      integer(4), parameter :: flag_spectralElem=1  ! 0 for Lagrange type, 1 for Chebyshev type

      !
      ! Solver params
      !
      integer(4), parameter :: maxIter=3
      real(8)   , parameter :: tol=0.00001d0

      !
      ! Other constants
      !
      real(8), parameter :: v_pi = 2.0d0*asin(1.0d0) ! Value of Pi
      real(8), parameter :: ce = 1.0d0   
      real(8), parameter :: cmax = 0.05d0 ! for FEM 0.5 for SEM 0.05/p
      real(8), parameter :: cglob = 1.0d0
      real(8), parameter :: c_rho = 1.0d0
      real(8), parameter :: c_ener = 1.0d0
      real(8), parameter :: c_sgs = 0.07d0

      real(8) :: flag_mu_factor=1.0d0

end module mod_constants
