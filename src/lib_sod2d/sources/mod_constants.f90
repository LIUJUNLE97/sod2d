module mod_constants

      implicit none

      integer(4), parameter::rp = 4

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
      integer(4), parameter :: flag_les_ilsa=0
      integer(4), parameter :: flag_solver_type=1    ! 1 = Lumped, 2 = APINV, 3 = CG
      integer(4), parameter :: flag_spectralElem=1  ! 0 for Lagrange type, 1 for Chebyshev type
      integer(4), parameter :: flag_normalise_entropy=1
      integer(4), parameter :: flag_activate_wall_model=1


      logical, parameter :: save_vtk = .false.
      logical, parameter :: save_hdf5 = .true.

      !
      ! Solver params
      !
      integer(4), parameter :: maxIter=3
      real(rp)   , parameter :: tol=0.00001_rp

      !
      ! Other constants
      !
      real(rp), parameter :: v_pi = 2.0_rp*asin(1.0_rp) ! Value of Pi
      real(rp), parameter :: ce = 1.0_rp   
      real(rp), parameter :: cmax = 0.5_rp 
      real(rp), parameter :: cglob =1.0_rp
      real(rp), parameter :: c_rho =1.0_rp
      real(rp), parameter :: c_ener = 1.0_rp
      real(rp), parameter :: c_sgs = 0.07_rp
      real(rp), parameter :: stau   = 0.022_rp
      real(rp), parameter :: T_ilsa = 20.0_rp
      real(rp), parameter :: T_wmles = 0.01_rp

      integer(4), parameter :: max_num_bou_codes = 10

      real(rp) :: flag_mu_factor=1.0_rp

      !
      ! NSCBC parameters
      !
      real(rp) :: nscbc_u_inf   = 1.0_rp
      real(rp) :: nscbc_p_inf = 1.0_rp
      real(rp) :: nscbc_gamma_inf = 1.0_rp
      real(rp) :: nscbc_c_inf = 1.0_rp
      real(rp) :: nscbc_rho_inf   = 1.0_rp
      real(rp) :: nscbc_Rgas_inf   = 1.0_rp
      real(rp) :: nscbc_T_H   = 293.0_rp
      real(rp) :: nscbc_T_C   = 293.0_rp

      !
      ! Boundary Conditions Types
      !

      integer(4), parameter :: bc_type_inlet                = 1
      integer(4), parameter :: bc_type_non_slip_adiabatic   = 2
      integer(4), parameter :: bc_type_non_slip_hot         = 3
      integer(4), parameter :: bc_type_non_slip_cold        = 4
      integer(4), parameter :: bc_type_slip_adiabatic       = 5
      integer(4), parameter :: bc_type_slip_wall_model      = 6
      integer(4), parameter :: bc_type_outlet               = 7

end module mod_constants
