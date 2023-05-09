module mod_constants

      implicit none

      integer(4), parameter :: rp = 4 !(4/8)
      integer(4), parameter :: rp_vtk = 4 !(4/8)

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

      !
      ! Other constants
      !
      real(rp), parameter :: v_pi = 2.0_rp*asin(1.0_rp) ! Value of Pi

      ! No of boundary codes
      integer(4), parameter :: max_num_bou_codes = 10

      !
      ! Boundary Conditions Types
      !

      integer(4), parameter :: bc_type_far_field            = 1
      integer(4), parameter :: bc_type_non_slip_adiabatic   = 2
      integer(4), parameter :: bc_type_non_slip_hot         = 3
      integer(4), parameter :: bc_type_non_slip_cold        = 4
      integer(4), parameter :: bc_type_slip_adiabatic       = 5
      integer(4), parameter :: bc_type_slip_wall_model      = 6


      !
      ! Types of implicit solvers
      !

      integer(4), parameter :: implicit_solver_esdirk      = 1
      integer(4), parameter :: implicit_solver_bdf2_rk10   = 2
      

end module mod_constants
