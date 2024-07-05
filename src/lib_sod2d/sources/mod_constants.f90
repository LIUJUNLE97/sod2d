module mod_constants

	implicit none

	integer(4), parameter :: rp = 8 !(4/8)
	integer(4), parameter :: rp_vtk = 4 !(4/8)
	integer(4), parameter :: rp_avg = 8 !(4/8)

	!
	! Dimensions
	!
	integer(4), parameter :: ndime=3

	!
	! Element characteristics
	!
	integer(4), parameter :: porder=4
	integer(4), parameter :: nnode=(porder+1)**3
	integer(4), parameter :: ngaus=nnode
	integer(4), parameter :: npbou=(porder+1)**2


	!
	! Other constants
	!
	real(rp), parameter :: v_pi = 2.0_rp*asin(1.0_rp) ! Value of Pi

	! No of boundary codes
	integer(4), parameter :: max_num_bou_codes = 15

	! No of max saved fields (size of pointer arrays)
	integer(4), parameter :: max_num_saved_fields = 50

	! No of max boundaries per elements (used in mesh_conversion_tool)
	integer(4),parameter :: maxBoundsPerElem = 4

	!
	! Boundary Conditions Types
	!

	integer(4), parameter :: bc_type_far_field            = 1
	integer(4), parameter :: bc_type_far_field_supersonic = 2
	integer(4), parameter :: bc_type_recirculation_inlet  = 3
	integer(4), parameter :: bc_type_outlet_incomp        = 4
	integer(4), parameter :: bc_type_outlet_supersonic    = 5
	integer(4), parameter :: bc_type_non_slip_adiabatic   = 6
	integer(4), parameter :: bc_type_non_slip_hot         = 7
	integer(4), parameter :: bc_type_non_slip_cold        = 8
	integer(4), parameter :: bc_type_symmetry			  = 9
	integer(4), parameter :: bc_type_slip_adiabatic       = 10
	integer(4), parameter :: bc_type_slip_isothermal      = 11
	integer(4), parameter :: bc_type_slip_wall_model      = 12
	integer(4), parameter :: bc_type_slip_atmosphere      = 13
	

	! Type of wall models
	integer(4), parameter :: wmles_type_reichardt = 1
	integer(4), parameter :: wmles_type_abl		  = 2
	!
	! Types of implicit solvers
	!

	integer(4), parameter :: implicit_solver_imex   	 = 1


end module mod_constants
