module mod_numerical_params

    use mod_constants
    
    implicit none

        ! Time integration algorithm
        integer(4) :: flag_rk_order=4 
        integer(4) :: flag_implicit=0 !Explicit: RK, Implicit: BDF2


        ! LES 
        integer(4) :: flag_les=0
        integer(4) :: flag_les_ilsa=0
        real(rp) :: c_sgs = 0.025_rp
        real(rp) :: stau   = 0.022_rp
        real(rp) :: T_ilsa = 1.0_rp

        ! Discretization parameters
        logical :: flag_total_enthalpy = .false.
        integer(4)  :: flag_solver_type=1    ! 1 = Lumped, 2 = APINV, 3 = CG
        integer(4)  :: flag_spectralElem=1  ! 0 for Lagrange type, 1 for Chebyshev type
        integer(4)  :: flag_normalise_entropy=1
        real(rp) :: ce = 0.1_rp   
        real(rp) :: cmax = 0.5_rp 
        real(rp) :: cglob =1.0_rp
        real(rp) :: c_rho =1.0_rp
        real(rp) :: c_ener = 1.0_rp
        real(rp) :: flag_mu_factor=1.0_rp

        ! Implicit solver
        integer(4) :: maxIter=20
        real(8)    :: tol=1e-4
        integer(4) :: flag_use_constant_dt = 0
        integer(4) :: implicit_solver = implicit_solver_imex

        !
        ! Reference conditions
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
        ! Penalisation buffer zone
        !

        logical :: flag_buffer_on = .false.
        logical :: flag_buffer_on_east = .false.
        logical :: flag_buffer_on_west = .false.
        logical :: flag_buffer_on_north = .false.
        logical :: flag_buffer_on_south = .false.
        logical :: flag_buffer_on_top = .false.
        logical :: flag_buffer_on_bottom = .false.

        real(rp) :: flag_buffer_e_min = 0.0_rp
        real(rp) :: flag_buffer_e_size= 0.0_rp
        real(rp) :: flag_buffer_w_min = 0.0_rp
        real(rp) :: flag_buffer_w_size = 0.0_rp

        real(rp) :: flag_buffer_n_min = 0.0_rp
        real(rp) :: flag_buffer_n_size = 0.0_rp
        real(rp) :: flag_buffer_s_min = 0.0_rp
        real(rp) :: flag_buffer_s_size = 0.0_rp

        real(rp) :: flag_buffer_t_min = 0.0_rp
        real(rp) :: flag_buffer_t_size = 0.0_rp
        real(rp) :: flag_buffer_b_min = 0.0_rp
        real(rp) :: flag_buffer_b_size = 0.0_rp

        !
        ! Wall model averaging
        !
        real(rp)    :: period_walave   = 1.0_rp
        logical     :: flag_walave     = .false.
        integer(4)  :: flag_walex      = 3
        integer(4)  :: flag_type_wmles = wmles_type_reichardt

        !
	    ! viscosity
	    !
	    integer(4) :: flag_real_diff=1
	    integer(4) :: flag_diff_suth=1
        real(rp)   :: incomp_viscosity=0.0_rp
        !
        ! FS
        !
        logical :: flag_fs_fix_pressure = .true.

        !
        ! CG
        !

        logical :: flag_cg_prec_bdc = .false. !Block Diagonal Cholesky
end module mod_numerical_params
