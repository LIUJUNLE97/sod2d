module mod_solver

	use mod_numerical_params
	use mod_comms
	use mod_mpi
	use mod_nvtx
	use mod_time_ops
	use mod_bc_routines
	use elem_diffu

	implicit none
      
	real(rp)  , allocatable, dimension(:,:) :: x_vars, r0_vars, p0_vars, qn_vars, b_vars,z0_vars,z1_vars,M_vars
	real(rp)  , allocatable, dimension(:,:) :: aux_u_vars
   real(rp)  , allocatable, dimension(:) 	:: aux_Tem_vars,tau_stab
   real(rp)  , allocatable, dimension(:,:) :: ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ
   logical  :: flag_cg_mem_alloc_vars=.true.
	integer(4) , parameter :: nvars = 5


	contains

		subroutine lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,R)

			implicit none

			integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
			real(rp),    intent(in)    :: Ml(npoin)
			real(rp),    intent(inout) :: R(npoin)
			integer(4)                :: ipoin

			!$acc parallel loop
			do ipoin = 1,npoin_w
				R(lpoin_w(ipoin)) = R(lpoin_w(ipoin))/Ml(lpoin_w(ipoin))
			end do
			!$acc end parallel

		end subroutine lumped_solver_scal

		subroutine lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,R)

			implicit none

			integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
			real(rp),    intent(in)    :: Ml(npoin)
			real(rp),    intent(inout) :: R(npoin,ndime)
			integer(4)                :: idime, ipoin

			! TODO: reverse this loop, might be faster on CPU
			!$acc parallel loop collapse(2)
			do ipoin = 1,npoin_w
				do idime = 1,ndime
					R(lpoin_w(ipoin),idime) = R(lpoin_w(ipoin),idime)/Ml(lpoin_w(ipoin))
				end do
			end do
			!$acc end  parallel loop

		end subroutine lumped_solver_vect

end module mod_solver