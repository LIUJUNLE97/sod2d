module elem_qua

use mod_numerical_params
use mod_maths

	implicit none
	!TODO: make thhis array allocatbale
	integer(4), parameter :: quad_order_edges(4,2) = transpose(reshape([1,2,2,3,3,4,4,1],(/2,4/)))

	contains

		subroutine quad_highorder(mporder,mnpbou,xi,eta,atoIJ,N,dN) ! QUA16 element
			implicit none
			integer(4),intent(in) :: mporder,mnpbou
			real(rp),intent(in)   :: xi,eta
			integer(4),intent(in) :: atoIJ(mnpbou)
			real(rp),intent(out)  :: N(mnpbou), dN(2,mnpbou)
			real(rp)              :: xi_grid(mporder+1)
			call getGaussLobattoLegendre_roots(mporder,xi_grid)
			call DoubleTensorProduct(mporder,mnpbou,xi_grid,xi,eta,atoIJ,N,dN)
		end subroutine quad_highorder

end module