module elem_qua

	use mod_constants
	use mod_maths

	implicit none
	integer(4), allocatable :: quad_order_edges(:,:)

	contains

		subroutine init_quad_info()
			implicit none
			allocate(quad_order_edges(4,2))
			quad_order_edges = transpose(reshape([1,2,2,3,3,4,4,1],(/2,4/)))
			!$acc enter data copyin(quad_order_edges)
		end subroutine init_quad_info

		subroutine quad_highorder(mporder,mnpbou,xi,eta,atoIJ,N,dN)
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