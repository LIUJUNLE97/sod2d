module quadrature_rules

	use mod_constants
	use mod_maths

	contains

		!> @brief Computes the integral of a function over a QUA_XX
		!> @details Closed rule quadrature for a function f(x) over a
		!> QUA_X element type. The integral is computed using the
		!> given grid and weights associated with a particular element order.
		!> Notice that the order of the nodes follows that of the element itself.
		!> For boundary elements only!
		!> @param[in] atoIJ Node a to IJ rellationship
		!> @param[out] xgp Quadrature points
		!> @param[out] wgp Quadrature weights
		subroutine GaussLobattoLegendre_qua(mporder,mnpbou,atoIJ,xgp,wgp)

			implicit none
			integer(4),intent(in) :: mporder,mnpbou,atoIJ(mnpbou)
			real(8),intent(out)  :: xgp(mnpbou,ndime-1),wgp(mnpbou)
			integer(4)            :: inode,i,j,lorder(mporder+1)
			real(8)               :: xi(mporder+1),w1d(mporder+1)

			call getGaussLobattoLegendre_weights_and_roots(mporder,w1d,xi)

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
				lorder(i) = i-1
			end do

			inode = 0
			do i = 1,mporder+1
				do j = 1,mporder+1
					inode = inode + 1
					xgp(atoIJ(inode),1:2) = [xi(lorder(i)), xi(lorder(j))]
					wgp(atoIJ(inode)) = w1d(lorder(i))*w1d(lorder(j))
				end do
			end do
		end subroutine GaussLobattoLegendre_qua

		!> @brief Computes the integral of a function over a HEX_XX
		!> @details Closed rule quadrature for a function f(x) over a
		!> HEX_XX element type. The integral is computed using the
		!> given grid and weights associated with a particular element order.
		!> Notice that the order of the nodes follows that of the element itself.
		!> @param[in] atoIJK Node a to IJK rellationship
		!> @param[out] xgp Quadrature points
		!> @param[out] wgp Quadrature weights
		subroutine GaussLobattoLegendre_hex(mporder,mngaus,atoIJK,xgp,wgp)

			implicit none
			integer(4),intent(in) :: mporder,mngaus,atoIJK(mngaus)
			real(8),intent(out)  :: xgp(mngaus,ndime),wgp(mngaus)
			integer(4)            :: inode,i,j,k,lorder(mporder+1)
			real(8)               :: xi(mporder+1),w1d(mporder+1)

			call getGaussLobattoLegendre_weights_and_roots(mporder,w1d,xi)

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
			lorder(i) = i-1
			end do

			inode = 0
			do k = 1,mporder+1
				do i = 1,mporder+1
					do j = 1,mporder+1
						inode = inode + 1
						xgp(atoIJK(inode),1:3) = [xi(lorder(i)), xi(lorder(j)), xi(lorder(k))]
						wgp(atoIJK(inode)) = w1d(lorder(i))*w1d(lorder(j))*w1d(lorder(k))
					end do
				end do
			end do

		end subroutine GaussLobattoLegendre_hex

		subroutine GaussLobattoLegendre_hex_2(mporder,mngaus,xgp,wgp)

			implicit none
			integer(4),intent(in) :: mporder,mngaus
			real(8),intent(out)  :: xgp(mngaus,ndime),wgp(mngaus)
			integer(4)            :: inode,i,j,k,lorder(mporder+1)
			real(8)              :: xi(mporder+1),w1d(mporder+1),w0,w1,w2,w3
			
			call getGaussLobattoLegendre_weights_and_roots(mporder,w1d,xi)

			inode = 0
			do k = 1,mporder+1
				do i = 1,mporder+1
					do j = 1,mporder+1
						inode = inode + 1
						xgp(inode,1:3) = [xi(i), xi(j), xi(k)]
						wgp(inode) = w1d(i)*w1d(j)*w1d(k)
					end do
				end do
			end do

		end subroutine GaussLobattoLegendre_hex_2

		subroutine getGaussPoints_equispaced_hex(mporder,mngaus,atoIJK,xgp)

			implicit none
			integer(4),intent(in) :: mporder,mngaus,atoIJK(mngaus)
			real(8),intent(out)   :: xgp(mngaus,ndime)
			integer(4)            :: inode,i,j,k,lorder(mporder+1)
			real(8)               :: xi_equi(mporder+1)

			call getEquispaced_roots(mporder,xi_equi)

			!write(*,*) 'xi_equi',xi_equi

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
			lorder(i) = i-1
			end do

			inode = 0
			do k = 1,mporder+1
				do i = 1,mporder+1
					do j = 1,mporder+1
						inode = inode + 1
						xgp(atoIJK(inode),1:3) = [xi_equi(lorder(i)), xi_equi(lorder(j)), xi_equi(lorder(k))]
					end do
				end do
			end do

			!write(*,*) 'xgp',xgp(:,:)

		end subroutine getGaussPoints_equispaced_hex

end module quadrature_rules