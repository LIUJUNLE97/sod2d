module mod_postpro

   use mod_constants

      contains

      subroutine compute_fieldDerivs(nelem,npoin,connec,lelpn,He,dNgp,rho,gradRho)

         implicit none

			integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode), lelpn(npoin)
			real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem), dNgp(ndime,nnode,ngaus)
			real(8),    intent(in)  :: rho(npoin)
			real(8),    intent(out) :: gradRho(npoin,ndime)
			integer(4)              :: ielem, idime, inode, igaus, ipoin
			real(8)                 :: gpcar(ndime,nnode), rho_e(nnode), aux1

			!$acc kernels
			gradRho(:,:) = 0.0d0
			!$acc end kernels
			!$acc parallel loop gang
			do ielem = 1,nelem
				!$acc loop vector
				do inode = 1,nnode
					rho_e(inode) = rho(connec(ielem,inode))
				end do
				!$acc loop worker
				do igaus = 1,ngaus
					!$acc loop vector collapse(2)
					do idime = 1,ndime
						do inode = 1,nnode
							gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
						end do
					end do
					!$acc loop seq
					do  idime = 1,ndime
						aux1 = 0.0d0
						!$acc loop vector
						do inode = 1,nnode
							aux1 = aux1+gpcar(idime,inode)*rho_e(inode)
						end do
						gradRho(connec(ielem,igaus),idime) = gradRho(connec(ielem,igaus),idime) + aux1
					end do
				end do
			end do
			!$acc end parallel loop gang

			!$acc parallel loop collapse(2)
			do idime = 1,ndime
				do ipoin = 1,npoin
					gradRho(ipoin,idime) = gradRho(ipoin,idime)/dble(lelpn(ipoin))
				end do
			end do
			!$acc end parallel loop

      end subroutine compute_fieldDerivs

end module mod_postpro