module mass_matrix

	use mod_numerical_params
	use mod_nvtx
	use mod_mpi
	use mod_mpi_mesh
	use mod_comms
	
	contains

		subroutine lumped_mass_spectral(nelem,npoin,connec,gpvol,Ml)
		
			implicit none
			
			integer(4), intent(in)  :: nelem,npoin, connec(nelem,nnode)
			real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
			real(rp),    intent(out) :: Ml(npoin)
			integer(4)              :: ielem, inode
			
			!$acc kernels
			Ml(:) = 0.0_rp
			!$acc end kernels
			!$acc parallel loop gang
			do ielem = 1,nelem
				!$acc loop vector
				do inode = 1,nnode
					!$acc atomic update
					Ml(connec(ielem,inode)) = Ml(connec(ielem,inode))+gpvol(1,inode,ielem)
					!$acc end atomic
				end do
			end do
			!$acc end parallel loop
			
			if(mpi_size.ge.2) then
				call nvtxStartRange("MPI_comms_mass")
				call mpi_halo_atomic_update_real(Ml)
				call nvtxEndRange
			end if
			
			!$acc update host(Ml(:))
		
		end subroutine lumped_mass_spectral
		
end module mass_matrix