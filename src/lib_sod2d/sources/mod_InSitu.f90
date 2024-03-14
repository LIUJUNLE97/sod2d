module mod_InSitu
   use mod_mpi_mesh
!  use mod_arrays
   implicit none

contains

   subroutine init_InSitu()
      implicit none
      integer(4) :: ii

!     do ii = 1,numElemsRankPar*nnode            This is no longer needed
!        connecVTK(ii) = connecVTK(ii)-1
!     end do
      !write(mpi_rank+100000,*) 'init_InSitu'

   end subroutine init_InSitu

   subroutine end_InSitu()
      implicit none
      !write(mpi_rank+100000,*) 'end_InSitu'

   end subroutine end_InSitu

   subroutine run_InSitu(u,istep)
      implicit none
      integer(4),intent(in) :: istep
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u
      integer(4) :: ielem,inode,inodeL

      do ielem = 1,numElemsRankPar
        !write(mpi_rank+100*istep+10000,'(64(i8,1x))') (connecVTK(nnode*(ielem-1)+inode),inode = 1,nnode)
      end do

      do iNodeL = 1,numNodesRankPar
        !write(mpi_rank+100*istep+20000,'(3(e14.7,1x))') coordPar(iNodeL,1:3)
      end do

      do iNodeL = 1,numNodesRankPar
        !write(mpi_rank+100*istep+30000,'(3(e14.7,1x))') u(iNodeL,1:3)
      end do

   end subroutine run_InSitu

end module mod_InSitu
