module mod_InSitu
   use mod_mpi_mesh
   use mod_constants
   implicit none

contains

   subroutine init_InSitu()
      implicit none
      integer rank, size, err
      integer(4) :: ii
      include "mpif.h"
      call print_sensei("init sensei")

      call mpi_comm_size(MPI_COMM_WORLD,size,err)
      call mpi_comm_rank(MPI_COMM_WORLD,rank,err)
      call init_sensei(rank, size)

!     do ii = 1,numElemsRankPar*nnode
!        connecVTK(ii) = connecVTK(ii)-1
!     end do

   end subroutine init_InSitu

   subroutine end_InSitu()
      call finalize_sensei()
   end subroutine end_InSitu

   subroutine run_InSitu(u,istep)
      implicit none
      integer(4),intent(in) :: istep
      real(rp),intent(inout),dimension(numNodesRankPar,ndime) :: u
      integer(4) :: ielem,inode,inodeL
      logical, save :: first_time=.true.
      if(rp==4) then
         ! element list
         ! do ielem = 1,numElemsRankPar
         !    write(mpi_rank+100*istep+10000,'(64(i8,1x))') (connecVTK(nnode*(ielem-1)+inode),inode = 1,nnode)
         ! end do
         if(first_time) then
            call print_sensei("numElemsRankPar")
            call print_sensei_int(numElemsRankPar*(porder**3))
            call print_sensei("nnode")
            call print_sensei_int(8)
            call creategrid(coordPar, numNodesRankPar, connecVTK, numElemsRankPar*8*(porder**3))
         endif
         ! coords as array of structs
         ! do iNodeL = 1,numNodesRankPar
         !    write(mpi_rank+100*istep+20000,'(3(e14.7,1x))') coordPar(iNodeL,1:3)
         ! end do

         ! velocity per verticy
         ! do iNodeL = 1,numNodesRankPar
         !    write(mpi_rank+100*istep+30000,'(3(e14.7,1x))') u(iNodeL,1:3)
         ! end do
         call add_vector_field(u, "velocity")
         call process_sensei(istep)

      end if
      first_time = .false.
   end subroutine run_InSitu

end module mod_InSitu