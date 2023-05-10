module mod_test_funcs
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms
   use mod_hdf5
   use mod_nvtx
#ifndef NOACC   
   use openacc
#endif
   !use cudafor

   implicit none

contains

   subroutine test_mpi_cudaware(n,numIters)
      implicit none
      integer,intent(in) :: n,numIters
      integer :: i,iter
      real(4), allocatable :: x(:),y(:)
      integer :: tag_send,tag_recv,ngb_rank

      allocate(x(n))
      allocate(y(n))

      tag_send = 0
      tag_recv = tag_send
      ngb_rank = merge(1, 0, mpi_rank.eq.0)

      if(mpi_rank.eq.0) write(*,*) 'Running test_mpi_cudaware'

      call nvtxStartRange("full_loop")
      !$acc data create(x(:)) copyout(y(:))
      do iter=1,numIters

         if(mpi_rank.eq.0) write(*,*) ' iter',iter

         call nvtxStartRange("loop_data")
         !$acc parallel loop
         do i=1,n
            x(i) = mpi_rank + 0.5
            y(i) = 1.5
         end do

         !$acc parallel loop
         do i=1,n
            y(i) = 2.0*x(i)**2.+y(i)**2. - 2.0*x(i)**2.+y(i)**2. + (mpi_rank+1)*1.
         end do
         call nvtxEndRange

         call nvtxStartRange("data_transfer")
         !$acc host_data use_device (y,x)
         call MPI_Sendrecv(y, n, mpi_datatype_real, ngb_rank, tag_send, &
                           x, n, mpi_datatype_real, ngb_rank, tag_recv, &
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
         !$acc end host_data
         call nvtxEndRange

      end do
      call nvtxEndRange
      !$acc end data

      deallocate(x)
      deallocate(y)

   end subroutine test_mpi_cudaware


end module mod_test_funcs