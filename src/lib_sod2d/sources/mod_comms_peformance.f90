module mod_comms_peformance
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms
   use mod_parTimer
   use openacc
   use mod_nvtx
   implicit none

   real(8), dimension(:), allocatable :: res_dfield
   real(4), dimension(:), allocatable :: res_ffield,aux_ffield

   real(4), dimension(:), allocatable :: xVec,yVec

contains

    subroutine init_comms_peformance(numNodesSrl,numNodesB_1r)
      implicit none
      integer,intent(in) :: numNodesSrl,numNodesB_1r
      integer :: i,j,iRank,iNodeL,iNodeGSrl,iNodeGPar,iNS,iNE
      integer :: numNodesBound_1rank,numNodesBoundRank,numNodesRankSrl
      integer, dimension(0:mpi_size-1) :: iNodeStartSrl, iNodeEndSrl
      integer, dimension(0:mpi_size-1) :: iNodeStartPar, iNodeEndPar
      integer, allocatable :: boundaryNodes(:),vecSharedBN_full(:),nodesInRank(:)
      logical :: useIntInComms=.false.,useFloatInComms=.false.,useDoubleInComms=.false.

        !1.fer la meva propia 'malla dummy'
        totalNumNodesSrl = numNodesSrl
        numNodesBound_1rank = numNodesB_1r
        numNodesBoundRank = numNodesBound_1rank*2

        call get_serialNodePartitioning(numNodesRankSrl,iNodeStartSrl,iNodeEndSrl)

        !write(*,*) 'rank[',mpi_rank,'] rankSrl ',numNodesRankSrl
        !do iRank=0,mpi_size-1
            !write(*,*) ' ##rank ', iRank , ' iNS ', iNodeStartSrl(iRank), ' iNE ', iNodeEndSrl(iRank)
        !end do

        !2. crear artificialment els nodes boundaries d'aquesta malla
        numNodesRankPar = numNodesRankSrl + numNodesBound_1rank
        allocate(boundaryNodes(numNodesBoundRank))
        boundaryNodes(:)=-1

        !write(*,*) 'rank[',mpi_rank,'] rankPar ',numNodesRankPar

        i=1
        iNS=iNodeStartSrl(mpi_rank)
        iNE=iNodeStartSrl(mpi_rank)+numNodesBound_1rank-1
        do iNodeGSrl=iNS,iNE
            !write(*,*) 'rank[',mpi_rank,'] iNodeGSrl ',iNodeGSrl,' i ',i
            boundaryNodes(i)=iNodeGSrl
            i=i+1
        end do

        if(mpi_rank.ne.(mpi_size-1)) then
            iNS=iNodeEndSrl(mpi_rank)+1
            iNE=iNodeEndSrl(mpi_rank)+numNodesBound_1rank
        else
            iNS=iNodeStartSrl(0)
            iNE=iNodeStartSrl(0)+numNodesBound_1rank-1
        end if
        do iNodeGSrl=iNS,iNE
            !write(*,*) 'rank[',mpi_rank,'] iNodeGSrl ',iNodeGSrl,' i ',i
            boundaryNodes(i)=iNodeGSrl
            i=i+1
        end do

        !3. crear llista amb id node global d'aquest rank
        allocate(nodesInRank(numNodesRankPar))
        nodesInRank(:)=-1

        i=1
        iNS=iNodeStartSrl(mpi_rank)
        iNE=iNodeStartSrl(mpi_rank)+numNodesRankSrl-1
        do iNodeGSrl=iNS,iNE
            !write(*,*) 'rank[',mpi_rank,'] iNodeGSrl ',iNodeGSrl,' i ',i
            nodesInRank(i)=iNodeGSrl
            i=i+1
        end do
        iNS=numNodesBoundRank-numNodesBound_1rank+1
        iNE=numNodesBoundRank
        do j=iNS,iNE
            iNodeGSrl=boundaryNodes(j)
            nodesInRank(i)=iNodeGSrl
            i=i+1
        end do

        !write(*,*) '#rank ', mpi_rank, ' nodesInRank ', nodesInRank(:)

        !3. fer el particionament dels nodes
        call define_parallelNodePartitioning(iNodeStartPar,iNodeEndPar)

        write(*,*) '#rank ', mpi_rank, ' nodesInRank ', numNodesRankPar, ' iNodeS ', rankNodeStart, ' iNodeE ', rankNodeEnd
        !write(*,*) 'rank[',mpi_rank,'] -> start ',iNodeStartPar(:),' end ', iNodeEndPar(:), 'tNNP ', totalNumNodesPar

        call define_mpi_boundaries_inPar(boundaryNodes,vecSharedBN_full)

        !4. fer un pseudo-reordering

        allocate(globalIdSrl(numNodesRankPar))
        allocate(globalIdPar(numNodesRankPar))

        do iNodeL=1,numNodesRankPar
               !iNodeL = iPos
               !isNodeAdded(iNodeGSrl)  = iNodeL
            iNodeGSrl = nodesInRank(iNodeL)
            iNodeGPar = iNodeL + iNodeStartPar(mpi_rank) - 1

            globalIdSrl(iNodeL) = iNodeGsrl
            globalIdPar(iNodeL) = iNodeGPar

        end do

        !5. generar esquema de comunicacio
        call generate_mpi_comm_scheme(vecSharedBN_full)

        deallocate(boundaryNodes)
        deallocate(vecSharedBN_full)

        !6. aloquetejem el vector que farem servir per testejar comms
        allocate(res_ffield(numNodesRankPar))
        allocate(aux_ffield(numNodesRankPar))
        allocate(res_dfield(numNodesRankPar))

        useIntInComms=.true.
        useFloatInComms=.true.
        useDoubleInComms=.true.
        call init_comms(useIntInComms,useFloatInComms,useDoubleInComms)

   end subroutine init_comms_peformance

   subroutine saxpy(n,y,alpha,x)
      !saxpy: y <- y+alpha*x
      implicit none
      integer, intent(in) :: n
      real(4), dimension(n), intent(inout) :: y
      real(4), intent(in) :: alpha
      real(4), dimension(n), intent(in) :: x
      integer :: i

      !$acc kernels
      do i=1,n
          y(i) = y(i) + alpha*x(i)
      end do
      !$acc end kernels
   end subroutine saxpy

   subroutine sgemv(n,y,alpha,beta,A,x)
      !sgemv: y <- beta*y+alpha*[A]*x
      implicit none
      integer, intent(in) :: n
      real(4), dimension(n), intent(inout) :: y
      real(4), intent(in) :: alpha,beta
      real(4), dimension(n,n), intent(in) :: A
      real(4), dimension(n), intent(in) :: x
      integer :: i,j

      !test order 
#if 1
      !$acc kernels
      do j=1,n !cols
        y(j) = 0 
        do i=1,n !rows
           y(i) = beta*y(i) + alpha*A(i,j)*x(j) 
        end do 
      end do
      !$acc end kernels
#else
      !$acc kernels
      do i=1,n !rows
         y(j) = 0 
         do j=1,n !cols
            y(i) = beta*y(i) + alpha*A(i,j)*x(j) 
         end do
      end do
      !$acc end kernels
#endif
   end subroutine sgemv

   subroutine do_saxpy_loop(numIters)
      implicit none
      integer, intent(in) :: numIters
      integer :: iter
      real(4) :: alpha

      res_ffield(:) = 1.
      aux_ffield(:) = 1.
      alpha = 0.

      write(*,*) 'rank',mpi_rank,' calling SAXPY LOOP!'


      do iter=1,numIters

         call nvtxStartRange("saxpy op")
         call saxpy(numNodesRankPar,res_ffield,alpha,aux_ffield)
         call nvtxEndRange

         call nvtxStartRange("comms")
         call update_and_comm_floatField(res_ffield)
         call nvtxEndRange

      end do

   end subroutine do_saxpy_loop

   subroutine do_comms(numIters)
       implicit none
       integer, intent(in) :: numIters
       integer :: iter
       type(parTimer) :: timer1,timer2
       real(8) :: wc_start,wc_end

       res_dfield(:) = 1.d0
       res_ffield(:) = 1.

#if 1
        !call timer1%init_timer()
        !call timer2%init_timer()

        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        !call timer1%start_timer()
        !wc_start = MPI_Wtime()
        !call MPI_Barrier(MPI_COMM_WORLD, mpi_err)


        do iter=1,numIters
            res_dfield(:) = 0.1
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
            call update_and_comm_doubleField(res_dfield)
            !write(*,*) 'res_dfield[',mpi_rank,'] ', res_dfield(:)
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        end do

        !call timer1%stop_timer()
        !wc_end = MPI_Wtime()
        !!WRITE(*, '(A,I0,A,F4.2,A)') '[MPI process ', mpi_rank, '] time elapsed during the job: ', wc_end - wc_start, 's.'
        !write(*, *) '[MPI process ', mpi_rank, '] time elapsed during the job: ', wc_end - wc_start, 's.'
        !wc_end = timer1%get_totalTime()
        !write(*, *) '[MPI process ', mpi_rank, '] parTimer: ', wc_end, 's.'
        !call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        call print_dtimers()

#endif
#if 1
        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

        do iter=1,numIters
            res_ffield(:) = 0.1
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
            call update_and_comm_floatField(res_ffield)
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        end do

        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        !write(*,*) 'res_ffield [',mpi_rank,'] ', res_ffield(:)
        call print_ftimers()
#endif

#if 0
        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        call init_shared_mem_window()

        do iter=1,numIters
            res_ffield(:) = 0.1
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
            call update_and_comm_shared_mem_floatField(res_ffield)

            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        end do

        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        !write(*,*) 'res_sm[',mpi_rank,'] ', res_ffield(:)
        call print_smtimers()


        call close_shared_mem_windows()
#endif

    end subroutine do_comms

    subroutine do_comms_sharedMem()
        implicit none
        integer :: iter,numIters

        call init_shared_mem_window()

        numIters = 10000

        do iter=1,numIters
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
            call update_and_comm_shared_mem_floatField(res_ffield)
            call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        end do

        call print_smtimers()


        call close_shared_mem_windows()

    end subroutine do_comms_sharedMem

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


      call nvtxStartRange("full_loop")
      !$acc data create(x(:)) copyout(y(:))
      do iter=1,numIters

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
         call MPI_Sendrecv(y, n, MPI_FLOAT, ngb_rank, tag_send, &
                           x, n, MPI_FLOAT, ngb_rank, tag_recv, &
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
         !$acc end host_data
         call nvtxEndRange

      end do
      call nvtxEndRange
      !$acc end data

      deallocate(x)
      deallocate(y)

   end subroutine test_mpi_cudaware


   subroutine do_crazy_mpi_test(N,numIters)
      implicit none
      integer,intent(in) :: N,numIters
      integer :: iter

      call alloc_vecs(N)

      call nvtxStartRange("crazy_loop")
      do iter=1,numIters

         call nvtxStartRange("loop_data")
         call do_crazy_loops(N)
         call nvtxEndRange

         call nvtxStartRange("data_transfer")
         call do_crazy_comms(N)
         call nvtxEndRange

      end do
      call nvtxEndRange

      call dealloc_vecs()

   end subroutine

   subroutine do_crazy_loops(N)
      implicit none
      integer, intent(in) :: N
      integer :: i

      !$acc parallel loop present(xVec(:),yVec(:))
      do i=1,N
         xVec(i) = mpi_rank + 0.5
         yVec(i) = 1.5
      end do

      !$acc parallel loop present(xVec(:),yVec(:))
      do i=1,N
         yVec(i) = 2.0*xVec(i)**2.+yVec(i)**2. - 2.0*xVec(i)**2.+yVec(i)**2. + (mpi_rank+1)*1.
      end do
   end subroutine do_crazy_loops

   subroutine do_crazy_comms(N)
      implicit none
      integer, intent(in) :: N
      integer :: tag_send,tag_recv,ngb_rank

      tag_send = 0
      tag_recv = tag_send
      ngb_rank = merge(1, 0, mpi_rank.eq.0)

      !$acc host_data use_device (yVec,xVec)
      call MPI_Sendrecv(yVec, N, MPI_FLOAT, ngb_rank, tag_send, &
                        xVec, N, MPI_FLOAT, ngb_rank, tag_recv, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
      !$acc end host_data
   end subroutine do_crazy_comms

   subroutine alloc_vecs(N)
      implicit none
      integer, intent(in) :: N

      allocate(xVec(N))
      allocate(yVec(N))
      !$acc enter data create(xVec(:))
      !$acc enter data create(yVec(:))
   end subroutine alloc_vecs

   subroutine dealloc_vecs()
      implicit none

      !$acc exit data delete(xVec)
      !$acc exit data delete(yVec)
      deallocate(xVec)
      deallocate(yVec)

   end subroutine dealloc_vecs






end module mod_comms_peformance
