module mod_solver_incomp

      use mod_numerical_params
      use mod_comms
      use mod_mpi
      use mod_nvtx
      use mod_time_ops
      use mod_bc_routines
      use mod_operators

      implicit none

	   real(rp)  , allocatable, dimension(:) :: x, r0, p0, qn, v, b,z0,z1,M,x0
	   logical  :: flag_cg_mem_alloc=.true.

      contains

        subroutine conjGrad_scalar_incomp(igtime,save_logFile_next,nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,Rp0,R)

           implicit none

           integer(4),           intent(in)    :: igtime,save_logFile_next
           integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
           real(rp)   , intent(in)    :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
           real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin),Rp0(npoin)
           integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
           real(rp)   , intent(inout) :: R(npoin)
           integer(4)                :: ipoin, iter,ialpha
           real(rp)                   :: T1, alphaCG, betaCG,Q1(2)
           real(8)                     :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,alpha(5),alpha2(5),aux_alpha
          
           call nvtxStartRange("CG solver scalar")
          if (flag_cg_mem_alloc .eqv. .true.) then
				allocate(x(npoin), r0(npoin), p0(npoin), qn(npoin), v(npoin), b(npoin),z0(npoin),z1(npoin),M(npoin),x0(npoin))
            !$acc enter data create(x(:), r0(:), p0(:), qn(:), v(:), b(:),z0(:),z1(:),M(:),x0(:))
				flag_cg_mem_alloc = .false.
			 end if

           !
           ! Initialize solver
           !
           !$acc kernels
           x(:) = 0.0_rp
           r0(:) = 0.0_rp
           p0(:) = 0.0_rp
           qn(:) = 0.0_rp
           v(:) = 0.0_rp
           b(:) = 0.0_rp
           z0(:) = 0.0_rp
           z1(:) = 0.0_rp
           M(:) = Ml(:)
           !$acc end kernels
            call eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Rp0,qn)! A*x0
            !$acc parallel loop
            do ipoin = 1,npoin_w
               x0(lpoin_w(ipoin)) = Rp0(lpoin_w(ipoin))
               b(lpoin_w(ipoin)) = R(lpoin_w(ipoin)) - qn(lpoin_w(ipoin))
               x(lpoin_w(ipoin)) = 0.0_rp
            end do
            !$acc end parallel loop
               
            ! Real solver form here

            if((mpi_rank.eq.0) .and. (flag_fs_fix_pressure .eqv. .true.)) then
               b(lpoin_w(1)) = 0.0_rp
               x(lpoin_w(1)) = 0.0_rp
            end if
            call eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,x,qn)! A*x0
           !$acc parallel loop
           do ipoin = 1,npoin_w
              r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-qn(lpoin_w(ipoin)) ! b-A*x0
              z0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))/M(lpoin_w(ipoin))
              p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))
           end do
            !$acc end parallel loop


            !auxT1 = 0.0d0
            !!$acc parallel loop reduction(+:auxT1)
            !do ipoin = 1,npoin
            !   auxT1 = auxT1+real(b(ipoin)*b(ipoin),8)
            !end do
            !!$acc end parallel loop
            !call MPI_Allreduce(auxT1,auxB,1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
            auxB = 1.0_rp !sqrt(auxB)


           !
           ! Start iterations
           !
           do iter = 1,maxIter
              call nvtxStartRange("Iteration")
              call eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,p0,qn) ! A*s_k-1
              auxQ1 = 0.0d0
              auxQ2 = 0.0d0
              !$acc parallel loop reduction(+:auxQ1,auxQ2) 
              do ipoin = 1,npoin_w
                 auxQ1 = auxQ1+real(r0(lpoin_w(ipoin))*z0(lpoin_w(ipoin)),8) ! <s_k-1,r_k-1>
                 auxQ2 = auxQ2+real(p0(lpoin_w(ipoin))*qn(lpoin_w(ipoin)),8) ! <s_k-1,A*s_k-1>
              end do
              !$acc end parallel loop
              auxQ(1) = auxQ1
              auxQ(2) = auxQ2
              call MPI_Allreduce(auxQ,Q1,2,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
              if(Q1(2) .lt. 1e-10) Q1(2) = 1.0_rp
              alphaCG = real(Q1(1)/Q1(2),rp)
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+alphaCG*p0(lpoin_w(ipoin)) ! x_k = x_k-1 + alpha*s_k-1
              end do
              !$acc end parallel loop
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 r0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))-alphaCG*qn(lpoin_w(ipoin)) ! b-A*p0
                 z1(lpoin_w(ipoin)) = z0(lpoin_w(ipoin)) 
                 z0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))/M(lpoin_w(ipoin)) 
              end do
              !$acc end parallel loop
              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin
                 auxT1 = auxT1+real(r0(ipoin)*r0(ipoin),8)
              end do

               call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

               T1 = real(auxT2,rp)
              !
              ! Stop cond
              !
              if (sqrt(T1) .lt. (tol*auxB)) then
                 call nvtxEndRange
                 exit
              end if
              !
              ! Update p
              !
              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin
                 auxT1 = auxT1+real(r0(ipoin)*(z0(ipoin)-z1(ipoin)),8) ! <r_k,A*s_k-1>
              end do
              !$acc end parallel loop
              call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
              if(Q1(1) .lt. 1e-10) Q1(1) = 1.0_rp
              betaCG = real(auxT2/Q1(1),rp)
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))+betaCG*p0(lpoin_w(ipoin)) ! s_k = r_k+beta*s_k-1
              end do
              !$acc end parallel loop
              call nvtxEndRange
              !if(mpi_rank.eq.0) write(111,*) "--|[in] CG, iters: ",iter," tol ",sqrt(T1)
           end do
           if (iter == maxIter) then
              !if(mpi_rank.eq.0) write(111,*) "--| TOO MANY ITERATIONS!"
              !call nvtxEndRange
              !stop 1
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[out] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           else
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[out] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           endif
            
            !$acc kernels
            R(:) = x0(:)+x(:)
            !$acc end kernels

           call nvtxEndRange

        end subroutine conjGrad_scalar_incomp
end module mod_solver_incomp
