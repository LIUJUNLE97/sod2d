module mod_solver_meshElasticity

  use mod_numerical_params
  use mod_comms
  use mod_mpi
  use mod_nvtx
  use mod_time_ops
  use mod_bc_routines_meshElasticity
  use mod_operators
  use elem_diffu_meshElasticity

  implicit none

  real(rp)  , allocatable, dimension(:) :: x, r0, p0, qn, v, b,z0,z1,M,x0,diag
  real(rp)  , allocatable, dimension(:,:) :: x_u, r0_u, p0_u, qn_u, v_u, b_u,z0_u,z1_u,M_u
  real(rp)  , allocatable, dimension(:,:,:) :: L,Lt
  real(rp)  , allocatable, dimension(:,:) :: TauPX,TauPY,TauPZ
  real(rp)  , allocatable, dimension(:) :: tau
  logical  :: flag_cg_mem_alloc_pres=.true.
  logical  :: flag_cg_mem_alloc_veloc=.true.


contains

subroutine conjGrad_meshElasticity(igtime,save_logFile_next,noBoundaries,nelem,npoin,npoin_w,nboun,connec,lpoin_w,invAtoIJK,&
  gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,helem_k,nu,E,Rp0,R, &
  bou_codes_nodes,normalsAtNodes,u_buffer,metric) ! Optional args

  implicit none

  logical,    intent(in) :: noBoundaries
  integer(4), intent(in) :: igtime,save_logFile_next
  integer(4), intent(in) :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w),nboun
  real(rp),   intent(in) :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
  real(rp),   intent(in) :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin),Rp0(npoin,ndime)
  integer(4), intent(in) :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
  real(rp),   intent(in) :: nu, E,helem_k(nelem)
  real(rp), optional, intent(in) :: metric(ndime,ndime,npoin)
  integer(4),optional, intent(in) :: bou_codes_nodes(npoin)
  real(rp),optional,   intent(in) :: normalsAtNodes(npoin,ndime)
  real(rp),optional,   intent(in) :: u_buffer(npoin,ndime)
  real(rp), intent(inout) :: R(npoin,ndime)
  integer(4) :: ipoin,iter,ialpha,idime
  real(rp)   :: alphaCG,betaCG
  real(8)    :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,alpha(5),alpha2(5),aux_alpha,Q1(2)
  !
  ! Allocate
  !
  call nvtxStartRange("CG solver veloc")
  if (flag_cg_mem_alloc_veloc .eqv. .true.) then
    allocate(x_u(npoin,ndime), r0_u(npoin,ndime), p0_u(npoin,ndime), qn_u(npoin,ndime), v_u(npoin,ndime), &
      b_u(npoin,ndime),z0_u(npoin,ndime),z1_u(npoin,ndime),M_u(npoin,ndime))
    !$acc enter data create(x_u(:,:), r0_u(:,:), p0_u(:,:), qn_u(:,:), v_u(:,:), b_u(:,:),z0_u(:,:),z1_u(:,:),M_u(:,:))

    allocate(tau(nelem),TauPX(npoin,ndime),TauPY(npoin,ndime),TauPZ(npoin,ndime))
    !$acc enter data create(tau(:), TauPX(:,:), TauPY(:,:), TauPZ(:,:))

    flag_cg_mem_alloc_veloc = .false.
  end if

 !
 ! Initialize solver
 !
 call nvtxStartRange("CG_u init")
  !$acc parallel loop
  do ipoin = 1,npoin
     !$acc loop seq
     do idime = 1,ndime
        r0_u(ipoin,idime) = 0.0_rp
        p0_u(ipoin,idime) = 0.0_rp
        qn_u(ipoin,idime) = 0.0_rp
        v_u(ipoin,idime) = 0.0_rp
        b_u(ipoin,idime) = 0.0_rp
        z0_u(ipoin,idime) = 0.0_rp
        z1_u(ipoin,idime) = 0.0_rp
        M_u(ipoin,idime) = 1.0_rp!Ml(ipoin)
     end do
  end do
  !$acc end parallel loop

  !$acc parallel loop
  do ipoin = 1,npoin_w
     !$acc loop seq
     do idime = 1,ndime
        b_u(lpoin_w(ipoin),idime) = R(lpoin_w(ipoin),idime)
        x_u(lpoin_w(ipoin),idime) = Rp0(lpoin_w(ipoin),idime)
     end do
  end do
  !$acc end parallel loop

  ! Real solver form here

  call full_diffusion_ijk_meshElasticity(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,&
    invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
    x_u,nu,E,Ml,qn_u,metric=metric)

  if(mpi_size.ge.2) then
     call nvtxStartRange("CG_u halo")
     call mpi_halo_atomic_update_real_arrays(ndime,qn_u(:,:))
     call nvtxEndRange
  end if
  !$acc parallel loop
  do ipoin = 1,npoin_w
     !$acc loop seq
     do idime = 1,ndime
        r0_u(lpoin_w(ipoin),idime) = b_u(lpoin_w(ipoin),idime)-qn_u(lpoin_w(ipoin),idime) ! b-A*x0
    end do
  end do
  !$acc end parallel loop            
  if (noBoundaries .eqv. .false.) then
     call temporary_bc_routine_dirichlet_prim_residual_meshElasticity(npoin,nboun,bou_codes_nodes,normalsAtNodes,r0_u,u_buffer)
  end if            
  !$acc parallel loop
  do ipoin = 1,npoin_w
     !$acc loop seq
     do idime = 1,ndime
        z0_u(lpoin_w(ipoin),idime) = r0_u(lpoin_w(ipoin),idime)/M_u(lpoin_w(ipoin),idime)
        p0_u(lpoin_w(ipoin),idime) = z0_u(lpoin_w(ipoin),idime)
    end do
  end do
  !$acc end parallel loop


  auxT1 = 0.0d0
  !$acc parallel loop reduction(+:auxT1)
  do ipoin = 1,npoin_w
     !$acc loop seq
    do idime = 1,ndime 
     auxT1 = auxT1+real(r0_u(lpoin_w(ipoin),idime)*r0_u(lpoin_w(ipoin),idime),8)
    end do
  end do

  call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

  auxB = sqrt(auxT2)
  call nvtxEndRange

  !
  ! Start iterations
  !
  call nvtxStartRange("CG_u iters")
  do iter = 1,maxIter
    call nvtxStartRange("Iter_u")
    call full_diffusion_ijk_meshElasticity(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,&
      invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
      p0_u,nu,E,Ml,qn_u,metric=metric)
    if(mpi_size.ge.2) then
        call mpi_halo_atomic_update_real_arrays(ndime,qn_u(:,:))
     end if
 
     auxQ1 = 0.0d0
     auxQ2 = 0.0d0
     !$acc parallel loop reduction(+:auxQ1,auxQ2)
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime
           auxQ1 = auxQ1+real(r0_u(lpoin_w(ipoin),idime)*z0_u(lpoin_w(ipoin),idime),8) ! <s_k-1,r_k-1>
           auxQ2 = auxQ2+real(p0_u(lpoin_w(ipoin),idime)*qn_u(lpoin_w(ipoin),idime),8) ! <s_k-1,A*s_k-1>
        end do
     end do
     !$acc end parallel loop
     auxQ(1) = auxQ1
     auxQ(2) = auxQ2
     call MPI_Allreduce(auxQ,Q1,2,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
     alphaCG = real(Q1(1)/Q1(2),rp)
     !$acc parallel loop
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime
           x_u(lpoin_w(ipoin),idime) = x_u(lpoin_w(ipoin),idime)+alphaCG*p0_u(lpoin_w(ipoin),idime) ! x_k = x_k-1 + alpha*s_k-1
        end do
     end do
     !$acc end parallel loop
     !$acc parallel loop
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime 
           r0_u(lpoin_w(ipoin),idime) = r0_u(lpoin_w(ipoin),idime)-alphaCG*qn_u(lpoin_w(ipoin),idime) ! b-A*p0
        end do
     end do
     !$acc end parallel loop
     if (noBoundaries .eqv. .false.) then
        call temporary_bc_routine_dirichlet_prim_residual_meshElasticity(npoin,nboun,bou_codes_nodes,normalsAtNodes,r0_u,u_buffer)
     end if
     !$acc parallel loop
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime 
           z1_u(lpoin_w(ipoin),idime) = z0_u(lpoin_w(ipoin),idime) 
           z0_u(lpoin_w(ipoin),idime) = r0_u(lpoin_w(ipoin),idime)/M_u(lpoin_w(ipoin),idime) 
        end do
     end do
     !$acc end parallel loop
     auxT1 = 0.0d0
     !$acc parallel loop reduction(+:auxT1)
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime 
           auxT1 = auxT1+real(r0_u(lpoin_w(ipoin),idime)*r0_u(lpoin_w(ipoin),idime),8)
        end do
     end do

     call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

     !
     ! Stop cond
     !
     if (sqrt(auxT2) .lt. (tol*auxB)) then
        call nvtxEndRange
        exit
     end if
     !
     ! Update p
     !
     auxT1 = 0.0d0
     !$acc parallel loop reduction(+:auxT1)
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime 
           auxT1 = auxT1+real(r0_u(lpoin_w(ipoin),idime)*(z0_u(lpoin_w(ipoin),idime)-z1_u(lpoin_w(ipoin),idime)),8) ! <r_k,A*s_k-1>
        end do
     end do
     !$acc end parallel loop
     call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
     betaCG = real(auxT2/Q1(1),rp)
     !$acc parallel loop
     do ipoin = 1,npoin_w
        !$acc loop seq
        do idime = 1,ndime
           p0_u(lpoin_w(ipoin),idime) = z0_u(lpoin_w(ipoin),idime)+betaCG*p0_u(lpoin_w(ipoin),idime) ! s_k = r_k+beta*s_k-1
        end do
     end do
     !$acc end parallel loop
     call nvtxEndRange
     if(mpi_rank.eq.0) write(111,*) "--|[veloc] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB
  end do
  call nvtxEndRange

  if (iter == maxIter) then
     if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[veloc] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB
  else
     if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[veloc] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB
  endif

  !$acc kernels
  R(:,:) = x_u(:,:)
  !$acc end kernels
 
!   print*,'displacement:  min    /     max'
!   print*,'     ',minval(x_u(:,1)),' / ',maxval(x_u(:,1))
!   print*,'     ',minval(x_u(:,2)),' / ',maxval(x_u(:,2))
!   print*,'     ',minval(x_u(:,3)),' / ',maxval(x_u(:,3))
 
  call nvtxEndRange

end subroutine conjGrad_meshElasticity

end module mod_solver_meshElasticity
