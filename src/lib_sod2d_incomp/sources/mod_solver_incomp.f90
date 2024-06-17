module mod_solver_incomp

      use mod_numerical_params
      use mod_comms
      use mod_mpi
      use mod_nvtx
      use mod_time_ops
      use mod_bc_routines_incomp
      use mod_operators
      use elem_diffu_incomp




      implicit none

	   real(rp)  , allocatable, dimension(:) :: x, r0, p0, qn, v, b,z0,z1,M,x0,diag
      real(rp)  , allocatable, dimension(:,:) :: x_u, r0_u, p0_u, qn_u, v_u, b_u,z0_u,z1_u,M_u
      real(rp)  , allocatable, dimension(:,:,:) :: L,Lt
	   logical  :: flag_cg_mem_alloc_pres=.true.
      logical  :: flag_cg_mem_alloc_veloc=.true.


      contains

            subroutine conjGrad_veloc_incomp(igtime,fact,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,connec,lpoin_w,invAtoIJK,&
                                             gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,mu_fluid,mu_e,mu_sgs,Rp0,R, &
                                             bou_codes_nodes,normalsAtNodes,u_buffer) ! Optional args

            implicit none

            logical,    intent(in) :: noBoundaries
            integer(4), intent(in) :: igtime,save_logFile_next
            integer(4), intent(in) :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w),nboun
            real(rp),   intent(in) :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode),dt,fact
            real(rp),   intent(in) :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin),Rp0(npoin,ndime)
            integer(4), intent(in) :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),   intent(inout) :: mu_fluid(npoin)
            real(rp),   intent(inout) :: mu_e(nelem,ngaus)
            real(rp),   intent(inout) :: mu_sgs(nelem,ngaus)
            integer(4),optional, intent(in) :: bou_codes_nodes(npoin)
            real(rp),optional,   intent(in) :: normalsAtNodes(npoin,ndime)
            real(rp),optional,   intent(in) :: u_buffer(npoin,ndime)
            real(rp), intent(inout) :: R(npoin,ndime)
            integer(4) :: ipoin,iter,ialpha,idime
            real(rp)   :: alphaCG,betaCG
            real(8)    :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,alpha(5),alpha2(5),aux_alpha,T1,Q1(2)
          
          call nvtxStartRange("CG solver veloc")
          if (flag_cg_mem_alloc_veloc .eqv. .true.) then
				allocate(x_u(npoin,ndime), r0_u(npoin,ndime), p0_u(npoin,ndime), qn_u(npoin,ndime), v_u(npoin,ndime), b_u(npoin,ndime),z0_u(npoin,ndime),z1_u(npoin,ndime),M_u(npoin,ndime))
            !$acc enter data create(x_u(:,:), r0_u(:,:), p0_u(:,:), qn_u(:,:), v_u(:,:), b_u(:,:),z0_u(:,:),z1_u(:,:),M_u(:,:))

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
                  M_u(ipoin,idime) = Ml(ipoin)/dt
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

            call full_diffusion_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,x_u,mu_fluid,mu_e,mu_sgs,Ml,qn_u)
            if(mpi_size.ge.2) then
               call nvtxStartRange("CG_u halo")
               do idime = 1,ndime
                  call mpi_halo_atomic_update_real(qn_u(:,idime))
               end do
               call nvtxEndRange
            end if
            !$acc parallel loop
            do ipoin = 1,npoin_w
               !$acc loop seq
               do idime = 1,ndime
                  qn_u(lpoin_w(ipoin),idime) = x_u(lpoin_w(ipoin),idime)*Ml(lpoin_w(ipoin))+qn_u(lpoin_w(ipoin),idime)*fact*dt
                  r0_u(lpoin_w(ipoin),idime) = b_u(lpoin_w(ipoin),idime)-qn_u(lpoin_w(ipoin),idime) ! b-A*x0
              end do
            end do
            !$acc end parallel loop            
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_prim_residual_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,r0_u,u_buffer)
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
              call full_diffusion_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,p0_u,mu_fluid,mu_e,mu_sgs,Ml,qn_u)
              if(mpi_size.ge.2) then
               do idime = 1,ndime
                  call mpi_halo_atomic_update_real(qn_u(:,idime))
                  end do              
               end if
             !$acc parallel loop
             do ipoin = 1,npoin_w
               !$acc loop seq
               do idime = 1,ndime  
                  qn_u(lpoin_w(ipoin),idime) = p0_u(lpoin_w(ipoin),idime)*Ml(lpoin_w(ipoin))+qn_u(lpoin_w(ipoin),idime)*fact*dt
              end do
             end do
            
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
                  call temporary_bc_routine_dirichlet_prim_residual_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,r0_u,u_buffer)
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

               T1 = auxT2
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
           end do
           call nvtxEndRange

           if (iter == maxIter) then
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[veloc] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           else
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[veloc] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           endif
            
            !$acc kernels
            R(:,:) = x_u(:,:)
            !$acc end kernels

           call nvtxEndRange

        end subroutine conjGrad_veloc_incomp

        subroutine conjGrad_pressure_incomp(igtime,save_logFile_next,noBoundaries,nelem,npoin,npoin_w,connec,lpoin_w,lelpn,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,dNgp,Ml,Rp0,R,nboun,bou_codes_nodes,normalsAtNodes)

           implicit none

           logical,              intent(in)   :: noBoundaries
           integer(4),           intent(in)    :: igtime,save_logFile_next
           integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w),lelpn(npoin)
           real(rp)   , intent(in)    :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
           real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin),Rp0(npoin)
           integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
           real(rp)   , intent(inout) :: R(npoin)
           integer(4), intent(in)     :: nboun,bou_codes_nodes(npoin)
           real(rp), intent(in)     :: normalsAtNodes(npoin,ndime)
           integer(4)                :: ipoin, iter,ialpha,ielem
           real(rp)                   :: T1, alphaCG, betaCG
           real(8)                     :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,Q1(2)
          
          call nvtxStartRange("CG solver press")
          if (flag_cg_mem_alloc_pres .eqv. .true.) then
				allocate(x(npoin), r0(npoin), p0(npoin), qn(npoin), v(npoin), b(npoin),z0(npoin),z1(npoin),M(npoin),x0(npoin),diag(npoin))
            !$acc enter data create(x(:), r0(:), p0(:), qn(:), v(:), b(:),z0(:),z1(:),M(:),x0(:),diag(:))

            call eval_laplacian_diag(nelem,npoin,connec,He,dNgp,gpvol,diag)

            if(flag_cg_prec_bdc .eqv. .true.) then
               allocate(L(nnode,nnode,nelem),Lt(nnode,nnode,nelem))
               !$acc enter data create(L(:,:,:),Lt(:,:,:))
               call eval_laplacian_BDL(nelem,npoin,connec,He,dNgp,invAtoIJK,gpvol,diag,L)

               !$acc parallel loop gang 
               do ielem = 1,nelem
                  Lt(:,:,ielem) = transpose(L(:,:,ielem))
               end do
               !$acc end parallel loop
            end if

				flag_cg_mem_alloc_pres = .false.
			 end if

           !
           ! Initialize solver
           !
           call nvtxStartRange("CG_p init")
           !$acc parallel loop
           do ipoin = 1,npoin           
               x(ipoin) = 0.0_rp
               r0(ipoin) = 0.0_rp
               p0(ipoin) = 0.0_rp
               qn(ipoin) = 0.0_rp
               v(ipoin) = 0.0_rp
               b(ipoin) = 0.0_rp
               z0(ipoin) = 0.0_rp
               z1(ipoin) = 0.0_rp
               M(ipoin) = diag(ipoin) !Ml(ipoin)
               x0(ipoin) = Rp0(ipoin)
            end do

            !$acc end parallel loop
            call eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,x0,qn)! A*x0
            !$acc parallel loop
            do ipoin = 1,npoin_w
               b(lpoin_w(ipoin)) = R(lpoin_w(ipoin)) - qn(lpoin_w(ipoin))
            end do
            !$acc end parallel loop
               
            ! Real solver form here

            if(flag_fs_fix_pressure) then
               !$acc kernels
               b(inode_fix_press) = 0.0_rp
               x(inode_fix_press) = nscbc_p_inf               
               !$acc end kernels
               
            end if

            call nvtxStartRange("CG_p precond")
            call eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,x,qn)! A*x0
           !$acc parallel loop
           do ipoin = 1,npoin_w
              r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-qn(lpoin_w(ipoin)) ! b-A*x0
           end do
            !$acc end parallel loop            
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_pressure_residual_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,r0)              
            end if            
           !$acc parallel loop
           do ipoin = 1,npoin_w
              z0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))/M(lpoin_w(ipoin))
              p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))
           end do
            !$acc end parallel loop

            if(flag_cg_prec_bdc .eqv. .true.) then
               call smoother_cholesky(nelem,npoin,npoin_w,lpoin_w,lelpn,connec,r0,z0)   
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))
               end do
               !$acc end parallel loop
            end if
            call nvtxEndRange

            auxT1 = 0.0d0
            !$acc parallel loop reduction(+:auxT1)
            do ipoin = 1,npoin_w
               auxT1 = auxT1+real(b(lpoin_w(ipoin))*b(lpoin_w(ipoin)),8)
            end do

            call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

            auxB = sqrt(auxT2) 
            call nvtxEndRange

           !
           ! Start iterations
           !
           call nvtxStartRange("CG_p iters")
           do iter = 1,maxIter
              call nvtxStartRange("Iter_p")
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
              call MPI_Allreduce(auxQ,Q1,2,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
              alphaCG = real(Q1(1)/Q1(2),rp)
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+alphaCG*p0(lpoin_w(ipoin)) ! x_k = x_k-1 + alpha*s_k-1
              end do
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 r0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))-alphaCG*qn(lpoin_w(ipoin)) ! b-A*p0
              end do
              !$acc end parallel loop
              if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_pressure_residual_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,r0)              
              end if
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 z1(lpoin_w(ipoin)) = z0(lpoin_w(ipoin)) 
                 z0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))/M(lpoin_w(ipoin)) 
              end do
              !$acc end parallel loop                            
              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin_w
                 auxT1 = auxT1+real(r0(lpoin_w(ipoin))*r0(lpoin_w(ipoin)),8)
              end do
              !$acc end parallel loop

               call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

               T1 = real(auxT2,rp)
              !
              ! Stop cond
              !
              if (sqrt(T1) .lt. (tol*auxB)) then
                 call nvtxEndRange
                 exit
              end if
               
               if(flag_cg_prec_bdc .eqv. .true.) then
                  call smoother_cholesky(nelem,npoin,npoin_w,lpoin_w,lelpn,connec,r0,z0)
               endif

              !
              ! Update p
              !
              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin_w
                 auxT1 = auxT1+real(r0(lpoin_w(ipoin))*(z0(lpoin_w(ipoin))-z1(lpoin_w(ipoin))),8) ! <r_k,A*s_k-1>
              end do
              !$acc end parallel loop
              call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
              betaCG = real(auxT2/Q1(1),rp)
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))+betaCG*p0(lpoin_w(ipoin)) ! s_k = r_k+beta*s_k-1
              end do
              !$acc end parallel loop
              call nvtxEndRange
           end do
           call nvtxEndRange

           if (iter == maxIter) then
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[pres] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           else
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[pres] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           endif
            
            !$acc kernels
            R(:) = x0(:)+x(:)
            !$acc end kernels

            if(flag_fs_fix_pressure) then
               !$acc kernels
               R(inode_fix_press) = nscbc_p_inf
               !$acc end kernels
            end if

           call nvtxEndRange

        end subroutine conjGrad_pressure_incomp

        subroutine smoother_cholesky(nelem,npoin,npoin_w,lpoin_w,lelpn,connec,b,x)

           implicit none

           integer(4), intent(in)    :: nelem, npoin,npoin_w,lpoin_w(npoin),lelpn(npoin),connec(nelem,nnode)
           real(rp)   , intent(in)    :: b(npoin)
           real(rp)   , intent(inout) :: x(npoin)
           integer(4)                :: inode,ielem
           integer(4)              :: ipoin(nnode),iNodeL,ipoin_w,jnode
           real(rp)                 :: bl(nnode),xl(nnode)

           call nvtxStartRange("smoother_cholesky")
           !$acc parallel loop gang private(bl,ipoin,xl)
           do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                  bl(inode)  = b(ipoin(inode))
               end do
               xl(1) = bl(1)/L(1,1,ielem)
               !$acc loop vector 
               do inode=2,nnode
                  jnode = inode-1
                  xl(inode) = (bl(inode) - dot_product(L(inode,1:jnode,ielem),xl(1:jnode)))/L(inode,inode,ielem)
               end do

               x(ipoin(nnode)) = xl(nnode)/Lt(nnode,nnode,ielem)

               !$acc loop vector 
               do inode=nnode-1,1,-1
                  jnode = inode+1
                  x(ipoin(inode)) = (xl(inode) - dot_product(Lt(inode,jnode:nnode,ielem),x(ipoin(jnode:nnode))))/Lt(inode,inode,ielem)
               end do              

           end do

         if(mpi_size.ge.2) then
            call mpi_halo_atomic_update_real(x)
         end if

         !$acc parallel loop
         do ipoin_w = 1,npoin_w
            iNodeL=lpoin_w(ipoin_w)
            x(iNodeL) = x(iNodeL)/real(lelpn(iNodeL),rp)
         end do
         !$acc end parallel loop
         call nvtxEndRange
          
        end subroutine smoother_cholesky   
end module mod_solver_incomp
