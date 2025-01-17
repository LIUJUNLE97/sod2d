module mod_solver_species

      use mod_numerical_params
      use mod_comms
      use mod_mpi
      use mod_nvtx
      use mod_time_ops
      use mod_bc_routines_species
      use elem_diffu_species
      use mod_solver
      use mod_operators
      use elem_stab_species
      use mod_bc_routines_species


      implicit none

	   real(rp)  , allocatable, dimension(:) :: x, r0, p0, qn,b,z0,z1,M
      real(rp)  , allocatable, dimension(:,:) :: grad
	   logical  :: flag_cg_mem_alloc_species=.true.


      contains

        subroutine conjGrad_species(ispc,igtime,fact,dt,save_logFile_next,noBoundaries,nelem,npoin,npoin_w,nboun,connec,lpoin_w,invAtoIJK,&
                                    gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ngp,Ml,mu_fluid,mu_e,mu_sgs,tau,Cp,Prt,rho,u,Rp0,R, &
                                    bou_codes,bound,nbnodes,lbnodes,lnbn_nodes,bou_codes_nodes,normalsAtNodes,Yk_buffer) ! Optional args
           implicit none

           logical,              intent(in)   :: noBoundaries
           integer(4),           intent(in)    :: ispc,igtime,save_logFile_next
           integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
           real(rp)   , intent(in)    :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode),fact,dt,u(npoin,ndime)
           real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin),Rp0(npoin),tau(nelem)
           integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
           real(rp)   , intent(inout) :: R(npoin)
           integer(4), intent(in)     :: nboun, nbnodes, lbnodes(*),lnbn_nodes(npoin)
           integer(4), intent(in)     :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
           real(rp), intent(in)     :: normalsAtNodes(npoin,ndime)
           real(rp), intent(in)    :: Yk_buffer(npoin,nspecies),Cp,Prt,rho(npoin),mu_e(nelem,nnode),mu_fluid(npoin),mu_sgs(nelem,nnode)
           integer(4)                :: ipoin, iter,ialpha,ielem
           real(rp)                   :: alphaCG, betaCG
           real(8)                     :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,Q1(2)

          call nvtxStartRange("CG solver press")
          if (flag_cg_mem_alloc_species .eqv. .true.) then
				allocate(x(npoin), r0(npoin), p0(npoin), qn(npoin), b(npoin),z0(npoin),z1(npoin),M(npoin),grad(npoin,ndime))
            !$acc enter data create(x(:), r0(:), p0(:), qn(:), b(:),z0(:),z1(:),M(:),grad(:,:))
				flag_cg_mem_alloc_species = .false.
			 end if

           !
           ! Initialize solver
           !
           call nvtxStartRange("CG_p init")
           !$acc parallel loop
           do ipoin = 1,npoin
               x(ipoin) = Rp0(ipoin)
               r0(ipoin) = 0.0_rp
               p0(ipoin) = 0.0_rp
               qn(ipoin) = 0.0_rp
               b(ipoin) = R(ipoin)
               z0(ipoin) = 0.0_rp
               z1(ipoin) = 0.0_rp
               M(ipoin) = Ml(ipoin)/dt
            end do
            !$acc end parallel loop

            ! Real solver form here

            call nvtxStartRange("CG_Yk precond")
            call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,x,grad,.true.)
            call species_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,x,grad,Cp,Prt,rho,tau,Ml,qn,.true.,-1.0_rp)
            call species_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho,x,mu_fluid,mu_e,mu_sgs,Ml,qn,.false.,1.0_rp)
            if(mpi_size.ge.2) then
               call nvtxStartRange("CG_Yk halo")
               call mpi_halo_atomic_update_real(qn)
               call nvtxEndRange
            end if

            !$acc parallel loop
            do ipoin = 1,npoin_w
              qn(lpoin_w(ipoin)) = x(lpoin_w(ipoin))*Ml(lpoin_w(ipoin))+qn(lpoin_w(ipoin))*fact*dt
              r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-qn(lpoin_w(ipoin)) ! b-A*x0
            end do
            !$acc end parallel loop
            if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_residual_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,r0,Yk_buffer)
            end if
            auxT1 = 0.0d0
            !$acc parallel loop reduction(+:auxT1)
            do ipoin = 1,npoin_w
              z0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))/M(lpoin_w(ipoin))
              p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))
              auxT1 = auxT1+real(r0(lpoin_w(ipoin))*r0(lpoin_w(ipoin)),8)
            end do
            !$acc end parallel loop

            call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
            auxB = sqrt(auxT2)

           !
           ! Start iterations
           !
           call nvtxStartRange("CG_p iters")
           do iter = 1,maxIter
              call nvtxStartRange("Iter_p")
              call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,p0,grad,.true.)
              !call eval_u_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,u,p0,grad,.true.)
              call species_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,p0,grad,Cp,Prt,rho,tau,Ml,qn,.true.,-1.0_rp)
              !call species_stab_2_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,p0,grad,Cp,Prt,rho,u,tau,Ml,qn,.true.,-1.0_rp)
              call species_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho,p0,mu_fluid,mu_e,mu_sgs,Ml,qn,.false.,1.0_rp)          
              if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_real(qn)
              end if
              auxQ1 = 0.0d0
              auxQ2 = 0.0d0
              !$acc parallel loop reduction(+:auxQ1,auxQ2)
              do ipoin = 1,npoin_w  
                  qn(lpoin_w(ipoin)) = p0(lpoin_w(ipoin))*Ml(lpoin_w(ipoin))+qn(lpoin_w(ipoin))*fact*dt
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
                 x(lpoin_w(ipoin))  = x(lpoin_w(ipoin))+real(alphaCG,rp)*p0(lpoin_w(ipoin)) ! x_k = x_k-1 + alpha*s_k-1
                 r0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))-real(alphaCG,rp)*qn(lpoin_w(ipoin)) ! b-A*p0
              end do
              !$acc end parallel loop   
              if (noBoundaries .eqv. .false.) then
               call temporary_bc_routine_dirichlet_residual_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,r0,Yk_buffer)
              end if
              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin_w
                 z1(lpoin_w(ipoin)) = z0(lpoin_w(ipoin)) 
                 z0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))/M(lpoin_w(ipoin)) 
                 auxT1 = auxT1+real(r0(lpoin_w(ipoin))*r0(lpoin_w(ipoin)),8)
              end do
              !$acc end parallel loop     

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
                 auxT1 = auxT1+real(r0(lpoin_w(ipoin))*(z0(lpoin_w(ipoin))-z1(lpoin_w(ipoin))),8) ! <r_k,A*s_k-1>
              end do
              !$acc end parallel loop
              call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
              betaCG = real(auxT2/Q1(1),rp)
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 p0(lpoin_w(ipoin)) = z0(lpoin_w(ipoin))+real(betaCG,rp)*p0(lpoin_w(ipoin)) ! s_k = r_k+beta*s_k-1
              end do
              !$acc end parallel loop
              call nvtxEndRange
           end do
           call nvtxEndRange

           if (iter == maxIter) then
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[species:",ispc,"] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB 
           else
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[species:",ispc,"] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB 
           endif

            !$acc kernels
            R(:) = x(:)
            !$acc end kernels

        end subroutine conjGrad_species

        end module mod_solver_species

