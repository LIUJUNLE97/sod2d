module mod_solver_imex

	use mod_numerical_params
	use mod_comms
	use mod_mpi
	use mod_nvtx
	use mod_time_ops
	use mod_bc_routines
	use elem_diffu
   use elem_stab, only : full_stab_ijk, full_proj_ijk
   use elem_stab, only : comp_tau

	implicit none
      
	real(rp)  , allocatable, dimension(:,:) :: x_vars, r0_vars, p0_vars, qn_vars, b_vars,z0_vars,z1_vars,M_vars
	real(rp)  , allocatable, dimension(:,:) :: aux_u_vars,b_stab
   real(rp)  , allocatable, dimension(:) 	:: aux_Tem_vars,tau_stab
   real(rp)  , allocatable, dimension(:,:) :: ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ
   logical  :: flag_cg_mem_alloc_vars=.true.
	integer(4) , parameter :: nvars = 5


	contains
           
		subroutine conjGrad_imex(fact,igtime,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,numBoundsWM,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                                             dlxigp_ip,He,gpvol,Ngp,Ml,helem_k,gamma_gas,Rgas,Cp,Prt,csound,mu_fluid,mu_e,mu_sgs,Rp0_mass,R_mass,Rp0_ener,R_ener,Rp0_mom,R_mom, &
                                             ndof,nbnodes,ldof,lbnodes,lnbn_nodes,bound,bou_codes,bou_codes_nodes,&                   ! Optional args
                                             listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer) ! Optional args

           implicit none

           logical,              intent(in)   :: noBoundaries
           integer(4),           intent(in)    :: igtime,save_logFile_next
           integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w),nboun,lnbn_nodes(npoin)
           real(rp)   , intent(in)    :: fact,gpvol(1,ngaus,nelem), Ngp(ngaus,nnode),dt,gamma_gas,Rgas,Cp,Prt
           real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin),helem_k(nelem)
           integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),             intent(inout) :: mu_fluid(npoin),csound(npoin)
            real(rp),             intent(inout) :: mu_e(nelem,ngaus)
            real(rp),             intent(inout) :: mu_sgs(nelem,ngaus)
            integer(4),            intent(in)    :: numBoundsWM
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: u_buffer(npoin,ndime) 
			real(rp),   intent(in)    :: Rp0_mass(npoin),Rp0_ener(npoin)
            real(rp)   , intent(inout) :: R_mass(npoin),R_ener(npoin)
		   	real(rp),   intent(in)    :: Rp0_mom(npoin,ndime)
            real(rp)   , intent(inout) :: R_mom(npoin,ndime)
           integer(4)                :: ipoin, iter,ialpha,idime,ivars, ipoinl
           real(rp)                   :: alphaCG, betaCG,umag,rhol,El,e_int_l
           real(8)                     :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,alpha(5),alpha2(5),aux_alpha,Q1(2)

           !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG begin"
          
          if (flag_cg_mem_alloc_vars .eqv. .true.) then
				allocate(x_vars(npoin,nvars), r0_vars(npoin,nvars), p0_vars(npoin,nvars), qn_vars(npoin,nvars), b_vars(npoin,nvars),z0_vars(npoin,nvars),z1_vars(npoin,nvars),M_vars(npoin,nvars))
            !$acc enter data create(x_vars(:,:), r0_vars(:,:), p0_vars(:,:), qn_vars(:,:), b_vars(:,:),z0_vars(:,:),z1_vars(:,:),M_vars(:,:))
				allocate(aux_u_vars(npoin,ndime), aux_Tem_vars(npoin),b_stab(npoin,nvars))
				!$acc enter data create(aux_u_vars(:,:), aux_Tem_vars(:),b_stab(:,:))
            allocate(ProjMass(npoin,ndime),ProjEner(npoin,ndime),ProjMX(npoin,ndime),ProjMY(npoin,ndime),ProjMZ(npoin,ndime),tau_stab(nelem))
            !$acc enter data create(ProjMass(:,:),ProjEner(:,:),ProjMX(:,:),ProjMY(:,:),ProjMZ(:,:),tau_stab(:))
				flag_cg_mem_alloc_vars = .false.
			 end if

           !
           ! Initialize solver
           !
           call nvtxStartRange("PCG init")

            !$acc parallel loop
            do ipoin = 1,npoin_w
               ipoinl = lpoin_w(ipoin)

               b_vars(ipoinl,1) = R_mass(ipoinl)
               b_vars(ipoinl,2) = R_ener(ipoinl)
               x_vars(ipoinl,1) = Rp0_mass(ipoinl)
               x_vars(ipoinl,2) = Rp0_ener(ipoinl)
               rhol = x_vars(ipoinl,1)
               El = x_vars(ipoinl,2)
               umag = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime   
                  b_vars(ipoinl,2+idime) = R_mom(ipoinl  ,idime)
                  x_vars(ipoinl,2+idime) = Rp0_mom(ipoinl,idime)
                  aux_u_vars(ipoinl,idime) = x_vars(ipoinl,2+idime)/rhol
                  umag = umag + (aux_u_vars(ipoinl,idime)*aux_u_vars(ipoinl,idime))
               end do
               e_int_l = (El/rhol)-0.5_rp*umag
               aux_Tem_vars(ipoinl) = rhol*(gamma_gas-1.0_rp)*e_int_l/(rhol*Rgas)
               !$acc loop seq
               do ivars = 1,nvars         
                  M_vars(ipoinl ,ivars) = Ml(ipoinl)/dt
               end do
            end do
            !$acc end parallel loop

            call comp_tau(nelem,npoin,connec,csound,aux_u_vars,helem_k,dt,tau_stab)
            call full_proj_ijk(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),aux_u_vars,aux_Tem_vars,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ)
            call full_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),x_vars(:,1),aux_u_vars,aux_Tem_vars,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ,tau_stab,b_stab(:,1),b_stab(:,3:nvars),b_stab(:,2),.true.,1.0_rp)

            call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),x_vars(:,1),aux_u_vars,&
                 aux_Tem_vars,mu_fluid,mu_e,mu_sgs,Ml,qn_vars(:,1),qn_vars(:,3:nvars),qn_vars(:,2))

            if(mpi_size.ge.2) then
               call nvtxStartRange("PCG halo")
               call mpi_halo_atomic_update_real_arrays_iSendiRcv(nvars,b_stab(:,:))
               call mpi_halo_atomic_update_real_arrays_iSendiRcv(nvars,qn_vars(:,:))
               call nvtxEndRange
            end if
            
             !$acc parallel loop 
            do ipoin = 1,npoin_w
               ipoinl = lpoin_w(ipoin)
               !$acc loop seq
               do ivars = 1,nvars
                  qn_vars(ipoinl,ivars) = x_vars(ipoinl ,ivars)*Ml(ipoinl)+qn_vars(ipoinl,ivars)*fact*dt
                  r0_vars(ipoinl,ivars) = b_vars(ipoinl ,ivars)-qn_vars(ipoinl,ivars)+b_stab(ipoinl,ivars)*fact*dt ! b-A*x0
              end do
            end do
            !$acc end parallel loop

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,r0_vars(:,1),r0_vars(:,3:5),r0_vars(:,2))
               call nvtxEndRange
            end if
            if(flag_force_2D) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  r0_vars(ipoin,5) =  0.0_rp
               end do
               !$acc end parallel loop
              end if

            auxT1 = 0.0d0
            !$acc parallel loop reduction(+:auxT1) 
            do ipoin = 1,npoin_w
               ipoinl = lpoin_w(ipoin)
               !$acc loop seq
               do ivars = 1,nvars
                  z0_vars(ipoinl,ivars) = r0_vars(ipoinl,ivars)/M_vars(ipoinl,ivars)
                  p0_vars(ipoinl,ivars) = z0_vars(ipoinl,ivars)
                  auxT1 = auxT1+real(r0_vars(ipoinl,ivars)*r0_vars(ipoinl,ivars),8)
              end do
              rhol = x_vars(ipoinl,1)
              El = p0_vars(ipoinl,2)
              umag = 0.0_rp
              !$acc loop seq
              do idime = 1,ndime
                 aux_u_vars(ipoinl,idime) = p0_vars(ipoinl,2+idime)/rhol
                 umag = umag + (aux_u_vars(ipoinl,idime)*aux_u_vars(ipoinl,idime))
              end do
              e_int_l = (El/rhol)-0.5_rp*umag
              aux_Tem_vars(ipoinl) = rhol*(gamma_gas-1.0_rp)*e_int_l/(rhol*Rgas)
            end do
            !$acc end parallel loop

            call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
            auxB = sqrt(auxT2)

            call nvtxEndRange

           !
           ! Start iterations
           !
           call nvtxStartRange("PCG iters")
           do iter = 1,maxIter
              !call full_proj_ijk(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,p0_vars(:,1),aux_u_vars,aux_Tem_vars,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ)
              !call full_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),p0_vars(:,1),aux_u_vars,aux_Tem_vars,Ml,ProjMass,ProjEner,ProjMX,ProjMY,ProjMZ,tau_stab,qn_vars(:,1),qn_vars(:,3:nvars),qn_vars(:,2),.true.,-1.0_rp)
              call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),p0_vars(:,1),aux_u_vars,&
                 aux_Tem_vars,mu_fluid,mu_e,mu_sgs,Ml,qn_vars(:,1),qn_vars(:,3:nvars),qn_vars(:,2))              
			      if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real_arrays_iSendiRcv(nvars,qn_vars(:,:))
               end if
               
               call nvtxStartRange("PCG qn_vars + PCG alpha")
               
               auxQ1 = 0.0d0
               auxQ2 = 0.0d0
               !$acc parallel loop reduction(+:auxQ1,auxQ2) 
               do ipoin = 1,npoin_w
                  ipoinl = lpoin_w(ipoin)
                 !$acc loop seq
                 do ivars = 1,nvars   
                    qn_vars(ipoinl,ivars) = p0_vars(ipoinl,ivars)*Ml(ipoinl)+qn_vars(ipoinl,ivars)*fact*dt
                    auxQ1 = auxQ1+real(r0_vars(ipoinl,ivars)*z0_vars(ipoinl,ivars),8) ! <s_k-1,r_k-1>
                    auxQ2 = auxQ2+real(p0_vars(ipoinl,ivars)*qn_vars(ipoinl,ivars),8) ! <s_k-1,A*s_k-1>
                end do
               end do
               !$acc end parallel loop
                        
              auxQ(1) = auxQ1
              auxQ(2) = auxQ2
              call MPI_Allreduce(auxQ,Q1,2,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
              alphaCG = Q1(1)/Q1(2)
              call nvtxEndRange

              call nvtxStartRange("PCG x^[n+1] + PCG r^[n+1]")
              !$acc parallel loop
              do ipoin = 1,npoin_w
                  ipoinl = lpoin_w(ipoin)
                 !$acc loop seq
                  do ivars = 1,nvars 
                     x_vars(ipoinl,ivars) = x_vars(ipoinl,ivars)+real(alphaCG,rp)*p0_vars(ipoinl,ivars) ! x_k = x_k-1 + alpha*s_k-1
                     r0_vars(ipoinl,ivars) = r0_vars(ipoinl,ivars)-real(alphaCG,rp)*qn_vars(ipoinl,ivars) ! b-A*p0
                  end do
              end do
              !$acc end parallel loop

              if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,r0_vars(:,1),r0_vars(:,3:5),r0_vars(:,2))
               call nvtxEndRange
              end if            
              
              if(flag_force_2D) then
               !$acc parallel loop
               do ipoin = 1,npoin
                  r0_vars(ipoin,5) =  0.0_rp
               end do
               !$acc end parallel loop
              end if

              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin_w
               ipoinl = lpoin_w(ipoin)
               !$acc loop seq
               do ivars = 1,nvars 
                  z1_vars(ipoinl,ivars) = z0_vars(ipoinl,ivars) 
                  z0_vars(ipoinl,ivars) = r0_vars(ipoinl,ivars)/M_vars(ipoinl,ivars)
                  auxT1 = auxT1+real(r0_vars(ipoinl,ivars)*r0_vars(ipoinl,ivars),8)
               end do
              end do
              !$acc end parallel loop
              call nvtxEndRange

               call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

              !
              ! Stop cond
              !
              if (sqrt(auxT2) .lt. (tol*auxB)) then
                 exit
              end if
              !
              ! Update p
              !
              call nvtxStartRange("PCG p^[k+1]")
              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin
                 !$acc loop seq
                  do ivars = 1,nvars  
                     auxT1 = auxT1+real(r0_vars(ipoin,ivars)*(z0_vars(ipoin,ivars)-z1_vars(ipoin,ivars)),8) ! <r_k,A*s_k-1>
                  end do
              end do
              !$acc end parallel loop
              call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
              betaCG = auxT2/Q1(1)
              !$acc parallel loop 
              do ipoin = 1,npoin_w
                  ipoinl = lpoin_w(ipoin)
                  !$acc loop seq
                  do ivars = 1,nvars 
                     p0_vars(ipoinl,ivars) = z0_vars(ipoinl,ivars)+real(betaCG,rp)*p0_vars(ipoinl,ivars) ! s_k = r_k+beta*s_k-1
                  end do
			
                  rhol = x_vars(ipoinl,1)
                  El = p0_vars(ipoinl,2)
                  umag = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     aux_u_vars(ipoinl,idime) = p0_vars(ipoinl,2+idime)/rhol
                     umag = umag + (aux_u_vars(ipoinl,idime)*aux_u_vars(ipoinl,idime))
                  end do
                  e_int_l = (El/rhol)-0.5_rp*umag
                  aux_Tem_vars(ipoinl) = rhol*(gamma_gas-1.0_rp)*e_int_l/(rhol*Rgas)
              end do
              !$acc end parallel loop
      
              call nvtxEndRange
           end do
           call nvtxEndRange

           if (iter == maxIter) then
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[IMEX] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB
           else
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[IMEX] CG, iters: ",iter," tol ",sqrt(auxT2)/auxB
           endif
            
			   !$acc parallel loop
            do ipoin = 1,npoin
			      R_mass(ipoin) = x_vars(ipoin,1) 
				   R_ener(ipoin) = x_vars(ipoin,2)
               !$acc loop seq
               do idime = 1,ndime   
                  R_mom(ipoin,idime) = x_vars(ipoin,2+idime)
               end do
            end do
            !$acc end parallel loop
            
            !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG end loop"
           

        end subroutine conjGrad_imex


end module mod_solver_imex