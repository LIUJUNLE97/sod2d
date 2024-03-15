module mod_solver

	use mod_numerical_params
	use mod_comms
	use mod_mpi
	use mod_nvtx
	use mod_time_ops
	use mod_bc_routines
	use elem_diffu


	implicit none
      
	real(rp)  , allocatable, dimension(:,:) :: x_vars, r0_vars, p0_vars, qn_vars, v_vars, b_vars,z0_vars,z1_vars,M_vars
	real(rp)  , allocatable, dimension(:,:) :: aux_u_vars
	real(rp)  , allocatable, dimension(:) 	:: aux_Tem_vars
    logical  :: flag_cg_mem_alloc_vars=.true.
	integer(4) , parameter :: nvars = 5


	contains

		subroutine lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,R)

			implicit none

			integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
			real(rp),    intent(in)    :: Ml(npoin)
			real(rp),    intent(inout) :: R(npoin)
			integer(4)                :: ipoin

			!$acc parallel loop
			do ipoin = 1,npoin_w
				R(lpoin_w(ipoin)) = R(lpoin_w(ipoin))/Ml(lpoin_w(ipoin))
			end do
			!$acc end parallel

		end subroutine lumped_solver_scal

		subroutine lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,R)

			implicit none

			integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
			real(rp),    intent(in)    :: Ml(npoin)
			real(rp),    intent(inout) :: R(npoin,ndime)
			integer(4)                :: idime, ipoin

			! TODO: reverse this loop, might be faster on CPU
			!$acc parallel loop collapse(2)
			do ipoin = 1,npoin_w
				do idime = 1,ndime
					R(lpoin_w(ipoin),idime) = R(lpoin_w(ipoin),idime)/Ml(lpoin_w(ipoin))
				end do
			end do
			!$acc end  parallel loop

		end subroutine lumped_solver_vect
           
		subroutine conjGrad_imex(fact,igtime,save_logFile_next,noBoundaries,dt,nelem,npoin,npoin_w,nboun,numBoundsWM,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                                             dlxigp_ip,He,gpvol,Ngp,Ml,gamma_gas,Rgas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Rp0_mass,R_mass,Rp0_ener,R_ener,Rp0_mom,R_mom, &
                                             ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&                   ! Optional args
                                             listBoundsWM,wgp_b,bounorm,normalsAtNodes,u_buffer) ! Optional args

           implicit none

           logical,              intent(in)   :: noBoundaries
           integer(4),           intent(in)    :: igtime,save_logFile_next
           integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w),nboun
           real(rp)   , intent(in)    :: fact,gpvol(1,ngaus,nelem), Ngp(ngaus,nnode),dt,gamma_gas,Rgas,Cp,Prt
           real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),Ml(npoin)
           integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),             intent(inout) :: mu_fluid(npoin)
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
           integer(4)                :: ipoin, iter,ialpha,idime,ivars
           real(rp)                   :: alphaCG, betaCG,umag,rhol,El,e_int_l
           real(8)                     :: auxT1,auxT2,auxQ(2),auxQ1,auxQ2,auxB,alpha(5),alpha2(5),aux_alpha,T1,Q1(2)

           !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG begin"
          
          if (flag_cg_mem_alloc_vars .eqv. .true.) then
				allocate(x_vars(npoin,nvars), r0_vars(npoin,nvars), p0_vars(npoin,nvars), qn_vars(npoin,nvars), v_vars(npoin,nvars), b_vars(npoin,nvars),z0_vars(npoin,nvars),z1_vars(npoin,nvars),M_vars(npoin,nvars))
            	!$acc enter data create(x_vars(:,:), r0_vars(:,:), p0_vars(:,:), qn_vars(:,:), v_vars(:,:), b_vars(:,:),z0_vars(:,:),z1_vars(:,:),M_vars(:,:))
				allocate(aux_u_vars(npoin,ndime), aux_Tem_vars(npoin))
				!$acc enter data create(aux_u_vars(:,:), aux_Tem_vars(:))	
				flag_cg_mem_alloc_vars = .false.
			 end if

           !
           ! Initialize solver
           !
           call nvtxStartRange("PCG init")
            !$acc parallel loop
            do ipoin = 1,npoin
               !$acc loop seq
               do ivars = 1,nvars           
                  r0_vars(ipoin,ivars) = 0.0_rp
                  p0_vars(ipoin,ivars) = 0.0_rp
                  qn_vars(ipoin,ivars) = 0.0_rp
                  v_vars(ipoin ,ivars) = 0.0_rp
                  b_vars(ipoin ,ivars) = 0.0_rp
                  z0_vars(ipoin,ivars) = 0.0_rp
                  z1_vars(ipoin,ivars) = 0.0_rp
                  M_vars(ipoin ,ivars) = Ml(ipoin)/dt
               end do
            end do 
            !$acc end parallel loop

            !$acc parallel loop
            do ipoin = 1,npoin_w
               b_vars(lpoin_w(ipoin),1) = R_mass(lpoin_w(ipoin))
               x_vars(lpoin_w(ipoin),1) = Rp0_mass(lpoin_w(ipoin))
               b_vars(lpoin_w(ipoin),2) = R_ener(lpoin_w(ipoin))
               x_vars(lpoin_w(ipoin),2) = Rp0_ener(lpoin_w(ipoin))
               !$acc loop seq
               do idime = 1,ndime   
                  b_vars(lpoin_w(ipoin),2+idime) = R_mom(lpoin_w(ipoin)  ,idime)
                  x_vars(lpoin_w(ipoin),2+idime) = Rp0_mom(lpoin_w(ipoin),idime)
               end do
               rhol = x_vars(lpoin_w(ipoin),1)
               El = x_vars(lpoin_w(ipoin),2)
               umag = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  aux_u_vars(lpoin_w(ipoin),idime) = x_vars(lpoin_w(ipoin),2+idime)/rhol
                  umag = umag + (aux_u_vars(lpoin_w(ipoin),idime)*aux_u_vars(lpoin_w(ipoin),idime))
               end do
               e_int_l = (El/rhol)-0.5_rp*umag
               aux_Tem_vars(lpoin_w(ipoin)) = rhol*(gamma_gas-1.0_rp)*e_int_l/(rhol*Rgas)
            end do
            !$acc end parallel loop
               
            ! Real solver form here

            call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),x_vars(:,1),aux_u_vars,&
                 aux_Tem_vars,mu_fluid,mu_e,mu_sgs,Ml,qn_vars(:,1),qn_vars(:,3:nvars),qn_vars(:,2))
			
            !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG before atomic"
			   if(mpi_size.ge.2) then
               call nvtxStartRange("PCG halo")
               do ivars = 1,nvars 
                  call mpi_halo_atomic_update_real(qn_vars(:,ivars))
               end do
               call nvtxEndRange
            end if
            
            !$acc parallel loop 
            do ipoin = 1,npoin_w
               !$acc loop seq
               do ivars = 1,nvars
                  qn_vars(lpoin_w(ipoin),ivars) = x_vars(lpoin_w(ipoin) ,ivars)*Ml(lpoin_w(ipoin))+qn_vars(lpoin_w(ipoin),ivars)*fact*dt
                  r0_vars(lpoin_w(ipoin),ivars) = b_vars(lpoin_w(ipoin) ,ivars)-qn_vars(lpoin_w(ipoin),ivars) ! b-A*x0
                  z0_vars(lpoin_w(ipoin),ivars) = r0_vars(lpoin_w(ipoin),ivars)/M_vars(lpoin_w(ipoin),ivars)
                  p0_vars(lpoin_w(ipoin),ivars) = z0_vars(lpoin_w(ipoin),ivars)
              end do
              rhol = x_vars(lpoin_w(ipoin),1)
              El = p0_vars(lpoin_w(ipoin),2)
              umag = 0.0_rp
              !$acc loop seq
              do idime = 1,ndime
                 aux_u_vars(lpoin_w(ipoin),idime) = p0_vars(lpoin_w(ipoin),2+idime)/rhol
                 umag = umag + (aux_u_vars(lpoin_w(ipoin),idime)*aux_u_vars(lpoin_w(ipoin),idime))
              end do
              e_int_l = (El/rhol)-0.5_rp*umag
              aux_Tem_vars(lpoin_w(ipoin)) = rhol*(gamma_gas-1.0_rp)*e_int_l/(rhol*Rgas)
            end do
            !$acc end parallel loop

            auxT1 = 0.0d0
            !$acc parallel loop reduction(+:auxT1)
            do ipoin = 1,npoin
               !$acc loop seq
               do ivars = 1,nvars 
               auxT1 = auxT1+real(r0_vars(ipoin,ivars)*r0_vars(ipoin,ivars),8)
              end do
            end do
            !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG before allreduce"
            call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
            auxB = sqrt(auxT2)

            call nvtxEndRange

            !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG before loop"
           !
           ! Start iterations
           !
           call nvtxStartRange("PCG iters")
           do iter = 1,maxIter
              call full_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,x_vars(:,1),p0_vars(:,1),aux_u_vars,&
                 aux_Tem_vars,mu_fluid,mu_e,mu_sgs,Ml,qn_vars(:,1),qn_vars(:,3:nvars),qn_vars(:,2))              
			      if(mpi_size.ge.2) then
               	do ivars = 1,nvars 
                  call mpi_halo_atomic_update_real(qn_vars(:,ivars))
                end do              
               end if
               call nvtxStartRange("PCG qn_vars")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                 !$acc loop seq
                 do ivars = 1,nvars   
                    qn_vars(lpoin_w(ipoin),ivars) = p0_vars(lpoin_w(ipoin),ivars)*Ml(lpoin_w(ipoin))+qn_vars(lpoin_w(ipoin),ivars)*fact*dt
                end do
               end do
               !$acc end parallel loop
               call nvtxEndRange
            
              call nvtxStartRange("PCG alpha")
              auxQ1 = 0.0d0
              auxQ2 = 0.0d0
              !$acc parallel loop reduction(+:auxQ1,auxQ2) 
              do ipoin = 1,npoin_w
                  !$acc loop seq
                  do ivars = 1,nvars 
                   auxQ1 = auxQ1+real(r0_vars(lpoin_w(ipoin),ivars)*z0_vars(lpoin_w(ipoin),ivars),8) ! <s_k-1,r_k-1>
                   auxQ2 = auxQ2+real(p0_vars(lpoin_w(ipoin),ivars)*qn_vars(lpoin_w(ipoin),ivars),8) ! <s_k-1,A*s_k-1>
                 end do
              end do
              !$acc end parallel loop
              auxQ(1) = auxQ1
              auxQ(2) = auxQ2
              call MPI_Allreduce(auxQ,Q1,2,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
              alphaCG = Q1(1)/Q1(2)
              call nvtxEndRange

              call nvtxStartRange("PCG x^[n+1]")
              !$acc parallel loop
              do ipoin = 1,npoin_w
                 !$acc loop seq
                  do ivars = 1,nvars 
                     x_vars(lpoin_w(ipoin),ivars) = x_vars(lpoin_w(ipoin),ivars)+alphaCG*p0_vars(lpoin_w(ipoin),ivars) ! x_k = x_k-1 + alpha*s_k-1
                 end do 
              end do
              !$acc end parallel loop
              call nvtxEndRange

              call nvtxStartRange("PCG r^[n+1]")
              !$acc parallel loop
              do ipoin = 1,npoin_w
                  !$acc loop seq
                  do ivars = 1,nvars 
                     r0_vars(lpoin_w(ipoin),ivars) = real(real(r0_vars(lpoin_w(ipoin),ivars),8)-alphaCG*real(qn_vars(lpoin_w(ipoin),ivars),8),rp) ! b-A*p0
                     z1_vars(lpoin_w(ipoin),ivars) = z0_vars(lpoin_w(ipoin),ivars) 
                     z0_vars(lpoin_w(ipoin),ivars) = r0_vars(lpoin_w(ipoin),ivars)/M_vars(lpoin_w(ipoin),ivars) 
                  end do
              end do
              !$acc end parallel loop
              call nvtxEndRange

              auxT1 = 0.0d0
              !$acc parallel loop reduction(+:auxT1)
              do ipoin = 1,npoin
                  !$acc loop seq
                 do ivars = 1,nvars  
                  auxT1 = auxT1+real(r0_vars(ipoin,ivars)*r0_vars(ipoin,ivars),8)
                 end do
              end do

               call MPI_Allreduce(auxT1,auxT2,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

               T1 = auxT2
              !
              ! Stop cond
              !
              if (sqrt(T1) .lt. (tol*auxB)) then
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
               !$acc loop seq
               do ivars = 1,nvars 
                  p0_vars(lpoin_w(ipoin),ivars) = real(real(z0_vars(lpoin_w(ipoin),ivars),8)+betaCG*real(p0_vars(lpoin_w(ipoin),ivars),8) , rp) ! s_k = r_k+beta*s_k-1
               end do
			
			      rhol = x_vars(lpoin_w(ipoin),1)
			      El = p0_vars(lpoin_w(ipoin),2)
			      umag = 0.0_rp
			      !$acc loop seq
			      do idime = 1,ndime
			      	aux_u_vars(lpoin_w(ipoin),idime) = p0_vars(lpoin_w(ipoin),2+idime)/rhol
			      	umag = umag + (aux_u_vars(lpoin_w(ipoin),idime)*aux_u_vars(lpoin_w(ipoin),idime))
			      end do
			      e_int_l = (El/rhol)-0.5_rp*umag
			      aux_Tem_vars(lpoin_w(ipoin)) = rhol*(gamma_gas-1.0_rp)*e_int_l/(rhol*Rgas)
              end do
              !$acc end parallel loop
            
               !if(mpi_rank.eq.0) write(111,*) "--|[IMEX] CG, iters: ",iter," tol ",sqrt(T1)/auxB

              call nvtxEndRange
           end do
           call nvtxEndRange

           if (iter == maxIter) then
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[IMEX] CG, iters: ",iter," tol ",sqrt(T1)/auxB
           else
               if(igtime==save_logFile_next.and.mpi_rank.eq.0) write(111,*) "--|[IMEX] CG, iters: ",iter," tol ",sqrt(T1)/auxB
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


end module mod_solver