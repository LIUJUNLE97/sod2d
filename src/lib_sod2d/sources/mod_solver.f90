! IMPORTANT: The GMRES rroutines are nnot to be cleaned! We havvee to come back and veriffy them.
module mod_solver

	use mod_numerical_params
	use mod_comms
	use mod_mpi
	use mod_nvtx
	use mod_time_ops
	use mod_bc_routines

	implicit none
	!-----------------------------------------------------------------------
	! Creating variables for the GMRES solver
	!-----------------------------------------------------------------------
	real(rp)                                  :: norm_bmass, norm_bener, norm_bmom(ndime),epsQ(3),epsR,epsE
	real(rp)                                  :: err_mass, err_ener, err_mom(ndime)

	real(rp)  ,allocatable, dimension(:)      :: e1_mass, e1_ener, beta_mass, beta_ener,cs_mass, cs_ener, sn_mass, sn_ener,updMass, updEner
	real(rp)  ,allocatable, dimension(:,:)    :: e1_mom,beta_mom, cs_mom, sn_mom, updMom

	real(rp)  , allocatable, dimension(:)     :: Jy_mass, Jy_ener, ymass, yener
	real(rp)  , allocatable, dimension(:)     :: Rmass_fix, Rener_fix, SDmass, SDener, SRmass, SRener,pEner,pMass
	real(rp)  , allocatable, dimension(:,:)   :: Jy_mom, ymom, Rmom_fix, SDmom, SRmom,pMom
	real(rp)  , allocatable, dimension(:,:)   :: Q_Mass, Q_Ener, H_mass, H_ener
	real(rp)  , allocatable, dimension(:,:,:) :: Q_Mom, H_mom
	logical                                   :: flag_gmres_mem_alloc=.true.
	logical                                   :: flag_gmres_mem_free=.false.
	real(rp)                                  :: eps=1e-16

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

		subroutine gmres_full(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
							atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
							noBoundaries,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,&
							rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
							gammaRK,dt,pt,bmass,bmom,bener,mass_sol,mom_sol,ener_sol,istep,gmresTol)
			implicit none
			integer(4), intent(in)    :: nelem, npoin, npoin_w, lpoin_w(npoin_w), connec(nelem,nnode)
			integer(4), intent(in)    :: atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
			real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,ngaus,nnode), He(ndime,ndime,ngaus,nelem)
			real(rp)  , intent(in)    :: gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
			real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin),gmresTol
			real(rp)  , intent(in)    :: Rgas, gamma_gas, Cp, Prt, gammaRK, dt,pt, Ml(npoin)
			real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,nnode), mu_sgs(nelem,nnode)
			real(rp)  , intent(in)    :: bmass(npoin), bmom(npoin,ndime), bener(npoin)
			integer(4), intent(in)     :: nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
			integer(4), intent(in)     :: nbnodes, lbnodes(nbnodes),lnbn_nodes(npoin)
			real(rp), intent(in)     :: normalsAtNodes(npoin,ndime)
			logical,    intent(in)   :: noBoundaries
			real(rp)  , intent(inout) :: mass_sol(npoin), mom_sol(npoin,ndime), ener_sol(npoin)
			integer(4), intent(in)    :: istep
			integer(4)                :: ik, jk, kk, ipoin, idime,npoin_w_g
			real(rp)                  :: errMax, errTol, aux,b, auxDot
			real(8)                  :: auxN(5),auxN2(5)

			b = sqrt(eps)
			errTol = gmresTol

			! Allocate the memory for the gmres solver if not yet allocated
			if (flag_gmres_mem_alloc .eqv. .true.) then
				allocate(Jy_mass(npoin), Jy_mom(npoin,ndime), Jy_ener(npoin))
				allocate(ymass(npoin), ymom(npoin,ndime), yener(npoin))
				allocate(Rmass_fix(npoin), Rmom_fix(npoin,ndime), Rener_fix(npoin))
				allocate(SRmass(npoin), SRmom(npoin,ndime), SRener(npoin))
				allocate(SDmass(npoin), SDmom(npoin,ndime),SDener(npoin),pEner(npoin),pMass(npoin),pMom(npoin,ndime))
				allocate(Q_Mass(npoin,maxIter+1), Q_Mom(npoin,ndime,maxIter+1), Q_Ener(npoin,maxIter+1))
				allocate(H_mass(maxIter+1,maxIter), H_mom(maxIter+1,maxIter,ndime), H_ener(maxIter+1,maxIter))
				allocate(e1_mass(maxIter+1), e1_ener(maxIter+1), beta_mass(maxIter+1), beta_ener(maxIter+1),cs_mass(maxIter), &
						cs_ener(maxIter), sn_mass(maxIter), sn_ener(maxIter),updMass(maxIter), updEner(maxIter))
				allocate(e1_mom(maxIter+1,ndime),beta_mom(maxIter+1,ndime), cs_mom(maxIter,ndime), sn_mom(maxIter,ndime), updMom(maxIter,ndime))
				flag_gmres_mem_alloc = .false.
			end if

			!$acc kernels
			ymass(:) = mass_sol(:)
			yener(:) = ener_sol(:)
			ymom(:,:) = mom_sol(:,:)
			!$acc end kernels

			auxN(:) = 0.0_rp
			call nvtxStartRange("GMRES: dot(bmass,bmass)")
			aux = 0.0_rp
			!$acc parallel loop reduction(+:aux)
			do ipoin = 1,npoin_w
				aux = aux + real(bmass(lpoin_w(ipoin))**2,8)
			end do
			!$acc end parallel loop
			auxN(1) = aux
			call nvtxEndRange()

			call nvtxStartRange("GMRES: dot(bener,bener)")
			aux = 0.0_rp
			!$acc parallel loop reduction(+:aux)
			do ipoin = 1,npoin_w
				aux = aux + real(bener(lpoin_w(ipoin))**2,8)
			end do
			!$acc end parallel loop
			auxN(2) = aux
			call nvtxEndRange()

			call nvtxStartRange("GMRES: dot(bmom,bmom)")
			do idime = 1,ndime
				aux = 0.0_rp
				!$acc parallel loop reduction(+:aux)
				do ipoin = 1,npoin_w
				aux = aux + real(bmom(lpoin_w(ipoin),idime)**2,8)
				end do
				!$acc end parallel loop
				auxN(idime+2) = aux
			end do
			call nvtxEndRange()

			call nvtxStartRange("GMRES: Update auxN2")
			call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
			call nvtxEndRange()

			if(auxN2(1)<epsilon(errTol)) auxN2(1) = 1.0_rp
			if(auxN2(2)<epsilon(errTol)) auxN2(2) = 1.0_rp
			if(auxN2(3)<epsilon(errTol)) auxN2(3) = 1.0_rp
			if(auxN2(4)<epsilon(errTol)) auxN2(4) = 1.0_rp
			if(auxN2(5)<epsilon(errTol)) auxN2(5) = 1.0_rp

			norm_bmass   = sqrt(real(auxN2(1),rp))
			norm_bener   = sqrt(real(auxN2(2),rp))
			norm_bmom(1) = sqrt(real(auxN2(3),rp))
			norm_bmom(2) = sqrt(real(auxN2(4),rp))
			norm_bmom(3) = sqrt(real(auxN2(5),rp))

			epsR = b
			epsE = b
			epsQ(1) = b
			epsQ(2) = b
			epsQ(3) = b

			if(istep == 1) then
				call nvtxStartRange("GMRES: form_approx_Jy initial")
				call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
				invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
				noBoundaries,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,&
				rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
				mass_sol,mom_sol,ener_sol,.true.)
				call nvtxEndRange()
			end if

			call nvtxStartRange("GMRES: callling init_gmres")
			call init_gmres(npoin,npoin_w,lpoin_w,bmass,bmom,bener,dt,pt,gammaRK)
			call nvtxEndRange()

			call nvtxStartRange("GMRES: iterations")
			outer:do ik = 1,maxIter

				! TODO: put this inside Jy and make gpu friendly
				! Verify: if this is properly done, ||Q|| = 1, so step should be

				!$acc kernels
				pMass(:) = Q_Mass(:,ik)
				pEner(:) = Q_Ener(:,ik)
				pMom(:,:) = Q_Mom(:,:,ik)
				!$acc end kernels

				!call  smooth_gmres(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
				!                  atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
				!                  rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
				!                  gammaRK,dt,pt,1)

				call nvtxStartRange("GMRES: arnoldi_iter")
				call arnoldi_iter(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
				atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
				noBoundaries,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,&
				rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
				gammaRK,dt,pt)
				call nvtxEndRange()

				call nvtxStartRange("GMRES: apply_givens_rotation")
				call apply_givens_rotation(ik)
				call nvtxEndRange()

				call nvtxStartRange("GMRES: update betas and errors")
				beta_mass(ik+1) = -sn_mass(ik)*beta_mass(ik)
				beta_mass(ik) = cs_mass(ik)*beta_mass(ik)
				err_mass = abs(beta_mass(ik+1))/norm_bmass
				beta_ener(ik+1) = -sn_ener(ik)*beta_ener(ik)
				beta_ener(ik) = cs_ener(ik)*beta_ener(ik)
				err_ener = abs(beta_ener(ik+1))/norm_bener
				do idime = 1,ndime
				beta_mom(ik+1,idime) = -sn_mom(ik,idime)*beta_mom(ik,idime)
				beta_mom(ik,idime) = cs_mom(ik,idime)*beta_mom(ik,idime)
				err_mom(idime) = abs(beta_mom(ik+1,idime))/norm_bmom(idime)
				end do

				errmax = 0.0_rp
				errmax = max(errmax,abs(err_mass))
				errmax = max(errmax,abs(err_ener))
				do idime = 1,ndime
					errmax = max(errmax,abs(err_mom(idime)))
				end do
				call nvtxEndRange()

				if (errMax .lt. errTol) then
					! Set upd* to beta_*
					call nvtxStartRange("GMRES: compute updates")
					!$acc kernels
					updMass(1:maxIter) = beta_mass(1:maxIter)
					updEner(1:maxIter) = beta_ener(1:maxIter)
					updMom(1:maxIter,1:ndime) = beta_mom(1:maxIter,1:ndime)
					!$acc end kernels
					do jk=ik,1,-1
						updMass(jk) = updMass(jk)/H_mass(jk,jk)
						updEner(jk) = updEner(jk)/H_ener(jk,jk)
						do kk = jk-1,1,-1
							updMass(kk) = updMass(kk) - H_mass(kk,jk)*updMass(jk)
							updEner(kk) = updEner(kk) - H_ener(kk,jk)*updEner(jk)
						end do
					end do
					do idime = 1,ndime
						do jk=ik,1,-1
							updMom(jk,idime) = updMom(jk,idime)/H_mom(jk,jk,idime)
							do kk = jk-1,1,-1
								updMom(kk,idime) = updMom(kk,idime) - H_mom(kk,jk,idime)*updMom(jk,idime)
							end do
						end do
					end do
					call nvtxEndRange()

					call nvtxStartRange("GMRES: update solution")
					!$acc parallel loop
					do ipoin = 1,npoin_w
						aux = 0.0_rp
						!$acc loop seq
						do kk = 1,ik
							aux = aux + Q_Mass(lpoin_w(ipoin),kk)*updMass(kk)
						end do
						mass_sol(lpoin_w(ipoin)) = mass_sol(lpoin_w(ipoin)) + aux
						aux = 0.0_rp
						!$acc loop seq
						do kk = 1,ik
							aux = aux + Q_Ener(lpoin_w(ipoin),kk)*updEner(kk)
						end do
						ener_sol(lpoin_w(ipoin)) = ener_sol(lpoin_w(ipoin)) + aux
						!$acc loop seq
						do idime = 1,ndime
							aux = 0.0_rp
							!$acc loop seq
							do kk = 1,ik
								aux = aux + Q_Mom(lpoin_w(ipoin),idime,kk)*updMom(kk,idime)
							end do
							mom_sol(lpoin_w(ipoin),idime) = mom_sol(lpoin_w(ipoin),idime) + aux
						end do
					end do
					!$acc end parallel loop
					call nvtxEndRange()
					exit outer
				end if
				!if(mpi_rank.eq.0)print*, " err ",errMax," it ",ik,' emass ',err_mass," eener ",err_ener," emom ",err_mom
			end do outer
			call nvtxEndRange()

			!if(mpi_rank.eq.0)print*, "(gmres) err ",errMax," it ",ik,' emass ',err_mass," eener ",err_ener," emom ",err_mom

			if(ik == (maxIter+1)) then
				call nvtxStartRange("GMRES: non-converged update")
				!$acc kernels
				updMass(1:maxIter) = beta_mass(1:maxIter)
				updEner(1:maxIter) = beta_ener(1:maxIter)
				updMom(1:maxIter,1:ndime) = beta_mom(1:maxIter,1:ndime)
				!$acc end kernels
				do jk=ik,1,-1
					updMass(jk) = updMass(jk)/H_mass(jk,jk)
					updEner(jk) = updEner(jk)/H_ener(jk,jk)
					do kk = jk-1,1,-1
						updMass(kk) = updMass(kk) - H_mass(kk,jk)*updMass(jk)
						updEner(kk) = updEner(kk) - H_ener(kk,jk)*updEner(jk)
					end do
				end do
				do idime = 1,ndime
					do jk=ik,1,-1
						updMom(jk,idime) = updMom(jk,idime)/H_mom(jk,jk,idime)
						do kk = jk-1,1,-1
							updMom(kk,idime) = updMom(kk,idime) - H_mom(kk,jk,idime)*updMom(jk,idime)
						end do
					end do
				end do

				!$acc parallel loop
				do ipoin = 1,npoin_w
					aux = 0.0_rp
					!$acc loop seq
					do kk = 1,maxIter
						aux = aux + Q_Mass(lpoin_w(ipoin),kk)*updMass(kk)
					end do
					mass_sol(lpoin_w(ipoin)) = mass_sol(lpoin_w(ipoin)) + aux
					aux = 0.0_rp
					!$acc loop seq
					do kk = 1,maxIter
						aux = aux + Q_Ener(lpoin_w(ipoin),kk)*updEner(kk)
					end do
					ener_sol(lpoin_w(ipoin)) = ener_sol(lpoin_w(ipoin)) + aux
					!$acc loop seq
					do idime = 1,ndime
						aux = 0.0_rp
						!$acc loop seq
						do kk = 1,maxIter
							aux = aux + Q_Mom(lpoin_w(ipoin),idime,kk)*updMom(kk,idime)
						end do
						mom_sol(lpoin_w(ipoin),idime) = mom_sol(lpoin_w(ipoin),idime) + aux
					end do
				end do
				!$acc end parallel loop
				call nvtxEndRange()
			end if

			if (flag_gmres_mem_free .eqv. .true.) then
				deallocate(Jy_mass, Jy_mom, Jy_ener)
				deallocate(ymass, ymom, yener)
				deallocate(Rmass_fix, Rmom_fix, Rener_fix)
				deallocate(SRmass, SRmom, SRener)
				deallocate(SDmass, SDmom, SDener)
				deallocate(Q_Mass, Q_Mom, Q_Ener)
			end if

		end subroutine gmres_full

		!-----------------------------------------------------------------------
		! Auxiliary routines for the gmres solver
		!-----------------------------------------------------------------------

		subroutine arnoldi_iter(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
								atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
								noBoundaries,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,&
								rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
								gammaRK,dt,pt)

			implicit none
			integer(4), intent(in)    :: ik
			integer(4), intent(in)    :: nelem, npoin, npoin_w, atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1)
			integer(4), intent(in)    :: connec(nelem,nnode), lpoin_w(npoin_w), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
			real(rp)  , intent(in)    :: gammaRK, dt,pt, gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
			real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus), He(ndime,ndime,ngaus,nelem), Ml(npoin)
			real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin),Rgas,gamma_gas, Cp, Prt
			integer(4), intent(in)     :: nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
			integer(4), intent(in)     :: nbnodes, lbnodes(nbnodes),lnbn_nodes(npoin)
			real(rp), intent(in)     :: normalsAtNodes(npoin,ndime)
			logical,    intent(in)   :: noBoundaries
			real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
			integer(4)                :: idime, jk, ipoin
			real(rp)                  :: zmass(npoin), zmom(npoin,ndime), zener(npoin)
			real(8)                  :: aux(5),aux2(5), auxDot

			! Compute the new J*Q(:,ik) arrays
			!$acc kernels
			zmass(:) = pMass(:)
			zmom (:,:) = pMom(:,:)
			zener(:) = pEner(:)
			!$acc end kernels

			call nvtxStartRange("Arnoldi: J*Q")
			call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
								invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
								noBoundaries,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,&
								rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
								zmass,zmom,zener,.false.)
			call nvtxEndRange()

			! Update Jy with the complement
			call nvtxStartRange("Arnoldi: Compute L*Q")
			!$acc parallel loop
			do ipoin = 1,npoin_w
				Q_Mass(lpoin_w(ipoin),ik+1) = zmass(lpoin_w(ipoin))/(pt) + zmass(lpoin_w(ipoin))/(gammaRK*dt) + Jy_mass(lpoin_w(ipoin))
				Q_Ener(lpoin_w(ipoin),ik+1) = zener(lpoin_w(ipoin))/(pt) + zener(lpoin_w(ipoin))/(gammaRK*dt) + Jy_ener(lpoin_w(ipoin))
				!$acc loop seq
				do idime = 1,ndime
					Q_Mom(lpoin_w(ipoin),idime,ik+1) =zmom(lpoin_w(ipoin),idime)/(pt) +  zmom(lpoin_w(ipoin),idime)/(gammaRK*dt) + Jy_mom(lpoin_w(ipoin),idime)
				end do
			end do
			!$acc end parallel loop
			call nvtxEndRange()

			! Compute the new H matrix
			! TODO: make it gpu friendly
			call nvtxStartRange("Arnoldi: iterations")
			do jk = 1,ik
				call nvtxStartRange("Arnoldi: Qjk*Qik+1")
				! New version
				aux(:) = 0.0_rp
				auxDot = 0.0_rp
				!$acc parallel loop reduction(+:auxDot)
				do ipoin = 1,npoin_w
					auxDot = auxDot + real(Q_Mass(lpoin_w(ipoin),jk)*Q_Mass(lpoin_w(ipoin),ik+1),8)
				end do
				!$acc end parallel loop
				aux(1) = auxDot
				auxDot = 0.0_rp
				!$acc parallel loop reduction(+:auxDot)
				do ipoin = 1,npoin_w
					auxDot = auxDot + real(Q_Ener(lpoin_w(ipoin),jk)*Q_Ener(lpoin_w(ipoin),ik+1),8)
				end do
				!$acc end parallel loop
				aux(2) = auxDot
				do idime = 1,ndime
					auxDot = 0.0_rp
					!$acc parallel loop reduction(+:auxDot)
					do ipoin = 1,npoin_w
						auxDot = auxDot + real(Q_Mom(lpoin_w(ipoin),idime,jk)*Q_Mom(lpoin_w(ipoin),idime,ik+1),8)
					end do
					!$acc end parallel loop
					aux(idime+2) = auxDot
				end do

				call MPI_Allreduce(aux,aux2,5,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
				call nvtxEndRange()

				call nvtxStartRange("Arnoldi: Update H and Q")
				H_mass(jk,ik) = real(aux2(1),rp)
				H_ener(jk,ik) = real(aux2(2),rp)
				!$acc parallel loop
				do ipoin = 1,npoin_w
					Q_Mass(lpoin_w(ipoin),ik+1) = Q_Mass(lpoin_w(ipoin),ik+1) - H_mass(jk,ik)*Q_Mass(lpoin_w(ipoin),jk)
					Q_Ener(lpoin_w(ipoin),ik+1) = Q_Ener(lpoin_w(ipoin),ik+1) - H_Ener(jk,ik)*Q_Ener(lpoin_w(ipoin),jk)
				end do
				!$acc end parallel loop
				do idime = 1,ndime
					H_mom(jk,ik,idime) = real(aux2(idime+2),rp)
					!$acc parallel loop
					do ipoin = 1,npoin_w
						Q_Mom(lpoin_w(ipoin),idime,ik+1) = Q_Mom(lpoin_w(ipoin),idime,ik+1) - H_mom(jk,ik,idime)*Q_Mom(lpoin_w(ipoin),idime,jk)
					end do
					!$acc end parallel loop
				end do
				call nvtxEndRange()
			end do
			call nvtxEndRange()

			! Fill H(ik+1,ik) with the norms of Q(:,ik+1)
			! New version
			call nvtxStartRange("Arnoldi: norm(Q)")
			aux(:) = 0.0_rp
			auxDot = 0.0_rp
			!$acc parallel loop reduction(+:auxDot)
			do ipoin = 1,npoin_w
				auxDot = auxDot + real(Q_Mass(lpoin_w(ipoin),ik+1)*Q_Mass(lpoin_w(ipoin),ik+1),8)
			end do
			!$acc end parallel loop
			aux(1) = auxDot
			auxDot = 0.0_rp
			!$acc parallel loop reduction(+:auxDot)
			do ipoin = 1,npoin_w
				auxDot = auxDot + real(Q_Ener(lpoin_w(ipoin),ik+1)*Q_Ener(lpoin_w(ipoin),ik+1),8)
			end do
			!$acc end parallel loop
			aux(2) = auxDot
			do idime = 1,ndime
				auxDot = 0.0_rp
				!$acc parallel loop reduction(+:auxDot)
				do ipoin = 1,npoin_w
					auxDot = auxDot + real(Q_Mom(lpoin_w(ipoin),idime,ik+1)*Q_Mom(lpoin_w(ipoin),idime,ik+1),8)
				end do
				!$acc end parallel loop
				aux(idime+2) = auxDot
			end do

			call MPI_Allreduce(aux,aux2,5,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
			call nvtxEndRange()

			call nvtxStartRange("Arnoldi: Update H")
			H_mass(ik+1,ik) = sqrt(real(aux2(1),rp))
			H_ener(ik+1,ik) = sqrt(real(aux2(2),rp))
			do idime = 1,ndime
				H_mom(ik+1,ik,idime) = sqrt(real(aux2(idime+2),rp))
			end do
			call nvtxEndRange()

			! Normalize every Q(:,ik+1)
			call nvtxStartRange("Arnoldi: Normalize Q(ik+1)")
			!$acc parallel loop
			do ipoin = 1,npoin_w
				Q_mass(lpoin_w(ipoin),ik+1) = Q_mass(lpoin_w(ipoin),ik+1)/H_mass(ik+1,ik)
				Q_ener(lpoin_w(ipoin),ik+1) = Q_ener(lpoin_w(ipoin),ik+1)/H_ener(ik+1,ik)
				!$acc loop seq
				do idime = 1,ndime
				Q_mom(lpoin_w(ipoin),idime,ik+1) = Q_mom(lpoin_w(ipoin),idime,ik+1)/H_mom(ik+1,ik,idime)
				end do
			end do
			!$acc end parallel loop
			call nvtxEndRange()

		end subroutine arnoldi_iter

		subroutine apply_givens_rotation(ik)

			implicit none
			integer(4), intent(in) :: ik
			integer(4)             :: jk, idime
			real(rp)                :: aux

			! Modify the ik_th column of H_* matrices
			do jk = 1,ik-1
				aux = cs_mass(jk)*H_mass(jk,ik) + sn_mass(jk)*H_mass(jk+1,ik)
				H_mass(jk+1,ik) = -sn_mass(jk)*H_mass(jk,ik) + cs_mass(jk)*H_mass(jk+1,ik)
				H_mass(jk,ik) = aux
				aux = cs_ener(jk)*H_ener(jk,ik) + sn_ener(jk)*H_ener(jk+1,ik)
				H_ener(jk+1,ik) = -sn_ener(jk)*H_ener(jk,ik) + cs_ener(jk)*H_ener(jk+1,ik)
				H_ener(jk,ik) = aux
				do idime = 1,ndime
					aux = cs_mom(jk,idime)*H_mom(jk,ik,idime) + sn_mom(jk,idime)*H_mom(jk+1,ik,idime)
					H_mom(jk+1,ik,idime) = -sn_mom(jk,idime)*H_mom(jk,ik,idime) + cs_mom(jk,idime)*H_mom(jk+1,ik,idime)
					H_mom(jk,ik,idime) = aux
				end do
			end do

			! Apply the Givens rotation to the ik_th column of H_* matrices
			call nvtxStartRange("Arnoldi: Givens rotation")
			call givens_rotation_full(H_mass(ik,ik),H_mass(ik+1,ik), &
									H_ener(ik,ik),H_ener(ik+1,ik), &
									H_mom(ik,ik,:),H_mom(ik+1,ik,:),ik)
			call nvtxEndRange()

			! Eliminate the ik+1_th row of H_* matrices
			aux = cs_mass(ik)*H_mass(ik,ik) + sn_mass(ik)*H_mass(ik+1,ik)
			H_mass(ik+1,ik) = -sn_mass(ik)*H_mass(ik,ik) + cs_mass(ik)*H_mass(ik+1,ik)
			H_mass(ik,ik) = aux
			aux = cs_ener(ik)*H_ener(ik,ik) + sn_ener(ik)*H_ener(ik+1,ik)
			H_ener(ik+1,ik) = -sn_ener(ik)*H_ener(ik,ik) + cs_ener(ik)*H_ener(ik+1,ik)
			H_ener(ik,ik) = aux
			do idime = 1,ndime
				aux = cs_mom(ik,idime)*H_mom(ik,ik,idime) + sn_mom(ik,idime)*H_mom(ik+1,ik,idime)
				H_mom(ik+1,ik,idime) = -sn_mom(ik,idime)*H_mom(ik,ik,idime) + cs_mom(ik,idime)*H_mom(ik+1,ik,idime)
				H_mom(ik,ik,idime) = aux
			end do

			H_mass(ik+1,ik) = 0.0_rp
			H_Ener(ik+1,ik) = 0.0_rp
			do idime = 1,ndime
				H_mom(ik+1,ik,idime) = 0.0_rp
			end do

		end subroutine apply_givens_rotation

		subroutine givens_rotation_full(v1mass,v2mass,v1ener,v2ener,v1mom,v2mom,ik)

			implicit none
			integer(4), intent(in) :: ik
			real(rp)   , intent(in) :: v1mass,v2mass,v1ener,v2ener,v1mom(ndime),v2mom(ndime)
			integer(4)             :: idime
			real(rp)                :: tmass,tener,tmom

			tmass = sqrt(v2mass**2 + v1mass**2)
			cs_mass(ik) = v1mass/tmass
			sn_mass(ik) = v2mass/tmass

			tener = sqrt(v2ener**2 + v1ener**2)
			cs_ener(ik) = v1ener/tener
			sn_ener(ik) = v2ener/tener

			do idime = 1,ndime
				tmom= sqrt(v2mom(idime)**2 + v1mom(idime)**2)
				cs_mom(ik,idime) = v1mom(idime)/tmom
				sn_mom(ik,idime) = v2mom(idime)/tmom
			end do

		end subroutine givens_rotation_full

		subroutine form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
								invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
								noBoundaries,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,&
								rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
								zmass,zmom,zener,flag_gmres_form_fix)

			use mod_comms
			use elem_convec
			use elem_diffu

			implicit none
			logical   , intent(in) :: flag_gmres_form_fix
			integer(4), intent(in) :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
			integer(4), intent(in) :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
			real(rp)  , intent(in) :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus), gpvol(1,ngaus,nelem), xgp(ngaus,ndime)
			real(rp)  , intent(in) :: He(ndime,ndime,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), Ml(npoin)
			real(rp)  , intent(in) :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin), Rgas,gamma_gas,Cp, Prt
			real(rp)  , intent(in) :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
			real(rp)  , intent(in) :: zmass(npoin), zener(npoin), zmom(npoin,ndime)
			integer(4), intent(in) :: nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
			integer(4), intent(in) :: nbnodes, lbnodes(nbnodes),lnbn_nodes(npoin)
			logical   , intent(in) :: noBoundaries
			real(rp)  , intent(in) :: normalsAtNodes(npoin,ndime)
			real(rp)               :: zu(npoin,ndime)
			real(rp)               :: auxQ(npoin,ndime),auxE(npoin),auxRho(npoin)
			integer(4)             :: idime,ipoin
			real(rp)               :: aux

			! Form the R(u^n) arrays if not formed already
			if (flag_gmres_form_fix .eqv. .true.) then
				call full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, invAtoIJK, &
				gmshAtoI, gmshAtoJ, gmshAtoK, u, q, rho, pr, E, Rmass_fix, Rmom_fix, Rener_fix)
			end if

			!$acc parallel loop
			do ipoin = 1,npoin_w
				auxRho(lpoin_w(ipoin)) = rho(lpoin_w(ipoin))+epsR*zmass(lpoin_w(ipoin))
				auxE(lpoin_w(ipoin)) = E(lpoin_w(ipoin))+epsE*zener(lpoin_w(ipoin))
				!$acc loop seq
				do idime = 1,ndime
					auxQ(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime)+epsQ(idime)*zmom(lpoin_w(ipoin),idime)
				end do
			end do
			!$acc end parallel loop
			
			call full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, invAtoIJK, &
				gmshAtoI, gmshAtoJ, gmshAtoK, u, auxQ, auxRho, pr, auxE, SRmass, SRmom, SRener)

			! Form the J*y arrays
			!$acc parallel loop
			do ipoin = 1,npoin_w
				Jy_mass(lpoin_w(ipoin)) = (SRmass(lpoin_w(ipoin)) - Rmass_fix(lpoin_w(ipoin)))/epsR
				Jy_ener(lpoin_w(ipoin)) = (SRener(lpoin_w(ipoin)) - Rener_fix(lpoin_w(ipoin)))/epsE
				!$acc loop seq
				do idime = 1,ndime
					Jy_mom(lpoin_w(ipoin),idime) = (SRmom(lpoin_w(ipoin),idime) - Rmom_fix(lpoin_w(ipoin),idime))/epsQ(idime)
				end do
			end do
			!$acc end parallel loop
			! Communicate before applying ML
			if(mpi_size.ge.2) then
				call nvtxStartRange("MPI_comms_tI")
				call mpi_halo_atomic_update_real(Jy_mass)
				call mpi_halo_atomic_update_real(Jy_ener)
				do idime = 1,ndime
				call mpi_halo_atomic_update_real(Jy_mom(:,idime)) ! TODO: need to pass the fulll aray?
				end do
				call nvtxEndRange
			end if

			! Multiply residuals by inverse Ml
			call lumped_solver_scal(npoin, npoin_w, lpoin_w, Ml, Jy_mass)
			call lumped_solver_scal(npoin, npoin_w, lpoin_w, Ml, Jy_ener)
			call lumped_solver_vect(npoin, npoin_w, lpoin_w, Ml, Jy_mom)

			if (noBoundaries .eqv. .false.) then
				call bc_fix_dirichlet_Jacobian(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Jy_mass,Jy_mom,Jy_ener)
			end if

		end subroutine form_approx_Jy

		subroutine init_gmres(npoin,npoin_w,lpoin_w,bmass,bmom,bener,dt,pt,gammaRK)

			implicit none
			integer(4), intent(in) :: npoin, npoin_w, lpoin_w(npoin_w)
			real(rp)  , intent(in) :: bmass(npoin), bmom(npoin,ndime), bener(npoin), dt,pt, gammaRK
			integer(4)             :: ipoin, idime
			real(rp)               :: amass(npoin), aener(npoin), amom(npoin,ndime)
			real(8)                :: aux(5), aux2(5), aux3

			! Zero Q_* arrays
			!$acc kernels
			Q_Mass(:,:) = 0.0_rp
			Q_Ener(:,:) = 0.0_rp
			Q_Mom(:,:,:) = 0.0_rp
			!$acc end kernels

			! Zero cs and sn arrays
			!$acc kernels
			cs_mass(:) = 0.0_rp
			cs_ener(:) = 0.0_rp
			cs_mom(:,:) = 0.0_rp
			sn_mass(:) = 0.0_rp
			sn_ener(:) = 0.0_rp
			sn_mom(:,:) = 0.0_rp
			!$acc end kernels

			! Zero H_* mmatrices
			!$acc kernels
			H_Mass(:,:) = 0.0_rp
			H_Ener(:,:) = 0.0_rp
			H_Mom(:,:,:) = 0.0_rp
			!$acc end kernels

			! Initialize e1_* arrays
			!$acc kernels
			e1_mass(:) = 0.0_rp
			e1_ener(:) = 0.0_rp
			e1_mom(:,:) = 0.0_rp
			!$acc end kernels
			e1_mass(1) = 1.0_rp
			e1_ener(1) = 1.0_rp
			e1_mom(1,:) = 1.0_rp

			! Add the remaining terms too form thhe L*y arrays
			!$acc parallel loop
			do ipoin = 1,npoin_w
				amass(lpoin_w(ipoin)) = ymass(lpoin_w(ipoin))/(pt) + ymass(lpoin_w(ipoin))/(gammaRK*dt) + Jy_mass(lpoin_w(ipoin))
				aener(lpoin_w(ipoin)) = yener(lpoin_w(ipoin))/(pt) + yener(lpoin_w(ipoin))/(gammaRK*dt) + Jy_ener(lpoin_w(ipoin))
				!$acc loop seq
				do idime = 1,ndime
					amom(lpoin_w(ipoin),idime) = ymom(lpoin_w(ipoin),idime)/(pt) + ymom(lpoin_w(ipoin),idime)/(gammaRK*dt) + Jy_mom(lpoin_w(ipoin),idime)
				end do

				Q_Mass(lpoin_w(ipoin),1) = bmass(lpoin_w(ipoin)) - amass(lpoin_w(ipoin))
				Q_Ener(lpoin_w(ipoin),1) = bener(lpoin_w(ipoin)) - aener(lpoin_w(ipoin))
				!$acc loop seq
				do idime = 1,ndime
					Q_Mom(lpoin_w(ipoin),idime,1) = bmom(lpoin_w(ipoin),idime) - amom(lpoin_w(ipoin),idime)
				end do
			end do
			!$acc end parallel loop

			! Normalize each residual
			! TODO: make it GPU compatible without using atomic update
			!$acc kernels
			aux(:) = 0.0_rp
			aux2(:) = 0.0_rp
			!$acc end kernels
			aux3 = 0.0_rp
			!$acc parallel loop reduction(+:aux3)
			do ipoin = 1,npoin_w
				aux3 = aux3 + real(Q_Mass(lpoin_w(ipoin),1)**2,8)
			end do
			!$acc end parallel loop
			aux(1) = aux3

			aux3 = 0.0_rp
			!$acc parallel loop reduction(+:aux3)
			do ipoin = 1,npoin_w
				aux3 = aux3 + real(Q_Ener(lpoin_w(ipoin),1)**2,8)
			end do
			!$acc end parallel loop
			aux(2) = aux3

			do idime = 1,ndime
				aux3 = 0.0_rp
				!$acc parallel loop reduction(+:aux3)
				do ipoin = 1,npoin_w
					aux3 = aux3 + real(Q_Mom(lpoin_w(ipoin),idime,1)**2,8)
				end do
				!$acc end parallel loop
				aux(idime+2) = aux3
			end do

			call MPI_Allreduce(aux,aux2,5,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

			if(aux2(1)<1e-10) aux2(1) = 1.0_rp
			if(aux2(2)<1e-10) aux2(2) = 1.0_rp
			if(aux2(3)<1e-10) aux2(3) = 1.0_rp
			if(aux2(4)<1e-10) aux2(4) = 1.0_rp
			if(aux2(5)<1e-10) aux2(5) = 1.0_rp

			!$acc kernels
			beta_mass(:) = sqrt(real(aux2(1),rp))*e1_mass(:)
			beta_ener(:) = sqrt(real(aux2(2),rp))*e1_ener(:)
			!$acc end kernels

			!$acc parallel loop
			do ipoin = 1,npoin_w
				Q_Mass(lpoin_w(ipoin),1) = Q_Mass(lpoin_w(ipoin),1)/sqrt(real(aux2(1),rp))
				Q_Ener(lpoin_w(ipoin),1) = Q_Ener(lpoin_w(ipoin),1)/sqrt(real(aux2(2),rp))
			end do
			!$acc end parallel loop

			do idime = 1,ndime
				!$acc kernels
				beta_mom(:,idime) = sqrt(real(aux2(idime+2),rp))*e1_mom(:,idime)
				!$acc end kernels
				!$acc parallel loop
				do ipoin = 1,npoin_w
					Q_Mom(lpoin_w(ipoin),idime,1) = Q_Mom(lpoin_w(ipoin),idime,1)/sqrt(real(aux2(idime+2),rp))
				end do
				!$acc end parallel loop
			end do
		end subroutine init_gmres

		subroutine smooth_gmres(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
								atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
								rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
								gammaRK,dt,pt,iters)

			implicit none
			integer(4), intent(in)    :: ik,iters
			integer(4), intent(in)    :: nelem, npoin, npoin_w, atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1)
			integer(4), intent(in)    :: connec(nelem,nnode), lpoin_w(npoin_w), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
			real(rp)  , intent(in)    :: gammaRK, dt, gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
			real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus), He(ndime,ndime,ngaus,nelem), Ml(npoin)
			real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin),Rgas,gamma_gas, Cp, Prt,pt
			real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
			integer(4)                :: idime, jk, ipoin,i
			real(rp)                  :: bmass(npoin), bmom(npoin,ndime), bener(npoin),relax=1.0_rp,diag
			real(8)                  :: auxN(5),auxN2(5)

			!$acc kernels
			bmass(:) = pMass(:)
			bener(:) = pEner(:)
			bmom(:,:) = pMom(:,:)
			!$acc end kernels

			diag = (1.0_rp/pt) !+ (1.0_rp/(gammaRK*dt))
			diag = 1.0_rp/diag

			!$acc parallel loop
			do ipoin = 1,npoin_w
				pMass(lpoin_w(ipoin)) = diag*(bmass(lpoin_w(ipoin)))
				pEner(lpoin_w(ipoin)) = diag*(bener(lpoin_w(ipoin)))
				!$acc loop seq
				do idime = 1,ndime
					pMom(lpoin_w(ipoin),idime) = diag*(bmom(lpoin_w(ipoin),idime))
				end do
			end do
			!$acc end parallel loop

		end subroutine smooth_gmres

end module mod_solver