module mod_analysis

   use mod_numerical_params
   use mod_veclen
   use mod_nvtx
   use mod_mpi

   contains

      subroutine volAvg_EK(nelem,npoin,connec,gpvol,Ngp,rho0,rho,u,EK)

         implicit none

         integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
         real(rp),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
         real(rp),    intent(in)  :: rho0, rho(npoin), u(npoin,ndime)
         real(rp),    intent(out) :: EK
         integer(4)               :: ielem, igaus, inode
         real(rp)                 :: R1
         real(8)                  :: EK_l, EK_d

         call nvtxStartRange("Compute EK")
         EK_l = 0.0
         !$acc parallel loop gang reduction(+:EK_l) vector_length(vecLength)
         do ielem = 1,nelem
            R1 = 0.0_rp
            !$acc loop vector collapse(2) reduction(+:R1)
            do igaus = 1,ngaus
               do inode = 1,nnode
                  R1 = R1 + gpvol(1,igaus,ielem)*Ngp(igaus,inode)*0.5_rp*rho(connec(ielem,inode))* &
                     dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))
               end do
            end do
            EK_l = EK_l+real(R1,8)
         end do
         !$acc end parallel loop
         EK_l = EK_l/real(rho0*((2.0_rp*3.14159_rp)**3.0_rp),8)

         call MPI_Allreduce(EK_l,EK_d,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
         EK = real(EK_d, rp)

         call nvtxEndRange

      end subroutine volAvg_EK

      subroutine visc_dissipationRate(nelem,npoin,connec,leviCivi,rho0,mu_fluid,mu_e,u,volT,gpvol,He,dNgp,eps_S,eps_D,eps_T)

         implicit none

         integer(4), intent(in)   :: nelem, npoin, connec(nelem,nnode)
         real(rp),    intent(in)  :: leviCivi(3,3,3), rho0, mu_fluid(npoin), u(npoin,ndime), volT
         real(rp),    intent(in)  :: gpvol(1,ngaus,nelem), He(ndime,ndime,ngaus,nelem), dNgp(ndime,ngaus,nelem),mu_e(nelem,ngaus)
         real(rp),    intent(out) :: eps_S, eps_D, eps_T
         integer(4)               :: ielem, igaus, inode, idime, jdime, kdime
         real(rp)                 :: R1, R2, div2U, curl2U, alpha, aux,aux2
         real(rp)                 :: gpcar(ndime,nnode), gradU(ndime,ndime)
         real(8)                  :: eps_S_l, eps_D_l, eps_S_d, eps_D_d

         if (flag_spectralElem .ne. 1) then
            write(1,*) "--| THIS ONLY WORKS WITH SPECTRAL ELEMENTS FOR NOW!"
            error stop
         end if

         eps_S_l = 0.0
         eps_D_l = 0.0
         eps_T = 0.0_rp

         call nvtxStartRange("Compute visc_dissipationRate")
         !$acc parallel loop gang private(gradU,gpcar) reduction(+:eps_S_l,eps_D_l) vector_length(vecLength)
         do ielem = 1,nelem
            R1 = 0.0_rp
            R2 = 0.0_rp
            !$acc loop seq
            do igaus = 1,ngaus
               !
               ! Compute gpcar
               !
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                  end do
               end do
               !
               ! Compute gradient of u
               !
               !$acc loop seq
               do jdime = 1,ndime
                  !$acc loop seq
                  do kdime = 1,ndime
                     aux = 0.0_rp
                     !$acc loop vector reduction(+:aux)
                     do inode = 1,nnode
                        aux = aux + gpcar(jdime,inode)*u(connec(ielem,inode),kdime)
                     end do
                     gradU(jdime,kdime) = aux
                  end do
               end do
               !
               ! Compute divergence of u
               !
               div2U = gradU(1,1)+gradU(2,2)+gradU(3,3)
               !
               ! Squared divergence of u
               !
               div2U = div2U**2
               !
               !
               ! Compute dot vproduct of curl u and curl u
               !
               curl2U = 0.0_rp
               !!$acc loop vector collapse(3) reduction(+:curl2U)
               !$acc loop seq
               do idime = 1,ndime
                  aux2= 0.0_rp
                  !$acc loop seq
                  do jdime = 1,ndime
                     !$acc loop seq
                     do kdime = 1,ndime
                        aux2 = aux2+(leviCivi(idime,jdime,kdime)*gradU(jdime,kdime))
                     end do
                  end do
                  curl2U = curl2U + aux2*aux2
               end do
               !
               ! Compute enstrophy
               !
               R1 = R1+gpvol(1,igaus,ielem)*(mu_fluid(connec(ielem,igaus))+mu_e(ielem,igaus))*curl2U
               !
               ! Compute dilational component
               !
               R2 = R2+gpvol(1,igaus,ielem)*mu_fluid(connec(ielem,igaus))*div2U
            end do
            eps_S_l = eps_S_l+real(R1,8)
            eps_D_l = eps_D_l+real(R2,8)
         end do
         !!$acc end parallel loop
         call nvtxEndRange

         call MPI_Allreduce(eps_S_l,eps_S_d,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)
         call MPI_Allreduce(eps_D_l,eps_D_d,1,mpi_datatype_real8,MPI_SUM,app_comm,mpi_err)

         alpha = 1.0_rp/(rho0*volT)
         eps_S_d = eps_S_d*real(alpha,8)
         eps_D_d = (4.0/3.0)*eps_D_d*real(alpha,8)
         eps_T = real(eps_S_d+eps_D_d, rp)
         eps_S = real(eps_S_d, rp)
         eps_D = real(eps_D_d, rp)

      end subroutine visc_dissipationRate

      subroutine maxMach(npoin,npoin_w,lpoin_w,machno,maxmachno)

         implicit none

         integer(4), intent(in)  :: npoin
         integer(4), intent(in)  :: npoin_w
         integer(4), intent(in)  :: lpoin_w(npoin_w)
         real(rp),   intent(in)  :: machno(npoin)
         real(rp),   intent(out) :: maxmachno
         real(rp)                :: maxmachno_l
         integer(4)              :: ipoin

         maxmachno = 0.0_rp
         maxmachno_l = 0.0_rp
         !$acc parallel loop reduction(max:maxmachno_l)
         do ipoin = 1,npoin_w
            maxmachno_l = max(maxmachno_l,machno(lpoin_w(ipoin)))
         end do
         !$acc end parallel loop

         call MPI_Allreduce(maxmachno_l,maxmachno,1,mpi_datatype_real,MPI_MAX,app_comm,mpi_err)

      end subroutine maxMach

      subroutine write_EK(time,EK,eps_S,eps_D,eps_T,maxmachno)

         implicit none

         real(rp), intent(in) :: time, EK, eps_S, eps_D, eps_T, maxmachno

         if(mpi_rank.eq.0) then
            write(666,10) time, EK, eps_S, eps_D, eps_T, maxmachno
            10 format(6(F12.8,2X))
         end if

      end subroutine write_EK

		!> @brife Compute the area of a particular physical boundary
		!> @details Using the boundary elements that compose the selected surface, the
		!> routine computes the area as the sum of all int(|n|,dxi,deta), where the normal
		!> of each element has been pre-computed and is stored into bounorm.
		!> @param[in] nbound The number of boundary elements
		!> @param[in] surfCode The surface code
		!> @param[in] bou_code Matching of surface codes to boundary elements
		!> @param[in] bounorm The normal of each boundary element at each reference node
		!> @param[out] surfArea The area of the selected surface
      subroutine surfInfo(iter,time,nelem,npoin,nbound,surfCode,connec,bound,point2elem,bou_code, &
                          bounorm,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coord, &
                          mu_fluid,mu_e,mu_sgs,rho,u,pr,surfArea,Fpr,Ftau)

         implicit none

			integer(4), intent(in)  :: iter, npoin, nbound, surfCode, bound(nbound,npbou), bou_code(nbound)
			integer(4), intent(in)  :: nelem, connec(nelem,nnode), point2elem(npoin)
			integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
         real(rp),    intent(in)  :: wgp_b(npbou), bounorm(nbound,ndime*npbou), time
         real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime), pr(npoin)
         real(rp),    intent(in)  :: mu_e(nelem,ngaus), mu_sgs(nelem,ngaus), mu_fluid(npoin)
         real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1),coord(npoin,ndime)
			real(rp),    intent(out) :: surfArea, Fpr(ndime), Ftau(ndime)
			real(rp)                :: surfArea_l,surfArea_s,Fpr_l(ndime),Ftau_l(ndime)
         integer(4)              :: ibound, idime, igaus, ipbou, ielem, jgaus
         integer(4)              :: numBelem, counter, isoI, isoJ, isoK, ii, jdime, kdime
         integer(4), allocatable :: lelbo(:)
         real(rp)                 :: bnorm(npbou*ndime), nmag, prl(npbou), ul(nnode,ndime), rhol(nnode)
         real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
         real(rp)                :: mu_fgp, mu_egp, mufluidl(nnode),face2centoid(ndime),sig,aux(ndime)

			! Create lelbo for the surface, where lelbo is a list of boundary elements belonging to that surface
			numBelem = 0
         !$acc parallel loop reduction(+:numBelem)
			do ibound = 1, nbound
				if (bou_code(ibound) == surfCode) then
					numBelem = numBelem + 1
				end if
			end do
         !$acc end parallel loop
			counter = 0
			allocate(lelbo(numBelem))
			do ibound = 1, nbound
				if (bou_code(ibound) == surfCode) then
					counter = counter + 1
					lelbo(counter) = ibound!bou_code(ibound,1)
				end if
			end do

			! Compute surface information through sum of Gaussian quadratures over all boundary elements
			surfArea_l = 0.0_rp
         surfArea = 0.0_rp
         !$acc kernels
         Fpr_l(:) = 0.0_rp
         Ftau_l(:) = 0.0_rp
         Fpr(:) = 0.0_rp
         Ftau(:) = 0.0_rp
         !$acc end kernels
         !$acc parallel loop gang private(bnorm,prl) reduction(+:surfArea_l)
			do ibound = 1, numBelem
				bnorm(1:npbou*ndime) = bounorm(lelbo(ibound),1:npbou*ndime)
            prl(1:npbou) = pr(bound(lelbo(ibound),1:npbou))
            surfArea_s = 0.0_rp
				! Element area
            !$acc loop vector private(tau,ul,rhol,mufluidl,gradIsoU,gradU,face2centoid) reduction(+:surfArea_s)
				do igaus = 1,npbou
               ielem = point2elem(bound(lelbo(ibound),igaus))
               jgaus = minloc(abs(connec(ielem,:)-bound(lelbo(ibound),igaus)),1)
              ! !$acc loop seq
              ! do ii=1, ngaus
              !    if(connec(ielem,ii) .eq. bound(lelbo(ibound),igaus)) jgaus=ii
              ! end do
               ul(1:nnode,1:ndime) = u(connec(ielem,:),1:ndime)
               rhol(1:nnode) = rho(connec(ielem,:))
               mufluidl(1:nnode) = mu_fluid(connec(ielem,1:nnode))
               mu_fgp = mufluidl(jgaus)+rhol(jgaus)*mu_sgs(ielem,jgaus)
               mu_egp = mu_e(ielem,jgaus)
               isoI = gmshAtoI(jgaus)
               isoJ = gmshAtoJ(jgaus)
               isoK = gmshAtoK(jgaus)

               sig=1.0_rp
               aux(1) = bnorm((igaus-1)*ndime+1)
               aux(2) = bnorm((igaus-1)*ndime+2)
               aux(3) = bnorm((igaus-1)*ndime+3)
               if(dot_product(coord(connec(ielem,nnode),:)-coord(connec(ielem,jgaus),:), aux(:)) .lt. 0.0_rp ) then
                  sig=-1.0_rp
               end if

               gradIsoU(:,:) = 0.0_rp
               !$acc loop seq
               do ii = 1,porder+1
                  !$acc loop seq
                  do idime = 1,ndime
                     gradIsoU(idime,1) = gradIsoU(idime,1)+dlxigp_ip(jgaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                     gradIsoU(idime,2) = gradIsoU(idime,2)+dlxigp_ip(jgaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                     gradIsoU(idime,3) = gradIsoU(idime,3)+dlxigp_ip(jgaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
                  end do
               end do
               gradU(:,:) = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     !$acc loop seq
                     do kdime = 1,ndime
                        gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,jgaus,ielem)*gradIsoU(idime,kdime)
                     end do
                  end do
               end do
               divU = gradU(1,1)+gradU(2,2)+gradU(3,3)
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
                  end do
                  tau(idime,idime) = tau(idime,idime) - (mu_fgp)*(2.0_rp/3.0_rp)*divU
               end do
					nmag = 0.0_rp
               !$acc loop seq
					do idime = 1,ndime
						nmag = nmag + bnorm((igaus-1)*ndime+idime)*bnorm((igaus-1)*ndime+idime)
					end do
               nmag = sqrt(nmag)
               !$acc loop seq
					do idime = 1,ndime
                  !$acc atomic update
                  Fpr_l(idime) = Fpr_l(idime)-wgp_b(igaus)*(prl(igaus)-nscbc_p_inf)*bnorm((igaus-1)*ndime+idime)*sig
                  !$acc end atomic
					end do
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     !$acc atomic update
                     Ftau_l(idime) = Ftau_l(idime)+wgp_b(igaus)*tau(idime,jdime)*bnorm((igaus-1)*ndime+jdime)*sig
                     !$acc end atomic
                  end do
               end do
					surfArea_s = surfArea_s + nmag*wgp_b(igaus)
				end do
            surfArea_l = surfArea_l + surfArea_s
			end do
         !$acc end parallel loop
			deallocate(lelbo)

         call MPI_Allreduce(surfArea_l,surfArea,1,mpi_datatype_real,MPI_SUM,app_comm,mpi_err)
         call MPI_Allreduce(Fpr_l,Fpr,ndime,mpi_datatype_real,MPI_SUM,app_comm,mpi_err)
         call MPI_Allreduce(Ftau_l,Ftau,ndime,mpi_datatype_real,MPI_SUM,app_comm,mpi_err)

         !write(*,*) '(',mpi_rank,')surfCode',surfCode,'t',time, ",", dble(surfArea), ",", dble(Fpr(1)), ",", dble(Fpr(2)), ",", dble(Fpr(3)), ",", dble(Ftau(1)), ",", dble(Ftau(2)), ",", dble(Ftau(3))

         if(mpi_rank.eq.0) then
            write(888+surfCode,"(I8,1X,A)",ADVANCE="NO") iter, ","
            write(888+surfCode,50) time, ",", dble(surfArea), ",", dble(Fpr(1)), ",", dble(Fpr(2)), ",", dble(Fpr(3)), ",", dble(Ftau(1)), ",", dble(Ftau(2)), ",", dble(Ftau(3)), ""
            !50 format(16(1X,E10.4,1X,A))
            50 format(16(1X,F16.8,1X,A))
         end if
      end subroutine surfInfo

      subroutine twInfo(nelem,npoin,nbound,surfCode,connec,bound,point2elem,bou_code, &
                          bounorm,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,dlxigp_ip,He,coord, &
                          mu_fluid,mu_e,mu_sgs,rho,u,Ftau_neg,Ftau_pos)

         implicit none

			integer(4), intent(in)  :: npoin, nbound, surfCode, bound(nbound,npbou), bou_code(nbound)
			integer(4), intent(in)  :: nelem, connec(nelem,nnode), point2elem(npoin)
			integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
         real(rp),    intent(in)  :: wgp_b(npbou), bounorm(nbound,ndime*npbou)
         real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime)
         real(rp),    intent(in)  :: mu_e(nelem,ngaus), mu_sgs(nelem,ngaus), mu_fluid(npoin)
         real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1),coord(npoin,ndime)
			real(rp),    intent(out) :: Ftau_neg(ndime), Ftau_pos(ndime)
			real(rp)                :: Ftau_n_l(ndime),Ftau_p_l(ndime)
         integer(4)              :: ibound, idime, igaus, ipbou, ielem, jgaus
         integer(4)              :: numBelem, counter, isoI, isoJ, isoK, ii, jdime, kdime
         integer(4), allocatable :: lelbo(:)
         real(rp)                 :: bnorm(npbou*ndime), nmag, ul(nnode,ndime), rhol(nnode)
         real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
         real(rp)                :: mu_fgp, mu_egp, mufluidl(nnode),face2centoid(ndime),sig,aux(ndime),tau_aux

			! Create lelbo for the surface, where lelbo is a list of boundary elements belonging to that surface
			numBelem = 0
         !$acc parallel loop reduction(+:numBelem)
			do ibound = 1, nbound
				if (bou_code(ibound) == surfCode) then
					numBelem = numBelem + 1
				end if
			end do
         !$acc end parallel loop
			counter = 0
			allocate(lelbo(numBelem))
			do ibound = 1, nbound
				if (bou_code(ibound) == surfCode) then
					counter = counter + 1
					lelbo(counter) = ibound!bou_code(ibound,1)
				end if
			end do

         !$acc kernels
         Ftau_n_l(:) = 0.0_rp
         Ftau_p_l(:) = 0.0_rp
         Ftau_neg(:) = 0.0_rp
         Ftau_pos(:) = 0.0_rp
         !$acc end kernels
         !$acc parallel loop gang private(bnorm)
			do ibound = 1, numBelem
				bnorm(1:npbou*ndime) = bounorm(lelbo(ibound),1:npbou*ndime)
            !$acc loop vector private(tau,ul,rhol,mufluidl,gradIsoU,gradU,face2centoid)
				do igaus = 1,npbou
               ielem = point2elem(bound(lelbo(ibound),igaus))
               jgaus = minloc(abs(connec(ielem,:)-bound(lelbo(ibound),igaus)),1)
               ul(1:nnode,1:ndime) = u(connec(ielem,:),1:ndime)
               rhol(1:nnode) = rho(connec(ielem,:))
               mufluidl(1:nnode) = mu_fluid(connec(ielem,1:nnode))
               mu_fgp = mufluidl(jgaus)+rhol(jgaus)*mu_sgs(ielem,jgaus)
               mu_egp = mu_e(ielem,jgaus)
               isoI = gmshAtoI(jgaus)
               isoJ = gmshAtoJ(jgaus)
               isoK = gmshAtoK(jgaus)

               sig=1.0_rp
               aux(1) = bnorm((igaus-1)*ndime+1)
               aux(2) = bnorm((igaus-1)*ndime+2)
               aux(3) = bnorm((igaus-1)*ndime+3)
               if(dot_product(coord(connec(ielem,nnode),:)-coord(connec(ielem,jgaus),:), aux(:)) .lt. 0.0_rp ) then
                  sig=-1.0_rp
               end if

               gradIsoU(:,:) = 0.0_rp
               !$acc loop seq
               do ii = 1,porder+1
                  !$acc loop seq
                  do idime = 1,ndime
                     gradIsoU(idime,1) = gradIsoU(idime,1)+dlxigp_ip(jgaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                     gradIsoU(idime,2) = gradIsoU(idime,2)+dlxigp_ip(jgaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                     gradIsoU(idime,3) = gradIsoU(idime,3)+dlxigp_ip(jgaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
                  end do
               end do
               gradU(:,:) = 0.0_rp
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     !$acc loop seq
                     do kdime = 1,ndime
                        gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,jgaus,ielem)*gradIsoU(idime,kdime)
                     end do
                  end do
               end do
               divU = gradU(1,1)+gradU(2,2)+gradU(3,3)
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
                  end do
                  tau(idime,idime) = tau(idime,idime) - (mu_fgp)*(2.0_rp/3.0_rp)*divU
               end do
					nmag = 0.0_rp
               !$acc loop seq
					do idime = 1,ndime
						nmag = nmag + bnorm((igaus-1)*ndime+idime)*bnorm((igaus-1)*ndime+idime)
					end do
               nmag = sqrt(nmag)
               !$acc loop seq
               do idime = 1,ndime
                  !$acc loop seq
                  do jdime = 1,ndime
                     tau_aux =  wgp_b(igaus)*tau(idime,jdime)*bnorm((igaus-1)*ndime+jdime)*sig
                     if(tau_aux .lt. 0) then
                        !$acc atomic update
                        Ftau_n_l(idime) = Ftau_n_l(idime) + tau_aux
                        !$acc end atomic
                     else
                        !$acc atomic update
                        Ftau_p_l(idime) = Ftau_p_l(idime) + tau_aux
                        !$acc end atomic
                     end if
                  end do
               end do
				end do
			end do
         !$acc end parallel loop
			deallocate(lelbo)

         call MPI_Allreduce(Ftau_n_l,Ftau_neg,ndime,mpi_datatype_real,MPI_SUM,app_comm,mpi_err)
         call MPI_Allreduce(Ftau_p_l,Ftau_pos,ndime,mpi_datatype_real,MPI_SUM,app_comm,mpi_err)

      end subroutine twInfo
end module mod_analysis
