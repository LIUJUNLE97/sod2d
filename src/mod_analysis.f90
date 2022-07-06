module mod_analysis

   use mod_constants
   use mod_veclen
   use mod_nvtx

   contains

      subroutine volAvg_EK(nelem,npoin,connec,gpvol,Ngp,rho0,rho,u,EK)

         implicit none

         integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
         real(rp),    intent(in)  :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
         real(rp),    intent(in)  :: rho0, rho(npoin), u(npoin,ndime)
         real(rp),    intent(out) :: EK
         integer(4)              :: ielem, igaus, inode
         real(rp)                 :: R1

         call nvtxStartRange("Compute EK")
         EK = 0.0_rp
         !$acc parallel loop gang reduction(+:EK) vector_length(vecLength)
         do ielem = 1,nelem
            R1 = 0.0_rp
            !$acc loop vector collapse(2) reduction(+:R1)
            do igaus = 1,ngaus
               do inode = 1,nnode
                  R1 = R1 + gpvol(1,igaus,ielem)*Ngp(igaus,inode)*0.5_rp*rho(connec(ielem,inode))* &
                     dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))
               end do
            end do
            EK = EK+R1
         end do
         !$acc end parallel loop
         EK = EK/(rho0*((2.0_rp*3.14159_rp)**3.0_rp))
         call nvtxEndRange

      end subroutine volAvg_EK

      subroutine visc_dissipationRate(nelem,npoin,connec,leviCivi,rho0,mu_fluid,mu_e,u,volT,gpvol,He,dNgp,eps_S,eps_D,eps_T)

         implicit none

         integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
         real(rp),    intent(in)  :: leviCivi(3,3,3), rho0, mu_fluid(npoin), u(npoin,ndime), volT
         real(rp),    intent(in)  :: gpvol(1,ngaus,nelem), He(ndime,ndime,ngaus,nelem), dNgp(ndime,ngaus,nelem),mu_e(nelem,ngaus)
         real(rp),    intent(out) :: eps_S, eps_D, eps_T
         integer(4)              :: ielem, igaus, inode, idime, jdime, kdime
         real(rp)                 :: R1, R2, div2U, curl2U, alpha, aux,aux2
         real(rp)                 :: gpcar(ndime,nnode), gradU(ndime,ndime)

         if (flag_spectralElem .ne. 1) then
            write(1,*) "--| THIS ONLY WORKS WITH SPECTRAL ELEMENTS FOR NOW!"
            error stop
         end if

         eps_S = 0.0_rp
         eps_D = 0.0_rp
         eps_T = 0.0_rp

         call nvtxStartRange("Compute visc_dissipationRate")
         !$acc parallel loop gang private(gradU,gpcar) reduction(+:eps_S,eps_D) vector_length(vecLength)
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
            eps_S = eps_S+R1
            eps_D = eps_D+R2
         end do
         !!$acc end parallel loop
         call nvtxEndRange

         alpha = 1.0_rp/(rho0*volT)
         eps_S = eps_S*alpha
         eps_D = (4.0_rp/3.0_rp)*eps_D*alpha
         eps_T = eps_S+eps_D

      end subroutine visc_dissipationRate

      subroutine maxMach(npoin,npoin_w,lpoin_w,machno,maxmachno)

         implicit none

         integer(4), intent(in)  :: npoin
         integer(4), intent(in)  :: npoin_w
         integer(4), intent(in)  :: lpoin_w(npoin)
         real(rp),    intent(in)  :: machno(npoin)
         real(rp),    intent(out) :: maxmachno
         integer(4)              :: ipoin

         maxmachno = 0.0_rp
         !$acc parallel loop reduction(max:maxmachno)
         do ipoin = 1,npoin_w
            maxmachno = max(maxmachno,machno(lpoin_w(ipoin)))
         end do
         !$acc end parallel loop

      end subroutine maxMach

      subroutine write_EK(time,EK,eps_S,eps_D,eps_T,maxmachno)

         implicit none

         real(rp), intent(in) :: time, EK, eps_S, eps_D, eps_T, maxmachno

         write(666,10) time, EK, eps_S, eps_D, eps_T, maxmachno
10       format(6(F12.8,2X))

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
      subroutine surfInfo(nelem,npoin,nbound,surfCode,connec,bound,point2elem,bou_code, &
                          bounorm,wgp_b,dlxigp_ip,He,u,pr,surfArea,Fpr,Ftau)

         implicit none

			integer(4), intent(in)  :: npoin, nbound, surfCode, bound(nbound,npbou), bou_code(nbound,2)
			integer(4), intent(in)  :: nelem, connec(nelem,nnode), point2elem(npoin)
         real(rp),    intent(in)  :: wgp_b(npbou), bounorm(nbound,ndime*npbou)
         real(rp),    intent(in)  :: u(npoin,ndime), pr(npoin)
         real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1)
			real(rp),    intent(out) :: surfArea, Fpr(ndime), Ftau(ndime)
         integer(4)              :: ibound, idime, igaus, ipbou
         integer(4)              :: numBelem, counter
         integer(4), allocatable :: lelbo(:)
         real(rp)                 :: bnorm(npbou*ndime), nmag, prl(npbou)

			! Create lelbo for the surface, where lelbo is a list of boundary elements belonging to that surface
			numBelem = 0
         !$acc parallel loop reduction(+:numBelem)
			do ibound = 1, nbound
				if (bou_code(ibound,2) == surfCode) then
					numBelem = numBelem + 1
				end if
			end do
			counter = 0
			allocate(lelbo(numBelem))
			do ibound = 1, nbound
				if (bou_code(ibound,2) == surfCode) then
					counter = counter + 1
					lelbo(counter) = bou_code(ibound,1)
				end if
			end do

			! Compute surface information through sum of Gaussian quadratures over all boundary elements
			surfArea = 0.0_rp
         !$acc kernels
         Fpr(:) = 0.0_rp
         Ftau(:) = 0.0_rp
         !$acc end kernels
         !$acc parallel loop gang private(bnorm,prl) reduction(+:surfArea)
			do ibound = 1, numBelem
				bnorm(1:npbou*ndime) = bounorm(lelbo(ibound),1:npbou*ndime)
            prl(1:npbou) = pr(bound(lelbo(ibound),1:npbou))
				! Element area
            !$acc loop vector
				do igaus = 1,npbou
					nmag = 0.0_rp
               !$acc loop seq
					do idime = 1,ndime
						nmag = nmag + bnorm((igaus-1)*ndime+idime)*bnorm((igaus-1)*ndime+idime)
                  !$acc atomic update
                  Fpr(idime) = Fpr(idime)+wgp_b(igaus)*prl(igaus)*bnorm((igaus-1)*ndime+idime)
                  !$acc end atomic
					end do
               nmag = sqrt(nmag)
					surfArea = surfArea + nmag*wgp_b(igaus)
				end do
			end do
         !$acc end parallel loop
			deallocate(lelbo)
      end subroutine surfInfo

end module mod_analysis