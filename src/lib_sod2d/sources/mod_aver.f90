module mod_aver

	use mod_constants

	contains

	subroutine  favre_average(nelem,npoin,npoin_w,lpoin_w,connec,dt,rho,u,pr, &
									  mu_fluid,mu_e,mu_sgs,acutim,acurho,acupre,acuvel,acuve2,acuvex,acumueff)

		implicit none

		integer(4), intent(in)                             :: nelem, npoin, npoin_w, lpoin_w(npoin_w), connec(nelem,nnode)
		real(rp),    intent(in)                             :: dt
		real(rp),    intent(in),    dimension(npoin)        :: rho, pr, mu_fluid
		real(rp),    intent(in),    dimension(npoin,ndime)  :: u
		real(rp),    intent(in),    dimension(nelem,ngaus)  :: mu_e, mu_sgs
		real(rp),    intent(inout)                          :: acutim
		real(rp),    intent(inout), dimension(npoin)        :: acurho, acupre, acumueff
		real(rp),    intent(inout), dimension(npoin,ndime)  :: acuvel, acuve2, acuvex
		integer(4)                                         :: ipoin, idime, ielem, inode
		real(rp)                                            :: envit(npoin), mut(npoin)

		! Compute accumulated time
		acutim = acutim+dt
		! Compute accumulated tally for density times current dt and other variables times density times current dt
		!$acc parallel loop
		do ipoin = 1,npoin_w
			!$acc atomic update
			acurho(lpoin_w(ipoin)) = acurho(lpoin_w(ipoin))+rho(lpoin_w(ipoin))*dt
			!$acc end atomic
		end do
		!$acc end parallel loop
		!$acc parallel loop
		do ipoin = 1,npoin_w
			!$acc atomic update
			acupre(lpoin_w(ipoin)) = acupre(lpoin_w(ipoin))+pr(lpoin_w(ipoin))*dt
			!$acc end atomic
		end do
		!$acc end parallel loop
		!$acc parallel loop collapse(2)
		do ipoin = 1,npoin_w
			do idime = 1,ndime
				!$acc atomic update
				acuvel(lpoin_w(ipoin),idime) = acuvel(lpoin_w(ipoin),idime)+rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),idime)*dt
				!$acc end atomic
			end do
		end do
		!$acc end parallel loop
		!$acc parallel loop collapse(2)
		do ipoin = 1,npoin_w
			do idime = 1,ndime
				!$acc atomic update
				acuve2(lpoin_w(ipoin),idime) = acuve2(lpoin_w(ipoin),idime)+ &
					rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),idime)*u(lpoin_w(ipoin),idime)*dt
				!$acc end atomic
			end do
         acuvex(lpoin_w(ipoin),1) = acuve2(lpoin_w(ipoin),1)+ &
            rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),1)*u(lpoin_w(ipoin),2)*dt
         acuvex(lpoin_w(ipoin),2) = acuve2(lpoin_w(ipoin),2)+ &
            rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),1)*u(lpoin_w(ipoin),3)*dt
         acuvex(lpoin_w(ipoin),3) = acuve2(lpoin_w(ipoin),3)+ &
            rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),2)*u(lpoin_w(ipoin),3)*dt
		end do
		!$acc end parallel loop

		! Compute accumulated tally for effective viscosity times current dt
		!$acc parallel loop collapse(2)
		do ielem = 1,nelem
         do inode = 1, nnode
            envit(connec(ielem,inode)) =  mu_e(ielem,inode)
            mut(connec(ielem,inode))   =  mu_sgs(ielem,inode)
         end do
      end do
		!$acc end parallel loop
		!$acc parallel loop
		do ipoin = 1,npoin_w
			!$acc atomic update
			acumueff(lpoin_w(ipoin)) = acumueff(lpoin_w(ipoin))+ &
				(mu_fluid(lpoin_w(ipoin))+envit(lpoin_w(ipoin))+mut(lpoin_w(ipoin)))*dt
			!$acc end atomic
		end do
		!$acc end parallel loop

	end subroutine favre_average

   subroutine eval_average_window(isPeriodic,npoin,nelem,acuvel,acuve2,acuvex,acurho,acupre,acumueff,acutim,&
											avvel,avve2,avvex,avrho,avpre,avmueff,nper,masSla)
      implicit none
      logical, intent(in)                             :: isPeriodic
      integer(4), intent(in)                          :: npoin,nelem
      !integer(4), intent(in), optional						:: nper  
		integer(4), intent(in)          						:: nper   ! Oriol not good this is just for Itel compilation
      integer(4), intent(in), optional                :: masSla(nper,2)
      real(rp), intent(inout)                         :: acutim
      real(rp), intent(inout), dimension(npoin)       :: acurho, acupre, acumueff
      real(rp), intent(inout), dimension(npoin,ndime) :: acuvel, acuve2, acuvex
      real(rp), intent(inout), dimension(npoin)       :: avrho, avpre, avmueff
      real(rp), intent(inout), dimension(npoin,ndime) :: avvel, avve2, avvex
      integer(4) :: iper, idime, ipoin

      !
      ! If case is periodic, adjust slave nodes
      !
      if (isPeriodic .and. present(masSla)) then
      !if (isPeriodic .and. present(masSla) .and. present(nper)) then
         !$acc parallel loop
         do iper = 1,nper
            acuvel(masSla(iper,2),1) = acuvel(masSla(iper,1),1)
            acuvel(masSla(iper,2),2) = acuvel(masSla(iper,1),2)
            acuvel(masSla(iper,2),3) = acuvel(masSla(iper,1),3)
            acuve2(masSla(iper,2),1) = acuve2(masSla(iper,1),1)
            acuve2(masSla(iper,2),2) = acuve2(masSla(iper,1),2)
            acuve2(masSla(iper,2),3) = acuve2(masSla(iper,1),3)
            acuvex(masSla(iper,2),1) = acuvex(masSla(iper,1),1)
            acuvex(masSla(iper,2),2) = acuvex(masSla(iper,1),2)
            acuvex(masSla(iper,2),3) = acuvex(masSla(iper,1),3)
            acurho(masSla(iper,2)) = acurho(masSla(iper,1))
            acupre(masSla(iper,2)) = acupre(masSla(iper,1))
            acumueff(masSla(iper,2)) = acumueff(masSla(iper,1))
         end do
         !$acc end parallel loop
      end if

      !
      ! Divide accumulated vars by the accumulated time
      !
      !$acc kernels
      avrho(:) = acurho(:) / acutim
      avpre(:) = acupre(:) / acutim
      avmueff(:) = acumueff(:) / acutim
      avvel(:,:) = acuvel(:,:) / acutim
      avve2(:,:) = acuve2(:,:) / acutim
      avvex(:,:) = acuvex(:,:) / acutim
      !$acc end kernels

     ! !
     ! ! Favre average the rho*phi reynolds-averagedd variables
     ! !
     ! !$acc parallel loop
     ! do ipoin = 1,npoin
     !    !$acc loop seq
     !    do idime = 1,ndime
     !       avvel(ipoin,idime) = avvel(ipoin,idime)!/avrho(ipoin)
     !       avve2(ipoin,idime) = avve2(ipoin,idime)!/avrho(ipoin)
     !    end do
     ! end do
     ! !$acc end parallel loop

      !
      ! Reset the accumulated variables
      !
      !$acc kernels
      acurho(:) = 0.0_rp
      acupre(:) = 0.0_rp
      acumueff(:) = 0.0_rp
      acuvel(:,:) = 0.0_rp
      acuve2(:,:) = 0.0_rp
      acuvex(:,:) = 0.0_rp
      !$acc end kernels
      acutim = 0.0_rp

   end subroutine eval_average_window


end module mod_aver
