module mod_aver
	contains
	subroutine  favre_average()
		implicit none
		! Compute accumulated time
		avtim = avtim+dt
		! Compute accumulated tally for density times current dt and other variables times density times current dt
		do ipoin = 1,npoin_w
			avden(lpoin_w(ipoin)) = avden(lpoin_w(ipoin)) + rho(lpoin_w(ipoin))*dt ! Sum(rho*dt)
			avrhou(lpoin_w(ipoin),1) = avrhou(lpoin_w(ipoin),1) + rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),1)*dt ! Sum(rho*u1*dt)
			avrhou(lpoin_w(ipoin),2) = avrhou(lpoin_w(ipoin),2) + rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),2)*dt ! Sum(rho*u2*dt)
			avrhou(lpoin_w(ipoin),3) = avrhou(lpoin_w(ipoin),3) + rho(lpoin_w(ipoin))*u(lpoin_w(ipoin),3)*dt ! Sum(rho*u3*dt)
		end do
	end subroutine favre_average
end module mod_aver