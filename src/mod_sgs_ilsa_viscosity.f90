module mod_sgs_ilsa_viscosity

   use mod_constants
   use mod_nvtx
   use mod_veclen

   ! TODO: Finish module and create unit tests

   implicit none

   real(8), save :: time_ilsa = 0.0d0

contains

   ! it implents the ILSA SGS

   subroutine sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dt, &
         rho,u,mu_sgs,mu_fluid,mue,kres,etot,au,ax1,ax2,ax3)

      implicit none

      integer(4), intent(in)  :: nelem, npoin, npoin_w,lpoin_w(npoin_w),connec(nelem,nnode)
      real(8),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(8),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dt
      real(8),    intent(in)  :: rho(npoin), u(npoin,ndime)
      real(8),    intent(out) :: mu_sgs(nelem,ngaus)
      real(8),    intent(in) :: mu_fluid(npoin),mue(nelem,ngaus)
      real(8),    intent(inout) :: kres(npoin),etot(npoin),au(npoin,ndime),ax1(npoin),ax2(npoin),ax3(npoin)
      integer(4)              :: ielem, inode, igaus, kdime, idime, jdime,ipoin
      real(8)                 :: gpcar(ndime,nnode), aux,aux2,mueff(nnode),ul(nnode,ndime),kresl(nnode),etotl(nnode)
      real(8)                 :: gradU(ndime,ndime), gradUf(ndime,ndime),eliti,ave,strain(ndime,ndime),strain_m,strainf(ndime,ndime),uf(nnode,ndime)
      real(8)                 :: gpkres,ax1l(nnode),ax2l(nnode),ax3l(nnode),gpax1,gpax2,gpax3,c_k,a,b,c,d,gprij(ndime,ndime),gplest,gpepst,aul(nnode,ndime)

      if(time_ilsa>T_ilsa) then
         time_ilsa = 0.0
      end if
      ave = dt/(time_ilsa+dt)
      eliti = time_ilsa/(time_ilsa+dt)
      time_ilsa = time_ilsa+dt

      !$acc parallel loop collapse(2)
      do ipoin = 1,npoin_w
         do idime = 1,ndime
            au(lpoin_w(ipoin),idime) =  ave*u(lpoin_w(ipoin),idime) + eliti*au(lpoin_w(ipoin),idime) 
         end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang private(ul,mueff,kresl,etotl,ax1l,ax2l,ax3l,uf,aul) 
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            mueff(inode) = mu_fluid(connec(ielem,inode))/rho(connec(ielem,inode)) + mue(ielem,inode)/rho(connec(ielem,inode))+mu_sgs(ielem,inode)
            kresl(inode) = kres(connec(ielem,inode))
            etotl(inode) = etot(connec(ielem,inode)) 
            ax1l(inode)  = ax1(connec(ielem,inode)) 
            ax2l(inode)  = ax2(connec(ielem,inode)) 
            ax3l(inode)  = ax3(connec(ielem,inode)) 
         end do
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               ul(inode,idime) = u(connec(ielem,inode),idime)
               uf(inode,idime) = u(connec(ielem,inode),idime) - au(connec(ielem,inode),idime)
               aul(inode,idime) = au(connec(ielem,inode),idime)
            end do
         end do
         !$acc loop worker private(gpcar,gradU,strain,strainf,gradUf,gprij,aux,aux2,strain_m,gpkres,gpepst,gplest,gpax1,gpax2,gpax3,a,b,c,d)
         do igaus = 1,ngaus
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop vector
               do inode = 1,nnode
                  aux =  dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                  gpcar(idime,inode) = aux
               end do
            end do

            ! Compute strain
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  aux  = 0.0d0
                  aux2 = 0.0d0
                  !$acc loop vector reduction(+:aux,aux2)
                  do inode = 1,nnode
                     aux  = aux +gpcar(jdime,inode)*ul(inode,idime)
                     aux2 = aux2+gpcar(jdime,inode)*aul(inode,idime)
                  end do
                  gradU(idime,jdime) = aux
                  gradUf(idime,jdime) = aux-aux2
               end do
            end do
            strain_m = 0.0d0
            gpkres = eliti*kresl(igaus)+ave*0.5d0*(uf(igaus,1)**2+uf(igaus,2)**2+uf(igaus,3)**2)
            aux = 0.0d0
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  strain(idime,jdime)  = 0.5d0*(gradU(idime,jdime)+gradU(jdime,idime))
                  strainf(idime,jdime) = 0.5d0*(gradUf(idime,jdime)+gradUf(jdime,idime))
                  strain_m = strain_m + strain(idime,jdime)*strain(idime,jdime)
                  aux = aux + strainf(idime,jdime)*strainf(idime,jdime)
                  gprij(idime,jdime) =  uf(igaus,idime)*uf(igaus,jdime)
               end do
               gprij(idime,idime) = gprij(idime,idime) - (1.0d0/3.0d0)*uf(igaus,idime)*uf(igaus,idime)
            end do
            strain_m = sqrt(2.0d0*strain_m)
            gpepst = eliti*etotl(igaus) + ave*2.0d0*mueff(igaus)*aux

            ! Isla evaluations for ck
            gplest = sqrt(gpkres**3)/(1e-10+gpepst)
            gpax1 = eliti*ax1l(igaus)+ave*2.0d0*(gplest**4)*(strain_m**4)
            aux = 0.0d0
            aux2= 0.0d0
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  aux  = aux + strain(idime,jdime)*gprij(idime,jdime)
                  aux2 = aux2 + gprij(idime,jdime)*gprij(idime,jdime)
               end do
            end do
            gpax2 = eliti*ax2l(igaus) + ave*4.0d0*gplest*gplest*strain_m*aux
            gpax3 = eliti*ax3l(igaus) + ave*aux2
            a = gpax1*(1.0d0 - (1.0d0/stau)**2)
            b = -gpax2
            c = gpax3
            if(b> 0.0d0) then
               d = -0.5d0*(b + sqrt(b*b - 4.0d0*a*c))
            else 
               d = -0.5d0*(b - sqrt(b*b - 4.0d0*a*c))
            end if
            if(abs(d)> 1e-10) then
               if(abs(a)>1e-10) then 
                  aux = d / a
                  aux2 = c / d   
               else
                  aux = 0.0d0
                  aux2 = c / d
               end if
            else
               aux = 0.0d0
               aux2 = 0.0d0 
            end if
            c_k = max(max(aux,aux2), 0.0d0)


            ! final viscosity
            if(etotl(igaus)>1.0e-10) then
               mu_sgs(ielem,igaus) =  (c_k)*(gplest*gplest)*strain_m
            else
               mu_sgs(ielem,igaus) = 0.0d0
            end if

            kres(connec(ielem,igaus)) = gpkres
            etot(connec(ielem,igaus)) = gpepst
            ax1(connec(ielem,igaus))  = gpax1
            ax2(connec(ielem,igaus))  = gpax2
            ax3(connec(ielem,igaus))  = gpax3
         end do
      end do
      !$acc end parallel loop

   end subroutine sgs_ilsa_visc
end module mod_sgs_ilsa_viscosity
