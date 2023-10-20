module mod_sgs_ilsa_viscosity

   use mod_numerical_params
   use mod_nvtx
   use mod_filters,only:convertIJK,al_weights,am_weights,an_weights

   ! TODO: Finish module and create unit tests

   implicit none

   real(rp), save :: time_ilsa = 0.0_rp

contains

   ! it implents the ILSA SGS

   subroutine sgs_ilsa_visc(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dt, &
         rho,u,mu_sgs,mu_fluid,mue,kres,etot,au,ax1,ax2,ax3,mue_l)

      implicit none

      integer(4),intent(in)    :: nelem, npoin, npoin_w,lpoin_w(npoin_w),connec(nelem,nnode)
      real(rp),  intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),  intent(in)    :: He(ndime,ndime,ngaus,nelem),dt
      real(rp),  intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1)
      integer(4),intent(in)    :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
      real(rp),  intent(in)    :: rho(npoin), u(npoin,ndime)
      real(rp),  intent(out)   :: mu_sgs(nelem,ngaus)
      real(rp),  intent(in)    :: mu_fluid(npoin),mue(nelem,ngaus)
      real(rp),  intent(inout) :: kres(npoin),etot(npoin),au(npoin,ndime),ax1(npoin),ax2(npoin),ax3(npoin),mue_l(nelem,nnode)
      integer(4)               :: ielem, inode, igaus, idime, jdime,ipoin,kdime,isoI, isoJ, isoK 
      real(rp)                 :: gpcar(ndime,nnode),aux,aux1,aux2,mueff(nnode),ul(nnode,ndime),kresl(nnode),etotl(nnode)
      real(rp)                 :: gradU(ndime,ndime), gradUf(ndime,ndime),eliti,ave,strain(ndime,ndime),strain_m,strainf(ndime,ndime),uf(nnode,ndime)
      real(rp)                 :: gpkres,ax1l(nnode),ax2l(nnode),ax3l(nnode),gpax1,gpax2,gpax3,c_k,a,b,c,d,gprij(ndime,ndime),gplest,gpepst,aul(nnode,ndime),aux3(nnode)
      real(rp)                 :: gradIsoUf(ndime,ndime),gradIsoU(ndime,ndime)
      integer(4)               :: ii,jj,kk,mm,nn,ll

      !$acc kernels
      mue_l(:,:) = mu_sgs(:,:)
      !$acc end kernels

      if(time_ilsa>T_ilsa) then
         time_ilsa = 0.0_rp
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

      !$acc parallel loop gang private(ul,mueff,kresl,etotl,ax1l,ax2l,ax3l,uf,aul,aux3) 
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            !mueff(inode) = mu_fluid(connec(ielem,inode))/rho(connec(ielem,inode)) + mue(ielem,inode)/rho(connec(ielem,inode))+mu_sgs(ielem,inode)
            mueff(inode) = mu_fluid(connec(ielem,inode))/rho(connec(ielem,inode)) + mu_sgs(ielem,inode)
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
         !$acc loop vector private(gpcar,gradU,strain,strainf,gradUf,gprij,aux,aux2,strain_m,gpkres,gpepst,gplest,gpax1,gpax2,gpax3,a,b,c,d,gradIsoUf, gradIsoU)
         do igaus = 1,ngaus
            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoU(:,:) = 0.0_rp
            gradIsoUf(:,:) = 0.0_rp
            !$acc loop seq
            do ii=1,porder+1
               !$acc loop seq
               do idime=1,ndime
                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                  gradIsoUf(idime,1) = gradIsoUf(idime,1) + dlxigp_ip(igaus,1,ii)*aul(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoUf(idime,2) = gradIsoUf(idime,2) + dlxigp_ip(igaus,2,ii)*aul(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoUf(idime,3) = gradIsoUf(idime,3) + dlxigp_ip(igaus,3,ii)*aul(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            gradU(:,:) = 0.0_rp
            gradUf(:,:) = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  !$acc loop seq
                  do kdime=1,ndime
                     gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                     gradUf(idime,jdime) = gradUf(idime,jdime) + He(jdime,kdime,igaus,ielem) * (gradIsoU(idime,kdime)-gradIsoUf(idime,kdime))
                  end do
               end do
            end do

            strain_m = 0.0_rp
            gpkres = eliti*kresl(igaus)+ave*0.5_rp*(uf(igaus,1)**2+uf(igaus,2)**2+uf(igaus,3)**2)
            aux = 0.0_rp
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  strain(idime,jdime)  = 0.5_rp*(gradU(idime,jdime)+gradU(jdime,idime))
                  strainf(idime,jdime) = 0.5_rp*(gradUf(idime,jdime)+gradUf(jdime,idime))
                  strain_m = strain_m + strain(idime,jdime)*strain(idime,jdime)
                  aux = aux + strainf(idime,jdime)*strainf(idime,jdime)
                  gprij(idime,jdime) =  uf(igaus,idime)*uf(igaus,jdime)
               end do
               gprij(idime,idime) = gprij(idime,idime) - (1.0_rp/3.0_rp)*uf(igaus,idime)*uf(igaus,idime)
            end do
            strain_m = sqrt(2.0_rp*strain_m)
            gpepst = eliti*etotl(igaus) + ave*2.0_rp*mueff(igaus)*aux

            ! Isla evaluations for ck
            gplest = sqrt(gpkres**3)/(1e-10+gpepst)
            gpax1 = eliti*ax1l(igaus)+ave*2.0_rp*(gplest**4)*(strain_m**4)
            aux = 0.0_rp
            aux2= 0.0_rp
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  aux  = aux + strain(idime,jdime)*gprij(idime,jdime)
                  aux2 = aux2 + gprij(idime,jdime)*gprij(idime,jdime)
               end do
            end do
            gpax2 = eliti*ax2l(igaus) + ave*4.0_rp*gplest*gplest*strain_m*aux
            gpax3 = eliti*ax3l(igaus) + ave*aux2
            a = gpax1*(1.0_rp - (1.0_rp/stau)**2)
            b = -gpax2
            c = gpax3
            if(b> 0.0_rp) then
               d = -0.5_rp*(b + sqrt(b*b - 4.0_rp*a*c))
            else 
               d = -0.5_rp*(b - sqrt(b*b - 4.0_rp*a*c))
            end if
            if(abs(d)> 1e-10) then
               if(abs(a)>1e-10) then 
                  aux = d / a
                  aux2 = c / d   
               else
                  aux = 0.0_rp
                  aux2 = c / d
               end if
            else
               aux = 0.0_rp
               aux2 = 0.0_rp 
            end if
            c_k = max(max(aux,aux2), 0.0_rp)

            mue_l(ielem,igaus)  = (c_k)*(gplest*gplest)*strain_m

            kres(connec(ielem,igaus)) = gpkres
            etot(connec(ielem,igaus)) = gpepst
            ax1(connec(ielem,igaus))  = gpax1
            ax2(connec(ielem,igaus))  = gpax2
            ax3(connec(ielem,igaus))  = gpax3
         end do

         !$acc loop vector collapse(3)
         do ii=1,porder+1
            do jj=1,porder+1
               do kk=1,porder+1           
                  aux1 = 0.00_rp
                  !$acc loop seq
                  do ll=-1,1
                     !$acc loop seq
                     do mm=-1,1
                        !$acc loop seq
                        do nn=-1,1           
                           aux1 =   aux1 +  al_weights(ll)*am_weights(mm)*an_weights(nn)*mue_l(ielem,invAtoIJK(convertIJK(ii+ll),convertIJK(jj+mm),convertIJK(kk+nn)))
                        end do
                     end do 
                  end do
                  mu_sgs(ielem,invAtoIJK(convertIJK(ii),convertIJK(jj),convertIJK(kk))) = aux1
               end do
            end do 
         end do

      end do
      !$acc end parallel loop

   end subroutine sgs_ilsa_visc
end module mod_sgs_ilsa_viscosity
