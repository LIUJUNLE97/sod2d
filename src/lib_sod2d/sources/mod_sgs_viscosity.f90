module mod_sgs_viscosity

   use mod_numerical_params
   use mod_nvtx
   use mod_veclen

   ! TODO: Finish module and create unit tests

contains

   ! it implents the Vreman SGS

   subroutine sgs_visc(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
         rho,u,Ml,mu_sgs)

      implicit none

      integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
      real(rp),   intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: rho(npoin), u(npoin,ndime),Ml(npoin)
      real(rp),   intent(out) :: mu_sgs(nelem,ngaus)
      integer(4)              :: ielem, inode, igaus, kdime, idime, jdime,isoI, isoJ, isoK 
      real(rp)                :: hLES, evol,gpcar(ndime,nnode), aux, mue(npoin), ave(npoin)
      real(rp)                :: gradU(ndime,ndime),gradV2(ndime,ndime),Bbeta, alpha, tau, aux2
      real(rp)                :: gradIsoU(ndime,ndime)
      real(rp)                :: ul(nnode,ndime)
      integer(4)              :: convertIJK(0:porder+2),ii,jj,kk,mm,nn,ll
      real(rp)                :: mue_l(nelem,nnode),al(-1:1),am(-1:1),an(-1:1),aux1

      do ii=3,porder+1
         convertIJK(ii-1) = ii
      end do 
      convertIJK(0) = 3
      convertIJK(1) = 1
      convertIJK(porder+1) = 2
      convertIJK(porder+2) = porder

      al(-1) = 1.0_rp/4.0_rp
      al(0)  = 2.0_rp/4.0_rp
      al(1)  = 1.0_rp/4.0_rp

      am(-1) = 1.0_rp/4.0_rp
      am(0)  = 2.0_rp/4.0_rp
      am(1)  = 1.0_rp/4.0_rp
      
      an(-1) = 1.0_rp/4.0_rp
      an(0)  = 2.0_rp/4.0_rp
      an(1)  = 1.0_rp/4.0_rp

      !$acc kernels
      mue_l(:,:) = mu_sgs(:,:)
      !$acc end kernels

      !$acc parallel loop gang private(ul) 
      do ielem = 1,nelem
         aux2 = 0.0_rp
         !$acc loop vector collapse(2)
         do inode = 1,nnode
            do idime = 1,ndime
               ul(inode,idime) = u(connec(ielem,inode),idime)
            end do
         end do
         !$acc loop vector private(gradU,gradV2,gradIsoU) 
         do igaus = 1,ngaus
            hLES = (gpvol(1,igaus,ielem)**(1.0_rp/3.0_rp))

            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoU(:,:) = 0.0_rp
            !$acc loop seq
            do ii=1,porder+1
               !$acc loop seq
               do idime=1,ndime
                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            gradU(:,:) = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  !$acc loop seq
                  do kdime=1,ndime
                     gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                  end do
               end do
            end do

            alpha = 0.0_rp
            !$acc loop seq
            do idime = 1, ndime
               !$acc loop seq
               do jdime = 1, ndime
                  gradV2(idime,jdime) = 0.0_rp
               end do
            end do

            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  !$acc loop seq
                  do kdime = 1,ndime
                     gradV2(idime,jdime) = gradV2(idime,jdime) + gradU(kdime,idime)*gradU(kdime,jdime)*hLES*hLES
                  end do
                  alpha = alpha + gradU(idime,jdime)*gradU(idime,jdime)
               end do
            end do
            Bbeta = gradV2(1,1)*gradV2(2,2) + gradV2(1,1)*gradV2(3,3) + gradV2(2,2)*gradV2(3,3) &
               - gradV2(1,2)*gradV2(1,2) - gradV2(1,3)*gradV2(1,3) - gradV2(2,3)*gradV2(2,3)
            if(alpha > 1.0e-10) then
               mue_l(ielem,igaus) = c_sgs*sqrt(max(Bbeta,1.0e-10)/alpha)
            end if
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
                           aux1 =   aux1 +  al(ll)*am(mm)*an(nn)*mue_l(ielem,invAtoIJK(convertIJK(ii+ll),convertIJK(jj+mm),convertIJK(kk+nn)))
                        end do
                     end do 
                  end do
                  mu_sgs(ielem,invAtoIJK(convertIJK(ii),convertIJK(jj),convertIJK(kk))) = aux1
               end do
            end do 
         end do

      end do
      !$acc end parallel loop


   end subroutine sgs_visc
end module mod_sgs_viscosity
