module mod_postpro

   use mod_constants

contains

   subroutine compute_fieldDerivs(nelem,npoin,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho,u,gradRho,curlU,divU,Qcrit)

      implicit none

      integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode), lelpn(npoin)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime)
      real(rp),    intent(in)  :: leviCivi(ndime,ndime,ndime)
      real(rp),    intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1)
      integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),    intent(out) :: gradRho(npoin,ndime), curlU(npoin,ndime), divU(npoin), Qcrit(npoin)
      integer(4)              :: ielem, idime, inode, igaus, ipoin, jdime,kdime,isoI, isoJ, isoK , ii
      real(rp)                 :: gpcar(ndime,nnode), rho_e(nnode), u_e(nnode,ndime), aux1, aux2
      real(rp)                 :: gradu_e(ndime,ndime)
      real(rp)                 :: gradIsoRho(ndime),gradIsoU(ndime,ndime)

      !$acc kernels
      divU(:) = 0.0_rp
      gradRho(:,:) = 0.0_rp
      curlU(:,:) = 0.0_rp
      Qcrit(:) = 0.0_rp
      !$acc end kernels
      !$acc parallel loop gang private(rho_e,u_e)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            rho_e(inode) = rho(connec(ielem,inode))
            u_e(inode,1) = u(connec(ielem,inode),1)
            u_e(inode,2) = u(connec(ielem,inode),2)
            u_e(inode,3) = u(connec(ielem,inode),3)
         end do
         !$acc loop vector private(gpcar,gradu_e,gradIsoRho,gradIsoU)
         do igaus = 1,ngaus
            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoRho(:) = 0.0_rp
            gradIsoU(:,:) = 0.0_rp
            !$acc loop seq
            do ii=1,porder+1
               gradIsoRho(1) = gradIsoRho(1) + dlxigp_ip(igaus,1,ii)*rho_e(invAtoIJK(ii,isoJ,isoK))
               gradIsoRho(2) = gradIsoRho(2) + dlxigp_ip(igaus,2,ii)*rho_e(invAtoIJK(isoI,ii,isoK))
               gradIsoRho(3) = gradIsoRho(3) + dlxigp_ip(igaus,3,ii)*rho_e(invAtoIJK(isoI,isoJ,ii))

               !$acc loop seq
               do idime=1,ndime
                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*u_e(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*u_e(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*u_e(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            gradRho(connec(ielem,igaus),:) = 0.0_rp
            gradu_e(:,:) = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  gradRho(connec(ielem,igaus),idime) = gradRho(connec(ielem,igaus),idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                  !$acc loop seq
                  do kdime=1,ndime
                     gradu_e(idime,jdime) = gradu_e(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                  end do
               end do
            end do
            divU(connec(ielem,igaus)) = gradu_e(1,1)+gradu_e(2,2)+gradu_e(3,3)

            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  !$acc loop seq
                  do kdime = 1,ndime
                     !$acc atomic update
                     curlU(connec(ielem,igaus),idime) = curlU(connec(ielem,igaus),idime) + leviCivi(idime,jdime,kdime)*gradu_e(jdime,kdime)
                     !$acc end atomic
                  end do
                  !$acc atomic update
                  Qcrit(connec(ielem,igaus)) = Qcrit(connec(ielem,igaus))-0.5_rp*(gradu_e(idime,jdime)*gradu_e(jdime,idime))
                  !$acc end atomic
               end do
            end do
         end do
      end do
      !$acc end parallel loop

      !$acc parallel loop
      do ipoin = 1,npoin
         !$acc loop seq
         do idime = 1,ndime
            gradRho(ipoin,idime) = gradRho(ipoin,idime)/real(lelpn(ipoin),rp)
            curlU(ipoin,idime) = curlU(ipoin,idime)/real(lelpn(ipoin),rp)
         end do
         divU(ipoin) = divU(ipoin)/real(lelpn(ipoin),rp)
         Qcrit(ipoin) = Qcrit(ipoin)/real(lelpn(ipoin),rp)
      end do
      !$acc end parallel loop

   end subroutine compute_fieldDerivs

end module mod_postpro
