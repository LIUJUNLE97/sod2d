module mod_postpro

   use mod_constants

contains

   subroutine compute_fieldDerivs(nelem,npoin,connec,lelpn,He,dNgp,leviCivi,rho,u,gradRho,curlU,divU,Qcrit)

      implicit none

      integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode), lelpn(npoin)
      real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime)
      real(rp),    intent(in)  :: leviCivi(ndime,ndime,ndime)
      real(rp),    intent(out) :: gradRho(npoin,ndime), curlU(npoin,ndime), divU(npoin), Qcrit(npoin)
      integer(4)              :: ielem, idime, inode, igaus, ipoin, jdime, kdime
      real(rp)                 :: gpcar(ndime,nnode), rho_e(nnode), u_e(nnode,ndime), aux1, aux2
      real(rp)                 :: gradu_e(ndime,ndime)

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
         !$acc loop worker private(gpcar,gradu_e)
         do igaus = 1,ngaus
            !$acc loop vector collapse(2)
            do idime = 1,ndime
               do inode = 1,nnode
                  gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
               end do
            end do
            !$acc loop seq
            do idime = 1,ndime
               aux1 = 0.0_rp
               !$acc loop vector reduction(+:aux1)
               do inode = 1,nnode
                  aux1 = aux1+gpcar(idime,inode)*rho_e(inode)
               end do
               !$acc loop seq
               do jdime = 1,ndime
                  aux2 = 0.0_rp
                  !$acc loop vector reduction(+:aux2)
                  do inode = 1,nnode
                     aux2 = aux2+gpcar(jdime,inode)*u_e(inode,idime)
                  end do
                  gradu_e(idime,jdime) = aux2
               end do
               !
               ! Gradient of density
               !
               !$acc atomic update
               gradRho(connec(ielem,igaus),idime) = gradRho(connec(ielem,igaus),idime) + aux1
               !$acc end atomic
               !
               ! Divergence of velocity field
               !
               !$acc atomic update
               divU(connec(ielem,igaus)) = divU(connec(ielem,igaus)) + gradu_e(idime,idime)
               !$acc end atomic
            end do
            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                  !$acc loop seq
                  do kdime = 1,ndime
                     !$acc atomic update
                     curlU(connec(ielem,igaus),idime) = curlU(connec(ielem,igaus),idime) + leviCivi(idime,jdime,kdime)*gradu_e(jdime,kdime)
                     !$$acc end atomic update
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
