module mod_postpro

   use mod_constants
   use mod_nvtx
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

contains

   subroutine compute_fieldDerivs(nelem,npoin,npoin_w,lpoin_w,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho,u,gradRho,curlU,divU,Qcrit)

      implicit none

      integer(4), intent(in)  :: nelem, npoin, npoin_w, lpoin_w(npoin), connec(nelem,nnode), lelpn(npoin)
      real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem), dNgp(ndime,nnode,ngaus)
      real(rp),   intent(in)  :: rho(npoin), u(npoin,ndime)
      real(rp),   intent(in)  :: leviCivi(ndime,ndime,ndime)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1)
      integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(out) :: gradRho(npoin,ndime), curlU(npoin,ndime), divU(npoin), Qcrit(npoin)
      integer(4)              :: iNodeL, ielem, idime, inode, igaus, ipoin, jdime,kdime,isoI, isoJ, isoK , ii
      real(rp)                :: gpcar(ndime,nnode), rho_e(nnode), u_e(nnode,ndime), ru_e(nnode,ndime), aux1, aux2
      real(rp)                :: gradu_e(ndime,ndime),gradru_e(ndime,ndime),divul,divql
      real(rp)                :: gradIsoRho(ndime),gradIsoU(ndime,ndime),gradIsoRU(ndime,ndime)

      !$acc kernels
      divU(:) = 0.0_rp
      gradRho(:,:) = 0.0_rp
      curlU(:,:) = 0.0_rp
      Qcrit(:) = 0.0_rp
      !$acc end kernels
      !$acc parallel loop gang private(rho_e,u_e,ru_e)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            rho_e(inode) = rho(connec(ielem,inode))
            u_e(inode,1) = u(connec(ielem,inode),1)
            u_e(inode,2) = u(connec(ielem,inode),2)
            u_e(inode,3) = u(connec(ielem,inode),3)
            ru_e(inode,1) = rho(connec(ielem,inode))*u(connec(ielem,inode),1)
            ru_e(inode,2) = rho(connec(ielem,inode))*u(connec(ielem,inode),2)
            ru_e(inode,3) = rho(connec(ielem,inode))*u(connec(ielem,inode),3)
         end do
         !$acc loop vector private(gpcar,gradu_e,gradIsoRho,gradIsoU,gradIsoRU,gradru_e)
         do igaus = 1,ngaus
            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoRho(:) = 0.0_rp
            gradIsoU(:,:) = 0.0_rp
            gradIsoRU(:,:) = 0.0_rp
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

                  gradIsoRU(idime,1) = gradIsoRU(idime,1) + dlxigp_ip(igaus,1,ii)*ru_e(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoRU(idime,2) = gradIsoRU(idime,2) + dlxigp_ip(igaus,2,ii)*ru_e(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoRU(idime,3) = gradIsoRU(idime,3) + dlxigp_ip(igaus,3,ii)*ru_e(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            !gradRho(connec(ielem,igaus),:) = 0.0_rp
            gradu_e(:,:) = 0.0_rp
            divql = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  !$acc atomic update
                  gradRho(connec(ielem,igaus),idime) = gradRho(connec(ielem,igaus),idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                  !$acc end atomic
                  divql = divql + He(idime,jdime,igaus,ielem) * gradIsoRU(idime,jdime)
                  !$acc loop seq
                  do kdime=1,ndime
                     gradu_e(idime,jdime) = gradu_e(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                  end do
               end do
            end do
            divul  = gradu_e(1,1)  + gradu_e(2,2)  + gradu_e(3,3) 

            !$acc atomic update
            divU(connec(ielem,igaus)) = divU(connec(ielem,igaus))+0.5_rp*(divql+rho(connec(ielem,igaus))*divul)
            !$acc end atomic

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

      if(mpi_size.ge.2) then
         call nvtxStartRange("MPI_comms_post")
         call mpi_halo_atomic_update_real(Qcrit)
         call mpi_halo_atomic_update_real(divU)
         do idime = 1,ndime
            call mpi_halo_atomic_update_real(curlU(:,idime))
            call mpi_halo_atomic_update_real(gradRho(:,idime))
         end do
         call nvtxEndRange
      end if

      !$acc parallel loop
      do ipoin = 1,npoin_w
         iNodeL=lpoin_w(ipoin)
         !$acc loop seq
         do idime = 1,ndime
            gradRho(iNodeL,idime) = gradRho(iNodeL,idime)/real(lelpn(iNodeL),rp)
            curlU(iNodeL,idime) = curlU(iNodeL,idime)/real(lelpn(iNodeL),rp)
         end do
         divU(iNodeL) = divU(iNodeL)/real(lelpn(iNodeL),rp)
         Qcrit(iNodeL) = Qcrit(iNodeL)/real(lelpn(iNodeL),rp)
      end do
      !$acc end parallel loop

   end subroutine compute_fieldDerivs

end module mod_postpro
