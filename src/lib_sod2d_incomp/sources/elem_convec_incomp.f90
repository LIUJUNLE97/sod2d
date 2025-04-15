module elem_convec_incomp

   use mod_nvtx
   use mod_numerical_params
   
   use mod_maths
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Computes convective term for Euler/NS equation system, as well             !
      ! as for any generic scalar transport that might occur. Based                !
      ! on Ljunkvist matrix-free implementation (assembles only rhs vector).       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains
         subroutine full_convec_ijk_incomp(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,Rmom)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)   :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin)
            real(rp),    intent(inout) :: Rmom(npoin,ndime)
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode,ndime)
            real(rp)                 :: gradIsoRho(ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradQ(ndime,ndime),divF(ndime),divU
            real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Full convection")
            !$acc kernels
            Rmom(:,:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(ipoin,Re_mom,ul,ql,rhol)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                  !$acc loop seq
                  do idime = 1,ndime
                     ul(inode,idime)      = u(ipoin(inode),idime)
                     ql(inode,idime)      = q(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho, gradIsoF, gradIsoQ,gradRho,gradQ,divF,divU)
               do igaus = 1,ngaus
                  !$acc loop seq
                  do ii=1,porder+1
                     dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                     dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                     dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                  end do
                  isoI = gmshAtoI(igaus) 
                  isoJ = gmshAtoJ(igaus) 
                  isoK = gmshAtoK(igaus) 

                  gradIsoF(:,:,:) = 0.0_rp
                  gradIsoQ(:,:) = 0.0_rp
                  !$acc loop seq
                  do ii=1,porder+1
                     
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)
                        !$acc loop seq
                        do jdime=1, ndime
                            gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)*ul(invAtoIJK(ii,isoJ,isoK),jdime)
                            gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)*ul(invAtoIJK(isoI,ii,isoK),jdime)
                            gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)*ul(invAtoIJK(isoI,isoJ,ii),jdime)
                        end do
                     end do
                  end do

                  gradQ(:,:) = 0.0_rp
                  divF(:) = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         !$acc loop seq
                         do kdime=1,ndime
                            gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                            divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                         end do
                     end do
                  end do
                  divU = gradQ(1,1) + gradQ(2,2) + gradQ(3,3)
                  !$acc loop seq
                  do idime=1, ndime
                     Re_mom(igaus,idime) = 0.5_rp*(divF(idime)+ul(igaus,idime)*divU)                     
                     !$acc loop seq
                     do jdime=1, ndime
                        Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.5_rp*(ul(igaus,jdime)*gradQ(idime,jdime))
                     end do
                     Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                  end do
               end do
               !
               ! Final assembly
               !
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     !$acc atomic update
                     Rmom(ipoin(inode),idime) = Rmom(ipoin(inode),idime)+Re_mom(inode,idime)
                     !$acc end atomic
                  end do
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

         end subroutine full_convec_ijk_incomp     
end module elem_convec_incomp
