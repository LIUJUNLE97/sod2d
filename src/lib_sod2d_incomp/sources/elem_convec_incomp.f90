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
         subroutine full_convec_ijk_incomp(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,Rmom)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)   :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin)
            real(rp),    intent(out) :: Rmom(npoin,ndime)
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode,ndime)
            real(rp)                 :: gradIsoRho(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradQ(ndime,ndime),divF(ndime),divU,divQ_star, gradIsoQ_star(ndime,ndime)
            real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), fl(nnode,ndime,ndime),ql_star(nnode,ndime)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Full convection")
            !$acc kernels
            Rmom(:,:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(ipoin,Re_mom,ul,ql,rhol,fl,ql_star)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
               end do
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ul(inode,idime)      = u(ipoin(inode),idime)
                     ql(inode,idime)      = q(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector collapse(3)
               do idime = 1,ndime
                  do jdime = 1,ndime
                     do inode = 1,nnode
                        fl(inode,idime,jdime)  = q(ipoin(inode),idime)*u(ipoin(inode),jdime)
                     end do
                  end do
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoU, gradIsoF, gradIsoQ,gradIsoQ_star,gradRho,gradQ,divF,divU)
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

                  gradIsoU(:,:) = 0.0_rp
                  gradIsoF(:,:,:) = 0.0_rp
                  gradIsoQ(:,:) = 0.0_rp
                  !$acc loop seq
                  do ii=1,porder+1
                     
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)
                        
                        !$acc loop seq
                        do jdime=1, ndime
                            gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*fl(invAtoIJK(ii,isoJ,isoK),idime,jdime)
                            gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*fl(invAtoIJK(isoI,ii,isoK),idime,jdime)
                            gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*fl(invAtoIJK(isoI,isoJ,ii),idime,jdime)
                        end do
                     end do
                  end do

                  gradQ(:,:) = 0.0_rp
                  divF(:) = 0.0_rp
                  divU = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         divU = divU + He(idime,jdime,igaus,ielem) * gradIsoU(idime,jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                            divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                         end do
                     end do
                  end do
                  !$acc loop seq
                  do idime=1, ndime
                     Re_mom(igaus,idime) = 0.5_rp*(divF(idime))
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

         subroutine generic_scalar_convec_ijk_incomp(nelem,npoin,connec,Ngp, &
               dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,q,eta,u,Rconvec,alpha)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime)
            real(rp),    intent(in)  :: eta(npoin)
            real(rp),    intent(in)  :: u(npoin,ndime)
            real(rp),    intent(in)  :: alpha(npoin)
            real(rp),    intent(out) :: Rconvec(npoin)
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ii
            integer(4)               :: ipoin(nnode)
            real(rp)                 :: gradIsoE(ndime),gradIsoU(ndime,ndime),gradIsoFe(ndime,ndime)
            real(rp)                 :: gradE(ndime),gradU(ndime,ndime),divU,divFe
            real(rp)                 :: ul(nnode,ndime), fel(nnode,ndime), etal(nnode), Re(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Generic Convection")
            !$acc kernels
            Rconvec(:) = 0.0_rp
            !$acc end kernels
            !$acc parallel loop gang  private(ipoin,Re,ul,fel,etal) !!vector_length(vecLength)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
               end do
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ul(inode,idime) = u(ipoin(inode),idime)
                     fel(inode,idime) = q(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector
               do inode = 1,nnode
                  etal(inode) = eta(ipoin(inode))
               end do
               !$acc loop vector private(gradIsoE,gradIsoU,gradIsoFe,gradE,gradU,divU,divFe, dlxi_ip, dleta_ip, dlzeta_ip)
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

                  gradIsoE(:) = 0.0_rp
                  gradIsoFe(:,:) = 0._rp
                  !$acc loop seq
                  do ii=1,porder+1
                     gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*etal(invAtoIJK(ii,isoJ,isoK))
                     gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*etal(invAtoIJK(isoI,ii,isoK))
                     gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*etal(invAtoIJK(isoI,isoJ,ii))
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                     end do
                  end do

                  gradE(:) = 0.0_rp
                  divFe = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradE(idime) = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                         divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                     end do
                  end do

                  Re(igaus) = 0.5_rp*(divFe)
                  !$acc loop seq
                  do idime=1, ndime
                     Re(igaus) = Re(igaus) + 0.5_rp*(ul(igaus,idime)*gradE(idime))
                  end do
                  Re(igaus) = gpvol(1,igaus,ielem)*Re(igaus)
               end do
               !
               ! Assembly
               !
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rconvec(ipoin(inode)) = Rconvec(ipoin(inode))+Re(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

         end subroutine generic_scalar_convec_ijk_incomp
end module elem_convec_incomp
