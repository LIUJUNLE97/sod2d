module elem_convec

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

        subroutine full_convec_ijk_H(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener,initialze,fact)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
            real(rp),    intent(inout) :: Rmass(npoin)
            real(rp),    intent(inout) :: Rmom(npoin,ndime)
            real(rp),    intent(inout) :: Rener(npoin)
            logical, optional, intent(in)    :: initialze
            real(rp), optional, intent(in)  :: fact
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode,ndime)
            real(rp)                 :: Re_mass(nnode), Re_ener(nnode)
            real(rp)                 :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ,gradQ(ndime,ndime),gradIsoFuu(ndime,ndime,ndime)
            real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime),fuul(nnode,ndime,ndime),divFuu(ndime)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
            real(rp)  :: aux_fact = 1.0_rp

            call nvtxStartRange("Full convection")
            if(present(initialze)) then
               if (initialze .eqv. .true.) then
               !$acc kernels
                  Rmom(:,:) = 0.0_rp
                  Rmass(:) = 0.0_rp
                  Rener(:) = 0.0_rp
                  !$acc end kernels
               end if
            else
               !$acc kernels
               Rmom(:,:) = 0.0_rp
               Rmass(:) = 0.0_rp
               Rener(:) = 0.0_rp
               !$acc end kernels
            end if
            if(present(fact)) then
               aux_fact = fact
            end if

            !$acc parallel loop gang private(ipoin,Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,fl,fel,fuul) !!vector_length(vecLength)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
               end do
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ul(inode,idime) = u(ipoin(inode),idime)
                     ql(inode,idime) = q(ipoin(inode),idime)
                     fel(inode,idime) = (E(ipoin(inode))+pr(ipoin(inode)))*u(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector collapse(3)
               do idime = 1,ndime
                  do jdime = 1,ndime
                     do inode = 1,nnode
                        fl(inode,idime,jdime)  = q(ipoin(inode),idime)*u(ipoin(inode),jdime)
                        fuul(inode,idime,jdime)  = u(ipoin(inode),idime)*u(ipoin(inode),jdime)
                     end do
                  end do
               end do
               !$acc loop vector
               do inode = 1,nnode
                  rhol(inode) = rho(ipoin(inode))
                  El(inode) = E(ipoin(inode))
                  prl(inode) = pr(ipoin(inode))
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoP, gradIsoE,gradIsoU, gradIsoF, gradIsoFuu, gradIsoQ, gradIsoFe,gradRho,gradP,gradE,gradU,divF,divU,divQ,divFe,gradQ,divFuu)
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

                  gradIsoRho(:) = 0.0_rp
                  gradIsoP(:) = 0.0_rp
                  gradIsoE(:) = 0.0_rp
                  gradIsoU(:,:) = 0.0_rp
                  gradIsoF(:,:,:) = 0.0_rp
                  gradIsoFuu(:,:,:) = 0.0_rp
                  gradIsoQ(:,:) = 0._rp
                  gradIsoFe(:,:) = 0._rp
                  !$acc loop seq
                  do ii=1,porder+1
                     gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                     gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                     gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                     gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*prl(invAtoIJK(ii,isoJ,isoK))
                     gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*prl(invAtoIJK(isoI,ii,isoK))
                     gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*prl(invAtoIJK(isoI,isoJ,ii))

                     gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*(prl(invAtoIJK(ii,isoJ,isoK))  + El(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*(prl(invAtoIJK(isoI,ii,isoK)) + El(invAtoIJK(isoI,ii,isoK)))
                     gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*(prl(invAtoIJK(isoI,isoJ,ii))+ El(invAtoIJK(isoI,isoJ,ii)))
                     
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                        
                        !$acc loop seq
                        do jdime=1, ndime
                            gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*fl(invAtoIJK(ii,isoJ,isoK),idime,jdime)
                            gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*fl(invAtoIJK(isoI,ii,isoK),idime,jdime)
                            gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*fl(invAtoIJK(isoI,isoJ,ii),idime,jdime)
                            gradIsoFuu(idime,jdime,1) = gradIsoFuu(idime,jdime,1) + dlxi_ip(ii)*fuul(invAtoIJK(ii,isoJ,isoK),idime,jdime)
                            gradIsoFuu(idime,jdime,2) = gradIsoFuu(idime,jdime,2) + dleta_ip(ii)*fuul(invAtoIJK(isoI,ii,isoK),idime,jdime)
                            gradIsoFuu(idime,jdime,3) = gradIsoFuu(idime,jdime,3) + dlzeta_ip(ii)*fuul(invAtoIJK(isoI,isoJ,ii),idime,jdime)
                        end do
                     end do
                  end do

                  gradRho(:) = 0.0_rp
                  gradP(:) = 0.0_rp
                  gradE(:) = 0.0_rp
                  gradU(:,:) = 0.0_rp
                  gradQ(:,:) = 0.0_rp
                  divF(:) = 0.0_rp
                  divFuu(:) = 0.0_rp
                  divQ = 0.0_rp
                  divFe = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                         gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                         gradE(idime)   = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                         divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                            gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                            divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                            divFuu(idime) = divFuu(idime) + He(jdime,kdime,igaus,ielem)*gradIsoFuu(idime,jdime,kdime)
                         end do
                     end do
                  end do
                  divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                  divQ  = gradQ(1,1)  + gradQ(2,2)  + gradQ(3,3) 
                  Re_mass(igaus) = 0.5_rp*(divQ+rhol(igaus)*divU) 
                  Re_ener(igaus) = 0.5_rp*(divFe+(El(igaus)+prl(igaus))*divU)
                  !$acc loop seq
                  do idime=1, ndime
                     Re_mom(igaus,idime) = 0.25_rp*(divF(idime)+ul(igaus,idime)*divQ+ql(igaus,idime)*divU+rhol(igaus)*divFuu(idime)) + gradP(idime)
                     Re_mass(igaus) = Re_mass(igaus) + 0.5_rp*(ul(igaus,idime)*gradRho(idime))
                     Re_ener(igaus) = Re_ener(igaus) + 0.5_rp*(ul(igaus,idime)*gradE(idime))
                     !$acc loop seq
                     do jdime=1, ndime
                        Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.25_rp*(fuul(igaus,idime,jdime)*gradRho(jdime)  &
                                                                  + ql(igaus,jdime)*gradU(idime,jdime)+ul(igaus,jdime)*gradQ(idime,jdime))
                     end do
                     Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                  end do
                  Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
                  Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
               end do
               !
               ! Final assembly
               !
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     !$acc atomic update
                     Rmom(ipoin(inode),idime) = Rmom(ipoin(inode),idime)+aux_fact*Re_mom(inode,idime)
                     !$acc end atomic
                  end do
               end do
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rmass(ipoin(inode)) = Rmass(ipoin(inode))+aux_fact*Re_mass(inode)
                  !$acc end atomic
                  !$acc atomic update
                  Rener(ipoin(inode)) = Rener(ipoin(inode))+aux_fact*Re_ener(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

         end subroutine full_convec_ijk_H


         subroutine full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener,initialze,fact)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
            real(rp),    intent(inout) :: Rmass(npoin)
            real(rp),    intent(inout) :: Rmom(npoin,ndime)
            real(rp),    intent(inout) :: Rener(npoin)
            logical, optional, intent(in)    :: initialze
            real(rp), optional, intent(in)  :: fact
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode,ndime)
            real(rp)                 :: Re_mass(nnode), Re_ener(nnode),aux
            real(rp)                 :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime),gradIsoRE(ndime),gradIsoFue(ndime,ndime), gradIsok(ndime),gradIsoRk(ndime),gradIsoFuk(ndime,ndime),gradIsoFk(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ,gradQ(ndime,ndime),gradIsoFuu(ndime,ndime,ndime),gradRE(ndime),divFue,gradk(ndime),divFk,gradRk(ndime),divFuk
            real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime),fuul(nnode,ndime,ndime),divFuu(ndime),REl(nnode),fuel(nnode,ndime),kl(nnode),fkl(nnode,ndime),fukl(nnode,ndime),Rkl(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
            real(rp)  :: aux_fact = 1.0_rp

            call nvtxStartRange("Full convection")
            if(present(initialze)) then
               if (initialze .eqv. .true.) then
               !$acc kernels
                  Rmom(:,:) = 0.0_rp
                  Rmass(:) = 0.0_rp
                  Rener(:) = 0.0_rp
                  !$acc end kernels
               end if
            else
               !$acc kernels
               Rmom(:,:) = 0.0_rp
               Rmass(:) = 0.0_rp
               Rener(:) = 0.0_rp
               !$acc end kernels
            end if
            if(present(fact)) then
               aux_fact = fact
            end if

            !$acc parallel loop gang private(ipoin,Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,REl,kl,Rkl) !!vector_length(vecLength)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                  kl(inode) = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     kl(inode) = kl(inode) + u(ipoin(inode),idime)*u(ipoin(inode),idime)*0.5_rp
                  end do
                  rhol(inode) = rho(ipoin(inode))
                  El(inode) = E(ipoin(inode))
                  REl(inode) = rho(ipoin(inode))*E(ipoin(inode))
                  Rkl(inode) = rho(ipoin(inode))*kl(inode)
                  prl(inode) = pr(ipoin(inode))
                  !$acc loop seq
                  do idime = 1,ndime
                     ul(inode,idime) = u(ipoin(inode),idime)
                     ql(inode,idime) = q(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoP,gradIsoRE,gradIsoRk, gradIsoE, gradIsok,gradIsoU, gradIsoF, gradIsoFuu, gradIsoQ, gradIsoFe,gradIsoFue,gradIsoFk,gradIsoFuk,gradRho,gradP,gradE,gradRE,gradk,gradRk,gradU,divF,divU,divQ,divFe,divFue,divFk,divFuk,gradQ,divFuu)
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

                  gradIsoRho(:) = 0.0_rp
                  gradIsoP(:) = 0.0_rp
                  gradIsoE(:) = 0.0_rp
                  gradIsoRE(:) = 0.0_rp
                  gradIsok(:) = 0.0_rp
                  gradIsoRk(:) = 0.0_rp
                  gradIsoU(:,:) = 0.0_rp
                  gradIsoF(:,:,:) = 0.0_rp
                  gradIsoFuu(:,:,:) = 0.0_rp
                  gradIsoQ(:,:) = 0.0_rp
                  gradIsoFe(:,:) = 0.0_rp
                  gradIsoFue(:,:) = 0.0_rp
                  gradIsoFk(:,:) = 0.0_rp
                  gradIsoFuk(:,:) = 0.0_rp
                  !$acc loop seq
                  do ii=1,porder+1
                     gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                     gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                     gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                     gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*prl(invAtoIJK(ii,isoJ,isoK))
                     gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*prl(invAtoIJK(isoI,ii,isoK))
                     gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*prl(invAtoIJK(isoI,isoJ,ii))

                     gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*(El(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*(El(invAtoIJK(isoI,ii,isoK)))
                     gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*(El(invAtoIJK(isoI,isoJ,ii)))

                     gradIsoRE(1) = gradIsoRE(1) + dlxi_ip(ii)*(REl(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoRE(2) = gradIsoRE(2) + dleta_ip(ii)*(REl(invAtoIJK(isoI,ii,isoK)))
                     gradIsoRE(3) = gradIsoRE(3) + dlzeta_ip(ii)*(REl(invAtoIJK(isoI,isoJ,ii)))

                     gradIsok(1) = gradIsok(1) + dlxi_ip(ii)*(kl(invAtoIJK(ii,isoJ,isoK)))
                     gradIsok(2) = gradIsok(2) + dleta_ip(ii)*(kl(invAtoIJK(isoI,ii,isoK)))
                     gradIsok(3) = gradIsok(3) + dlzeta_ip(ii)*(kl(invAtoIJK(isoI,isoJ,ii)))

                     gradIsoRk(1) = gradIsoRk(1) + dlxi_ip(ii)*(Rkl(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoRk(2) = gradIsoRk(2) + dleta_ip(ii)*(Rkl(invAtoIJK(isoI,ii,isoK)))
                     gradIsoRk(3) = gradIsoRk(3) + dlzeta_ip(ii)*(Rkl(invAtoIJK(isoI,isoJ,ii)))


                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))*El(invAtoIJK(ii,isoJ,isoK))*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))*El(invAtoIJK(isoI,ii,isoK))*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))*El(invAtoIJK(isoI,isoJ,ii))*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFue(idime,1) = gradIsoFue(idime,1) + dlxi_ip(ii)*El(invAtoIJK(ii,isoJ,isoK))*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFue(idime,2) = gradIsoFue(idime,2) + dleta_ip(ii)*El(invAtoIJK(isoI,ii,isoK))*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFue(idime,3) = gradIsoFue(idime,3) + dlzeta_ip(ii)*El(invAtoIJK(isoI,isoJ,ii))*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFk(idime,1) = gradIsoFk(idime,1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))*kl(invAtoIJK(ii,isoJ,isoK))*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFk(idime,2) = gradIsoFk(idime,2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))*kl(invAtoIJK(isoI,ii,isoK))*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFk(idime,3) = gradIsoFk(idime,3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))*kl(invAtoIJK(isoI,isoJ,ii))*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFuk(idime,1) = gradIsoFuk(idime,1) + dlxi_ip(ii)*kl(invAtoIJK(ii,isoJ,isoK))*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFuk(idime,2) = gradIsoFuk(idime,2) + dleta_ip(ii)*kl(invAtoIJK(isoI,ii,isoK))*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFuk(idime,3) = gradIsoFuk(idime,3) + dlzeta_ip(ii)*kl(invAtoIJK(isoI,isoJ,ii))*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        !$acc loop seq
                        do jdime=1, ndime
                           gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)*ul(invAtoIJK(ii,isoJ,isoK),jdime)
                           gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)*ul(invAtoIJK(isoI,ii,isoK),jdime)
                           gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)*ul(invAtoIJK(isoI,isoJ,ii),jdime)
                           gradIsoFuu(idime,jdime,1) = gradIsoFuu(idime,jdime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)*ul(invAtoIJK(ii,isoJ,isoK),jdime)
                           gradIsoFuu(idime,jdime,2) = gradIsoFuu(idime,jdime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)*ul(invAtoIJK(isoI,ii,isoK),jdime)
                           gradIsoFuu(idime,jdime,3) = gradIsoFuu(idime,jdime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)*ul(invAtoIJK(isoI,isoJ,ii),jdime)
                        end do
                     end do
                  end do

                  gradRho(:) = 0.0_rp
                  gradP(:) = 0.0_rp
                  gradE(:) = 0.0_rp
                  gradRE(:) = 0.0_rp
                  gradk(:) = 0.0_rp
                  gradRk(:) = 0.0_rp
                  gradU(:,:) = 0.0_rp
                  gradQ(:,:) = 0.0_rp
                  divF(:) = 0.0_rp
                  divFuu(:) = 0.0_rp
                  divQ = 0.0_rp
                  divFe = 0.0_rp
                  divFue = 0.0_rp
                  divFk = 0.0_rp
                  divFuk = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                         gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                         gradE(idime)   = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                         gradRE(idime)  = gradRE(idime) + He(idime,jdime,igaus,ielem) * gradIsoRE(jdime)
                         divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                         divFue = divFue + He(idime,jdime,igaus,ielem) * gradIsoFue(idime,jdime)
                         gradk(idime)   = gradk(idime) + He(idime,jdime,igaus,ielem) * gradIsok(jdime)
                         gradRk(idime)  = gradRk(idime) + He(idime,jdime,igaus,ielem) * gradIsoRk(jdime)
                         divFk = divFk + He(idime,jdime,igaus,ielem) * gradIsoFk(idime,jdime)
                         divFuk = divFuk + He(idime,jdime,igaus,ielem) * gradIsoFuk(idime,jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                            gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                            divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                            divFuu(idime) = divFuu(idime) + He(jdime,kdime,igaus,ielem)*gradIsoFuu(idime,jdime,kdime)
                         end do
                     end do
                  end do
                  divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                  divQ  = gradQ(1,1)  + gradQ(2,2)  + gradQ(3,3) 
                  Re_mass(igaus) = 0.5_rp*(divQ+rhol(igaus)*divU) 
                  Re_ener(igaus) =0.25_rp*(divFe+El(igaus)*divQ+REl(igaus)*divU+rhol(igaus)*divFue) + &
                                  0.25_rp*(divFk+kl(igaus)*divQ+Rkl(igaus)*divU+rhol(igaus)*divFuk)
                  aux = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     Re_mom(igaus,idime) = 0.25_rp*(divF(idime)+ul(igaus,idime)*divQ+ql(igaus,idime)*divU+rhol(igaus)*divFuu(idime)) + gradP(idime)
                     Re_mass(igaus) = Re_mass(igaus) + 0.5_rp*(ul(igaus,idime)*gradRho(idime))
                     Re_ener(igaus) = Re_ener(igaus) + 0.25_rp*(ql(igaus,idime)*gradE(idime)+ul(igaus,idime)*gradRE(idime)+&
                                                                El(igaus)*ul(igaus,idime)*gradRho(idime)) +&
                                                     0.25_rp*(ql(igaus,idime)*gradk(idime)+ul(igaus,idime)*gradRk(idime)+&
                                                                kl(igaus)*ul(igaus,idime)*gradRho(idime)) 
                     !$acc loop seq
                     do jdime=1, ndime
                        aux = ul(igaus,idime)*ul(igaus,jdime)
                        Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.25_rp*(aux*gradRho(jdime)  &
                                                                  + ql(igaus,jdime)*gradU(idime,jdime)+ul(igaus,jdime)*gradQ(idime,jdime))
                     end do
                     Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                  end do
                  Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
                  Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
               end do
               !
               ! Final assembly
               !
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rmass(ipoin(inode)) = Rmass(ipoin(inode))+aux_fact*Re_mass(inode)
                  !$acc end atomic
                  !$acc atomic update
                  Rener(ipoin(inode)) = Rener(ipoin(inode))+aux_fact*Re_ener(inode)
                  !$acc end atomic
                  !$acc loop seq
                  do idime = 1,ndime
                     !$acc atomic update
                     Rmom(ipoin(inode),idime) = Rmom(ipoin(inode),idime)+aux_fact*Re_mom(inode,idime)
                     !$acc end atomic
                  end do
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

         end subroutine full_convec_ijk

         subroutine generic_scalar_convec_ijk(nelem,npoin,connec,Ngp, &
               dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,q,eta,u,Rconvec,initialze,fact)

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
            real(rp),    intent(inout) :: Rconvec(npoin)
            logical, optional, intent(in)    :: initialze
            real(rp), optional, intent(in)  :: fact
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ii
            integer(4)               :: ipoin(nnode)
            real(rp)                 :: gradIsoE(ndime),gradIsoU(ndime,ndime),gradIsoFe(ndime,ndime)
            real(rp)                 :: gradE(ndime),gradU(ndime,ndime),divU,divFe
            real(rp)                 :: ul(nnode,ndime), fel(nnode,ndime), etal(nnode), Re(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
            real(rp)  :: aux_fact = 1.0_rp

            call nvtxStartRange("Generic Convection")
            if(present(initialze)) then
               if (initialze .eqv. .true.) then
                  !$acc kernels
                  Rconvec(:) = 0.0_rp
                  !$acc end kernels
               end if
            else
               !$acc kernels
               Rconvec(:) = 0.0_rp
               !$acc end kernels
            end if
            if(present(fact)) then
               aux_fact = fact
            end if
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
                  gradIsoU(:,:) = 0.0_rp
                  gradIsoFe(:,:) = 0._rp
                  !$acc loop seq
                  do ii=1,porder+1
                     gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*etal(invAtoIJK(ii,isoJ,isoK))
                     gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*etal(invAtoIJK(isoI,ii,isoK))
                     gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*etal(invAtoIJK(isoI,isoJ,ii))
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                     end do
                  end do

                  gradE(:) = 0.0_rp
                  divU = 0.0_rp
                  divFe = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradE(idime) = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                         divU = divU + He(idime,jdime,igaus,ielem) * gradIsoU(idime,jdime)
                         divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                     end do
                  end do

                  Re(igaus) = 0.5_rp*(divFe+etal(igaus)*divU)
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
                  Rconvec(ipoin(inode)) = Rconvec(ipoin(inode))+aux_fact*Re(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

         end subroutine generic_scalar_convec_ijk

         subroutine generic_scalar_convec_projection_residual_ijk(nelem,npoin,connec,Ngp, &
            dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,q,eta,u,Rproj,Rconvec,initialze,fact)

         implicit none

         integer(4), intent(in)  :: nelem, npoin
         integer(4), intent(in)  :: connec(nelem,nnode)
         real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
         real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
         real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
         integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
         real(rp),    intent(in)  :: q(npoin,ndime)
         real(rp),    intent(in)  :: eta(npoin),Rproj(npoin)
         real(rp),    intent(in)  :: u(npoin,ndime)
         real(rp),    intent(inout) :: Rconvec(npoin)
         logical, optional, intent(in)    :: initialze
         real(rp), optional, intent(in)  :: fact
         integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ii
         integer(4)               :: ipoin(nnode)
         real(rp)                 :: gradIsoE(ndime),gradIsoU(ndime,ndime),gradIsoFe(ndime,ndime)
         real(rp)                 :: gradE(ndime),gradU(ndime,ndime),divU,divFe
         real(rp)                 :: ul(nnode,ndime), fel(nnode,ndime), etal(nnode), Re(nnode),Rprojl(nnode)
         real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip
         real(rp)  :: aux_fact = 1.0_rp

         call nvtxStartRange("Generic Convection")
         if(present(initialze)) then
            if (initialze .eqv. .true.) then
               !$acc kernels
               Rconvec(:) = 0.0_rp
               !$acc end kernels
            end if
         else
            !$acc kernels
            Rconvec(:) = 0.0_rp
            !$acc end kernels
         end if
         if(present(fact)) then
            aux_fact = fact
         end if
         !$acc parallel loop gang  private(ipoin,Re,ul,fel,etal,Rprojl) !!vector_length(vecLength)
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
               Rprojl(inode) = Rproj(ipoin(inode))
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
               gradIsoU(:,:) = 0.0_rp
               gradIsoFe(:,:) = 0._rp
               !$acc loop seq
               do ii=1,porder+1
                  gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*etal(invAtoIJK(ii,isoJ,isoK))
                  gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*etal(invAtoIJK(isoI,ii,isoK))
                  gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*etal(invAtoIJK(isoI,isoJ,ii))
                  !$acc loop seq
                  do idime=1,ndime
                     gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                     gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                     gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                     gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                     gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                     gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                  end do
               end do

               gradE(:) = 0.0_rp
               divU = 0.0_rp
               divFe = 0.0_rp
               !$acc loop seq
               do idime=1, ndime
                  !$acc loop seq
                  do jdime=1, ndime
                      gradE(idime) = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                      divU = divU + He(idime,jdime,igaus,ielem) * gradIsoU(idime,jdime)
                      divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                  end do
               end do

               Re(igaus) = 0.5_rp*(divFe+etal(igaus)*divU)
               !$acc loop seq
               do idime=1, ndime
                  Re(igaus) = Re(igaus) + 0.5_rp*(ul(igaus,idime)*gradE(idime))
               end do
               Re(igaus) = gpvol(1,igaus,ielem)*(Re(igaus)-Rprojl(igaus))
            end do
            !
            ! Assembly
            !
            !$acc loop vector
            do inode = 1,nnode
               !$acc atomic update
               Rconvec(ipoin(inode)) = Rconvec(ipoin(inode))+aux_fact*Re(inode)
               !$acc end atomic
            end do
         end do
         !$acc end parallel loop
         call nvtxEndRange

      end subroutine generic_scalar_convec_projection_residual_ijk

         subroutine full_convec_ijk_jacobian(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),   intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),   intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
            real(rp),   intent(inout) :: Rmass(npoin)
            real(rp),   intent(inout) :: Rmom(npoin,ndime)
            real(rp),   intent(inout) :: Rener(npoin)
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            real(rp)                :: Re_mom(nnode,ndime)
            real(rp)                :: Re_mass(nnode), Re_ener(nnode)
            real(rp)                :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime)
            real(rp)                :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ
            real(rp)                :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Full convection")
            !$acc kernels
            Rmom(:,:) = 0.0_rp
            Rmass(:) = 0.0_rp
            Rener(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,fl,fel) !!vector_length(vecLength)
            do ielem = 1,nelem
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ql(inode,idime) = rho(connec(ielem,inode))*u(connec(ielem,inode),idime)
                     fel(inode,idime) = (E(connec(ielem,inode))+pr(connec(ielem,inode)))*u(connec(ielem,inode),idime)
                  end do
               end do
               !$acc loop vector collapse(3)
               do idime = 1,ndime
                  do jdime = 1,ndime
                     do inode = 1,nnode
                        fl(inode,idime,jdime)  = q(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime)
                     end do
                  end do
               end do
               !$acc loop vector
               do inode = 1,nnode
                  prl(inode) = pr(connec(ielem,inode))
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoP, gradIsoE,gradIsoU, gradIsoF, gradIsoQ, gradIsoFe,gradRho,gradP,gradE,gradU,divF,divU,divQ,divFe)
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

                  gradIsoRho(:) = 0.0_rp
                  gradIsoP(:) = 0.0_rp
                  gradIsoE(:) = 0.0_rp
                  gradIsoU(:,:) = 0.0_rp
                  gradIsoF(:,:,:) = 0._rp
                  gradIsoQ(:,:) = 0._rp
                  gradIsoFe(:,:) = 0._rp
                  !$acc loop seq
                  do ii=1,porder+1

                     gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*prl(invAtoIJK(ii,isoJ,isoK))
                     gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*prl(invAtoIJK(isoI,ii,isoK))
                     gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*prl(invAtoIJK(isoI,isoJ,ii))
                     
                     !$acc loop seq
                     do idime=1,ndime

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
                        
                        !$acc loop seq
                        do jdime=1, ndime
                            gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*fl(invAtoIJK(ii,isoJ,isoK),idime,jdime)
                            gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*fl(invAtoIJK(isoI,ii,isoK),idime,jdime)
                            gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*fl(invAtoIJK(isoI,isoJ,ii),idime,jdime)
                        end do
                     end do
                  end do

                  gradRho(:) = 0.0_rp
                  gradP(:) = 0.0_rp
                  gradE(:) = 0.0_rp
                  gradU(:,:) = 0.0_rp
                  divF(:) = 0.0_rp
                  divQ = 0.0_rp
                  divFe = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                         divQ = divQ + He(idime,jdime,igaus,ielem) * gradIsoQ(idime,jdime)
                         divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                         end do
                     end do
                  end do
                  Re_mass(igaus) = divQ 
                  Re_ener(igaus) = divFe
                  !$acc loop seq
                  do idime=1, ndime
                     Re_mom(igaus,idime) = divF(idime)+gradP(idime)
                     Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                  end do
                  Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
                  Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
               end do
               !
               ! Final assembly
               !
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     !$acc atomic update
                     Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                     !$acc end atomic
                  end do
               end do
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re_mass(inode)
                  !$acc end atomic
                  !$acc atomic update
                  Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re_ener(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

         end subroutine full_convec_ijk_jacobian
end module elem_convec
