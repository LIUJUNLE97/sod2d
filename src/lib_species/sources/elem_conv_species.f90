module elem_convec_species

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
 
       subroutine species_convec_ijk(nelem,npoin,connec,Ngp, &
        He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho,Yk,u,Rconvec,initialze,fact)

        implicit none

        integer(4), intent(in)  :: nelem, npoin
        integer(4), intent(in)  :: connec(nelem,nnode)
        real(rp),    intent(in)  :: Ngp(ngaus,nnode)
        real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
        real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
        integer(4), intent(in)  ::  invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
        real(rp),    intent(in)  :: Yk(npoin)
        real(rp),    intent(in)  :: rho(npoin)
        real(rp),    intent(in)  :: u(npoin,ndime)
        real(rp),    intent(out) :: Rconvec(npoin)
        logical, optional, intent(in)    :: initialze
        real(rp), optional, intent(in)  :: fact        
        integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ii
        real(rp)                 :: gradIsoY(ndime),gradIsoU(ndime,ndime),gradIsoFy(ndime,ndime),gradIsoRho(ndime)
        real(rp)                 :: gradY(ndime),gradU(ndime,ndime),divU,divFy,gradRho(ndime)
        real(rp)                 :: ul(nnode,ndime), fykl(nnode,ndime), ykl(nnode),rhol(nnode), Re(nnode)
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

        !$acc parallel loop gang  private(Re,ul,fykl,rhol,ykl) 
        do ielem = 1,nelem
            !$acc loop vector collapse(2)
            do idime = 1,ndime
                do inode = 1,nnode
                    ul(inode,idime) = u(connec(ielem,inode),idime)
                    fykl(inode,idime) = rho(connec(ielem,inode))*u(connec(ielem,inode),idime)*Yk(connec(ielem,inode))
                end do
            end do
            !$acc loop vector
            do inode = 1,nnode
                rhol(inode) = rho(connec(ielem,inode))
                ykl(inode) = Yk(connec(ielem,inode))
            end do
            !$acc loop vector private(gradIsoY,gradIsoU,gradIsoFy,gradY,gradU,divU,divFy,gradIsoRho,gradRho, dlxi_ip, dleta_ip, dlzeta_ip)
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

            gradIsoY(:) = 0.0_rp
            gradIsoRho(:) = 0.0_rp
            gradIsoU(:,:) = 0.0_rp
            gradIsoFy(:,:) = 0._rp
            !$acc loop seq
            do ii=1,porder+1
                gradIsoY(1) = gradIsoY(1) + dlxi_ip(ii)*ykl(invAtoIJK(ii,isoJ,isoK))
                gradIsoY(2) = gradIsoY(2) + dleta_ip(ii)*ykl(invAtoIJK(isoI,ii,isoK))
                gradIsoY(3) = gradIsoY(3) + dlzeta_ip(ii)*ykl(invAtoIJK(isoI,isoJ,ii))
                gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))                          
                !$acc loop seq
                do idime=1,ndime
                    gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                    gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                    gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                    gradIsoFy(idime,1) = gradIsoFy(idime,1) + dlxi_ip(ii)*fykl(invAtoIJK(ii,isoJ,isoK),idime)
                    gradIsoFy(idime,2) = gradIsoFy(idime,2) + dleta_ip(ii)*fykl(invAtoIJK(isoI,ii,isoK),idime)
                    gradIsoFy(idime,3) = gradIsoFy(idime,3) + dlzeta_ip(ii)*fykl(invAtoIJK(isoI,isoJ,ii),idime)
                end do
            end do

            gradY(:) = 0.0_rp
            gradRho(:) = 0.0_rp
            divU = 0.0_rp
            divFy = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
                !$acc loop seq
                do jdime=1, ndime
                    gradY(idime) = gradY(idime) + He(idime,jdime,igaus,ielem) * gradIsoY(jdime)
                    gradRho(idime) = gradY(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                    divU = divU + He(idime,jdime,igaus,ielem) * gradIsoU(idime,jdime)
                    divFy = divFy + He(idime,jdime,igaus,ielem) * gradIsoFy(idime,jdime)
                end do
            end do

            Re(igaus) = 0.5_rp*(divFy+ykl(igaus)*rhol(igaus)*divU)
            !Re(igaus) = 0.5_rp*(divFy)
            !$acc loop seq
            do idime=1, ndime
                Re(igaus) = Re(igaus) + 0.5_rp*(ykl(igaus)*ul(igaus,idime)*gradRho(idime)+rhol(igaus)*ul(igaus,idime)*gradY(idime))
                !Re(igaus) = Re(igaus) + 0.5_rp*(rhol(igaus)*ul(igaus,idime)*gradY(idime))
            end do
            Re(igaus) = gpvol(1,igaus,ielem)*Re(igaus)
            end do
            !
            ! Assembly
            !
            !$acc loop vector
            do inode = 1,nnode
            !$acc atomic update
            Rconvec(connec(ielem,inode)) = Rconvec(connec(ielem,inode))+aux_fact*Re(inode)
            !$acc end atomic
            end do
        end do
        !$acc end parallel loop
        call nvtxEndRange

    end subroutine species_convec_ijk

end module elem_convec_species
