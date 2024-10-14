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
        real(rp)                 :: gradIsoFy(ndime,ndime)
        real(rp)                 :: divFy
        real(rp)                 :: fykl(nnode,ndime), Re(nnode)
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

        !$acc parallel loop gang  private(Re,fykl) 
        do ielem = 1,nelem
            !$acc loop vector collapse(2)
            do idime = 1,ndime
                do inode = 1,nnode
                    fykl(inode,idime) = rho(connec(ielem,inode))*u(connec(ielem,inode),idime)*Yk(connec(ielem,inode))
                end do
            end do
            !$acc loop vector private(gradIsoFy,divFy, dlxi_ip, dleta_ip, dlzeta_ip)
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

                gradIsoFy(:,:) = 0._rp
                !$acc loop seq
                do ii=1,porder+1
                    !$acc loop seq
                    do idime=1,ndime
                        gradIsoFy(idime,1) = gradIsoFy(idime,1) + dlxi_ip(ii)*fykl(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFy(idime,2) = gradIsoFy(idime,2) + dleta_ip(ii)*fykl(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFy(idime,3) = gradIsoFy(idime,3) + dlzeta_ip(ii)*fykl(invAtoIJK(isoI,isoJ,ii),idime)
                    end do
                end do

                divFy = 0.0_rp
                !$acc loop seq
                do idime=1, ndime
                    !$acc loop seq
                    do jdime=1, ndime
                        divFy = divFy + He(idime,jdime,igaus,ielem) * gradIsoFy(idime,jdime)
                    end do
                end do                
                Re(igaus) = gpvol(1,igaus,ielem)*divFy
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
