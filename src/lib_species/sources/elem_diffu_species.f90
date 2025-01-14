module elem_diffu_species

    use mod_nvtx
    use mod_numerical_params
    
    use mod_mpi
    use mod_mpi_mesh
    use mod_hdf5
    use mod_comms
 
    contains

    subroutine species_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho,Yk,mu_fluid,mu_e,mu_sgs,Ml,RYk,initialze,fact)
            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  ::  invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: Cp,Pr,rho(npoin), Yk(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus),Ml(npoin)
            real(rp),    intent(in)  :: mu_fluid(npoin)
            real(rp),    intent(inout) :: RYk(npoin)
            logical, optional, intent(in)    :: initialze
            real(rp), optional, intent(in)  :: fact
            integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
            real(rp)                 :: kappa_y, mu_fgp
            real(rp)                 :: gradY(ndime),tmp1,vol,arho
            real(rp)                 :: gradIsoY(ndime)
            real(rp)                 :: divDy
            real(rp)                 :: rhol(nnode),ykl(nnode),mufluidl(nnode),muel(nnode)
            real(rp)                 :: gradYl(nnode,ndime)
            real(rp)  :: aux_fact = 1.0_rp

            call nvtxStartRange("Full diffusion")

            if(present(initialze)) then
                if (initialze .eqv. .true.) then
                    !$acc kernels
                    RYk(:) = 0.0_rp
                    !$acc end kernels
                end if
             else
                !$acc kernels
                RYk(:) = 0.0_rp
                !$acc end kernels
             end if
             if(present(fact)) then
                aux_fact = fact
             end if
 
            !$acc parallel loop gang  private(ykl,rhol,mufluidl,muel,gradYl)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                    rhol(inode) = rho(connec(ielem,inode))
                    ykl(inode) = Yk(connec(ielem,inode))
                    mufluidl(inode) = mu_fluid(connec(ielem,inode))/Pr + rhol(inode)*mu_sgs(ielem,inode)/0.9_rp
                    muel(inode) = mu_e(ielem,inode)/Pr
                end do
                
                gradYl(:,:) = 0.0_rp
                !$acc loop vector private(gradY,gradIsoY)
                do igaus = 1,ngaus
                    isoI = gmshAtoI(igaus) 
                    isoJ = gmshAtoJ(igaus) 
                    isoK = gmshAtoK(igaus) 

                    gradIsoY(:) = 0.0_rp
                    !$acc loop seq
                    do ii=1,porder+1
                        gradIsoY(1) = gradIsoY(1) + dlxigp_ip(igaus,1,ii)*ykl(invAtoIJK(ii,isoJ,isoK))
                        gradIsoY(2) = gradIsoY(2) + dlxigp_ip(igaus,2,ii)*ykl(invAtoIJK(isoI,ii,isoK))
                        gradIsoY(3) = gradIsoY(3) + dlxigp_ip(igaus,3,ii)*ykl(invAtoIJK(isoI,isoJ,ii))
                    end do

                    gradY(:) = 0.0_rp
                    !$acc loop seq
                    do idime=1, ndime
                        !$acc loop seq
                        do jdime=1, ndime
                            gradY(idime) = gradY(idime) + He(idime,jdime,igaus,ielem) * gradIsoY(jdime)
                        end do
                    end do

                    !$acc loop seq
                    do idime = 1,ndime
                        gradYl(igaus,idime) =  gradY(idime)
                    end do
                end do

                !$acc loop vector private(divDy) 
                do igaus = 1,ngaus
                    kappa_y = (mufluidl(igaus) + muel(igaus))

                    isoI = gmshAtoI(igaus) 
                    isoJ = gmshAtoJ(igaus) 
                    isoK = gmshAtoK(igaus) 

                    divDy = 0.0_rp
                    
                    !$acc loop seq
                    do ii=1,porder+1
                        !$acc loop seq
                        do idime=1,ndime
                            divDy = divDy + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*(kappa_y*gradYl(invAtoIJK(ii,isoJ,isoK),idime))
                            divDy = divDy + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*(kappa_y*gradYl(invAtoIJK(isoI,ii,isoK),idime))
                            divDy = divDy + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*(kappa_y*gradYl(invAtoIJK(isoI,isoJ,ii),idime))
                        end do
                    end do

                    !$acc atomic update
                    RYk(connec(ielem,igaus)) = RYk(connec(ielem,igaus))+aux_fact*divDy
                    !$acc end atomic
                end do
            end do
            !$acc end parallel loop
        call nvtxEndRange
    end subroutine species_diffusion_ijk

end module elem_diffu_species
