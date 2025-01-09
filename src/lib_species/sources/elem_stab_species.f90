module elem_stab_species

    use mod_nvtx
    use mod_numerical_params
    
    use mod_mpi
    use mod_mpi_mesh
    use mod_hdf5
    use mod_comms
 
    contains

    subroutine species_stab_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Yk,gradYk,Cp,Prt,rho,tau,Ml,RYk,initialze,fact)
            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  ::  invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: Yk(npoin), gradYk(npoin,ndime), tau(nelem),Ml(npoin),rho(npoin),Cp,Prt
            real(rp),    intent(inout) :: RYk(npoin)
            logical, optional, intent(in)    :: initialze
            real(rp), optional, intent(in)  :: fact
            integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
            real(rp)                 :: kappa_y, mu_fgp
            real(rp)                 :: gradY(ndime),tmp1,vol
            real(rp)                 :: gradIsoY(ndime)
            real(rp)                 :: divDy
            real(rp)                 :: ykl(nnode),gradykl(nnode,ndime),rhol(nnode)
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
 
            !$acc parallel loop gang  private(ykl,gradykl,gradYl,kappa_y,rhol)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                    ykl(inode)  = Yk(connec(ielem,inode))
                    rhol(inode) = rho(connec(ielem,inode))
                    !$acc loop seq
                    do idime=1,ndime
                        gradykl(inode,idime) = gradYk(connec(ielem,inode),idime)
                    end do
                end do
                
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
                        gradYl(igaus,idime) =  (gradykl(igaus,idime) - gradY(idime))
                    end do
                end do

                !$acc loop vector private(divDy) 
                do igaus = 1,ngaus
                    isoI = gmshAtoI(igaus) 
                    isoJ = gmshAtoJ(igaus) 
                    isoK = gmshAtoK(igaus) 

                    divDy = 0.0_rp
                    kappa_y = 0.1_rp*tau(ielem)/Prt
                    
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
    end subroutine species_stab_ijk


    subroutine species_tau(nelem,npoin,connec,u,helem_k,dt,tau)

        ! TODO: Compute element size h

        implicit none

        integer(4), intent(in)   :: nelem, npoin,connec(nelem,nnode)
        real(rp),    intent(in)  :: u(npoin,ndime),helem_k(nelem),dt
        real(rp),    intent(out) :: tau(nelem)
        integer(4)               :: ielem, inode
        real(rp)                 :: taul
        real(rp)                 :: aux1


        !$acc parallel loop gang 
        do ielem = 1,nelem
            taul = 0.0_rp
            !$acc loop vector reduction(max:taul)
            do inode = 1,nnode
                aux1 = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
                taul = max(taul,(helem_k(ielem))*(c_species_stab/real(porder,rp))*aux1)
            end do
            tau(ielem) = taul
        end do
        !$acc end parallel loop

    end subroutine species_tau


end module elem_stab_species
