module mod_entropy_viscosity_species

    use mod_numerical_params
    use mod_nvtx
    
    use mod_mpi
    use mod_mpi_mesh
    use mod_comms
    use mod_filters,only:convertIJK,al_weights,am_weights,an_weights
 
    implicit none
 
    contains
 
    subroutine species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Ngp,coord,dNgp,gpvol,wgp, &
        rho,u,eta,helem,helem_k,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)

        ! TODO: Compute element size h

        implicit none

        integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
        real(rp),    intent(in)  :: Reta(npoin), Ngp(ngaus,nnode)
        real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime), eta(npoin),helem(nelem,nnode),helem_k(nelem),Ml(npoin)
        real(rp),    intent(out) :: mu_e(nelem,ngaus)
        real(rp),   intent(in)  :: coord(npoin,ndime), dNgp(ndime,nnode,ngaus), wgp(ngaus)
        real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
        integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
        real(rp),intent(inout)  :: mue_l(nelem,nnode)
        integer(4)               :: ielem, inode,igaus,ipoin,npoin_w_g,idime,jdime
        real(rp)                 :: R1, R2, Ve
        real(rp)                 :: betae,mu,vol,vol2
        real(rp)                 :: L3, aux1
        real(rp)                 :: maxEta_r,maxEta, maxRho, norm_r,norm,maxV, maxC
        real(rp)                :: Je(ndime,ndime), maxJe, minJe,ced,magJe, M, ceM
        integer(4)              :: ii,jj,kk,mm,nn,ll

        !$acc kernels
        mue_l(:,:) = mu_e(:,:)
        !$acc end kernels

        if(flag_normalise_entropy .eq. 1) then
            maxEta_r = 0.0_rp
            !$acc parallel loop reduction(+:maxEta_r)
            do ipoin = 1,npoin_w
                maxEta_r = maxEta_r + eta(lpoin_w(ipoin))
            end do
            !$acc end parallel loop

            call MPI_Allreduce(maxEta_r,maxEta,1,mpi_datatype_real,MPI_SUM,MPI_COMM_WORLD,mpi_err)
            call MPI_Allreduce(npoin_w,npoin_w_g,1,mpi_datatype_int,MPI_SUM,MPI_COMM_WORLD,mpi_err)

            maxEta = maxEta/real(npoin_w_g,rp)

            norm_r = 0.0_rp
            !$acc parallel loop reduction(max:norm_r)
            do ipoin = 1,npoin_w
                norm_r = max(norm_r, abs(eta(lpoin_w(ipoin))-maxEta))
            end do
            !$acc end parallel loop

            call MPI_Allreduce(norm_r,norm,1,mpi_datatype_real,MPI_MAX,MPI_COMM_WORLD,mpi_err)
        else
            norm = 1.0_rp
        end if

        !$acc parallel loop gang 
        do ielem = 1,nelem
            maxJe=0.0_rp
            minJe=1000000.0_rp
            !$acc loop seq
            do igaus = 1,ngaus
                minJe = min(minJe,gpvol(1,igaus,ielem)/wgp(igaus))
                maxJe = max(maxJe,gpvol(1,igaus,ielem)/wgp(igaus))
            end do
            ced = max(1.0_rp-(minJe/maxJe)**2,ce)

            mu = 0.0_rp
            betae = 0.0_rp
            !$acc loop vector reduction(max:betae)
            do inode = 1,nnode
                aux1 = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
                betae = max(betae,(rho(connec(ielem,inode))*helem_k(ielem))*(cmax/real(porder,rp))*aux1)
            end do
            !$acc loop vector
            do inode = 1,nnode
                R1 = rho(connec(ielem,inode))*abs(Reta(connec(ielem,inode)))/norm
                Ve = ced*R1*(helem(ielem,inode))**2
                mue_l(ielem,inode) = cglob*min(Ve,betae)
            end do
            !$acc loop vector collapse(3)
            do ii=1,porder+1
                do jj=1,porder+1
                    do kk=1,porder+1
                        aux1 = 0.00_rp
                        !$acc loop seq
                        do ll=-1,1
                            !$acc loop seq
                            do mm=-1,1
                                !$acc loop seq
                                do nn=-1,1
                                    aux1 =   aux1 +  al_weights(ll)*am_weights(mm)*an_weights(nn)*mue_l(ielem,invAtoIJK(convertIJK(ii+ll),convertIJK(jj+mm),convertIJK(kk+nn)))
                                end do
                            end do
                        end do
                        mu_e(ielem,invAtoIJK(convertIJK(ii),convertIJK(jj),convertIJK(kk))) = aux1
                    end do
                end do
            end do
        end do
        !$acc end parallel loop

    end subroutine species_smart_visc_spectral


    
end module mod_entropy_viscosity_species
