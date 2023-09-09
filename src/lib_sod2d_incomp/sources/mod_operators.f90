module mod_operators

        use mod_nvtx
        use mod_numerical_params
   
        use mod_mpi
        use mod_mpi_mesh
        use mod_hdf5
        use mod_comms
        use mod_solver

    implicit none

       contains

        subroutine eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,p,ResP)

            implicit none

            integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),   intent(inout)  :: ResP(npoin)
            real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),p(npoin)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4)               :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode),ipoin2
            integer(4)              :: convertIJK(0:porder+2),ii,jj,kk,mm,nn,ll,iter
            real(rp)                :: p_l(npoin),al(-1:1),am(-1:1),an(-1:1),aux1,auxP(npoin)
            real(rp)                :: gradIsoP(ndime)
            real(rp)                :: gradP(ndime),divDp
            real(rp)                :: pl(nnode),gradPl(nnode,ndime)

            !$acc kernels
             ResP(:) = 0.0_rp
             auxP(:) = p(:)
            !$acc end kernels
            if((mpi_rank.eq.0) .and. (flag_fs_fix_pressure .eqv. .true.)) then
               auxP(lpoin_w(1)) = 0.0_rp
            end if
            !$acc parallel loop gang  private(ipoin,pl,gradPl)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                   ipoin(inode) = connec(ielem,inode)
                end do
                !$acc loop vector
                do inode = 1,nnode
                   pl(inode) = auxP(ipoin(inode))
                end do
                gradPl(:,:) = 0.0_rp

                !$acc loop vector private(gradIsoP,gradP)
                do igaus = 1,ngaus

                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   gradIsoP(:) = 0.0_rp
                   !$acc loop seq
                   do ii=1,porder+1
                      gradIsoP(1) = gradIsoP(1) + dlxigp_ip(igaus,1,ii)*pl(invAtoIJK(ii,isoJ,isoK))
                      gradIsoP(2) = gradIsoP(2) + dlxigp_ip(igaus,2,ii)*pl(invAtoIJK(isoI,ii,isoK))
                      gradIsoP(3) = gradIsoP(3) + dlxigp_ip(igaus,3,ii)*pl(invAtoIJK(isoI,isoJ,ii))
                   end do

                   gradP(:) = 0.0_rp
                   !$acc loop seq
                   do idime=1, ndime
                      !$acc loop seq
                      do jdime=1, ndime
                         gradP(idime) = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                      end do
                   end do

                   !$acc loop seq
                   do idime = 1,ndime
                      gradPl(igaus,idime) =  gradP(idime)
                   end do
                end do

                !$acc loop vector private(divDp) 
                do igaus = 1,ngaus

                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   divDp = 0.0_rp
                   
                   !$acc loop seq
                   do ii=1,porder+1
                      !$acc loop seq
                      do idime=1,ndime
                         divDp = divDp + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*gradPl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDp = divDp + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*gradPl(invAtoIJK(isoI,ii,isoK),idime)
                         divDp = divDp + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*gradPl(invAtoIJK(isoI,isoJ,ii),idime)
                      end do
                   end do

                   !$acc atomic update
                   ResP(ipoin(igaus)) = ResP(ipoin(igaus))+divDp
                   !$acc end atomic
                end do
             end do
             !$acc end parallel loop

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(ResP(:))
               call nvtxEndRange
            end if
            if((mpi_rank.eq.0) .and. (flag_fs_fix_pressure .eqv. .true.)) then
               ResP(lpoin_w(1)) = 0.0_rp
            end if
        end subroutine eval_laplacian_mult

        subroutine eval_laplacian_diag(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,ResP)

            implicit none

            integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),   intent(inout)  :: ResP(npoin)
            real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4)               :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode),ipoin2
            integer(4)              :: convertIJK(0:porder+2),ii,jj,kk,mm,nn,ll,iter
            real(rp)                :: p_l(npoin),al(-1:1),am(-1:1),an(-1:1),aux1
            real(rp)                :: gradIsoP(ndime)
            real(rp)                :: gradP(ndime),divDp

            !$acc kernels
             ResP(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang  private(ipoin)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                   ipoin(inode) = connec(ielem,inode)
                end do

                !$acc loop vector private(gradIsoP,gradP)
                do igaus = 1,ngaus
                  isoI = gmshAtoI(igaus) 
                  isoJ = gmshAtoJ(igaus) 
                  isoK = gmshAtoK(igaus) 

                  gradIsoP(1) = dlxigp_ip(igaus,1,isoI)
                  gradIsoP(2) = dlxigp_ip(igaus,2,isoJ)
                  gradIsoP(3) = dlxigp_ip(igaus,3,isoK)                  

                  gradP(:) = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                        gradP(idime) = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                     end do
                  end do

                  !write(*,*) "Grad Pres ",gradP

                  divDp = 0.0_rp
                  !$acc loop seq
                  do idime=1,ndime
                     divDp = divDp + gradP(idime)*gradP(idime)
                  end do
                  !write(*,*) "AA ",divDp*gpvol(1,igaus,ielem)
                  !$acc atomic update
                  ResP(ipoin(igaus)) = ResP(ipoin(igaus))+gpvol(1,igaus,ielem)*divDp
                  !$acc end atomic
                  !write(*,*) "ResP ",ResP(ipoin(igaus))
                end do
             end do
             !$acc end parallel loop

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(ResP(:))
               call nvtxEndRange
            end if

            if((mpi_rank.eq.0) .and. (flag_fs_fix_pressure .eqv. .true.)) then
               ResP(lpoin_w(1)) = 1.0_rp
            end if


        end subroutine eval_laplacian_diag

        subroutine eval_laplacian_mult2(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,p,ResP)

            implicit none

            integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),   intent(inout)  :: ResP(npoin)
            real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),p(npoin),Ml(npoin)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp)                :: auxP(npoin),auxGradP(npoin,ndime)
            !$acc kernels
             ResP(:) = 0.0_rp
             auxP(:) = p(:)
            !$acc end kernels
            if(mpi_rank.eq.0) then
               auxP(lpoin_w(1)) = nscbc_p_inf
            end if

            call eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,auxP,auxGradP,.true.)
            call eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,auxGradP,ResP)

            if(mpi_rank.eq.0) then
               ResP(lpoin_w(1)) = 0.0_rp
            end if
        end subroutine eval_laplacian_mult2


        subroutine eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,x,gradX,lump)

            implicit none

            integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),   intent(inout)  :: gradX(npoin,ndime)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),x(npoin)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            logical,    intent(in)  :: lump
            integer(4)               :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode),ipoin2
            integer(4)              :: convertIJK(0:porder+2),ii,jj,kk,mm,nn,ll,iter
            real(rp)                :: al(-1:1),am(-1:1),an(-1:1),aux1
            real(rp)                :: gradIsoP(ndime)
            real(rp)                :: gradP(ndime),divDp
            real(rp)                :: pl(nnode),gradPl(nnode,ndime)

            !$acc kernels
             gradX(:,:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang  private(ipoin,pl,gradPl)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                   ipoin(inode) = connec(ielem,inode)
                end do
                !$acc loop vector
                do inode = 1,nnode
                   pl(inode) = x(ipoin(inode))
                end do
                gradPl(:,:) = 0.0_rp
                !$acc loop vector private(gradIsoP,gradP)
                do igaus = 1,ngaus
                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 
                   gradIsoP(:) = 0.0_rp
                   !$acc loop seq
                   do ii=1,porder+1
                      gradIsoP(1) = gradIsoP(1) + dlxigp_ip(igaus,1,ii)*pl(invAtoIJK(ii,isoJ,isoK))
                      gradIsoP(2) = gradIsoP(2) + dlxigp_ip(igaus,2,ii)*pl(invAtoIJK(isoI,ii,isoK))
                      gradIsoP(3) = gradIsoP(3) + dlxigp_ip(igaus,3,ii)*pl(invAtoIJK(isoI,isoJ,ii))
                   end do

                   gradP(:) = 0.0_rp
                   !$acc loop seq
                   do idime=1, ndime
                      !$acc loop seq
                      do jdime=1, ndime
                         gradP(idime) = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                      end do
                   end do

                   !$acc loop seq
                   do idime = 1,ndime
                      gradPl(igaus,idime) =  gpvol(1,igaus,ielem)*gradP(idime)
                   end do
               end do
               !
               ! Final assembly
               !
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     !$acc atomic update
                     gradX(ipoin(inode),idime) = gradX(ipoin(inode),idime)+gradPl(inode,idime)
                     !$acc end atomic
                  end do
               end do
            end do
            !$acc end parallel loop

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               do idime = 1,ndime
                  call mpi_halo_atomic_update_real(gradX(:,idime))
               end do
               call nvtxEndRange
            end if
               
            !
            ! Call lumped mass matrix solver
            !
            if(lump .eqv. .true.) then
               call nvtxStartRange("Call solver")
               call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,gradX(:,:))
               call nvtxEndRange
            end if

    end subroutine eval_gradient

    subroutine eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,Rmom)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)   :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: u(npoin,ndime)
            real(rp),    intent(out) :: Rmom(npoin)
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode)
            real(rp)                 :: gradIsoU(ndime,ndime)
            real(rp)                 :: divU
            real(rp)                 :: ul(nnode,ndime)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Full convection")
            !$acc kernels
            Rmom(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(ipoin,Re_mom,ul)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
               end do
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ul(inode,idime)  = u(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoU,divU)
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
                  !$acc loop seq
                  do ii=1,porder+1
                     
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
  
                     end do
                  end do

                  divU = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         divU = divU + He(idime,jdime,igaus,ielem) * gradIsoU(idime,jdime)
                     end do
                  end do
                  Re_mom(igaus) = gpvol(1,igaus,ielem)*divU
               end do
               !
               ! Final assembly
               !
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rmom(ipoin(inode)) = Rmom(ipoin(inode))+Re_mom(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(Rmom(:))
               call nvtxEndRange
            end if
    end subroutine eval_divergence
end module mod_operators