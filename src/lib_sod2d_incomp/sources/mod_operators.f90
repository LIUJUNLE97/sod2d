module mod_operators

        use mod_nvtx
        use mod_numerical_params
   
        use mod_mpi
        use mod_mpi_mesh
        use mod_hdf5
        use mod_comms
        use mod_solver

    implicit none

      logical :: allocate_memory_mod_oprators = .true.
      !for eval_laplacian_mult
      real(rp),allocatable :: op_auxP(:)

       contains

        subroutine eval_laplacian_mult(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,p,ResP)

            implicit none

            integer(4),intent(in)    :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),  intent(inout) :: ResP(npoin)
            real(rp),  intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),p(npoin)
            integer(4),intent(in)    :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4)               :: ii, ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode)
            !real(rp)                :: p_l(npoin),aux1,auxP(npoin)
            real(rp)                :: gradIsoP(ndime)
            real(rp)                :: gradP(ndime),divDp
            real(rp)                :: pl(nnode),gradPl(nnode,ndime)

            if(allocate_memory_mod_oprators) then
               allocate_memory_mod_oprators = .false.
               
               allocate(op_auxP(npoin))
               !$acc enter data create(op_auxP(:))
            end if

            call nvtxStartRange("Eval laplacian mult")
            !$acc kernels
             ResP(:) = 0.0_rp
             op_auxP(:) = p(:)
            !$acc end kernels

            !$acc parallel loop gang  private(ipoin,pl,gradPl)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                   pl(inode) = op_auxP(ipoin(inode))
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
            call nvtxEndRange

            if(flag_fs_fix_pressure) then
               !$acc kernels
               ResP(inode_fix_press) = 0.0_rp              
               !$acc end kernels
            end if
        end subroutine eval_laplacian_mult

        subroutine eval_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,invMl,x,gradX,lump)

            implicit none

            integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
            real(rp),   intent(inout) :: gradX(npoin,ndime)
            real(rp),   intent(in)    :: invMl(npoin)
            real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),x(npoin)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            logical,    intent(in)  :: lump
            integer(4)              :: ii, ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode)
            real(rp)                :: gradIsoP(ndime)
            real(rp)                :: gradP(ndime),divDp
            real(rp)                :: pl(nnode)

            call nvtxStartRange("Eval grad")
            !$acc kernels
             gradX(:,:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang  private(ipoin,pl)
            do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                   pl(inode) = x(ipoin(inode))
                end do
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
                     !$acc atomic update
                     gradX(ipoin(igaus),idime) = gradX(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*gradP(idime)
                     !$acc end atomic
                   end do
               end do
            end do
            !$acc end parallel loop

            if(lump .eqv. .true.) then
               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real_arrays(ndime,gradX(:,:))
                  call nvtxEndRange
               end if
               
               !
               ! Call lumped mass matrix solver
               !
            
               call nvtxStartRange("Call solver")
               call lumped_solver_vect_opt(npoin,npoin_w,lpoin_w,invMl,gradX(:,:))
               call nvtxEndRange
            end if
            call nvtxEndRange

    end subroutine eval_gradient

    subroutine eval_u_gradient(nelem,npoin,npoin_w,connec,lpoin_w,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Ml,u,x,gradX,lump)

      implicit none

      integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
      real(rp),   intent(inout) :: gradX(npoin,ndime)
      real(rp),   intent(in)    :: Ml(npoin),u(npoin,ndime)
      real(rp),   intent(in)    :: dlxigp_ip(ngaus,ndime,porder+1),He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),x(npoin)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      logical,    intent(in)  :: lump
      integer(4)              :: ii, ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,ipoin(nnode)
      real(rp)                :: gradIsoP(ndime)
      real(rp)                :: gradP(ndime),divDp
      real(rp)                :: pl(nnode),ul(nnode,ndime)

      call nvtxStartRange("Eval grad")
      !$acc kernels
       gradX(:,:) = 0.0_rp
      !$acc end kernels

      !$acc parallel loop gang  private(ipoin,pl,ul)
      do ielem = 1,nelem
          !$acc loop vector
          do inode = 1,nnode
            ipoin(inode) = connec(ielem,inode)
             pl(inode) = x(ipoin(inode))
             !$acc loop seq
             do idime = 1, ndime
               ul(inode,idime) = u(ipoin(inode),idime)
             end do
          end do
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
               !$acc atomic update
               gradX(ipoin(igaus),idime) = gradX(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*gradP(idime)*ul(igaus,idime)
               !$acc end atomic
             end do
         end do
      end do
      !$acc end parallel loop

      if(lump .eqv. .true.) then
         if(mpi_size.ge.2) then
            call nvtxStartRange("MPI_comms_tI")
            call mpi_halo_atomic_update_real_arrays(ndime,gradX(:,:))
            call nvtxEndRange
         end if
         
         !
         ! Call lumped mass matrix solver
         !
      
         call nvtxStartRange("Call solver")
         call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,gradX(:,:))
         call nvtxEndRange
      end if
      call nvtxEndRange

   end subroutine eval_u_gradient


    subroutine eval_tau_veloc(nelem,npoin,npoin_w,connec,lpoin_w,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,Ml,GradX,GradY,GradZ)
      implicit none

      integer(4), intent(in)  :: nelem, npoin,npoin_w
      integer(4), intent(in)  :: connec(nelem,nnode),lpoin_w(npoin_w)
      real(rp),   intent(in)  :: Ngp(ngaus,nnode)
      real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
      real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: u(npoin,ndime),Ml(npoin)
      real(rp),   intent(out) :: GradX(npoin,ndime), GradY(npoin,ndime), GradZ(npoin,ndime)
      integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
      integer(4)              :: ipoin(nnode)
      real(rp)                :: gradU(ndime,ndime), tmp1,vol,arho
      real(rp)                :: gradIsoU(ndime,ndime),tau(ndime,ndime)
      real(rp)                :: divDm(ndime)
      real(rp)                :: ul(nnode,ndime)
      real(rp)                :: gradRhol(nnode,ndime)

      call nvtxStartRange("Full diffusion")
      !$acc kernels
      GradX(:,:) = 0.0_rp
      GradY(:,:) = 0.0_rp
      GradZ(:,:) = 0.0_rp
      !$acc end kernels

      !$acc parallel loop gang  private(ipoin,ul)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            ipoin(inode) = connec(ielem,inode)
            !$acc loop seq
            do idime = 1,ndime
               ul(inode,idime) = u(ipoin(inode),idime)
            end do
         end do

         !$acc loop vector private(gradU,gradIsoU,tau)
         do igaus = 1,ngaus

            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoU(:,:) = 0.0_rp
            !$acc loop seq
            do ii=1,porder+1
               !$acc loop seq
               do idime=1,ndime
                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            gradU(:,:) = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  !$acc loop seq
                  do kdime=1,ndime
                     gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                  end do
               end do
            end do

            !$acc loop seq
            do idime = 1,ndime
               !$acc loop seq
               do jdime = 1,ndime
                    tau(idime,jdime) = (gradU(idime,jdime)+gradU(jdime,idime))
               end do
            end do

            !$acc loop seq
            do idime = 1,ndime
               !$acc atomic update
               GradX(ipoin(igaus),idime) = GradX(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*tau(1,idime)
               !$acc end atomic
               !$acc atomic update
               GradY(ipoin(igaus),idime) = GradY(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*tau(2,idime)
               !$acc end atomic
               !$acc atomic update
               GradZ(ipoin(igaus),idime) = GradZ(ipoin(igaus),idime)+gpvol(1,igaus,ielem)*tau(3,idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop
     call nvtxEndRange

     if(mpi_size.ge.2) then
      call nvtxStartRange("MPI_comms_tI")
      call mpi_halo_atomic_update_real_arrays(ndime,GradX(:,:))
      call mpi_halo_atomic_update_real_arrays(ndime,GradY(:,:))
      call mpi_halo_atomic_update_real_arrays(ndime,GradZ(:,:))
      call nvtxEndRange
     end if
   
      !
      ! Call lumped mass matrix solver
      !

      call nvtxStartRange("Call solver")
      call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,gradX(:,:))
      call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,gradY(:,:))
      call lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,gradZ(:,:))
      call nvtxEndRange

   end subroutine eval_tau_veloc

    subroutine eval_divergence(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,Rmom)

            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),   intent(in)  :: u(npoin,ndime)
            real(rp),   intent(out) :: Rmom(npoin)
            integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                :: divU,Re_mom(nnode),gradIsoU(ndime,ndime),ul(nnode,ndime)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Eval div")
            !$acc kernels
            Rmom(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(ipoin,Re_mom,ul)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                  !$acc loop seq
                  do idime = 1,ndime                  
                     ul(inode,idime)  = u(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip,gradIsoU,divU)
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

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(Rmom(:))
               call nvtxEndRange
            end if
            call nvtxEndRange
    end subroutine eval_divergence

    subroutine eval_laplacian_BDL(nelem,npoin,connec,He,dNgp,invAtoIJK,gpvol,diag,L)

            implicit none

            integer(4), intent(in)   :: nelem, npoin, connec(nelem,nnode)
            real(rp),   intent(inout)  :: L(nnode,nnode,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1)
            real(rp),   intent(in)    :: He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),dNgp(ndime,nnode,ngaus),diag(npoin)
            integer(4)               :: ielem, igaus, idime, inode,jnode,knode
            real(rp)                  :: gpcar(ndime,nnode),tmp,A(nnode,nnode)

            call nvtxStartRange("Eval laplacian BDL")
            !$acc parallel loop gang  private(gpcar,A)
            do ielem = 1,nelem
               !$acc loop vector collapse(2)
               do inode = 1,nnode
                  do jnode = 1,nnode
                     A(inode,jnode) = 0.0_rp
                     L(inode,jnode,ielem) = 0.0_rp
                  end do
               end do 
               !$acc loop seq
               do igaus = 1,ngaus
                  !$acc loop vector collapse(2)
                  do idime = 1,ndime
                     do inode = 1,nnode
                        gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                     end do
                  end do
                  !$acc loop vector collapse(2)
                  do inode = 1,nnode
                     do jnode = 1,nnode
                        A(inode,jnode) = A(inode,jnode) + gpvol(1,igaus,ielem)*dot_product(gpcar(:,inode),gpcar(:,jnode))
                     end do
                  end do
               end do
               !$acc loop vector 
               do inode = 1,nnode
                  A(inode,inode) = diag(connec(ielem,inode))   
               end do

               L(1,1,ielem) = sqrt(A(1,1))
               !$acc loop vector 
               do inode = 2,nnode
                  L(inode,1,ielem) = A(inode,1)/L(1,1,ielem)
               end do
               !$acc loop vector 
               do jnode = 2,nnode
                  knode = (jnode-1)
                  tmp = A(jnode,jnode) - dot_product(A(jnode,1:knode),A(jnode,1:knode))
                  L(jnode,jnode,ielem) = sqrt(tmp)
               end do
               !$acc loop seq
               do jnode = 2,nnode-1
                  !$acc loop seq
                  do inode = jnode+1,nnode
                     knode = (jnode-1)
                     L(inode,jnode,ielem) = (A(inode,jnode) - dot_product(L(inode,1:knode,ielem),L(jnode,1:knode,ielem)))/L(jnode,jnode,ielem)
                  end do
               end do
            end do
            !$acc end parallel loop
            call nvtxEndRange

    end subroutine eval_laplacian_BDL

        subroutine eval_laplacian_diag(nelem,npoin,connec,He,dNgp,gpvol,d)

            implicit none

            integer(4), intent(in)   :: nelem,npoin, connec(nelem,nnode)
            real(rp),   intent(inout)  :: d(npoin)
            real(rp),   intent(in)    :: He(ndime,ndime,ngaus,nelem),gpvol(1,ngaus,nelem),dNgp(ndime,nnode,ngaus)
            integer(4)               :: ielem, igaus, idime, inode,jnode,knode
            real(rp)                  :: gpcar(ndime,nnode),tmp

            call nvtxStartRange("Eval lapl diag")
            !$acc kernels
            d(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang  private(gpcar)
            do ielem = 1,nelem
               !$acc loop seq
               do igaus = 1,ngaus
                  !$acc loop vector collapse(2)
                  do idime = 1,ndime
                     do inode = 1,nnode
                        gpcar(idime,inode) = dot_product(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                     end do
                  end do
                  !$acc loop vector collapse(2)
                  do idime = 1,ndime
                     do inode = 1,nnode
                        !$acc atomic update
                        d(connec(ielem,inode)) = d(connec(ielem,inode)) +  gpvol(1,igaus,ielem)*gpcar(idime,inode)*gpcar(idime,inode)
                        !$acc end atomic
                     end do
                  end do
               end do
            end do
            !$acc end parallel loop
            
            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(d(:))
               call nvtxEndRange
            end if
            call nvtxEndRange

    end subroutine eval_laplacian_diag

   subroutine compute_vorticity(nelem,npoin,npoin_w,lpoin_w,connec,lelpn,He,dNgp,leviCivi,dlxigp_ip,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,curlU,makeConti)

      implicit none

      integer(4), intent(in)  :: nelem, npoin, npoin_w, lpoin_w(npoin), connec(nelem,nnode), lelpn(npoin)
      real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem), dNgp(ndime,nnode,ngaus)
      real(rp),   intent(in)  :: u(npoin,ndime)
      real(rp),   intent(in)  :: leviCivi(ndime,ndime,ndime)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1)
      integer(4), intent(in)  :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(out) :: curlU(npoin,ndime)
      logical,    intent(in)  :: makeConti
      integer(4)              :: iNodeL, ielem, idime, inode, igaus, ipoin, jdime,kdime,isoI, isoJ, isoK , ii
      real(rp)                :: gpcar(ndime,nnode), u_e(nnode,ndime), aux1, aux2
      real(rp)                :: gradu_e(ndime,ndime)
      real(rp)                :: gradIsoU(ndime,ndime)

      !$acc kernels
      curlU(:,:) = 0.0_rp
      !$acc end kernels
      !$acc parallel loop gang private(u_e)
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            u_e(inode,1) = u(connec(ielem,inode),1)
            u_e(inode,2) = u(connec(ielem,inode),2)
            u_e(inode,3) = u(connec(ielem,inode),3)
         end do
         !$acc loop vector private(gradu_e,gradIsoU)
         do igaus = 1,ngaus
            isoI = gmshAtoI(igaus) 
            isoJ = gmshAtoJ(igaus) 
            isoK = gmshAtoK(igaus) 

            gradIsoU(:,:) = 0.0_rp
            !$acc loop seq
            do ii=1,porder+1
               !$acc loop seq
               do idime=1,ndime
                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*u_e(invAtoIJK(ii,isoJ,isoK),idime)
                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*u_e(invAtoIJK(isoI,ii,isoK),idime)
                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*u_e(invAtoIJK(isoI,isoJ,ii),idime)
               end do
            end do

            gradu_e(:,:) = 0.0_rp
            !$acc loop seq
            do idime=1, ndime
               !$acc loop seq
               do jdime=1, ndime
                  !$acc loop seq
                  do kdime=1,ndime
                     gradu_e(idime,jdime) = gradu_e(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                  end do
               end do
            end do

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
               end do
            end do
         end do
      end do
      !$acc end parallel loop

      if(mpi_size.ge.2) then
         call nvtxStartRange("MPI_comms_post")
         call mpi_halo_atomic_update_real_arrays(ndime,curlU(:,:))
         call nvtxEndRange
      end if

      if(makeConti .eqv. .true.) then
         !$acc parallel loop
         do ipoin = 1,npoin_w
            iNodeL=lpoin_w(ipoin)
            !$acc loop seq
            do idime = 1,ndime
               curlU(iNodeL,idime) = curlU(iNodeL,idime)/real(lelpn(iNodeL),rp)
            end do
         end do
         !$acc end parallel loop
      end if

   end subroutine compute_vorticity    
end module mod_operators