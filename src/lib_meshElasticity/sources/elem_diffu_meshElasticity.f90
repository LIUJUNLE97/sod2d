module elem_diffu_meshElasticity

   use mod_nvtx
   use mod_numerical_params
   
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

      contains
        subroutine full_diffusion_ijk_meshElasticity(&
          nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,nu,E,Ml,Rmom,metric)
             implicit none

             integer(4), intent(in)  :: nelem, npoin
             integer(4), intent(in)  :: connec(nelem,nnode)
             real(rp),   intent(in)  :: Ngp(ngaus,nnode)
             real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
             real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
             integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
             real(rp),   intent(in)  :: u(npoin,ndime), nu,E,Ml(npoin)
             real(rp), optional, intent(in) :: metric(ndime,ndime,npoin)
             real(rp),   intent(out) :: Rmom(npoin,ndime)
             integer(4)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
             integer(4)              :: ipoin(nnode)
             real(rp)                :: mu_fgp, mu_egp,divU,nu_e,tau(ndime,ndime)
             real(rp)                :: lambda, mu, Cx(ndime), Cy(ndime), Cz(ndime)
             real(rp)                :: gradU(ndime,ndime), tmp1,vol,arho
             real(rp)                :: gradIsoU(ndime,ndime)
             real(rp)                :: divDm(ndime)
             real(rp)                :: ul(nnode,ndime), mufluidl(nnode)
             real(rp)                :: tauXl(nnode,ndime), tauYl(nnode,ndime), tauZl(nnode,ndime)
             real(rp)                :: gradRhol(nnode,ndime),muel(nnode)
             real(rp) :: t_aux,eigen,hdesired,hz_maz,hz_min
!              real(rp) :: V(3,3),vaux(3,1)

             call nvtxStartRange("Full diffusion")
             !$acc kernels
             Rmom(:,:) = 0.0_rp
             !$acc end kernels
             
             lambda =  E*nu/((1.0_rp+nu)*(1.0_rp-2.0_rp*nu))
             mu = E/(2.0_rp*(1.0_rp+nu))
             
             Cx = (/lambda+2.0_rp*mu,mu,mu/)
             Cy = (/mu,lambda+2.0_rp*mu,mu/)
             Cz = (/mu,mu,lambda+2.0_rp*mu/)
             
             ! constant basis for illustrative example
!              V(1,:)=(/1.0_rp,-1.0_rp,0.0_rp/)/sqrt(2.0_rp)
!              V(2,:)=(/0.0_rp,0.0_rp,1.0_rp/)
!              V(3,:)=(/1.0_rp,1.0_rp,0.0_rp/)/sqrt(2.0_rp)
!              V(:,1)=(/1.0_rp,-1.0_rp,0.0_rp/)/sqrt(2.0_rp)
!              V(:,2)=(/0.0_rp,0.0_rp,1.0_rp/)
!              V(:,3)=(/1.0_rp,1.0_rp,0.0_rp/)/sqrt(2.0_rp)
             
             !$acc parallel loop gang  private(ipoin,ul,tauXl,tauYl,tauZl,muel)
             do ielem = 1,nelem
                !$acc loop vector
                do inode = 1,nnode
                   ipoin(inode) = connec(ielem,inode)
                end do
                
                !$acc loop vector collapse(2)
                do inode = 1,nnode
                   do idime = 1,ndime
                      ul(inode,idime) = u(ipoin(inode),idime)
                   end do
!                    vaux = reshape(ul(inode,:), (/3,1/) )
!                    vaux =  matmul(V,vaux )
!                    ul(inode,:) = reshape(vaux,(/3/)) ! vaux(:)!reshape(vaux,(/1,3/))
                end do
                tauXl(:,:) = 0.0_rp
                tauYl(:,:) = 0.0_rp
                tauZl(:,:) = 0.0_rp
                
                !$acc loop vector private(tau,gradU,gradIsoU,divU)
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
                   
                   !gradU = matmul(V,gradU)
!                    gradU = matmul(gradU,V)
                   
                   !$acc loop seq
                   do idime = 1,ndime
                      !$acc loop seq
                      do jdime = 1,ndime
                           tau(idime,jdime) = (gradU(idime,jdime)+gradU(jdime,idime))
                      end do
                   end do
                   !tau = matmul(V,tau)  ! test here change of basis!!!
                   !
                   if(present(metric)) then
                     do idime = 1,ndime
                       tau(idime,idime)= tau(idime,idime) + 1.0_rp  ! this one is from Id from reference metric
                       do jdime = 1,ndime
                         tau(idime,jdime)= tau(idime,jdime) - metric(idime,jdime,connec(ielem,igaus))
                       end do
                     end do
!                      tau = matmul(tau,transpose(V))
!                      tau = matmul(V,tau)
                   end if

!                    t_aux = coordPar(connec(ielem,igaus),3)/6.28319_rp
!                    hz_maz = 2.0_rp !riemanian
!                    hz_min = 0.1_rp !riemanian
!                    t_aux = t_aux**1 ! quadratic transition from min (aniso) to max (iso)
!                    hdesired = hz_maz*t_aux  + hz_min*(1-t_aux) ! elem size in between [0.1 1] according to taux (height, normalized)
!                    eigen = 1.0_rp/hdesired**2
!
!                    tau(3,3)= tau(3,3) + 1.0_rp  - eigen
! !                    tau(3,3)= tau(3,3) + eigen   - 1.0_rp
!
! !                    tau(1,1)= tau(1,1) + 1.0_rp
! !                    tau(2,2)= tau(2,2) + 1.0_rp
! !                    tau(3,3)= tau(3,3) + eigen
!
! !                    tau(:,3)=tau(:,3)*eigen
! !                    tau(1,1)= tau(1,1) - 1.0_rp
! !                    tau(2,2)= tau(2,2) - 1.0_rp
! !                    tau(3,3)= tau(3,3) - 1.0_rp
!
! !                    hz_maz = 1.0_rp !riemanian
! !                    hz_min = 0.05_rp !riemanian
! !                    t_aux = sqrt(sum((coordPar(connec(ielem,igaus),:)-3.14)**2)) ! radius on center of cube
! !                    t_aux = t_aux/6.28_rp ! t\in[0,1]
! !                    !t_aux = t_aux*t_aux
! !                    hdesired = hz_maz*t_aux  + hz_min*(1-t_aux)
! !                    eigen = 1.0_rp/hdesired**2
! !                    tau(1,1)= tau(1,1) + 1.0_rp  - eigen
! !                    tau(2,2)= tau(2,2) + 1.0_rp  - eigen
! !                    tau(3,3)= tau(3,3) + 1.0_rp  - eigen
!                    ! here cross M_ij
!                    !tau(i,j)= tau(i,j) -Mij

                   !$acc loop seq
                   do idime = 1,ndime
                      tauXl(igaus,idime) =  tau(1,idime)*Cx(idime)
                      tauYl(igaus,idime) =  tau(2,idime)*Cy(idime)
                      tauZl(igaus,idime) =  tau(3,idime)*Cz(idime)
                   end do
                   ! Do we want an acc parall here?
                   tauXl(igaus, 1) = tauXl(igaus, 1) + tau(2, 2) * lambda + tau(3, 3) * lambda
                   tauYl(igaus, 2) = tauYl(igaus, 2) + tau(1, 1) * lambda + tau(3, 3) * lambda
                   tauZl(igaus, 3) = tauZl(igaus, 3) + tau(1, 1) * lambda + tau(2, 2) * lambda
                end do

                !$acc loop vector private(divDm) 
                do igaus = 1,ngaus
                   isoI = gmshAtoI(igaus) 
                   isoJ = gmshAtoJ(igaus) 
                   isoK = gmshAtoK(igaus) 

                   divDm(:) = 0.0_rp
                   
                   !$acc loop seq
                   do ii=1,porder+1
                      !$acc loop seq
                      do idime=1,ndime
                         divDm(1) = divDm(1) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauXl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDm(1) = divDm(1) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauXl(invAtoIJK(isoI,ii,isoK),idime)
                         divDm(1) = divDm(1) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauXl(invAtoIJK(isoI,isoJ,ii),idime)

                         divDm(2) = divDm(2) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauYl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDm(2) = divDm(2) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauYl(invAtoIJK(isoI,ii,isoK),idime)
                         divDm(2) = divDm(2) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauYl(invAtoIJK(isoI,isoJ,ii),idime)

                         divDm(3) = divDm(3) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauZl(invAtoIJK(ii,isoJ,isoK),idime)
                         divDm(3) = divDm(3) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauZl(invAtoIJK(isoI,ii,isoK),idime)
                         divDm(3) = divDm(3) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauZl(invAtoIJK(isoI,isoJ,ii),idime)
                      end do
                   end do

                   do idime = 1,ndime
                      !$acc atomic update
                      Rmom(ipoin(igaus),idime) = Rmom(ipoin(igaus),idime)+divDm(idime)
                      !$acc end atomic
                   end do
                end do
             end do
             !$acc end parallel loop
            call nvtxEndRange
        end subroutine full_diffusion_ijk_meshElasticity
end module elem_diffu_meshElasticity
