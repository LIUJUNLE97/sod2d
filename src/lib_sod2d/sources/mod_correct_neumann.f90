module mod_correct_neumann

   use mod_numerical_params
   use mod_veclen
   use mod_nvtx
   use mod_mpi

contains

   subroutine evalCorrectNeumann(nelem,npoin,nboun,connec,bound,point2elem,atoIJK, bou_code, &
         bounorm,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,mu_fluid,ui,Rdiff)

      implicit none

      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin),atoIJK(nnode)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: ui(npoin,ndime),mu_fluid(npoin)
      real(rp),   intent(in)  :: coord(npoin,ndime)
      real(rp),   intent(inout) :: Rdiff(npoin,ndime)
      real(rp)                :: gradIsoV(ndime), gradV(ndime),vl(nnode),auxmag,mul(nnode)
      integer(4)              :: iBound,igaus,inode,ielem,ii, isoI, isoJ, isoK,jdime,idime,bcode
      real(rp)                :: bnorm(npbou*ndime),val
      real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

      !$acc parallel loop gang private(bnorm,vl)
      do iBound = 1,nboun
         bcode = bou_code(iBound) 
         if (bcode == bc_type_far_field_SB) then
            bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
            ielem = point2elem(bound(iBound,npbou))
            !$acc loop vector
            do inode = 1,nnode
                  vl(inode)  = ui(connec(ielem,inode),2)
                  mul(inode) = mu_fluid(connec(ielem,inode))
            end do

            val = 0.0_rp
            !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip,gradIsoV,gradV) reduction(+:val)
            do igaus = 1,npbou
               !$acc loop seq
               do ii=1,porder+1
                  dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                  dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                  dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
               end do
               isoI = gmshAtoI(igaus) 
               isoJ = gmshAtoJ(igaus) 
               isoK = gmshAtoK(igaus) 

               gradIsoV(:) = 0.0_rp
               !$acc loop seq
               do ii=1,porder+1
                  gradIsoV(1) = gradIsoV(1) + dlxi_ip(ii)*vl(invAtoIJK(ii,isoJ,isoK))
                  gradIsoV(2) = gradIsoV(2) + dleta_ip(ii)*vl(invAtoIJK(isoI,ii,isoK))
                  gradIsoV(3) = gradIsoV(3) + dlzeta_ip(ii)*vl(invAtoIJK(isoI,isoJ,ii))
               end do

               gradV(:) = 0.0_rp
               !$acc loop seq
               do idime=1, ndime
                  !$acc loop seq
                  do jdime=1, ndime
                     gradV(idime) = gradV(idime) + He(idime,jdime,igaus,ielem) * gradIsoV(jdime)
                  end do
               end do
               !val = val + gradV(1)
               val = gradV(1)*bnorm((igaus-1)*ndime+2)

               !auxmag = sqrt(bnorm((igaus-1)*ndime+1)**2 + bnorm((igaus-1)*ndime+2)**2 + bnorm((igaus-1)*ndime+3)**2)
               !$acc atomic update
               !Rdiff(bound(iBound,igaus),1) = Rdiff(bound(iBound,igaus),1)-mu_fluid(igaus)*auxmag*val*wgp_b(igaus)
               Rdiff(bound(iBound,igaus),1) = Rdiff(bound(iBound,igaus),1)-mu_fluid(igaus)*val*wgp_b(igaus)
               !$acc end atomic

            end do
            !val = val/real(npbou,rp)
            !!$acc loop vector 
            !do igaus = 1,npbou
            !   auxmag = sqrt(bnorm((igaus-1)*ndime+1)**2 + bnorm((igaus-1)*ndime+2)**2 + bnorm((igaus-1)*ndime+3)**2)
            !   !$acc atomic update
            !   Rdiff(bound(iBound,igaus),1) = Rdiff(bound(iBound,igaus),1)-mu_fluid(igaus)*auxmag*val*wgp_b(igaus)
            !   !$acc end atomic
            !end do
         end if
      end do
      !$acc end parallel loop

   end subroutine evalCorrectNeumann

end module mod_correct_neumann
