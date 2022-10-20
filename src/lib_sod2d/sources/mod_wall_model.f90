module mod_wall_model

   use mod_constants
   use mod_veclen
   use mod_nvtx
   use mod_mpi

contains

   subroutine evalWallModel(bouCodes2WallModel,nelem,npoin,nboun,numBoundCodes,connec,bound,point2elem,atoIJK, bou_code, &
         bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,ui,Rdiff)

      implicit none

      integer(4), intent(in)  :: bouCodes2WallModel(numBoundCodes),npoin, nboun,numBoundCodes, bound(nboun,npbou), bou_code(nboun)
      integer(4), intent(in)  :: nelem, connec(nelem,nnode),point2elem(npoin),atoIJK(nnode)
      real(rp),    intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),    intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)  :: rho(npoin), ui(npoin,ndime),mu_fluid(npoin)
      real(rp),    intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),    intent(inout) :: Rdiff(npoin,ndime)
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: ibound, idime, igaus, ipbou, ielem,jgaus,kgaus,jjgaus
      integer(4)              :: numBelem, counter, ii, jdime, kdime,iFace
      integer(4)              :: isoI, isoJ, isoK
      real(rp)                 :: bnorm(npbou*ndime), nmag,  rhol,tmag
      real(rp)                :: sig,aux(ndime),aux2(ndime)
      ! wall law stuff
      real(rp)  ::  y, ul, nul,uistar,tvelo(ndime),uiex(ndime),auxmag,auxvn,surf
      integer(4)              :: itera,ldime,icode
      real(rp)                 :: xmuit,fdvfr,devfr,point(ndime),pointF(ndime),normalF(ndime)
      real(rp)                 :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                 :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                 :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol

      !$acc parallel loop gang private(bnorm,uiex)
      do ibound = 1, nboun
         icode=bou_code(ibound)
         if (bouCodes2WallModel(icode) == 1) then
            bnorm(1:npbou*ndime) = bounorm(ibound,1:npbou*ndime)
            ielem = point2elem(bound(ibound,npbou)) ! I use an internal face node to be sure is the correct element
            
            point(1:ndime) = 0.0_rp
            uiex(1:ndime) = 0.0_rp
            !$acc loop vector reduction(+:point,uiex)
            do igaus = 1,nnode
               !acc loop seq
               do idime=1,ndime
                 point(idime) = point(idime) + coord(connec(ielem,igaus),idime)
                 uiex(idime) = uiex(idime) + ui(connec(ielem,igaus),idime)
               end do
            end do
            point(1) = point(1)/real(nnode,rp)
            point(2) = point(2)/real(nnode,rp)
            point(3) = point(3)/real(nnode,rp)
            uiex(1) = uiex(1)/real(nnode,rp)
            uiex(2) = uiex(2)/real(nnode,rp)
            uiex(3) = uiex(3)/real(nnode,rp)

            !$acc loop vector private(aux)
            do igaus = 1,npbou

               pointF(1:ndime) = 0.0_rp
               normalF(1:ndime) = 0.0_rp
               rhol = rho(bound(ibound,igaus))
               nul = mu_fluid(bound(ibound,igaus))/rhol
               !acc loop seq
               do idime=1,ndime
                  pointF(idime) = coord(bound(ibound,igaus),idime)
                  normalF(idime) = normalsAtNodes(bound(ibound,igaus),idime)
               end do

               y = abs(dot_product(normalF,point(:)-pointF(:)))

               auxvn = dot_product(normalF,uiex)
               !$acc loop seq
               do idime = 1,ndime     
                  tvelo(idime) = uiex(idime) - auxvn*normalF(idime)
               end do

               ul = sqrt(dot_product(tvelo(:),tvelo(:)))

               if( y > 0.0_rp .and. ul > 1.0e-10 ) then            
                  uistar = sqrt( ul * nul / y )
                  if( uistar * y / nul > 5.0_rp ) then
                     vkinv = 1.0_rp / 0.41_rp
                     onovu = 1.0_rp / ul
                     xmuit = y / nul
                     itera = 0
                     parco = 1.0_rp
                     oneoe = 1.0_rp/11.0_rp
                     do while( parco >= 1.0e-6_rp .and. itera < 100 )
                        itera = itera + 1
                        uistar = max(uistar,0.0_rp)
                        yplus = uistar * xmuit
                        ypele = yplus * oneoe
                        ypel2 = min(ypele,20.0_rp)
                        expye = exp(-ypel2)
                        yplu2 = min(yplus,70.0_rp)
                        expyt = exp(-yplu2 * 0.33_rp) 
                        firsl = vkinv*log(1.0_rp + 0.4_rp*yplus)
                        fdvfr = uistar*(firsl+7.8_rp*(1.0_rp-expye-ypele*expyt))-ul
                        diffd = firsl + vkinv*0.4_rp*yplus/(1.0_rp+0.4_rp*yplus)&
                           & + 7.8_rp*(1.0_rp-expye*(1.0_rp-ypele)&
                           &  - ypele*expyt*(2.0_rp-yplus*0.33_rp))               
                        devfr = -fdvfr/diffd
                        parco = abs(devfr*onovu)

                        uistar= uistar + devfr
                     end do
                  end if
               else
                  uistar = 0.0_rp
               end if

               if( y*uistar/nul < 5.0_rp ) then
                  tmag  = rhol*nul/y 
               else
                  tmag = rhol*uistar*uistar/ul
               end if

               aux(1) = bnorm((igaus-1)*ndime+1)
               aux(2) = bnorm((igaus-1)*ndime+2)
               aux(3) = bnorm((igaus-1)*ndime+3)
               auxmag = sqrt(dot_product(aux(:),aux(:)))
               !$acc loop seq
               do idime = 1,ndime
                  !$acc atomic update
                  Rdiff(bound(ibound,igaus),idime) = Rdiff(bound(ibound,igaus),idime)+auxmag*wgp_b(igaus)*tmag*tvelo(idime)
                  !$acc end atomic
               end do
            end do
            !!$acc loop vector 
            ! do igaus = 1,nnode
            !    !$acc loop seq
            !    do idime = 1,ndime
            !       !$acc atomic update
            !       Rdiff(connec(ielem,igaus),idime) = Rdiff(connec(ielem,igaus),idime)+gpvol(1,igaus,ielem)*tmag*tvelo(idime)*(surf/vol)
            !       !$acc end atomic
            !    end do
            ! end do
         end if
      end do
      !$acc end parallel loop


   end subroutine evalWallModel

end module mod_wall_model
