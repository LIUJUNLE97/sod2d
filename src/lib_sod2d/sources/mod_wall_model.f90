module mod_wall_model

   use mod_numerical_params
   use mod_veclen
   use mod_nvtx
   use mod_mpi

contains

   subroutine evalWallModel(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,atoIJK, bou_code, &
         bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,ui,tauw,Rdiff)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin),atoIJK(nnode)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: rho(npoin),ui(npoin,ndime),mu_fluid(npoin)
      real(rp),   intent(inout) :: tauw(npoin,ndime)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),   intent(inout) :: Rdiff(npoin,ndime)
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: iBound,iElem,idime,igaus,iAux
      real(rp)                :: bnorm(npbou*ndime),rhol,tmag
      real(rp)                :: aux(ndime)
      ! wall law stuff
      real(rp)                :: y,ul,nul,uistar,tvelo(ndime),uiex(ndime),auxmag,auxvn,surf
      integer(4)              :: itera
      real(rp)                :: xmuit,fdvfr,devfr,point(ndime),pointF(ndime),normalF(ndime)
      real(rp)                :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol
      real(rp)                :: ux,uy,uz,px,pz

      !$acc parallel loop gang private(bnorm,uiex,point)
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         point(1:ndime) = 0.0_rp
         uiex(1:ndime) = 0.0_rp
         iElem = point2elem(bound(iBound,npbou)) ! I use an internal face node to be sure is the correct element
#if 1
         px = 0.0_rp
         py = 0.0_rp
         pz = 0.0_rp
         ux = 0.0_rp
         uy = 0.0_rp
         uz = 0.0_rp
         !$acc loop vector reduction(+:px,py,pz,ux,uy,uz)
         do igaus = 1,nnode
            px = px + coord(connec(iElem,igaus),1)
            py = py + coord(connec(iElem,igaus),2)
            pz = pz + coord(connec(iElem,igaus),3)
            ux = ux + ui(connec(iElem,igaus),1)
            uy = uy + ui(connec(iElem,igaus),2)
            uz = uz + ui(connec(iElem,igaus),3)
         end do
         point(1) = px/real(nnode,rp)
         point(2) = py/real(nnode,rp)
         point(3) = pz/real(nnode,rp)
         uiex(1) = ux/real(nnode,rp)
         uiex(2) = uy/real(nnode,rp)
         uiex(3) = uz/real(nnode,rp)
#else
         point(1:ndime) =coord(connec(iElem,nnode),1:ndime)
         uiex(1:ndime) = ui(connec(iElem,nnode),1:ndime)
#endif         

         !$acc loop vector private(aux,pointF,normalF,tvelo)
         do igaus = 1,npbou

            pointF(1:ndime) = 0.0_rp
            normalF(1:ndime) = 0.0_rp
            rhol = rho(bound(iBound,igaus))
            nul = mu_fluid(bound(iBound,igaus))/rhol
            !acc loop seq
            do idime=1,ndime
               pointF(idime) = coord(bound(iBound,igaus),idime)
               normalF(idime) = normalsAtNodes(bound(iBound,igaus),idime)
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
               Rdiff(bound(iBound,igaus),idime) = Rdiff(bound(iBound,igaus),idime)+auxmag*wgp_b(igaus)*tmag*tvelo(idime)
               !$acc end atomic
               !$acc atomic write
               tauw(bound(iBound,igaus),idime) = tmag*tvelo(idime)
               !$acc end atomic
            end do
         end do
         !!$acc loop vector 
         ! do igaus = 1,nnode
         !    !$acc loop seq
         !    do idime = 1,ndime
         !       !$acc atomic update
         !       Rdiff(connec(iElem,igaus),idime) = Rdiff(connec(iElem,igaus),idime)+gpvol(1,igaus,iElem)*tmag*tvelo(idime)*(surf/vol)
         !       !$acc end atomic
         !    end do
         ! end do
      end do
      !$acc end parallel loop

   end subroutine evalWallModel

end module mod_wall_model
