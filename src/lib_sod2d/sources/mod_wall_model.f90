module mod_wall_model

   use mod_constants
   use mod_veclen
   use mod_nvtx
   use mod_mpi

contains

   subroutine evalWallModel(surfCode,nelem,npoin,nboun,connec,bound,point2elem,atoIJK, bou_code, &
         bounorm,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,mu_wallex,rho,u,ui,Rdiff)

      implicit none

      integer(4), intent(in)  :: surfCode,npoin, nboun, bound(nboun,npbou), bou_code(nboun)
      integer(4), intent(in)  :: nelem, connec(nelem,nnode),point2elem(npoin),atoIJK(nnode)
      real(rp),    intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),    intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime),ui(npoin,ndime),mu_fluid(npoin), mu_wallex(npoin)
      real(rp),    intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),    intent(inout) :: Rdiff(npoin,ndime)
      real(rp)                :: Ftau_l(ndime)
      real(rp)                :: Ftau_real(ndime)
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: ibound, idime, igaus, ipbou, ielem,jgaus,kgaus,jjgaus
      integer(4)              :: numBelem, counter, ii, jdime, kdime,iFace
      integer(4)              :: isoI, isoJ, isoK
      real(rp)                 :: bnorm(npbou*ndime), nmag,  rhol,tmag
      real(rp)                :: sig,aux(ndime),aux2(ndime)
      ! wall law stuff
      real(rp)  :: ulocal(nnode,ndime), y, ul, nul, k, ustar,uistar,tvelo(ndime),uex(ndime),uiex(ndime),auxmag,auxvn
      integer(4)              :: itera,ldime
      real(rp)                 :: xmuit,fdvfr,devfr
      real(rp)                 :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                 :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                 :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol,surf
      integer(4) :: atoIJ(16) ! no esta be fer ho aqui

      atoIJ = [1,4,12,11,2,3,7,8,5,10,13,16,6,9,14,15]


      !$acc parallel loop gang private(bnorm,uex,uiex)
      do ibound = 1, nboun
         if (bou_code(ibound) == surfCode) then
            bnorm(1:npbou*ndime) = bounorm(ibound,1:npbou*ndime)
            Ftau_l(1:ndime) = 0.0_rp
            ielem = point2elem(bound(ibound,npbou)) ! I use an internal face node to be sure is the correct element
            jgaus = connec(ielem,nnode)              ! exchange location
            uex(1:ndime) = u(jgaus,1:ndime)
            uiex(1:ndime) = ui(jgaus,1:ndime)

            !$acc loop vector private(Ftau_l,aux)
            do igaus = 1,npbou
               kgaus = bound(ibound,igaus) ! node at the boundary
               rhol = rho(kgaus)
               nul = mu_fluid(kgaus)/rhol
               !nul = mu_wallex(kgaus)/rhol

               sig=1.0_rp
               aux(1) = bnorm((igaus-1)*ndime+1)
               aux(2) = bnorm((igaus-1)*ndime+2)
               aux(3) = bnorm((igaus-1)*ndime+3)
               auxmag = sqrt(dot_product(aux(:),aux(:)))
               if(dot_product(coord(jgaus,:)-coord(kgaus,:), aux(:)) .lt. 0.0_rp ) then
                  sig=-1.0_rp
               end if
               !$acc loop seq
               do idime = 1,ndime     
                  aux(idime) = aux(idime)*sig/auxmag
               end do

               y = abs(dot_product(aux,coord(jgaus,:)-coord(kgaus,:)))

               auxvn = dot_product(aux,uiex)
               !$acc loop seq
               do idime = 1,ndime     
                  tvelo(idime) = uiex(idime) - auxvn*aux(idime)
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

#if 1
               auxvn = dot_product(aux,uex)
               !$acc loop seq
               do idime = 1,ndime     
                  tvelo(idime) = uex(idime) - auxvn*aux(idime)
               end do

               ul = sqrt(dot_product(tvelo(:),tvelo(:)))

               if( y > 0.0_rp .and. ul > 1.0e-10 ) then            
                  ustar = sqrt( ul * nul / y )
                  if( ustar * y / nul > 5.0_rp ) then
                     vkinv = 1.0_rp / 0.41_rp
                     onovu = 1.0_rp / ul
                     xmuit = y / nul
                     itera = 0
                     parco = 1.0_rp
                     oneoe = 1.0_rp/11.0_rp
                     do while( parco >= 1.0e-6_rp .and. itera < 100 )
                        itera = itera + 1
                        ustar = max(ustar,0.0_rp)
                        yplus = ustar * xmuit
                        ypele = yplus * oneoe
                        ypel2 = min(ypele,20.0_rp)
                        expye = exp(-ypel2)
                        yplu2 = min(yplus,70.0_rp)
                        expyt = exp(-yplu2 * 0.33_rp) 
                        firsl = vkinv*log(1.0_rp + 0.4_rp*yplus)
                        fdvfr = ustar*(firsl+7.8_rp*(1.0_rp-expye-ypele*expyt))-ul
                        diffd = firsl + vkinv*0.4_rp*yplus/(1.0_rp+0.4_rp*yplus)&
                           & + 7.8_rp*(1.0_rp-expye*(1.0_rp-ypele)&
                           &  - ypele*expyt*(2.0_rp-yplus*0.33_rp))               
                        devfr = -fdvfr/diffd
                        parco = abs(devfr*onovu)

                        ustar= ustar + devfr
                     end do
                  end if
               else
                  ustar = 0.0_rp
               end if               

               if( y*ustar/nul < 5.0_rp ) then
                  tmag  = rhol*nul/y 
               else
                  tmag = rhol*uistar*ustar/ul
               end if
#else 

               if( y*ustar/nul < 5.0_rp ) then
                  tmag  = rhol*nul/y 
               else
                  tmag = rhol*uistar*uistar/ul
               end if
#endif               

               !$acc loop seq 
               do idime = 1,ndime
                  Ftau_l(idime) = -tmag*tvelo(idime)
               end do

               !$acc loop seq
               do idime = 1,ndime
                  !$acc atomic update
                  Rdiff(bound(ibound,igaus),idime) = Rdiff(bound(ibound,igaus),idime)-auxmag*wgp_b(igaus)*Ftau_l(idime)
                  !$acc end atomic
               end do
            end do
         end if
      end do
      !$acc end parallel loop


   end subroutine evalWallModel

end module mod_wall_model
