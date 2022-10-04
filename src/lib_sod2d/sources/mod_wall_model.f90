module mod_wall_model

   use mod_constants
   use mod_veclen
   use mod_nvtx
   use mod_mpi

contains

   subroutine evalWallModel(surfCode,nelem,npoin,nboun,connec,bound,point2elem,atoIJK, bou_code, &
         bounorm,wgp_b,coord,gpvol,mu_fluid,rho,u,Rdiff)

      implicit none

      integer(4), intent(in)  :: surfCode,npoin, nboun, bound(nboun,npbou), bou_code(nboun)
      integer(4), intent(in)  :: nelem, connec(nelem,nnode),point2elem(npoin),atoIJK(nnode)
      real(rp),    intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou)
      real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime), mu_fluid(npoin)
      real(rp),    intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),    intent(inout) :: Rdiff(npoin,ndime)
      real(rp)                :: Ftau_l(ndime)
      integer(4)              :: ibound, idime, igaus, ipbou, ielem, jgaus,kgaus
      integer(4)              :: numBelem, counter, ii, jdime, kdime,iFace
      real(rp)                 :: bnorm(npbou*ndime), nmag,  rhol,tmag
      real(rp)                :: mu_fgp, mu_egp, mufluidl,sig,aux(ndime)
      ! wall law stuff
      real(rp)  :: y, ul, nul, k, ustar,tvelo(ndime),uex(ndime),auxmag,auxvn
      integer(4)              :: itera,ldime
      real(rp)                 :: xmuit,fdvfr,devfr
      real(rp)                 :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                 :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                 :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol,surf
      integer(4) :: atoIJ(16) ! no esta be fer ho aqui

      atoIJ = [1,4,12,11,2,3,7,8,5,10,13,16,6,9,14,15]


      !$acc parallel loop gang private(bnorm,Ftau_l,uex)
      do ibound = 1, nboun
         if (bou_code(ibound) == surfCode) then
            bnorm(1:npbou*ndime) = bounorm(ibound,1:npbou*ndime)
            Ftau_l(1:ndime) = 0.0_rp
            vol=0.0
            ielem = point2elem(bound(ibound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element
            jgaus = connec(ielem,atoIJK(nnode))              ! exchange location
            uex(1:ndime) = u(jgaus,1:ndime)
            rhol = rho(jgaus)
            nul = mu_fluid(jgaus)/rhol

            !$acc loop seq private(aux,tvelo) reduction(+:surf) 
            do igaus = 1,npbou
               kgaus = bound(ibound,igaus) ! node at the boundary

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

               y = dot_product(aux,coord(jgaus,:)-coord(kgaus,:))

               auxvn = dot_product(aux,uex)
               !$acc loop seq
               do idime = 1,ndime     
                  tvelo(idime) = uex(idime) - auxvn*aux(idime)
               end do

               ul = sqrt(dot_product(tvelo(:),tvelo(:)))

               if( y > 0.0_rp .and. ul /= 0.0_rp ) then            
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
               end if
               if( y*ustar/nul < 5.0_rp ) then
                  tmag  = rhol*nul/y 
               else
                  tmag = rhol*ustar*ustar/ul
               end if

               !$acc loop vector 
               do idime = 1,ndime
                  !$acc atomic update
                  Ftau_l(idime) = Ftau_l(idime)-auxmag*wgp_b(igaus)*tmag*tvelo(idime)
                  !$acc end atomic
               end do
               surf = surf + auxmag*wgp_b(igaus)
            end do
            vol = 0
            !!$acc loop vector reduction(+:vol)
            !do igaus = 1,ngaus
            !   vol = vol + gpvol(1,igaus,ielem)
            !end do
            !$acc loop vector collapse(2)
            do igaus = 1,ngaus
               do idime = 1,ndime
                  !$acc atomic update
                  Rdiff(connec(ielem,igaus),idime) = Rdiff(connec(ielem,igaus),idime)-Ftau_l(idime)*gpvol(1,igaus,ielem)/surf
                  !$acc end atomic
               end do
            end do
         end if
      end do
      !$acc end parallel loop

   end subroutine evalWallModel

end module mod_wall_model
