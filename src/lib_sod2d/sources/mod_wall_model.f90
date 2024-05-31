module mod_wall_model

   use mod_numerical_params
   
   use mod_nvtx
   use mod_mpi
   use mod_mpi_mesh,only:mesh_a2ij
   use elem_qua

contains

#if 1  

   subroutine evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
         bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,ui,tauw,Rdiff,fact)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: rho(npoin),ui(npoin,ndime),mu_fluid(npoin)
      real(rp),   intent(inout) :: tauw(npoin,ndime)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),   intent(inout) :: Rdiff(npoin,ndime)
      real(rp), optional, intent(in)  :: fact
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: iBound,iElem,idime,igaus,iAux
      real(rp)                :: bnorm(npbou*ndime),rhol,tmag
      real(rp)                :: aux(ndime)
      ! wall law stuff
      real(rp)                :: y,ul,nul,uistar,tvelo(ndime),uiex(ndime),auxmag,auxvn,surf
      integer(4)              :: itera,isoI,isoJ,isoK,jgaus,type_ijk,ii,isoII,isoJJ,isoKK
      real(rp)                :: xmuit,fdvfr,devfr,point(ndime),pointF(ndime),normalF(ndime)
      real(rp)                :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol
      real(rp)                :: ux,uy,uz,px,pz
      integer(4)              :: atoIJ(npbou)
      real(rp)  :: aux_fact = 1.0_rp

      atoIJ(:) = mesh_a2ij(:)

      if(present(fact)) then
         aux_fact = fact
      end if

      !$acc parallel loop gang private(bnorm,uiex,point)
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         point(1:ndime) = 0.0_rp
         uiex(1:ndime) = 0.0_rp
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector private(aux,pointF,normalF,tvelo)
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,atoIJ(npbou))),1)

            isoI = gmshAtoI(jgaus) 
            isoJ = gmshAtoJ(jgaus) 
            isoK = gmshAtoK(jgaus)
            type_ijk = 0
            if((isoI .eq. 1) .or. (isoI .eq. 2)) then
               type_ijk = 1
               if(isoI .eq. 1) then
                  isoI = 2
               else
                  isoI = 1
               end if
            else if ((isoJ .eq. 1) .or. (isoJ .eq. 2)) then
               type_ijk = 2
               if(isoJ .eq. 1) then
                  isoJ = 2
               else
                  isoJ = 1
               end if            
            else
               type_ijk = 3
               if(isoK .eq. 1) then
                  isoK = 2
               else
                  isoK = 1
               end if            
            end if      


            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,igaus)),1)

            isoII = gmshAtoI(jgaus) 
            isoJJ = gmshAtoJ(jgaus) 
            isoKK = gmshAtoK(jgaus)

            if(type_ijk .eq. 1) then
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
            else if(type_ijk .eq. 2) then
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
            else 
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
            end if          

            pointF(1:ndime) = 0.0_rp
            normalF(1:ndime) = 0.0_rp
            rhol = rho(bound(iBound,igaus))
            nul = mu_fluid(bound(iBound,igaus))/rhol
            !$acc loop seq
            do idime=1,ndime
               pointF(idime) = coord(bound(iBound,igaus),idime)
               normalF(idime) = normalsAtNodes(bound(iBound,igaus),idime)
            end do

            y = abs(normalF(1)*(point(1)-pointF(1)) + normalF(2)*(point(2)-pointF(2)) + normalF(3)*(point(3)-pointF(3)))

            auxvn = (normalF(1)*uiex(1) + normalF(2)*uiex(2) + normalF(3)*uiex(3))
            !$acc loop seq
            do idime = 1,ndime     
               tvelo(idime) = uiex(idime) - auxvn*normalF(idime)
            end do

            ul = sqrt(tvelo(1)*tvelo(1) + tvelo(2)*tvelo(2) +  tvelo(3)*tvelo(3))

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
            auxmag = sqrt(aux(1)*aux(1) + aux(2)*aux(2)  +  aux(3)*aux(3))
            !$acc loop seq
            do idime = 1,ndime
               !$acc atomic update
               Rdiff(bound(iBound,igaus),idime) = Rdiff(bound(iBound,igaus),idime)+aux_fact*auxmag*wgp_b(igaus)*tmag*tvelo(idime)
               !$acc end atomic
               !$acc atomic write
               tauw(bound(iBound,igaus),idime) = tmag*tvelo(idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop

   end subroutine evalWallModelReichardt

#else

   subroutine evalWallModelReichardt(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
         bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,ui,tauw,Rdiff)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: rho(npoin),ui(npoin,ndime),mu_fluid(npoin)
      real(rp),   intent(inout) :: tauw(npoin,ndime)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),   intent(inout) :: Rdiff(npoin,ndime)
      !
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: iBound,iElem,idime,igaus,iAux
      real(rp)                :: bnorm(npbou*ndime),rhol,tmag
      real(rp)                :: aux_1,aux_2,aux_3
      ! wall law stuff
      real(rp)                :: y,ul,nul,uistar,auxmag,auxvn,surf
      real(rp)                :: tvelo_1,tvelo_2,tvelo_3,uiex_1,uiex_2,uiex_3
      integer(4)              :: itera,isoI,isoJ,isoK,jgaus,type_ijk,ii,isoII,isoJJ,isoKK
      integer(4)              :: isoA, isoB, isoC
      real(rp)                :: xmuit,fdvfr,devfr
      real(rp)                :: point_1,point_2,point_3,pointF_1,pointF_2,pointF_3,normalF_1,normalF_2,normalF_3
      real(rp)                :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol
      real(rp)                :: ux,uy,uz,px,pz
      !integer(4)              :: atoIJ(npbou)

      ! atoIJ(:) = mesh_a2ij(:)

      !$acc parallel loop gang private(bnorm)
      do iAux = 1,numBoundsWM

         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)

         !point(1:ndime) = 0.0_rp
         !uiex(1:ndime) = 0.0_rp
         ! iElem = point2elem(bound(iBound,atoIJ(npbou)))
         iElem = point2elem(bound(iBound,mesh_a2ij(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector
         do igaus = 1,npbou
            
            !jgaus = minloc(abs(connec(iElem,:)-bound(iBound,atoIJ(npbou))),1)
            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,mesh_a2ij(npbou))),1)

            isoI = gmshAtoI(jgaus) 
            isoJ = gmshAtoJ(jgaus) 
            isoK = gmshAtoK(jgaus)
            ! type_ijk = 0
            if (isoI<3) then              ! ((isoI .eq. 1) .or. (isoI .eq. 2))
               type_ijk = 1
               isoI = MODULO((isoI+1),2)  ! Flip isoI value
            else if (isoJ<3) then         ! ((isoJ .eq. 1) .or. (isoJ .eq. 2))
               type_ijk = 2
               isoJ = MODULO((isoJ+1),2)  ! Flip isoJ value
            else                          ! ((isoK .eq. 1) .or. (isoK .eq. 2)) ???
               type_ijk = 3
               isoK = MODULO((isoK+1),2)  ! Flip isoK value
            end if  

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,igaus)),1)

            !isoII = gmshAtoI(jgaus) 
            !isoJJ = gmshAtoJ(jgaus) 
            !isoKK = gmshAtoK(jgaus)

            if(type_ijk .eq. 1) then
               !point(1:ndime) =coord(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               !uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               isoA = isoI
               isoB = gmshAtoJ(jgaus) 
               isoC = gmshAtoK(jgaus)
            else if(type_ijk .eq. 2) then
               !point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               !uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               isoA = gmshAtoI(jgaus)
               isoB = isoJ
               isoC = gmshAtoK(jgaus)
            else 
               !point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
               !uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
               isoA = gmshAtoI(jgaus)
               isoB = gmshAtoJ(jgaus) 
               isoC = isoK
            end if      

            point_1 =coord(connec(iElem,invAtoIJK(isoA,isoB,isoC)),1)
            point_2 =coord(connec(iElem,invAtoIJK(isoA,isoB,isoC)),2)
            point_3 =coord(connec(iElem,invAtoIJK(isoA,isoB,isoC)),3)

            uiex_1 = ui(connec(iElem,invAtoIJK(isoA,isoB,isoC)),1)
            uiex_2 = ui(connec(iElem,invAtoIJK(isoA,isoB,isoC)),2)
            uiex_3 = ui(connec(iElem,invAtoIJK(isoA,isoB,isoC)),3)

            !pointF(1:ndime) = 0.0_rp
            !normalF(1:ndime) = 0.0_rp
            rhol = rho(bound(iBound,igaus))
            nul = mu_fluid(bound(iBound,igaus))/rhol

            !do idime=1,ndime
            !   pointF(idime) = coord(bound(iBound,igaus),idime)
            !   normalF(idime) = normalsAtNodes(bound(iBound,igaus),idime)
            !end do

            pointF_1 = coord(bound(iBound,igaus),1)
            normalF_1 = normalsAtNodes(bound(iBound,igaus),1)

            pointF_2 = coord(bound(iBound,igaus),2)
            normalF_2 = normalsAtNodes(bound(iBound,igaus),2)

            pointF_3 = coord(bound(iBound,igaus),3)
            normalF_3 = normalsAtNodes(bound(iBound,igaus),3)

            y = abs( normalF_1*(point_1-pointF_1) + normalF_2*(point_2-pointF_2) + normalF_3*(point_3-pointF_3) )

            auxvn = (normalF_1*uiex_1) + (normalF_2*uiex_2) + (normalF_3*uiex_3)

            !do idime = 1,ndime     
            !   tvelo(idime) = uiex(idime) - auxvn*normalF(idime)
            !end do
            tvelo_1 = uiex_1 - auxvn*normalF_1
            tvelo_2 = uiex_2 - auxvn*normalF_2
            tvelo_3 = uiex_3 - auxvn*normalF_3

            ul = sqrt((tvelo_1*tvelo_1) + (tvelo_2*tvelo_2) +  (tvelo_3*tvelo_3))

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

            aux_1 = bnorm((igaus-1)*ndime+1)
            aux_2 = bnorm((igaus-1)*ndime+2)
            aux_3 = bnorm((igaus-1)*ndime+3)
            auxmag = sqrt((aux_1*aux_1) + (aux_2*aux_2)  +  (aux_3*aux_3))

            !do idime = 1,ndime
            !   !$acc atomic update
            !   Rdiff(bound(iBound,igaus),idime) = Rdiff(bound(iBound,igaus),idime)+auxmag*wgp_b(igaus)*tmag*tvelo(idime)
            !   !$acc end atomic
            !   !$acc atomic write
            !   tauw(bound(iBound,igaus),idime) = tmag*tvelo(idime)
            !   !$acc end atomic
            !end do

            ! (1)
            !$acc atomic update
            Rdiff(bound(iBound,igaus),1) = Rdiff(bound(iBound,igaus),1)+auxmag*wgp_b(igaus)*tmag*tvelo_1
            !$acc end atomic
            !$acc atomic write
            tauw(bound(iBound,igaus),1) = tmag*tvelo_1
            !$acc end atomic

            ! (2)
            !$acc atomic update
            Rdiff(bound(iBound,igaus),2) = Rdiff(bound(iBound,igaus),2)+auxmag*wgp_b(igaus)*tmag*tvelo_2
            !$acc end atomic
            !$acc atomic write
            tauw(bound(iBound,igaus),2) = tmag*tvelo_2
            !$acc end atomic

            ! (3)
            !$acc atomic update
            Rdiff(bound(iBound,igaus),3) = Rdiff(bound(iBound,igaus),3)+auxmag*wgp_b(igaus)*tmag*tvelo_3
            !$acc end atomic
            !$acc atomic write
            tauw(bound(iBound,igaus),3) = tmag*tvelo_3
            !$acc end atomic

         end do
      end do
      !$acc end parallel loop

   end subroutine evalWallModelReichardt

#endif

   subroutine evalWallModelABL(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
         bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,ui,zo,tauw,Rdiff,fact)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: rho(npoin),ui(npoin,ndime),mu_fluid(npoin),zo(npoin)
      real(rp),   intent(inout) :: tauw(npoin,ndime)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),   intent(inout) :: Rdiff(npoin,ndime)
      real(rp), optional, intent(in)  :: fact
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: iBound,iElem,idime,igaus,iAux
      real(rp)                :: bnorm(npbou*ndime),rhol,tmag
      real(rp)                :: aux(ndime)
      ! wall law stuff
      real(rp)                :: y,ul,nul,uistar,tvelo(ndime),uiex(ndime),auxmag,auxvn,surf
      integer(4)              :: isoI,isoJ,isoK,jgaus,type_ijk,ii,isoII,isoJJ,isoKK
      real(rp)                :: point(ndime),pointF(ndime),normalF(ndime)
      integer(4)              :: atoIJ(npbou)
      real(rp)  :: aux_fact = 1.0_rp
      

      atoIJ(:) = mesh_a2ij(:)

      if(present(fact)) then
         aux_fact = fact
      end if

      !$acc parallel loop gang private(bnorm,uiex,point)
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         point(1:ndime) = 0.0_rp
         uiex(1:ndime) = 0.0_rp
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element   

         !$acc loop vector private(aux,pointF,normalF,tvelo)
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,atoIJ(npbou))),1)

            isoI = gmshAtoI(jgaus) 
            isoJ = gmshAtoJ(jgaus) 
            isoK = gmshAtoK(jgaus)
            type_ijk = 0
            if((isoI .eq. 1) .or. (isoI .eq. 2)) then
               type_ijk = 1
               if(isoI .eq. 1) then
                  isoI = 2
               else
                  isoI = 1
               end if
            else if ((isoJ .eq. 1) .or. (isoJ .eq. 2)) then
               type_ijk = 2
               if(isoJ .eq. 1) then
                  isoJ = 2
               else
                  isoJ = 1
               end if            
            else
               type_ijk = 3
               if(isoK .eq. 1) then
                  isoK = 2
               else
                  isoK = 1
               end if            
            end if   

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,igaus)),1)

            isoII = gmshAtoI(jgaus) 
            isoJJ = gmshAtoJ(jgaus) 
            isoKK = gmshAtoK(jgaus)

            if(type_ijk .eq. 1) then
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
            else if(type_ijk .eq. 2) then
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
            else 
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
            end if          

            pointF(1:ndime) = 0.0_rp
            normalF(1:ndime) = 0.0_rp
            rhol = rho(bound(iBound,igaus))
            !$acc loop seq
            do idime=1,ndime
               pointF(idime) = coord(bound(iBound,igaus),idime)
               normalF(idime) = normalsAtNodes(bound(iBound,igaus),idime)
            end do

            y = abs(normalF(1)*(point(1)-pointF(1)) + normalF(2)*(point(2)-pointF(2)) + normalF(3)*(point(3)-pointF(3)))

            auxvn = (normalF(1)*uiex(1) + normalF(2)*uiex(2) + normalF(3)*uiex(3))
            !$acc loop seq
            do idime = 1,ndime     
               tvelo(idime) = uiex(idime) - auxvn*normalF(idime)
            end do

            ul = sqrt(tvelo(1)*tvelo(1) + tvelo(2)*tvelo(2) +  tvelo(3)*tvelo(3))
            uistar = ul*0.41_rp/log(1.0_rp+y/zo(bound(iBound,igaus)))
            tmag = rhol*uistar*uistar/ul

            aux(1) = bnorm((igaus-1)*ndime+1)
            aux(2) = bnorm((igaus-1)*ndime+2)
            aux(3) = bnorm((igaus-1)*ndime+3)
            auxmag = sqrt(aux(1)*aux(1) + aux(2)*aux(2)  +  aux(3)*aux(3))
            !$acc loop seq
            do idime = 1,ndime
               !$acc atomic update
               Rdiff(bound(iBound,igaus),idime) = Rdiff(bound(iBound,igaus),idime)+aux_fact*auxmag*wgp_b(igaus)*tmag*tvelo(idime)
               !$acc end atomic
               !$acc atomic write
               tauw(bound(iBound,igaus),idime) = tmag*tvelo(idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop

   end subroutine evalWallModelABL

end module mod_wall_model
