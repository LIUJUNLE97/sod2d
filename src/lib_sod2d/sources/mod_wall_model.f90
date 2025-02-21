module mod_wall_model

   use mod_numerical_params
   
   use mod_nvtx
   use mod_mpi
   use mod_mpi_mesh,only:mesh_a2ij
   use mod_comms_boundaries
   use elem_qua

   implicit none

   integer(rp), allocatable, dimension(:,:) :: iex, extype


contains

   subroutine evalEXAtFace(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
      bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
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

      allocate(iex(nelem,nnode),extype(nelem,nnode))
      !$acc enter data create(iex(:,:))
      !$acc enter data create(extype(:,:))

      !$acc parallel loop gang 
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector
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
               iex(ielem,igaus) = isoI
            else if ((isoJ .eq. 1) .or. (isoJ .eq. 2)) then
               type_ijk = 2
               if(isoJ .eq. 1) then
                  isoJ = 2
               else
                  isoJ = 1
               end if            
               iex(ielem,igaus) = isoJ
            else
               type_ijk = 3
               if(isoK .eq. 1) then
                  isoK = 2
               else
                  isoK = 1
               end if 
               iex(ielem,igaus) = isoK           
            end if     
            extype(ielem,igaus) = type_ijk
         end do
      end do
      !$acc end parallel loop

   end subroutine evalEXAtFace

   subroutine evalEXAt1OffNode(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
      bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
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

      allocate(iex(nelem,nnode),extype(nelem,nnode))
      !$acc enter data create(iex(:,:))
      !$acc enter data create(extype(:,:))

      !$acc parallel loop gang 
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,atoIJ(npbou))),1)

            isoI = gmshAtoI(jgaus) 
            isoJ = gmshAtoJ(jgaus) 
            isoK = gmshAtoK(jgaus)
            type_ijk = 0
            if((isoI .eq. 1) .or. (isoI .eq. 2)) then
               type_ijk = 1
               if(isoI .eq. 1) then
                  isoI = 3
               else
                  isoI = porder+1
               end if
               iex(ielem,igaus) = isoI
            else if ((isoJ .eq. 1) .or. (isoJ .eq. 2)) then
               type_ijk = 2
               if(isoJ .eq. 1) then
                  isoJ = 3
               else
                  isoJ = porder+1
               end if            
               iex(ielem,igaus) = isoJ
            else
               type_ijk = 3
               if(isoK .eq. 1) then
                  isoK = 3
               else
                  isoK = porder+1
               end if 
               iex(ielem,igaus) = isoK           
            end if     
            extype(ielem,igaus) = type_ijk
         end do
      end do
      !$acc end parallel loop

   end subroutine evalEXAt1OffNode

   subroutine evalEXAtHWM(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
      bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime), tau(ndime,ndime), divU
      integer(4)              :: iBound,iElem,idime,igaus,iAux
      real(rp)                :: bnorm(npbou*ndime),rhol,tmag
      real(rp)                :: aux(ndime)
      integer(4)              :: i, inode
      real(rp)                :: dist, auxVal, sig
      ! wall law stuff
      real(rp)                :: y,ul,nul,uistar,tvelo(ndime),uiex(ndime),auxmag,auxvn,surf
      integer(4)              :: itera,isoI,isoJ,isoK,jgaus,type_ijk,ii,isoII,isoJJ,isoKK
      integer(4)              :: nearestPointWM,isoNearestI,isoNearestJ,isoNearestK
      real(rp)                :: xmuit,fdvfr,devfr,point(ndime),pointF(ndime),normalF(ndime)
      real(rp)                :: normalBouWall(ndime),pointHWM(ndime)
      real(rp)                :: vkinv,diffd,parco,yplus,onovu,yplu2
      real(rp)                :: ypele,expye,expyt,oneoe,firsl,ypel2
      real(rp)                :: pplus,densi,gradp,grpr2,py,sq,inv,ln4,uplus,vol
      real(rp)                :: ux,uy,uz,px,pz
      integer(4)              :: atoIJ(npbou)
      real(rp)  :: aux_fact = 1.0_rp

      atoIJ(:) = mesh_a2ij(:)

      allocate(iex(nelem,nnode),extype(nelem,nnode))
      !$acc enter data create(iex(:,:))
      !$acc enter data create(extype(:,:))
   

      !$acc parallel loop gang 
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)         
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector 
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,atoIJ(npbou))),1)

            isoI = gmshAtoI(jgaus) 
            isoJ = gmshAtoJ(jgaus) 
            isoK = gmshAtoK(jgaus)
            type_ijk = 0
            
            if((isoI .eq. 1) .or. (isoI .eq. 2)) then
               type_ijk = 1
               auxVal = 1000000_rp
               !$acc loop seq
               do i = 1, porder+1 
                  dist = abs(sum((coord(connec(iElem, invAtoIJK(i,isoJ,isoK)), 1:ndime) - coord(connec(iElem, invAtoIJK(isoI,isoJ,isoK)), 1:ndime) )**2))
                  dist = abs(dist-wmles_walex)
                  if (dist .lt. auxVal) then
                     auxVal = dist
                     isoI = i
                  endif
               end do        
               iex(ielem,igaus) = isoI                 
            else if ((isoJ .eq. 1) .or. (isoJ .eq. 2)) then
               type_ijk = 2
               auxVal = 1000000_rp
               !$acc loop seq
               do i = 1, porder+1 
                  dist = sqrt(sum((coord(connec(iElem, invAtoIJK(isoI,i,isoK)), 1:ndime) - coord(connec(iElem, invAtoIJK(isoI,isoJ,isoK)), 1:ndime) )**2))
                  dist = abs(dist-wmles_walex)
                  if (dist .lt. auxVal) then
                     auxVal = dist
                     isoJ = i
                  endif
               end do             
               iex(ielem,igaus) = isoJ           
            else
               type_ijk = 3
               auxVal = 1000000_rp
               !$acc loop seq
               do i = 1, porder+1 
                  dist = sqrt(sum((coord(connec(iElem, invAtoIJK(isoI,isoJ,i)), 1:ndime) - coord(connec(iElem, invAtoIJK(isoI,isoJ,isoK)), 1:ndime) )**2))
                  dist = abs(dist-wmles_walex)
                  if (dist .lt. auxVal) then
                     auxVal = dist
                     isoK = i
                  endif
               end do  
               iex(ielem,igaus) = isoK                        
            end if  

            extype(ielem,igaus) = type_ijk
         end do
      end do
      !$acc end parallel loop
    
   end subroutine evalEXAtHWM   

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

      !$acc parallel loop gang private(bnorm)
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector private(aux,pointF,normalF,tvelo,uiex,point)
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,igaus)),1)

            isoII = gmshAtoI(jgaus) 
            isoJJ = gmshAtoJ(jgaus) 
            isoKK = gmshAtoK(jgaus)
            type_ijk = extype(ielem,igaus)
            if(type_ijk .eq. 1) then
               isoI = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
            else if(type_ijk .eq. 2) then
               isoJ = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
            else 
               isoK = iex(ielem,igaus)
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

      !$acc parallel loop gang private(bnorm)
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         point(1:ndime) = 0.0_rp
         uiex(1:ndime) = 0.0_rp
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element   

         !$acc loop vector private(aux,pointF,normalF,tvelo,uiex,point)
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,igaus)),1)

            isoII = gmshAtoI(jgaus) 
            isoJJ = gmshAtoJ(jgaus) 
            isoKK = gmshAtoK(jgaus)
            type_ijk = extype(ielem,igaus)

            if(type_ijk .eq. 1) then
               isoI = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
            else if(type_ijk .eq. 2) then
               isoJ = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
            else 
               isoK = iex(ielem,igaus)
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

   subroutine evalWallModelThinBLFit(numBoundsWM,listBoundsWM,nelem,npoin,nboun,connec,bound,point2elem,bou_code, &
      bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,ui,gradP,tauw,wmles_thinBL_fit_d,Rdiff,fact)

      implicit none

      integer(4), intent(in)  :: numBoundsWM,listBoundsWM(numBoundsWM)
      integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun)
      integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
      real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
      integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
      real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
      real(rp),   intent(in)  :: rho(npoin),ui(npoin,ndime),mu_fluid(npoin),gradP(npoin,ndime)
      real(rp),   intent(inout) :: tauw(npoin,ndime)
      real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
      real(rp),   intent(inout) :: Rdiff(npoin,ndime)
      integer(4), intent(inout) :: wmles_thinBL_fit_d(npoin)
      real(rp), optional, intent(in)  :: fact
      real(rp)                :: gradPex(ndime)
      integer(4)              :: iBound,iElem,idime,igaus,iAux
      real(rp)                :: bnorm(npbou*ndime),rhol,tmag,gradPmag
      real(rp)                :: aux(ndime)
      ! wall law stuff
      integer(4)              :: itera,isoI,isoJ,isoK,jgaus,type_ijk,ii,isoII,isoJJ,isoKK
      real(rp)                :: Re_ex, Re_out, N, phiP, gradPproj, Chi, Re_min,p_aux, Re_fit
      real(rp)                :: y,ul,nul,uistar,tvelo(ndime),uiex(ndime),auxmag,auxvn,surf,auxpn,tgradP(ndime)
      real(rp)                :: point(ndime),pointF(ndime),normalF(ndime)
      real(rp)                :: Beta1, Beta2, Alpha, Gamma, Theta, Sigma, Mu, stream(ndime)
      integer(4)              :: atoIJ(npbou)
      real(rp)  :: aux_fact = 1.0_rp

      atoIJ(:) = mesh_a2ij(:)
      !stream(1) = 1.0_rp
      !stream(2) = 0.0_rp
      !stream(3) = 0.0_rp

      if(present(fact)) then
         aux_fact = fact
      end if

      !$acc kernels
      wmles_thinBL_fit_d(1:npoin) = 0
      !$acc end kernels

      !$acc parallel loop gang private(bnorm)
      do iAux = 1,numBoundsWM
         iBound = listBoundsWM(iAux)
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         iElem = point2elem(bound(iBound,atoIJ(npbou))) ! I use an internal face node to be sure is the correct element

         !$acc loop vector private(aux,pointF,normalF,tvelo,uiex,point,gradPex,tgradP)
         do igaus = 1,npbou

            jgaus = minloc(abs(connec(iElem,:)-bound(iBound,igaus)),1)

            isoII = gmshAtoI(jgaus) 
            isoJJ = gmshAtoJ(jgaus) 
            isoKK = gmshAtoK(jgaus)
            type_ijk = extype(ielem,igaus)
            if(type_ijk .eq. 1) then
               isoI = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
               gradPex(1:ndime) = gradP(connec(iElem,invAtoIJK(isoI,isoJJ,isoKK)),1:ndime)
            else if(type_ijk .eq. 2) then
               isoJ = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               uiex(1:ndime) = ui(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
               gradPex(1:ndime) = gradP(connec(iElem,invAtoIJK(isoII,isoJ,isoKK)),1:ndime)
            else 
               isoK = iex(ielem,igaus)
               point(1:ndime) =coord(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
               gradPex(1:ndime) = gradP(connec(iElem,invAtoIJK(isoII,isoJJ,isoK)),1:ndime)
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
            auxpn = (normalF(1)*gradPex(1) + normalF(2)*gradPex(2) + normalF(3)*gradPex(3))
            !$acc loop seq
            do idime = 1,ndime     
               tvelo(idime)  = uiex(idime)    - auxvn*normalF(idime)
               tgradP(idime) = gradPex(idime) - auxpn*normalF(idime)
            end do

            ! Inputs at the exchange location
            ul = sqrt(tvelo(1)*tvelo(1) + tvelo(2)*tvelo(2) +  tvelo(3)*tvelo(3))
            tvelo(:) = tvelo(:)/ul
            Re_ex = ul*y/nul
            N = dot_product(tgradP,tvelo)/rhol
            phiP = N*(y**3)/(nul**2) 

            ! Equilibrium contribution  
            Beta1 = 1.0_rp/(1.0_rp+(0.155_rp/(Re_ex**0.03_rp)))
            Beta2 = 1.7_rp-(1.0_rp/(1.0_rp+(36.0_rp/(Re_ex**0.75_rp))))
            Re_fit = 0.005_rp**(Beta1-0.5_rp)*(Re_ex**Beta1)*((1.0_rp+(1.0_rp/((0.005_rp*Re_ex)**Beta2)))**((Beta1-0.5_rp)/Beta2))                        
            Chi=N*y*((Re_ex/(ul*Re_fit))**2)
            if(abs(Chi) .lt. 0.2_rp) then 
               ! Equilibrium contribution         
               Beta2 = 1.7_rp-(1.0_rp/(1.0_rp+(36.0_rp/(Re_ex**0.65_rp))))
               Re_fit = 0.005_rp**(Beta1-0.5_rp)*(Re_ex**Beta1)*((1.0_rp+(1.0_rp/((0.005_rp*Re_ex)**Beta2)))**((Beta1-0.5_rp)/Beta2))
               ! Mild pressure gradients      
               Alpha = 0.0296_rp+(0.15_rp*tanh(Chi-0.2_rp))
               Mu = 2.25_rp-(0.4_rp*tanh(0.9_rp*Chi))
               Sigma = 0.5_rp+0.1*tanh(20_rp*Chi)
               Gamma = Alpha/exp(0.5_rp*(((log(Re_ex)-Mu)/Sigma)**2.0_rp))
               Theta = 1.0_rp/(1.0_rp+(0.0025_rp*Re_ex))
               Re_out = Re_fit*((Theta/((1.0_rp+(0.5_rp*Chi))**0.5_rp))+1.0_rp-Theta+Gamma)
            else 
               ! Strong pressure gradients  
               if(phiP .lt. 0.0_rp) then
                  Re_min = 1.5_rp*((-phiP)**0.39_rp)*(1.0_rp/((1.0_rp+((1000.0_rp/phiP)**2.0_rp))**0.055_rp))
                  p_aux = 2.5_rp-0.6_rp*(1.0_rp+tanh((2.0_rp*log(-phiP))-6.0_rp))
                  Re_out = ((Re_min**p_aux)+(Re_fit**p_aux))**(1.0_rp/p_aux)               
               else
                  Re_min = 2.5_rp*(phiP**0.54_rp)*(1.0_rp/((1.0_rp +((30.0_rp/phiP)**0.5_rp))**0.88_rp))
                  if(Re_ex .gt. Re_min) then
                     Re_out = Re_fit*(1.0_rp-(1.0_rp/((1.0_rp+log(Re_ex/Re_min))**1.9_rp)))
                  else
                     Re_out = 0.0_rp
                     wmles_thinBL_fit_d(bound(iBound,igaus)) = 1
                  end if
               end if
            end if

            uistar = ul*(Re_out/Re_ex) 
            tmag = rhol*uistar*uistar

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

      
      !if(mpi_size.ge.2) then
      !   call mpi_halo_bnd_atomic_max_integer_iSendiRcv(wmles_thinBL_fit_d)
      !end if

end subroutine evalWallModelThinBLFit

end module mod_wall_model
