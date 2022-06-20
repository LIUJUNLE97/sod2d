module mod_geom

      use mod_constants
      use elem_qua
      use elem_hex

      contains

         subroutine char_length(ielem,nelem,npoin,connec,coord,he)

                 implicit none

                 integer(4), intent(in)  :: ielem, nelem, npoin, connec(nelem,nnode)
                 real(8),    intent(in)  :: coord(npoin,ndime)
                 real(8),    intent(out) :: he
                 integer(4)              :: iedge, ncorner, nedge
                 real(8)                 :: dist(12,ndime), dist2, aux

                 !
                 ! Compute r = x2-x1 for all element edges
                 !
                 if (ndime == 2) then
                    if (nnode == 3 .or. nnode == 6) then ! TRI_XX
                       write(*,*) "ELEMENT TRI_XX NOT CODED!"
                    else if (nnode == 4 .or. nnode == 9) then ! QUA_XX
                       call quad_edges(ielem,nelem,npoin,connec,coord,ncorner,nedge,dist(1:4,1:ndime))
                    else
                       write(*,*) "INCORRECT ELEMENT TYPE (NODES ERROR)!"
                    end if
                 else if (ndime == 3) then
                    if (nnode == 4 .or. nnode == 10) then ! TET_XX
                       write(*,*) "ELEMENT TET_XX NOT CODED!"
                    else if (nnode == 8 .or. nnode == 27 .or. nnode == 64) then ! HEX_XX
                       call hexa_edges(ielem,nelem,npoin,connec,coord,ncorner,nedge,dist(1:12,1:ndime))
                    else
                       write(*,*) "INCORRECT ELEMENT TYPE (NODES ERROR)!"
                    end if
                 else
                    write(*,*) "BY SIGMAR NO!"
                 end if

                 !
                 ! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
                 !
                 dist2 = 1000000000000000000000000000000000000000000000000000000000000000000000.0d0
                 do iedge = 1,nedge
                    aux = sqrt(dot_product(dist(iedge,:),dist(iedge,:)))
                    dist2 = min(dist2,aux)
                 end do
                 he = dist2

         end subroutine char_length

         subroutine char_length_spectral(ielem,nelem,npoin,connec,coord,Ml,he_l)

                 implicit none

                 integer(4), intent(in)    :: ielem, nelem, npoin, connec(nelem,nnode)
                 real(8),    intent(in)    :: coord(npoin,ndime), Ml(npoin)
                 real(8),    intent(inout) :: he_l(nelem,nnode)
                 integer(4)                :: inode, jnode,idime
                 real(8)                   :: aux,dist(ndime),dist2

                 !
                 ! Compute r = x2-x1 for all element nodes
                 !
                 !
                 ! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
                 !
                 !aux = 0.0d0
                 !aux = 1000000000000000000000000000000000000000000000000000000000000000000000.0d0
                 !do inode = 1,nnode
                 !   do jnode = 1,nnode
                 !      if(inode .ne. jnode) then 
                 !         dist = coord(connec(ielem,inode),:)-coord(connec(ielem,jnode),:)
                 !         dist2 = sqrt(dot_product(dist(:),dist(:)))
                 !         aux = min(dist2,aux)
                 !      end if
                 !   end do
                 !end do
                 do inode = 1,nnode
                    he_l(ielem,inode) = Ml(connec(ielem,inode))**(1.0d0/dble(ndime))
                 end do
         end subroutine char_length_spectral

         subroutine linearMeshOutput(nelem,connec,listHEX08,connecLINEAR)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Given a high order mesh, generates a linearized connectivity !
            ! table based on listHEX08 ordering for each elem. type.       !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            implicit none

            integer(4), intent(in)  :: nelem,connec(nelem,nnode), listHEX08((porder**ndime),2**ndime) ! TODO: make this more generic
            integer(4), intent(out) :: connecLinear(nelem*(porder**ndime),2**ndime)
            integer(4)              :: ind, iproxy,ielem

            ind = 0
            do ielem = 1,nelem
               do iproxy = 1,(porder**ndime)
                  ind = ind+1
                  connecLINEAR(ind,1:8) = connec(ielem,listHEX08(iproxy,1:8))
               end do
            end do
            write(*,*) " ind ",ind

         end subroutine linearMeshOutput

         subroutine create_connecVTK(nelem,connec,atoIJK,vtk_atoIJK,connecVTK)

            implicit none

            integer(4), intent(in)  :: nelem, connec(nelem,nnode), atoIJK(nnode), vtk_atoIJK(nnode)
            integer(4), intent(out) :: connecVTK(nelem,nnode)
            integer(4)              :: i, j, k, ielem, indGmsh, indVTK

            !$acc parallel loop gang
            do ielem = 1,nelem
               !$acc loop vector collapse(3)
               do k = 0,porder
                  do i = 0,porder
                     do j = 0,porder
                        indGmsh = atoIJK(((porder+1)**2)*k+(porder+1)*i+j+1)
                        indVTK = vtk_atoIJK(((porder+1)**2)*k+(porder+1)*i+j+1)
                        connecVTK(ielem,indVTK) = connec(ielem,indGmsh)
                     end do
                  end do
               end do
            end do
            !$acc end parallel loop

         end subroutine create_connecVTK

         subroutine elemPerNode(nelem,npoin,connec,lelpn,point2elem)

            implicit none

            integer(4), intent(in)  :: nelem, npoin, connec(nelem,nnode)
            integer(4), intent(out) :: lelpn(npoin),point2elem(npoin)
            integer(4)              :: aux, ipoin, inode, ielem

            !$acc kernels
            lelpn(:) = 0
            !$acc end kernels
            !$acc parallel loop gang
            do ielem = 1,nelem
               !$acc loop worker
               do inode = 1,nnode
                  point2elem(connec(ielem,inode)) = ielem
                  !$acc loop vector
                  do ipoin = 1,npoin
                     if (connec(ielem,inode) == ipoin) then
                        !$acc atomic update
                        lelpn(ipoin) = lelpn(ipoin)+1
                        !$acc end atomic
                     end if
                  end do
               end do
            end do
            !$acc end parallel loop

         end subroutine elemPerNode

         subroutine nearBoundaryNode(nelem,npoin,nboun,connec,coord,bound,point2elem,lnbn)

            implicit none

            integer(4), intent(in)  :: nelem,npoin,nboun,connec(nelem,nnode),bound(nboun,npbou),point2elem(npoin)
            real(8), intent(in) :: coord(npoin,ndime)
            integer(4), intent(out) :: lnbn(nboun,npbou)
            integer(4)              :: ipoin, inode, ielem,bnode,ipbou,iboun
            real(8)                 :: dist(ndime), dist2,aux,aux2

            !$acc kernels
            lnbn(:,:) = 0
            !$acc end kernels
            !$acc parallel loop gang
            do iboun = 1,nboun
               !$acc loop seq
               do ipbou = 1,npbou
                  bnode = bound(iboun,ipbou)
                  ielem = point2elem(bound(iboun,ipbou))
                  aux = 1000000000000000000000000000000000000000000000000000000000000000000000.0d0
                  !$acc loop seq
                  do inode = 1,nnode
                     if(inode .ne. bnode) then 
                        aux2 = coord(connec(ielem,inode),1)-coord(connec(ielem,bnode),1)
                        if(aux2>1.0e10) then
                           dist = coord(connec(ielem,inode),:)-coord(connec(ielem,bnode),:)
                           dist2 = sqrt(dot_product(dist(:),dist(:)))
                           if(dist2<aux) then
                              aux = dist2
                              lnbn(iboun,ipbou)=connec(ielem,inode)
                           end if
                        end if
                     end if
                  end do
               end do
            end do
            !$acc end parallel loop

         end subroutine nearBoundaryNode

end module mod_geom
