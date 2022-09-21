module mod_geom

   use mod_constants
   use elem_qua
   use elem_hex
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

      contains

         subroutine char_length(iElem,nelem,npoin,connec,coord,he)
            implicit none
            integer(4), intent(in)  :: iElem, nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),   intent(in)  :: coord(npoin,ndime)
            real(rp),   intent(out) :: he
            integer(4)              :: iedge, ncorner, nedge
            real(rp)                :: dist(12,ndime), dist2, aux

                 !
                 ! Compute r = x2-x1 for all element edges
                 !
#if 0
!OLD SHIT
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
                       !call hexa_edges(ielem,nelem,npoin,connec,coord,ncorner,nedge,dist(1:12,1:ndime))

                    else
                       write(*,*) "INCORRECT ELEMENT TYPE (NODES ERROR)!"
                    end if
                 else
                    write(*,*) "BY SIGMAR NO!"
                 end if
#endif
                 call hexa_edges(iElem,nelem,npoin,connec,coord,ncorner,nedge,dist(1:12,1:ndime))
                 !
                 ! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
                 !
                 dist2 = 1000000000000.0_rp
                 do iedge = 1,nedge
                    aux = sqrt(dot_product(dist(iedge,:),dist(iedge,:)))
                    dist2 = min(dist2,aux)
                 end do
                 he = dist2

         end subroutine char_length

         subroutine char_length_spectral(ielem,nelem,npoin,connec,coord,Ml,he_l)

                 implicit none

                 integer(4), intent(in)  :: ielem, nelem, npoin, connec(nelem,nnode)
                 real(rp), intent(in)    :: coord(npoin,ndime), Ml(npoin)
                 real(rp), intent(inout) :: he_l(nelem,nnode)
                 integer(4)              :: inode, jnode,idime
                 real(rp)                :: aux,dist(ndime),dist2

                 !
                 ! Compute r = x2-x1 for all element nodes
                 !
                 !
                 ! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
                 !
                 !aux = 0.0_rp
                 !aux = 1000000000000000000000000000000000000000000000000000000000000000000000.0_rp
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
                    he_l(ielem,inode) = Ml(connec(ielem,inode))**(1.0_rp/real(ndime,rp))
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
            !write(*,*) " ind ",ind

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
            integer(4)              ::  ipoin, inode, ielem, aux

            !$acc kernels
            lelpn(:) = 0
            point2elem(:) = 0
            !$acc end kernels
            !$acc parallel loop gang
            do ielem = 1,nelem
               !$acc loop seq
               do inode = 1,nnode
                  point2elem(connec(ielem,inode)) = ielem
               end do
               !$acc loop worker
               do inode = 1,nnode
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

            if(mpi_size.ge.2) then
               call mpi_halo_atomic_update_int(lelpn)
            end if

         end subroutine elemPerNode

         subroutine nearBoundaryNode(nelem,npoin,nboun,connec,coord,bound,point2elem,atoIJK,lnbn)

            implicit none

            integer(4), intent(in)  :: nelem,npoin,nboun,connec(nelem,nnode),bound(nboun,npbou),point2elem(npoin),atoIJK(nnode)
            real(rp), intent(in) :: coord(npoin,ndime)
            integer(4), intent(out) :: lnbn(nboun,npbou)
            integer(4)              :: ipoin, inode,ielem,bnode,ipbou,iboun,rnode,c,i,j,k,innode
            integer(4)              :: aux1, aux2

            !$acc parallel loop gang
            do iboun = 1,nboun
               !$acc loop seq
               do ipbou = 1,npbou
                  bnode = bound(iboun,ipbou)
                  ielem = point2elem(bnode)
                  !$acc loop seq
                  do inode = 1,nnode
                     if(connec(ielem,inode) .eq. bnode) then 
                        aux1=inode
                        exit
                     end if
                  end do
                  c = 0
                  !$acc loop seq
                  outer: do k = 1,porder+1
                     do i = 1,porder+1
                        do j = 1,porder+1
                           c = c+1
                           if(atoIJK(c) .eq. aux1) then
                              aux2=c
                              exit outer
                           end if
                        end do
                     end do
                  end do outer

                  !lnbn(iboun,ipbou) = connec(ielem,atoIJK(aux2+8))
                  lnbn(iboun,ipbou) = connec(ielem,atoIJK(64))
               end do
            end do
            !$acc end parallel loop

         end subroutine nearBoundaryNode

         subroutine atioIJKInverse(atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)

            implicit none

            integer(4), intent(in)  :: atoIJK(nnode)
            integer(4), intent(out) :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            integer(4)              :: i,j,k,c

            c=0
            !$acc loop seq
            do k = 1,porder+1
               do i = 1,porder+1
                  do j = 1,porder+1
                     c = c+1
                     invAtoIJK(i,j,k) = atoIJK(c)
                     gmshAtoI(atoIJK(c)) = i
                     gmshAtoJ(atoIJK(c)) = j
                     gmshAtoK(atoIJK(c)) = k
                  end do
               end do
            end do
         end subroutine atioIJKInverse

			subroutine boundary_normals(npoin,nboun,bound,leviCivi,coord,dNgp_b,bounorm)
				implicit none
				integer(4), intent(in)  :: npoin, nboun, bound(nboun,npbou)
				real(rp),    intent(in)  :: coord(npoin,ndime), dNgp_b(ndime-1,npbou,npbou), leviCivi(ndime,ndime,ndime)
				real(rp),    intent(out) :: bounorm(nboun,ndime*npbou)
				integer(4)              :: iboun, inode, jnode, idime, jdime, kdime
				real(rp)                 :: xyz(npbou,ndime), u(ndime), v(ndime), aux1, aux2
				do iboun = 1,nboun
					do idime = 1,ndime
						do inode = 1,npbou
							xyz(inode,idime) = coord(bound(iboun,inode),idime)
						end do
					end do
					do inode = 1,npbou
						do idime = 1,ndime
							aux1 = 0.0_rp
							aux2 = 0.0_rp
							do jnode = 1,npbou
								aux1 = aux1+dNgp_b(1,jnode,inode)*xyz(jnode,idime)
								aux2 = aux2+dNgp_b(2,jnode,inode)*xyz(jnode,idime)
							end do
							u(idime) = aux1
							v(idime) = aux2
						end do
						do idime = 1,ndime
							aux1 = 0.0_rp
							do jdime = 1,ndime
								do kdime = 1,ndime
									aux1 = aux1 + leviCivi(idime,jdime,kdime)*u(jdime)*v(kdime)
								end do
							end do
							bounorm(iboun,(inode-1)*ndime+idime) = aux1
						end do
					end do
				end do
			end subroutine boundary_normals

end module mod_geom
