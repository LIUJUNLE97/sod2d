module mod_geom

	use mod_numerical_params
	use elem_qua
	use elem_hex
	use mod_mpi
	use mod_mpi_mesh
	use mod_hdf5
	use mod_comms

	contains

		subroutine char_length(mnnode,iElem,nelem,npoin,connec,coord,he)
			implicit none
			integer(4),intent(in) :: mnnode,iElem,nelem,npoin
			integer(4),intent(in) :: connec(nelem,mnnode)
			real(rp),intent(in)   :: coord(npoin,ndime)
			real(rp),intent(out)  :: he
			integer(4)            :: iedge,ncorner,nedge
			integer(4)            :: inode,jnode,idime
			real(rp)              :: dist(12,ndime), dist2, aux

			!
			! Compute r = x2-x1 for all element edges
			!
			call hexa_edges(mnnode,iElem,nelem,npoin,connec,coord,ncorner,nedge,dist(1:12,1:ndime))
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

		subroutine char_length_spectral(mnnode,ielem,nelem,npoin,connec,coord,Ml,he_l)

			implicit none

			integer(4),intent(in)  :: mnnode,ielem,nelem,npoin,connec(nelem,mnnode)
			real(rp),intent(in)    :: coord(npoin,ndime),Ml(npoin)
			real(rp),intent(inout) :: he_l(nelem,mnnode)
			integer(4)             :: inode,jnode,idime
			real(rp)               :: aux,dist(ndime),dist2

			!
			! Compute r = x2-x1 for all element nodes
			!
			!
			! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
			!
			do inode = 1,nnode
				aux = 1000000000000.0_rp
				do jnode = 1,nnode
					if(inode .ne. jnode) then 
						dist = coord(connec(ielem,inode),:)-coord(connec(ielem,jnode),:)
						dist2 = sqrt(dot_product(dist(:),dist(:)))
						aux = min(dist2,aux)
					end if
				end do
				he_l(ielem,inode) = aux
			end do
		end subroutine char_length_spectral

		subroutine elemPerNode(mnnode,numElems,numNodes,connec,lelpn,point2elem)

			implicit none

			integer(4),intent(in)    :: mnnode,numElems,numNodes,connec(numElems,mnnode)
			integer(4),intent(inout) :: lelpn(numNodes),point2elem(numNodes)
			integer(4)               :: iNode,iNodeL,iElem

			!$acc kernels
			lelpn(:) = 0
			point2elem(:) = 0
			!$acc end kernels
			!$acc parallel loop gang
			do iElem = 1,numElems
				!$acc loop seq
				do iNode = 1,mnnode
					iNodeL = connec(iElem,iNode)
					point2elem(iNodeL) = iElem
					!$acc atomic update
					lelpn(iNodeL) = lelpn(iNodeL)+1
					!$acc end atomic
				end do
			end do
			!$acc end parallel loop

			if(mpi_size.ge.2) then
				call mpi_halo_atomic_update_int(lelpn)
			end if

		end subroutine elemPerNode

		subroutine nearBoundaryNode(mporder,mnnode,mnpbou,nelem,npoin,nboun,connec,coord,bound,bouCodesNodes,point2elem,atoIJK,lnbn,lnbnNodes)

			implicit none

			integer(4),intent(in)    :: mporder,mnnode,mnpbou
			integer(4),INTENT(IN)    :: nelem,npoin,nboun,connec(nelem,mnnode),bound(nboun,mnpbou),bouCodesNodes(npoin),point2elem(npoin),atoIJK(mnnode)
			real(rp),intent(in)      :: coord(npoin,ndime)
			integer(4),intent(inout) :: lnbn(nboun,mnpbou)
			integer(4),intent(inout) :: lnbnNodes(npoin)
			integer(4)               :: ipoin, inode,ielem,bnode,ipbou,iboun,rnode,c,i,j,k,innode
			integer(4)               :: aux1, aux2

			! TODO: @Jordi, apparently the "lnbn" array is not used anymore, check if thiis loop can be removed. At any rate, it has a bug :D
			!$acc parallel loop gang
			do iboun = 1,nboun
				!$acc loop seq
				do ipbou = 1,mnpbou
					bnode = bound(iboun,ipbou)
					ielem = point2elem(bnode)
					!$acc loop seq
					do inode = 1,mnnode
						if(connec(ielem,inode) .eq. bnode) then
							aux1=inode
							exit
						end if
					end do
					c = 0
					!$acc loop seq
					outer: do k = 1,mporder+1
						do i = 1,mporder+1
							do j = 1,mporder+1
								c = c+1
								if(atoIJK(c) .eq. aux1) then
									aux2=c
									exit outer
								end if
							end do
						end do
					end do outer
					!lnbn(iboun,ipbou) = connec(ielem,atoIJK(aux2+8))
					lnbn(iboun,ipbou) = connec(ielem,atoIJK(mnnode))
				end do
			end do
			!$acc end parallel loop

			!$acc parallel loop  
			do inode = 1,npoin
				if(bouCodesNodes(inode) .lt. max_num_bou_codes) then
					ielem = point2elem(inode)
					lnbnNodes(inode) = connec(ielem,atoIJK(mnnode)) !this can be done better... using the last node for the moment
				end if
			end do
			!$acc end parallel loop

		end subroutine nearBoundaryNode

		subroutine atoIJKInverse(mporder,mnnode,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)

			implicit none
			integer(4),intent(in)  :: mporder,mnnode,atoIJK(mnnode)
			integer(4),intent(out) :: invAtoIJK(mporder+1,mporder+1,mporder+1),gmshAtoI(mnnode),gmshAtoJ(mnnode),gmshAtoK(mnnode)
			integer(4)             :: i,j,k,c

			c=0
			!$acc loop seq
			do k = 1,mporder+1
				do i = 1,mporder+1
					do j = 1,mporder+1
						c = c+1
						invAtoIJK(i,j,k) = atoIJK(c)
						gmshAtoI(atoIJK(c)) = i
						gmshAtoJ(atoIJK(c)) = j
						gmshAtoK(atoIJK(c)) = k
					end do
				end do
			end do

		end subroutine atoIJKInverse

		subroutine boundary_normals(mnpbou,npoin,nboun,bound,leviCivi,coord,dNgp_b,bounorm)

			implicit none
			integer(4),intent(in) :: mnpbou,npoin,nboun,bound(nboun,mnpbou)
			real(rp),intent(in)   :: coord(npoin,ndime),dNgp_b(ndime-1,mnpbou,mnpbou),leviCivi(ndime,ndime,ndime)
			real(rp),intent(out)  :: bounorm(nboun,ndime*mnpbou)
			integer(4)            :: iboun,inode,jnode,idime,jdime,kdime
			real(rp)              :: xyz(mnpbou,ndime),u(ndime),v(ndime),aux1,aux2

			do iboun = 1,nboun
				do idime = 1,ndime
					do inode = 1,mnpbou
						xyz(inode,idime) = coord(bound(iboun,inode),idime)
					end do
				end do
				do inode = 1,mnpbou
					do idime = 1,ndime
						aux1 = 0.0_rp
						aux2 = 0.0_rp
						do jnode = 1,mnpbou
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
