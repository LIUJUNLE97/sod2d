module mod_geom

      use mod_constants
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

         subroutine char_length_spectral(ielem,nelem,npoin,connec,coord,he_l)

                 implicit none

                 integer(4), intent(in)    :: ielem, nelem, npoin, connec(nelem,nnode)
                 real(8),    intent(in)    :: coord(npoin,ndime)
                 real(8),    intent(inout) :: he_l(npoin)
                 integer(4)                :: inode, jnode,idime
                 real(8)                   :: aux,dist(ndime),dist2

                 !
                 ! Compute r = x2-x1 for all element nodes
                 !
                 !
                 ! Obtain ||dist||_2 for all edges and select minimum size as elem. characteristic size
                 !
                 aux = 1000000000000000000000000000000000000000000000000000000000000000000000.0d0
                 do inode = 1,nnode
                    do jnode = 1,nnode
                       if(inode .ne. jnode) then 
                          dist = coord(connec(ielem,inode),:)-coord(connec(ielem,jnode),:)
                          dist2 = sqrt(dot_product(dist(:),dist(:)))
                          aux = min(dist2,aux)
                       end if
                    end do
                    he_l(connec(ielem,inode)) = min(aux,he_l(connec(ielem,inode)))
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

         end subroutine linearMeshOutput

end module mod_geom
