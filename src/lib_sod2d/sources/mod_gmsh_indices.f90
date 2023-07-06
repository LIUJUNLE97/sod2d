module mod_gmsh_indices
   use mod_constants
   use elem_qua ! Using only the ordering table for edges
   use elem_hex ! Using only the ordering tables for edges and faces

   implicit none
   integer(4),parameter :: nnode_hexa_p3 = 64
   integer(4),parameter :: npbou_hexa_p3 = 16

   !--------------------------------------------------------------------------------------------
   ! GMSH Indices

   integer(4),parameter :: gmsh2ijk_p3(nnode_hexa_p3) = [1,4,11,12,2,3,15,16,9,20,33,34,10,19,36,35,&
               5,8,27,28,6,7,29,30,25,32,53,56,26,31,54,55,&
               13,23,41,44,17,21,45,46,37,50,57,60,38,49,58,59,&
               14,24,42,43,18,22,48,47,40,51,61,64,39,52,62,63]

   integer(4),parameter :: gmsh2ij_p3(npbou_hexa_p3) = [1,4,12,11,2,3,7,8,5,10,13,16,6,9,14,15]

   integer(4),parameter :: gmsh2ij_vertices_p3(4) = [1,2,3,4]
   integer(4),parameter :: gmsh2ijk_vertices_p3(8) = [1,2,3,4,5,6,7,8]
   integer(4),parameter :: gmsh2ij_vertInnerNodes_p3(4) = [11,12,15,16]


   integer, parameter :: facefront2ijk_p3(npbou_hexa_p3)  = [1,9,13,5,33,41,45,37,49,57,61,53,17,25,29,21]
   integer, parameter :: faceleft2ijk_p3(npbou_hexa_p3)   = [2,4,3,1,34,36,35,33,50,52,51,49,18,20,19,17]
   integer, parameter :: facetop2ijk_p3(npbou_hexa_p3)    = [17,25,29,21,19,27,31,23,20,28,32,24,18,26,30,22]
   integer, parameter :: faceback2ijk_p3(npbou_hexa_p3)   = [2,10,14,6,34,42,46,38,50,58,62,54,18,26,30,22]
   integer, parameter :: faceright2ijk_p3(npbou_hexa_p3)  = [6,8,7,5,38,40,39,37,54,56,55,53,22,24,23,21]
   integer, parameter :: facebottom2ijk_p3(npbou_hexa_p3) = [1,9,13,5,3,11,15,7,4,12,16,8,2,10,14,6]

   !--------------------------------------------------------------------------------------------
   ! vtk indices

   integer(4),parameter :: vtk2ijk_p3(nnode_hexa_p3) = [1,4,15,16,2,3,11,12,9,13,49,51,10,14,50,52, &
              5,8,23,24,6,7,19,20,17,21,53,55,18,22,54,56, &
              25,31,33,34,27,29,37,38,41,45,57,59,42,46,58,60,&
              26,32,35,36,28,30,39,40,43,47,61,63,44,48,62,64]

   integer(4),parameter :: vtk2ij_p3(npbou_hexa_p3) = [1,4,11,12,2,3,7,8,5,9,13,15,6,10,14,16]
 
   !--------------------------------------------------------------------------------------------
   ! Other Indices

   !integer(4),parameter :: dummy2ijk(nnode_hexa_p3)= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,&
   !              17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,&
   !              33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,&
   !              49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64]



!---------------------------------------------------------
!GARBAGE/DEPRECATED SECTION
#if 0
   integer(4),parameter :: tV=4,tE=3,tF=2,tI=1
   integer(4),parameter :: tn2ijk(nnode) = [tV,tV,tE,tE,tV,tV,tE,tE,tE,tE,tF,tF,tE,tE,tF,tF,&
                                   tV,tV,tE,tE,tV,tV,tE,tE,tE,tE,tF,tF,tE,tE,tF,tF,&
                                   tE,tE,tF,tF,tE,tE,tF,tF,tF,tF,tI,tI,tF,tF,tI,tI,&
                                   tE,tE,tF,tF,tE,tE,tF,tF,tF,tF,tI,tI,tF,tF,tI,tI]

   !the one from lucas old code
   !integer(4),parameter :: vtk2ijk(nnode) = [1,4,15,16,2,3,11,12,9,13,49,51,10,14,50,52, &
   !           5,8,23,24,6,7,19,20,17,21,53,55,18,22,54,56, &
   !           25,29,33,34,27,31,37,38,41,45,57,59,42,46,58,60,&
   !           26,30,35,36,28,32,39,40,43,47,61,63,44,48,62,64]

   integer(4),parameter :: cgns2ijk(nnode)= [1,4,16,15,2,3,11,12,9,14,33,36,10,13,34,35,&
                 5,8,32,31,6,7,27,28,25,30,53,56,26,29,54,55,&
                 17,21,50,49,19,23,41,42,37,46,57,60,38,45,58,59,&
                 18,22,51,52,20,24,44,43,40,47,61,64,39,48,62,63]

!  according to cgns documentation....-----------------------------------------
!  integer(4),parameter :: cgns2ijk(nnode)= [1,4,16,15,2,3,11,12,9,14,33,36,10,13,34,35,&
!                5,8,32,31,6,7,27,28,25,30,53,56,26,29,54,55,&
!                17,23,50,49,19,21,41,42,37,46,57,60,38,45,58,59,&
!                18,24,51,52,20,22,44,43,40,47,61,64,39,48,62,63]

   integer,parameter :: posFaceVertices(4) = [1,2,5,6]
   integer,parameter :: posElemVertices(8) = [1,2,5,6,17,18,21,22]

#endif
!---------------------------------------------------------

   contains

      recursive subroutine genHighOrderHex(p,indexTable)
         implicit none
         integer(4), intent(in)  :: p
         integer(4), intent(out) :: indexTable((p+1)**3,3) ! Set as I, J, K
         integer(4)              :: tableFace((p-1)**2,2), tableVolume((p-1)**3,3)
         integer(4)              :: inode, iedge, iface, i0, i1, i3, u(3), v(3), i

         ! Initialize corner node to 0th position, or generate element of order 0
         indexTable(1,:) = [0,0,0]

         ! Generate element of order 1 (corner nodes)
         if (p .gt. 0) then
            !                  i, j, k
            indexTable(2,:) = [p, 0, 0]
            indexTable(3,:) = [p, p, 0]
            indexTable(4,:) = [0, p, 0]
            indexTable(5,:) = [0, 0, p]
            indexTable(6,:) = [p, 0, p]
            indexTable(7,:) = [p, p, p]
            indexTable(8,:) = [0, p, p]
            if (p .gt. 1) then
               ! Generate high-order edges
               inode = 9
               do iedge = 1,12
                  i0 = hex_order_edges(iedge,1)
                  i1 = hex_order_edges(iedge,2)
                  u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
                  do i = 1,p-1
                     indexTable(inode,:) = indexTable(i0,:) + i*u(:)
                     inode = inode + 1
                  end do
               end do
               ! Generate a generic high-order face with p` = p-2
               call genHighOrderQuad(p-2,tableFace)
               tableFace = tableFace + 1
               ! Generate faces interior nodes
               do iface = 1,6
                  i0 = hex_order_faces(iface,1)
                  i1 = hex_order_faces(iface,2)
                  i3 = hex_order_faces(iface,4)
                  u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
                  v(:) = (indexTable(i3,:) - indexTable(i0,:))/p
                  do i = 1,((p-1)**2)
                     indexTable(inode,:) = indexTable(i0,:) + u(:)*tableFace(i,1) + v(:)*tableFace(i,2)
                     inode = inode + 1
                  end do
               end do
               ! Generate volume nodes
               call genHighOrderHex(p-2,tableVolume)
               tableVolume = tableVolume + 1
               print*, "tableVolume = ", tableVolume
               call joinTables([(p-1)**3,3],tableVolume,inode,[(p+1)**3,3],indexTable)
            end if
         end if
      end subroutine genHighOrderHex

      recursive subroutine genHighOrderQuad(p,indexTable)
         implicit none
         integer(4), intent(in)  :: p
         integer(4), intent(out) :: indexTable((p+1)**2,2) ! Set as I, J
         integer(4)              :: tableFace((p-1)**2,2)
         integer(4)              :: inode, iedge, iface, i0, i1, u(2), i

         indexTable(1,:) = [0,0]
         if (p .gt. 0) then
            indexTable(2,:) = [p,0]
            indexTable(3,:) = [p,p]
            indexTable(4,:) = [0,p]
            if (p .gt. 1) then
               inode = 5
               do iedge = 1,4
                  i0 = quad_order_edges(iedge,1)
                  i1 = quad_order_edges(iedge,2)
                  u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
                  do i = 1,p-1
                     indexTable(inode,:) = indexTable(i0,:) + i*u(:)
                     inode = inode + 1
                  end do
               end do
               call genHighOrderQuad(p-2,tableFace)
               tableFace = tableFace + 1
               call joinTables([(p-1)**2,2],tableFace,inode,[(p+1)**2,2],indexTable)
            end if
         end if
      end subroutine genHighOrderQuad

      subroutine joinTables(size1,table1,indexDesti,size2,table2)
         implicit none
         integer(4), intent(in)    :: indexDesti, size1(2), size2(2)
         integer(4), intent(in)    :: table1(size1(1),size1(2))
         integer(4), intent(inout) :: table2(size2(1),size2(2))
         integer(4)                :: i, j

         j = indexDesti
         do i = 1,size1(1)
            table2(j,:) = table1(i,:)
            j = j + 1
         end do
      end subroutine joinTables

end module mod_gmsh_indices