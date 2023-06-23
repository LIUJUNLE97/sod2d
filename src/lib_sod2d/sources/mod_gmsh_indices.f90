module mod_gmsh_indices
   use mod_constants
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

end module mod_gmsh_indices