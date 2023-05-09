module mod_gmsh_indices
   use mod_constants
   implicit none

   !-----------------------------------------
   integer(4),parameter :: tV=4,tE=3,tF=2,tI=1
   integer(4),parameter :: tn2ijk(nnode) = [tV,tV,tE,tE,tV,tV,tE,tE,tE,tE,tF,tF,tE,tE,tF,tF,&
                                   tV,tV,tE,tE,tV,tV,tE,tE,tE,tE,tF,tF,tE,tE,tF,tF,&
                                   tE,tE,tF,tF,tE,tE,tF,tF,tF,tF,tI,tI,tF,tF,tI,tI,&
                                   tE,tE,tF,tF,tE,tE,tF,tF,tF,tF,tI,tI,tF,tF,tI,tI]

   integer(4),parameter :: gmsh2ijk(nnode) = [1,4,11,12,2,3,15,16,9,20,33,34,10,19,36,35,&
               5,8,27,28,6,7,29,30,25,32,53,56,26,31,54,55,&
               13,23,41,44,17,21,45,46,37,50,57,60,38,49,58,59,&
               14,24,42,43,18,22,48,47,40,51,61,64,39,52,62,63]

   integer(4),parameter :: gmsh2ij(npbou) = [1,4,12,11,2,3,7,8,5,10,13,16,6,9,14,15]

   integer(4),parameter :: vtk2ijk(nnode) = [1,4,15,16,2,3,11,12,9,13,49,51,10,14,50,52, &
              5,8,23,24,6,7,19,20,17,21,53,55,18,22,54,56, &
              25,31,33,34,27,29,37,38,41,45,57,59,42,46,58,60,&
              26,32,35,36,28,30,39,40,43,47,61,63,44,48,62,64]

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

   integer(4),parameter :: dummy2ijk(nnode)= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,&
                 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,&
                 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,&
                 49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64]

   integer, parameter :: faceFront2ijk(npbou)  = [1,9,13,5,33,41,45,37,49,57,61,53,17,25,29,21]
   integer, parameter :: faceLeft2ijk(npbou)   = [2,4,3,1,34,36,35,33,50,52,51,49,18,20,19,17]
   integer, parameter :: faceTop2ijk(npbou)    = [17,25,29,21,19,27,31,23,20,28,32,24,18,26,30,22]
   integer, parameter :: faceBack2ijk(npbou)   = [2,10,14,6,34,42,46,38,50,58,62,54,18,26,30,22]
   integer, parameter :: faceRight2ijk(npbou)  = [6,8,7,5,38,40,39,37,54,56,55,53,22,24,23,21]
   integer, parameter :: faceBottom2ijk(npbou) = [1,9,13,5,3,11,15,7,4,12,16,8,2,10,14,6]

   integer,parameter :: maxBoundsPerElem = 4
   integer,parameter :: posFaceVertices(4) = [1,2,5,6]
   integer,parameter :: posFaceInnerNodes(4) = [11,12,15,16]
   integer,parameter :: posElemVertices(8) = [1,2,5,6,17,18,21,22]

end module mod_gmsh_indices