module mod_ijk_indices
   use mod_constants
   use mod_mpi
   use mod_gmsh_indices !TO DEL WHEN FULLY IMPLEMENTED
   implicit none

   !--------------------------------------------------------------------------------------------
   ! GMSH Indices
   integer(4) :: gmsh2ijk(nnode),gmsh2ij(npbou) 
   integer(4),dimension(npbou) :: faceFront2ijk,faceLeft2ijk,faceTop2ijk,faceBack2ijk,faceRight2ijk,faceBottom2ijk

   integer(4) :: gmsh2ij_vertices(4),gmsh2ijk_vertices(8),gmsh2ij_vertInnerNodes(4)

   !integer(4),parameter :: gmsh2ij_vertices(4) = [1,2,3,4]
   !integer(4),parameter :: gmsh2ijk_vertices(8) = [1,2,3,4,5,6,7,8]
   !integer(4),parameter :: gmsh2ij_vertInnerNodes(4) = [11,12,15,16]

   !--------------------------------------------------------------------------------------------
   ! VTK Indices
   integer(4):: vtk2ijk(nnode),vtk2ij(npbou) 

contains

   subroutine set_ijk_indices()
      implicit none
      integer(4) :: i,j,k,gmshCnt,vtkCnt,gmshIndex,vtkIndex

      if(porder.le.2) then
         write(*,*) 'SOD2D is not ready to work for porder <= 2... You know, #gobigorgohome and set proder >= 3'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      elseif(porder .eq. 3) then

         do i=1,nnode_hexa_p3
            gmsh2ijk(i) = gmsh2ijk_p3(i)
            vtk2ijk(i) = vtk2ijk_p3(i)
         end do

         do i=1,npbou_hexa_p3
            gmsh2ij(i)  = gmsh2ij_p3(i)
            vtk2ij(i)  = vtk2ij_p3(i)

            faceFront2ijk(i) =  faceFront2ijk_p3(i)
            faceLeft2ijk(i) =   faceLeft2ijk_p3(i)
            faceTop2ijk(i) =    faceTop2ijk_p3(i)
            faceBack2ijk(i) =   faceBack2ijk_p3(i)
            faceRight2ijk(i) =  faceRight2ijk_p3(i)
            faceBottom2ijk(i) = faceBottom2ijk_p3(i)
         end do

         gmsh2ij_vertices(:)       = gmsh2ij_vertices_p3(:)
         gmsh2ijk_vertices(:)      = gmsh2ijk_vertices_p3(:)
         gmsh2ij_vertInnerNodes(:) = gmsh2ij_vertInnerNodes_p3(:)

      else 
         write(*,*) 'porder >= 4 not yet implemented'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)

         !--------------------------------------------------------------------
         !Filling high order hexahedra 
         gmshCnt=0
         vtkCnt=0
         do i=0,porder
            do j=0,porder
               do k=0,porder

                  call vtkHigherOrderHexahedron_pointIndexFromIJK(i,j,k,vtkIndex)
                  call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)

                  gmshCnt = gmshCnt+1
                  gmsh2ijk(gmshCnt) = gmshIndex

                  vtkCnt  = vtkCnt+1
                  vtk2ijk(vtkCnt) = vtkIndex

               end do
            end do
         end do
         !--------------------------------------------------------------------

         !--------------------------------------------------------------------
         !Filling gmsh2ijk_vertices(:) 
         !--------------------------------------------------------------------
         i=0
         j=0
         k=0
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(1) = gmshIndex
         !--------------------------------------------------------------------
         i=porder
         j=0
         k=0
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(2) = gmshIndex
         !--------------------------------------------------------------------
         i=porder
         j=porder
         k=0
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(3) = gmshIndex
         !--------------------------------------------------------------------
         i=0
         j=porder
         k=0
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(4) = gmshIndex
         !--------------------------------------------------------------------
         i=0
         j=0
         k=porder
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(5) = gmshIndex
         !--------------------------------------------------------------------
         i=porder
         j=0
         k=porder
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(6) = gmshIndex
         !--------------------------------------------------------------------
         i=porder
         j=porder
         k=porder
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(7) = gmshIndex
         !--------------------------------------------------------------------
         i=0
         j=porder
         k=porder
         call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
         gmsh2ijk_vertices(8) = gmshIndex
         !--------------------------------------------------------------------

         !--------------------------------------------------------------------
         !Filling faceFront2ijk(:)
         gmshCnt=0
         j=0
         do i=0,porder
            do k=0,porder
               call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
               
               gmshCnt = gmshCnt+1
               faceFront2ijk(gmshCnt) = gmshIndex
            end do
         end do
         !--------------------------------------------------------------------
         !Filling faceBack2ijk(:)
         gmshCnt=0
         j=porder
         do i=0,porder
            do k=0,porder
               call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
               
               gmshCnt = gmshCnt+1
               faceBack2ijk(gmshCnt) = gmshIndex
            end do
         end do
         !--------------------------------------------------------------------
         !Filling faceBottom2ijk(:)
         gmshCnt=0
         k=0
         do i=0,porder
            do j=0,porder
               call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
               
               gmshCnt = gmshCnt+1
               faceBottom2ijk(gmshCnt) = gmshIndex
            end do
         end do
         !--------------------------------------------------------------------
         !Filling faceTop2ijk(:)
         gmshCnt=0
         k=porder
         do i=0,porder
            do j=0,porder
               call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
               
               gmshCnt = gmshCnt+1
               faceTop2ijk(gmshCnt) = gmshIndex
            end do
         end do
         !--------------------------------------------------------------------
         !Filling faceLeft2ijk(:)
         gmshCnt=0
         i=0
         do j=0,porder
            do k=0,porder
               call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
               
               gmshCnt = gmshCnt+1
               faceLeft2ijk(gmshCnt) = gmshIndex
            end do
         end do
         !--------------------------------------------------------------------
         !Filling faceRight2ijk(:)
         gmshCnt=0
         i=porder
         do j=0,porder
            do k=0,porder
               call gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,gmshIndex)
               
               gmshCnt = gmshCnt+1
               faceRight2ijk(gmshCnt) = gmshIndex
            end do
         end do
         !--------------------------------------------------------------------

         !--------------------------------------------------------------------
         !filling high order quads
         do i=0,porder
            do j=0,porder
               call vtkHigherOrderQuadrilateral_pointIndexFromIJ(i,j,vtkIndex)
               call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
            end do
         end do

         !--------------------------------------------------------------------
         !Filling gmsh2ij_vertices(:) 
         !--------------------------------------------------------------------
         i=0
         j=0
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertices(1) = gmshIndex
         !--------------------------------------------------------------------
         i=porder
         j=0
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertices(2) = gmshIndex
         !--------------------------------------------------------------------
         i=porder
         j=porder
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertices(3) = gmshIndex
         !--------------------------------------------------------------------
         i=0
         j=porder
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertices(4) = gmshIndex
         !--------------------------------------------------------------------

         !--------------------------------------------------------------------
         !Filling gmsh2ij_vertInnerNodes(:)
         !--------------------------------------------------------------------
         i=1
         j=1
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertInnerNodes(1) = gmshIndex
         !--------------------------------------------------------------------
         i=porder-1
         j=1
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertInnerNodes(2) = gmshIndex
         !--------------------------------------------------------------------
         i=porder-1
         j=porder-1
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertInnerNodes(3) = gmshIndex
         !--------------------------------------------------------------------
         i=1
         j=porder-1
         call gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,gmshIndex)
         gmsh2ij_vertInnerNodes(4) = gmshIndex
         !--------------------------------------------------------------------
         

      end if

   end subroutine set_ijk_indices

   subroutine vtkHigherOrderHexahedron_pointIndexFromIJK(i,j,k,pointIndex)
      implicit none
      integer(4),intent(in) :: i,j,k
      integer(4),intent(out) :: pointIndex
      integer(4) :: ibdy,jbdy,kbdy,nbdy
      integer(4) :: aux_pi

      !----------------------------------------------------------------------------
      
      if((i.eq.0).or.(i.eq.porder)) then
         ibdy = 1
      else 
         ibdy = 0
      endif

      if((j.eq.0).or.(j.eq.porder)) then
         jbdy = 1
      else 
         jbdy = 0
      endif

      if((k.eq.0).or.(k.eq.porder)) then
         kbdy = 1
      else 
         kbdy = 0
      endif

      nbdy = ibdy + jbdy + kbdy

      !----------------------------------------------------------------------------
      pointIndex = 1

      if(nbdy .eq. 3) then !Vertex
         !return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);

         if(i.ne.0) then
            if(j.ne.0) then
               pointIndex = pointIndex + 2
            else
               pointIndex = pointIndex + 1
            end if
         else 
            if(j.ne.0) then
               pointIndex = pointIndex + 3
            end if
         end if

         if(k.ne.0) then
            pointIndex = pointIndex + 4
         end if

         return
      end if

      pointIndex = pointIndex + 8
      if(nbdy .eq. 2) then !Edge

         if(ibdy.eq.0) then !iaxis
            !return (i - 1) + (j ? order[0] + order[1] - 2 : 0) + (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
            pointIndex = pointIndex + (i-1)
            if(j.ne.0) then
               pointIndex = pointIndex + (porder*2 - 2)
            end if
            if(k.ne.0) then
               pointIndex = pointIndex + 2*(porder*2 - 2)
            end if
         else if(jbdy.eq.0) then !jaxis
            !return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
            pointIndex = pointIndex + (j-1)
            if(i.ne.0) then
               pointIndex = pointIndex + (porder - 1)
            else 
                pointIndex = pointIndex + (2*(porder - 1) + porder - 1)
            end if
            if(k.ne.0) then
               pointIndex = pointIndex + 2*(2*porder - 2)
            end if
         else !kaxis
            !offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
            !return (k - 1) + (order[2] - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
            pointIndex = pointIndex + 4*(porder-1)+ 4*(porder-1)

            aux_pi = 0
            if(i.ne.0) then
               if(j.ne.0) then
                  aux_pi = 2
               else
                  aux_pi = 1
               end if
            else
               if(j.ne.0) then
                  aux_pi = 3
               end if
            end if

            pointIndex = pointIndex + (k-1) + (porder - 1)*aux_pi
         end if

         return
      end if

      pointIndex = pointIndex + 4*(3*porder - 3)
      if(nbdy .eq. 1) then !Face
         if(ibdy.ne.0) then ! on i-normal face
            pointIndex = pointIndex + (j-1) + ((porder-1)*(k-1))
            if(i.ne.0) then
               pointIndex = pointIndex + (porder-1)*(porder-1)
            end if
            return
         end if

         pointIndex = pointIndex + 2*(porder-1)*(porder-1)
         if(jbdy.ne.0) then! on j-normal face
            pointIndex = pointIndex + (i-1) + ((porder-1)*(k-1))
            if(j.ne.0) then
               pointIndex = pointIndex + (porder-1)*(porder-1)
            end if
            return
         end if

         ! on k-normal face
         pointIndex = pointIndex + 2*(porder-1)*(porder-1)
         pointIndex = pointIndex + (i-1) + ((porder-1)*(j-1))
         if(k.ne.0) then
            pointIndex = pointIndex + (porder-1)*(porder-1)
         end if
         return

      end if

      ! nbdy == 0: Body DOF
      pointIndex = pointIndex + 2* ((porder-1)*(porder-1)+(porder-1)*(porder-1)+(porder-1)*(porder-1))
      pointIndex = pointIndex + (i-1)+(porder-1)*((j-1)+(porder-1)*((k - 1)))

   end subroutine vtkHigherOrderHexahedron_pointIndexFromIJK

   subroutine vtkHigherOrderQuadrilateral_pointIndexFromIJ(i,j,pointIndex)
      implicit none
      integer(4),intent(in) :: i,j
      integer(4),intent(out) :: pointIndex
      integer(4) :: ibdy,jbdy,nbdy
      integer(4) :: aux_pi

      !----------------------------------------------------------------------------
      
      if((i.eq.0).or.(i.eq.porder)) then
         ibdy = 1
      else 
         ibdy = 0
      endif

      if((j.eq.0).or.(j.eq.porder)) then
         jbdy = 1
      else 
         jbdy = 0
      endif

      nbdy = ibdy + jbdy

      !----------------------------------------------------------------------------
      pointIndex = 1

      if(nbdy .eq. 2) then !Vertex
         !return (i ? (j ? 2 : 1) : (j ? 3 : 0));

         if(i.ne.0) then
            if(j.ne.0) then
               pointIndex = pointIndex + 2
            else
               pointIndex = pointIndex + 1
            end if
         else 
            if(j.ne.0) then
               pointIndex = pointIndex + 3
            end if
         end if

         return
      end if

      pointIndex = pointIndex + 4
      if(nbdy .eq. 1) then !Edge

         if(ibdy.eq.0) then !iaxis
            !return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + offset;

            pointIndex = pointIndex + (i-1)
            if(j.ne.0) then
              pointIndex = pointIndex + (porder*2 - 2)
            end if

         else !jaxis
            !return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + offset;

            pointIndex = pointIndex + (j-1)
            if(i.ne.0) then
               pointIndex = pointIndex + (porder - 1)
            else 
                pointIndex = pointIndex + (2*(porder - 1) + porder - 1)
            end if

         end if

         return
      end if

      ! nbdy == 0: Body DOF
      pointIndex = pointIndex + (4*porder-4)
      pointIndex = pointIndex + (i-1)+(porder-1)*(j-1)

   end subroutine vtkHigherOrderQuadrilateral_pointIndexFromIJ

   subroutine gmshHigherOrderHexahedron_pointIndexFromIJK(i,j,k,pointIndex)
      implicit none
      integer(4),intent(in) :: i,j,k
      integer(4),intent(out) :: pointIndex

      !YET TO BE IMPLEMENTED! COME ON BENET!
   end subroutine gmshHigherOrderHexahedron_pointIndexFromIJK

   subroutine gmshHigherOrderQuadrilateral_pointIndexFromIJ(i,j,pointIndex)
      implicit none
      integer(4),intent(in) :: i,j
      integer(4),intent(out) :: pointIndex

      !YET TO BE IMPLEMENTED! COME ON BENET!
   end subroutine gmshHigherOrderQuadrilateral_pointIndexFromIJ

#if 0

int vtkHigherOrderQuadrilateral::PointIndexFromIJK(int i, int j, const int* order)
{
  bool ibdy = (i == 0 || i == order[0]);
  bool jbdy = (j == 0 || j == order[1]);
  // How many boundaries do we lie on at once?
  int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

  if (nbdy == 2) // Vertex DOF
  {              // ijk is a corner node. Return the proper index (somewhere in [0,7]):
    return (i ? (j ? 2 : 1) : (j ? 3 : 0));
  }

  int offset = 4;
  if (nbdy == 1) // Edge DOF
  {
    if (!ibdy)
    { // On i axis
      return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + offset;
    }
    if (!jbdy)
    { // On j axis
      return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + offset;
    }
  }

  offset += 2 * (order[0] - 1 + order[1] - 1);
  // nbdy == 0: Face DOF
  return offset + (i - 1) + (order[0] - 1) * ((j - 1));
}



int vtkHigherOrderHexahedron::PointIndexFromIJK(int i, int j, int k, const int* order)
{
  bool ibdy = (i == 0 || i == order[0]);
  bool jbdy = (j == 0 || j == order[1]);
  bool kbdy = (k == 0 || k == order[2]);
  // How many boundaries do we lie on at once?
  int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

  if (nbdy == 3) // Vertex DOF
  {              // ijk is a corner node. Return the proper index (somewhere in [0,7]):
    return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
  }

  int offset = 8;
  if (nbdy == 2) // Edge DOF
  {
    if (!ibdy)
    { // On i axis
      return (i - 1) + (j ? order[0] + order[1] - 2 : 0) + (k ? 2 * (order[0] + order[1] - 2) : 0) +
        offset;
    }
    if (!jbdy)
    { // On j axis
      return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
        (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
    }
    // !kbdy, On k axis
    offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
    return (k - 1) + (order[2] - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
  }

  offset += 4 * (order[0] + order[1] + order[2] - 3);
  if (nbdy == 1) // Face DOF
  {
    if (ibdy) // On i-normal face
    {
      return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) +
        offset;
    }
    offset += 2 * (order[1] - 1) * (order[2] - 1);
    if (jbdy) // On j-normal face
    {
      return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) +
        offset;
    }
    offset += 2 * (order[2] - 1) * (order[0] - 1);
    // kbdy, On k-normal face
    return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) +
      offset;
  }

  // nbdy == 0: Body DOF
  offset += 2 *
    ((order[1] - 1) * (order[2] - 1) + (order[2] - 1) * (order[0] - 1) +
      (order[0] - 1) * (order[1] - 1));
  return offset + (i - 1) + (order[0] - 1) * ((j - 1) + (order[1] - 1) * ((k - 1)));
}
#endif


end module mod_ijk_indices