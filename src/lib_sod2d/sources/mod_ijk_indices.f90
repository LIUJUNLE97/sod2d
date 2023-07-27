module mod_ijk_indices
   use mod_constants
   use mod_mpi
   use mod_gmsh_indices !TO DEL WHEN FULLY IMPLEMENTED
   implicit none

   !--------------------------------------------------------------------------------------------
   ! GMSH Indices
   integer(4),allocatable :: gmshHexahedraHO_ijkTable(:,:,:),gmshQuadrilateralHO_ijTable(:,:)
   integer(4) :: gmsh_porder=0
   !logical :: gmshTablesGenerated=false
!   integer(4) :: gmsh2ijk(nnode),gmsh2ij(npbou) 
!   integer(4),dimension(npbou) :: faceFront2ijk,faceLeft2ijk,faceTop2ijk,faceBack2ijk,faceRight2ijk,faceBottom2ijk

   !integer(4) :: gmsh2ij_vertices(4),gmsh2ijk_vertices(8),gmsh2ij_vertInnerNodes(4)

   !integer(4),parameter :: gmsh2ij_vertices(4) = [1,2,3,4]
   !integer(4),parameter :: gmsh2ijk_vertices(8) = [1,2,3,4,5,6,7,8]
   !integer(4),parameter :: gmsh2ij_vertInnerNodes(4) = [11,12,15,16]

   !--------------------------------------------------------------------------------------------
   ! VTK Indices
   !integer(4):: vtk2ijk(nnode),vtk2ij(npbou) 

contains

   subroutine get_porder_values(mporder,mnnode,mngaus,mnpbou)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(out) :: mnnode,mngaus,mnpbou

      mnnode = (mporder+1)**3
      mngaus = mnnode
      mnpbou = (mporder+1)**2

   end subroutine get_porder_values

   subroutine set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout),dimension(:),allocatable :: ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d
      integer(4) :: i,j

      !--------------------------------------------------------------------
      !defining criteria for ijk
      allocate(ijk_sod2d_to_gmsh(mporder+1))
      allocate(ijk_gmsh_to_sod2d(0:mporder))

      i=1
      ijk_sod2d_to_gmsh(i) = 0
      ijk_gmsh_to_sod2d(ijk_sod2d_to_gmsh(i)) = i

      i=i+1
      ijk_sod2d_to_gmsh(i) = mporder
      ijk_gmsh_to_sod2d(ijk_sod2d_to_gmsh(i)) = i

      do j=1,(mporder-1)
         i=i+1
         ijk_sod2d_to_gmsh(i)=j
         ijk_gmsh_to_sod2d(ijk_sod2d_to_gmsh(i)) = i
      end do
   end subroutine set_allocate_array_ijk_sod2d_criteria

   function get_indexIJK_sod2d(mporder,i,j,k) result(indexIJK)
      implicit none
      integer(4),intent(in) :: mporder,i,j,k
      integer :: indexIJK
      
      indexIJK = ((mporder+1)**2)*k+(mporder+1)*i+j+1
   end function get_indexIJK_sod2d

   subroutine set_allocate_hexahedronHO_ijk_indices(mporder,gmsh2ijk,vtk2ijk,ijk2gmsh,ijk2vtk)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),allocatable,dimension(:),intent(inout) :: gmsh2ijk,vtk2ijk
      integer(4),allocatable,intent(inout),optional :: ijk2gmsh(:,:,:),ijk2vtk(:,:,:)
      integer(4) :: mnnode,mngaus,mnpbou
      integer(4) :: i,j,k,ip,jp,kp,gmshCnt,vtkCnt,gmshIndex,vtkIndex,pIndex
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      !-----------------------------------------------------------------------------------------

      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'mnnode',mnnode,'mngaus',mngaus,'mnpbou',mnpbou

      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      !-----------------------------------------------------------------------------------------

      allocate(gmsh2ijk(mnnode))
      allocate(vtk2ijk(mnnode))

      if(present(ijk2gmsh)) allocate(ijk2gmsh(0:porder,0:porder,0:porder))
      if(present(ijk2vtk))  allocate(ijk2vtk(0:porder,0:porder,0:porder))

      if(mporder.le.2) then
         write(*,*) 'SOD2D is not ready to work for mporder <= 2... You know, #gobigorgohome and set mporder >= 3'
         call MPI_Abort(MPI_COMM_WORLD,-1,mpi_err)
      end if

      !--------------------------------------------------------------------
      !Filling high order hexahedra 
      pIndex=0
      do kp=1,mporder+1
         k=ijk_sod2d_to_gmsh(kp)
         do ip=1,mporder+1
            i=ijk_sod2d_to_gmsh(ip)
            do jp=1,mporder+1
               j=ijk_sod2d_to_gmsh(jp)

               call vtkHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,vtkIndex)
               call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)

               pIndex               = pIndex + 1
               gmsh2ijk(pIndex)     = gmshIndex
               vtk2ijk(pIndex)      = vtkIndex

               if(present(ijk2gmsh)) then

               end if

               if(present(ijk2vtk)) then

               end if

               !write(*,*) 'i,j,k',i,j,k,'ip,jp,kp',ip,jp,kp,'pI',pIndex,'gmsh',gmshIndex,'vtk',vtkIndex

               !a2ijk(pIndex)        = gmshIndex!pIndex
               !a2i(gmshIndex)       = i+1
               !a2j(gmshIndex)       = j+1
               !a2k(gmshIndex)       = k+1
               !ijk2a(i+1,j+1,k+1)   = gmshIndex

            end do
         end do
      end do
      !--------------------------------------------------------------------

      deallocate(ijk_sod2d_to_gmsh)
      deallocate(ijk_gmsh_to_sod2d)

   end subroutine set_allocate_hexahedronHO_ijk_indices

   subroutine set_allocate_quadrilateralHO_ij_indices(mporder,gmsh2ij,vtk2ij)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),allocatable,dimension(:),intent(inout) :: gmsh2ij,vtk2ij
      integer(4) :: mnnode,mngaus,mnpbou
      integer(4) :: i,j,ip,jp,gmshCnt,vtkCnt,gmshIndex,vtkIndex,pIndex
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      !-----------------------------------------------------------------------------------------
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'mnnode',mnnode,'mngaus',mngaus,'mnpbou',mnpbou

      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      !-----------------------------------------------------------------------------------------

      allocate(gmsh2ij(mnpbou))
      allocate(vtk2ij(mnpbou))

      !--------------------------------------------------------------------
      !filling high order quads
      pIndex=0
      do ip=1,mporder+1
         i=ijk_sod2d_to_gmsh(ip)
         do jp=1,mporder+1
            j=ijk_sod2d_to_gmsh(jp)
            call vtkHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,vtkIndex)
            call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)

            pIndex = pIndex + 1
            gmsh2ij(pIndex) = gmshIndex
            vtk2ij(pIndex)  = vtkIndex

         end do
      end do

      deallocate(ijk_sod2d_to_gmsh)
      deallocate(ijk_gmsh_to_sod2d)

   end subroutine set_allocate_quadrilateralHO_ij_indices

   subroutine get_gmshHexHOVertIndex(mporder,gmshHexVertInd)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout) :: gmshHexVertInd(8)
      integer(4) :: i,j,k,gmshIndex

      !--------------------------------------------------------------------
      !Filling gmshHexVertInd(:) 
      !--------------------------------------------------------------------
      i=0
      j=0
      k=0
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(1) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=0
      k=0
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(2) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=mporder
      k=0
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(3) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=mporder
      k=0
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(4) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=0
      k=mporder
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(5) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=0
      k=mporder
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(6) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=mporder
      k=mporder
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(7) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=mporder
      k=mporder
      call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
      gmshHexVertInd(8) = gmshIndex
      !--------------------------------------------------------------------

   end subroutine get_gmshHexHOVertIndex

   subroutine get_gmshHexHOFacesIndex(mporder,mnpbou,gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd)
      implicit none
      integer(4),intent(in) :: mporder,mnpbou
      integer(4),intent(inout),dimension(mnpbou) :: gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd
      integer(4) :: i,j,k,gmshIndex,gmshCnt

      !--------------------------------------------------------------------
      !Filling faceFront2ijk(:)
      gmshCnt=0
      j=0
      do i=0,mporder
         do k=0,mporder
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            
            gmshCnt = gmshCnt+1
            gmshHexFaceFrontInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceBack2ijk(:)
      gmshCnt=0
      j=mporder
      do i=0,mporder
         do k=0,mporder
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            
            gmshCnt = gmshCnt+1
            gmshHexFaceBackInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceBottom2ijk(:)
      gmshCnt=0
      k=0
      do i=0,mporder
         do j=0,mporder
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            
            gmshCnt = gmshCnt+1
            gmshHexFaceBottomInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceTop2ijk(:)
      gmshCnt=0
      k=mporder
      do i=0,mporder
         do j=0,mporder
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            
            gmshCnt = gmshCnt+1
            gmshHexFaceTopInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceLeft2ijk(:)
      gmshCnt=0
      i=0
      do j=0,mporder
         do k=0,mporder
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            
            gmshCnt = gmshCnt+1
            gmshHexFaceLeftInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------
      !Filling faceRight2ijk(:)
      gmshCnt=0
      i=mporder
      do j=0,mporder
         do k=0,mporder
            call gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,gmshIndex)
            
            gmshCnt = gmshCnt+1
            gmshHexFaceRightInd(gmshCnt) = gmshIndex
         end do
      end do
      !--------------------------------------------------------------------

   end subroutine get_gmshHexHOFacesIndex

   subroutine get_gmshQuadHOVertIndex(mporder,gmshQuadVertInd)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout) :: gmshQuadVertInd(4)
      integer(4) :: i,j,gmshIndex

      i=0
      j=0
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadVertInd(1) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=0
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadVertInd(2) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder
      j=mporder
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadVertInd(3) = gmshIndex
      !--------------------------------------------------------------------
      i=0
      j=mporder
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadVertInd(4) = gmshIndex
      !--------------------------------------------------------------------

   end subroutine get_gmshQuadHOVertIndex

   subroutine get_gmshQuadHOInnerVertIndex(mporder,gmshQuadInnerVertInd)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4),intent(inout) :: gmshQuadInnerVertInd(4)
      integer(4) :: i,j,gmshIndex


      !--------------------------------------------------------------------
      !Filling gmshQuadInnerVertInd(:)
      !--------------------------------------------------------------------
      i=1
      j=1
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadInnerVertInd(1) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder-1
      j=1
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadInnerVertInd(2) = gmshIndex
      !--------------------------------------------------------------------
      i=mporder-1
      j=mporder-1
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadInnerVertInd(3) = gmshIndex
      !--------------------------------------------------------------------
      i=1
      j=mporder-1
      call gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,gmshIndex)
      gmshQuadInnerVertInd(4) = gmshIndex
      !--------------------------------------------------------------------

   end subroutine get_gmshQuadHOInnerVertIndex




   subroutine vtkHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j,k
      integer(4),intent(out) :: pointIndex
      integer(4) :: ibdy,jbdy,kbdy,nbdy
      integer(4) :: aux_pi

      !----------------------------------------------------------------------------
      
      if((i.eq.0).or.(i.eq.mporder)) then
         ibdy = 1
      else 
         ibdy = 0
      endif

      if((j.eq.0).or.(j.eq.mporder)) then
         jbdy = 1
      else 
         jbdy = 0
      endif

      if((k.eq.0).or.(k.eq.mporder)) then
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
               pointIndex = pointIndex + (mporder*2 - 2)
            end if
            if(k.ne.0) then
               pointIndex = pointIndex + 2*(mporder*2 - 2)
            end if
         else if(jbdy.eq.0) then !jaxis
            !return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
            pointIndex = pointIndex + (j-1)
            if(i.ne.0) then
               pointIndex = pointIndex + (mporder - 1)
            else 
                pointIndex = pointIndex + (2*(mporder - 1) + mporder - 1)
            end if
            if(k.ne.0) then
               pointIndex = pointIndex + 2*(2*mporder - 2)
            end if
         else !kaxis
            !offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
            !return (k - 1) + (order[2] - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
            pointIndex = pointIndex + 4*(mporder-1)+ 4*(mporder-1)

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

            pointIndex = pointIndex + (k-1) + (mporder - 1)*aux_pi
         end if

         return
      end if

      pointIndex = pointIndex + 4*(3*mporder - 3)
      if(nbdy .eq. 1) then !Face
         if(ibdy.ne.0) then ! on i-normal face
            pointIndex = pointIndex + (j-1) + ((mporder-1)*(k-1))
            if(i.ne.0) then
               pointIndex = pointIndex + (mporder-1)*(mporder-1)
            end if
            return
         end if

         pointIndex = pointIndex + 2*(mporder-1)*(mporder-1)
         if(jbdy.ne.0) then! on j-normal face
            pointIndex = pointIndex + (i-1) + ((mporder-1)*(k-1))
            if(j.ne.0) then
               pointIndex = pointIndex + (mporder-1)*(mporder-1)
            end if
            return
         end if

         ! on k-normal face
         pointIndex = pointIndex + 2*(mporder-1)*(mporder-1)
         pointIndex = pointIndex + (i-1) + ((mporder-1)*(j-1))
         if(k.ne.0) then
            pointIndex = pointIndex + (mporder-1)*(mporder-1)
         end if
         return

      end if

      ! nbdy == 0: Body DOF
      pointIndex = pointIndex + 2* ((mporder-1)*(mporder-1)+(mporder-1)*(mporder-1)+(mporder-1)*(mporder-1))
      pointIndex = pointIndex + (i-1)+(mporder-1)*((j-1)+(mporder-1)*((k - 1)))

   end subroutine vtkHigherOrderHexahedron_pointIndexFromIJK

   subroutine vtkHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j
      integer(4),intent(out) :: pointIndex
      integer(4) :: ibdy,jbdy,nbdy
      integer(4) :: aux_pi

      !----------------------------------------------------------------------------
      
      if((i.eq.0).or.(i.eq.mporder)) then
         ibdy = 1
      else 
         ibdy = 0
      endif

      if((j.eq.0).or.(j.eq.mporder)) then
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
              pointIndex = pointIndex + (mporder*2 - 2)
            end if

         else !jaxis
            !return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + offset;

            pointIndex = pointIndex + (j-1)
            if(i.ne.0) then
               pointIndex = pointIndex + (mporder - 1)
            else 
                pointIndex = pointIndex + (2*(mporder - 1) + mporder - 1)
            end if

         end if

         return
      end if

      ! nbdy == 0: Body DOF
      pointIndex = pointIndex + (4*mporder-4)
      pointIndex = pointIndex + (i-1)+(mporder-1)*(j-1)

   end subroutine vtkHigherOrderQuadrilateral_pointIndexFromIJ

   subroutine initGmshIJKTables(mporder)
      implicit none
      integer(4),intent(in) :: mporder
      integer(4) :: mnnode,mngaus,mnpbou,pIndex,i,j,k
      integer(4),allocatable :: auxHexHOtable(:,:),auxQuadHOtable(:,:)

      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      !if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'mnnode',mnnode,'mngaus',mngaus,'mnpbou',mnpbou

      if(gmsh_porder.eq.0) then !arrays not initialized
         if(mpi_rank.eq.0) write(*,*) 'Initialising GMSH IJK Tables to order',mporder

      else if(gmsh_porder.eq.mporder) then !arrays already initalized to current order, do nothing and exit!
         if(mpi_rank.eq.0) write(*,*) 'GMSH IJK Tables already initialised to order',mporder,'doing nothing! :)'
         return
      else !arrays initalized to another other
         if(mpi_rank.eq.0) write(*,*) 'GMSH IJK Tables initialised to order',gmsh_porder,'! -> changing to order',mporder

         deallocate(gmshHexahedraHO_ijkTable)
         deallocate(gmshQuadrilateralHO_ijTable)
      end if

      allocate(auxHexHOtable(mnnode,3))
      allocate(auxQuadHOtable(mnpbou,2))

      allocate(gmshHexahedraHO_ijkTable(0:mporder,0:mporder,0:mporder))
      allocate(gmshQuadrilateralHO_ijTable(0:mporder,0:mporder) )

      call genHighOrderHex(mporder,auxHexHOtable)
      call genHighOrderQuad(mporder,auxQuadHOtable)

      gmsh_porder = mporder

      !write(*,*) 'generating HexHO_ijktable'

      do pIndex=1,mnnode

         i = auxHexHOtable(pIndex,1)
         j = auxHexHOtable(pIndex,2)
         k = auxHexHOtable(pIndex,3)

         !write(*,*) 'ijk',i,j,k,'pI',pIndex

         gmshHexahedraHO_ijkTable(i,j,k) = pIndex
      end do

      !write(*,*) 'generating QuadHO_ijtable'

      do pIndex=1,mnpbou
         i = auxQuadHOtable(pIndex,1)
         j = auxQuadHOtable(pIndex,2)

         !write(*,*) 'ij',i,j,'pI',pIndex
         gmshQuadrilateralHO_ijTable(i,j) = pIndex
      end do

      deallocate(auxHexHOtable)
      deallocate(auxQuadHOtable)

   end subroutine initGmshIJKTables

   subroutine gmshHigherOrderHexahedron_pointIndexFromIJK(mporder,i,j,k,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j,k
      integer(4),intent(out) :: pointIndex

      if(gmsh_porder.ne.mporder) then !arrays not initialized
         if(mpi_rank.eq.0) write(*,*) 'ERROR! GMSH IJK TABLES NOT PROPERLY INITALISED!! gmsh_porder',gmsh_porder,'!= mporder',mporder
         pointIndex = 0
         return
      end if

      pointIndex = gmshHexahedraHO_ijkTable(i,j,k)

   end subroutine gmshHigherOrderHexahedron_pointIndexFromIJK

   subroutine gmshHigherOrderQuadrilateral_pointIndexFromIJ(mporder,i,j,pointIndex)
      implicit none
      integer(4),intent(in) :: mporder,i,j
      integer(4),intent(out) :: pointIndex

      if(gmsh_porder.ne.mporder) then !arrays not initialized
         if(mpi_rank.eq.0) write(*,*) 'ERROR! GMSH IJK TABLES NOT PROPERLY INITALISED!! gmsh_porder',gmsh_porder,'!= mporder',mporder
         pointIndex = 0
         return
      end if

      pointIndex = gmshQuadrilateralHO_ijTable(i,j)

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