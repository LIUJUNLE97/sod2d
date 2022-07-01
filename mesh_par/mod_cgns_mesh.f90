module mod_cgns_mesh
    use CGNS
    use mod_mpi_mesh
    use iso_c_binding
    implicit none
!-----------------------------------
#include "cgnstypes_f03.h"
#define int_cgns 8 
#define int_size 4   
!-----------------------------------
! ################################################################################################
! ----------------- VARS for CGNS mesh FORMAT ----------------------------------------------------
! ################################################################################################

character(500) :: cgns_pathfile_name, basename, zonename
integer(int_size) :: index_file, index_base, index_zone, index_elem, cg_err, Cx, Cy, Cz
integer(int_size) :: index_sol_V, index_sol_CC

integer(int_size), parameter :: num_max_fields = 25
integer(int_size) :: num_fields_V, num_fields_CC
character(128),dimension(num_max_fields) :: fieldNamesV
character(128),dimension(num_max_fields) :: fieldNamesCC

! ################################################################################################
! --------------------------------- END VARS  ----------------------------------------------------
! ################################################################################################

contains

! ################################################################################################
! ------------------------------ CREATE CGNS MESH FILE -------------------------------------------
! ################################################################################################

   subroutine create_CGNSmesh_par(file_path,file_name)
      implicit none     
      character(500), intent(in) :: file_path,file_name

      cgns_pathfile_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'.cgns'

      if(mpi_rank.eq.0) write(*,*) ' # Creating cgns file: ', cgns_pathfile_name

      call cgp_open_f(cgns_pathfile_name, CG_MODE_WRITE,index_file,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !### Creating base ###
      basename = 'Base'
      call create_CGNSmesh_base(basename, index_base)

      !### Creating zone ###
      zonename = 'Zone'
      call create_CGNSmesh_zone(zonename, index_zone)

      !### Adding nodes and their coordinates ###
      call add_CGNSmesh_nodesCoordinates()

      !### Adding elements and defining their connectivity ###
      call add_CGNSmesh_elements()    

      !### Creating Vertex(nodes) solution
      call create_CGNSmesh_vertexSolution()

      !### Creating Cell-Centered solution
      call create_CGNSmesh_cellCenterSolution()

      !### ---------------------------------------------------------------------------------
      call add_write_integerField_CGNSmesh_vertexSolution('gidSrl',globalIdSrl)
      call add_write_integerField_CGNSmesh_vertexSolution('gidPar',globalIdPar)
      call add_write_integerField_CGNSmesh_cellCenterSolution('elemGid',elemGid)

      call save_CGNSmesh_mpiRankVector()

      call save_CGNSmesh_connecArray()

      call save_CGNSmesh_parallel_data()


   end subroutine create_CGNSmesh_par

   subroutine create_CGNSmesh_base(basename_in, baseId)
      implicit none
      character(500), intent(in)  :: basename_in
      integer, intent(out) :: baseId

      call cg_base_write_f(index_file,basename_in,ndime,ndime,baseId,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      if(mpi_rank.eq.0) write(*,*) ' # Created CGNS base ',trim(basename_in),' baseId ',baseId
   end subroutine

   subroutine create_CGNSmesh_zone(zonename_in, zoneId)
      implicit none
      character(500), intent(in)  :: zonename_in
      integer, intent(out) :: zoneId
      integer(cgsize_t) :: isize(ndime)

      isize(1) = totalNumNodesPar
      isize(2) = totalNumElements
      isize(3) = 0
   
      call cg_zone_write_f(index_file, index_base, zonename_in, isize, CGNS_ENUMV(Unstructured), zoneId, cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      if(mpi_rank.eq.0) write(*,*) ' # Created CGNS zone ',trim(zonename_in),' zoneId ',zoneId
   end subroutine


   subroutine add_CGNSmesh_nodesCoordinates()
      implicit none
      integer :: iNodeL,iNodeGsrl

      allocate(coord_x(numNodesRankPar))
      allocate(coord_y(numNodesRankPar))
      allocate(coord_z(numNodesRankPar))

      !- create the coordinate data for this process
      do iNodeL=1,numNodesRankPar
          iNodeGSrl=globalIdSrl(iNodeL)
          !write (*,*) 'iL ', iNodeL, ' iG ', iNodeG,' [', coord(iNodeG,1),']',' [', coord(iNodeG,2),']'
          coord_x(iNodeL) = coordGMSH(iNodeGSrl,1)
          coord_y(iNodeL) = coordGMSH(iNodeGSrl,2)
          coord_z(iNodeL) = coordGMSH(iNodeGSrl,3)
      end do

      !-- create data nodes for coordinates
      call cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateX',Cx,cg_err)
      call cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateY',Cy,cg_err)
      call cgp_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateZ',Cz,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !- write the coordinate data in parallel
      call cgp_coord_write_data_f(index_file, index_base, index_zone, Cx, rankNodeStart, rankNodeEnd, coord_x, cg_err)
      call cgp_coord_write_data_f(index_file, index_base, index_zone, Cy, rankNodeStart, rankNodeEnd, coord_y, cg_err)
      call cgp_coord_write_data_f(index_file, index_base, index_zone, Cz, rankNodeStart, rankNodeEnd, coord_z, cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine add_CGNSmesh_nodesCoordinates

   subroutine add_CGNSmesh_elements()
      implicit none
      integer(cgsize_t) :: firstElem

      !---- create data node for elements
      firstElem=1
      call cgp_section_write_f(index_file, index_base, index_zone, 'Hex', CGNS_ENUMV(HEXA_64),&
                   firstElem, totalNumElements, 0, index_elem, cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !---- write the element connectivity in parallel 
      call cgp_elements_write_data_f(index_file,index_base,index_zone,index_elem,rankElemStart,rankElemEnd,connecCGNS,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine add_CGNSmesh_elements

   !------------------------------------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------------------------------------

   subroutine save_CGNSmesh_mpiRankVector()
      integer,dimension(numElemsInRank) :: mpiRankVec

      mpiRankVec(:) = mpi_rank

      call add_write_integerField_CGNSmesh_cellCenterSolution('mpiRank',mpiRankVec)

   end subroutine save_CGNSmesh_mpiRankVector

   subroutine save_CGNSmesh_connecArray()
      implicit none
      integer(cgsize_t) :: array_size,arrayStart,arrayEnd
      integer,dimension(:),allocatable :: aux_array
      integer :: i,j,iNodeL,iElemL
      integer :: index_array

      call cg_goto_f(index_file,index_base,cg_err,zonename,0,'end')
      call cg_user_data_write_f('connec_data',cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cg_gorel_f(index_file,cg_err,'connec_data',0,'end')
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      array_size = totalNumElements*nnode

      call cgp_array_write_f('connecPar',CGNS_ENUMV(Integer),1,array_size,index_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f


      allocate(aux_array(numElemsInRank*nnode))

      i=1
      do iElemL=1,numElemsInRank
         do j=1,nnode
            iNodeL = connecPar(iElemL,j)

            aux_array(i) = iNodeL
            i=i+1
         end do
      end do

      arrayStart = (rankElemStart-1)*nnode + 1
      arrayEnd = (rankElemEnd)*nnode
      !write(*,*) 'arrayStart ',arrayStart, ' arrayEnd ',arrayEnd
      !write(*,*) 'connecPar(10) ', connecPar(10,:)

      call cgp_array_write_data_f(index_array,arrayStart,arrayEnd,aux_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      deallocate(aux_array)

   end subroutine save_CGNSmesh_connecArray

   subroutine load_CGNSmesh_coordinates()
   !-- read grid coordinates
      implicit none
   
      Cx=1
      Cy=2
      Cz=3

      allocate(coord_x(numNodesRankPar))
      allocate(coord_y(numNodesRankPar))
      allocate(coord_z(numNodesRankPar))

      call cgp_coord_read_data_f(index_file,index_base,index_zone,Cx,rankNodeStart,rankNodeEnd,coord_x,cg_err)
      call cgp_coord_read_data_f(index_file,index_base,index_zone,Cy,rankNodeStart,rankNodeEnd,coord_y,cg_err)
      call cgp_coord_read_data_f(index_file,index_base,index_zone,Cz,rankNodeStart,rankNodeEnd,coord_z,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !write(*,*) '[',mpi_rank,'] coordx(10) ',coord_x(10)
   end subroutine load_CGNSmesh_coordinates

   subroutine load_CGNSmesh_connecArray()
      implicit none
      integer(cgsize_t) :: array_size,arrayStart,arrayEnd
      integer,dimension(:),allocatable :: aux_array
      integer :: i,j,iNodeL,iElemL,iArray

      call cg_goto_f(index_file,index_base,cg_err,zonename,0,'connec_data',0,'end')

      arrayStart = mpi_rank+1
      arrayEnd   = mpi_rank+1

      allocate( aux_array(numElemsInRank*nnode) )
      allocate( connecPar(numElemsInRank,nnode) )

      iArray = 1
      arrayStart = (rankElemStart-1)*nnode + 1
      arrayEnd = (rankElemEnd)*nnode
      
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,aux_array,cg_err)

      i=1
      do iElemL=1,numElemsInRank
         do j=1,nnode
            iNodeL = aux_array(i)
            connecPar(iElemL,j)= iNodeL
            i=i+1
         end do
      end do

      deallocate(aux_array)

   end subroutine load_CGNSmesh_connecArray

   subroutine save_CGNSmesh_parallel_data()
      implicit none
      integer(cgsize_t) :: array_size,arrayStart,arrayEnd
      integer,dimension(:),allocatable :: aux_array
      integer :: i,iArray, accumVal

      call cg_goto_f(index_file,index_base,cg_err,zonename,0,'end')
      call cg_user_data_write_f('Parallel_data',cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cg_gorel_f(index_file,cg_err,'Parallel_data',0,'end')
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !write num of CPUs

      !write rankNodeStart,rankNodeEnd
      !--------------------------------------------------------------------------------------------------
      array_size = mpi_size
      arrayStart = mpi_rank+1
      arrayEnd   = mpi_rank+1

      !--- rankNodeStart -> iArray = 1
      call cgp_array_write_f('rankNodeStart',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankNodeStart,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- rankNodeEnd  -> iArray = 2
      call cgp_array_write_f('rankNodeEnd',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankNodeEnd,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- rankElemStart -> iArray = 3
      call cgp_array_write_f('rankElemStart',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankElemStart,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- rankElemEnd -> iArray = 4
      call cgp_array_write_f('rankElemEnd',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankElemEnd,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- numRanksWithComms -> iArray = 5
      call cgp_array_write_f('numRanksWithComms',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,numRanksWithComms,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !---- READING THE VALUES FROM ALL THE PROCS TO KNOW THE TOTAL VECTOR PLACES REQUIRED AND STORE ASSOCIATED STUFF
      arrayStart = 1
      arrayEnd   = mpi_size
      allocate( aux_array(mpi_size) )
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,aux_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      
      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      array_size = accumVal

      arrayStart=1
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         arrayStart=arrayStart+aux_array(i)
      end do
      arrayEnd=arrayStart+aux_array(mpi_rank+1)-1
      !write(*,*) '[',mpi_rank,']nrwc aux: ',aux_array(:),' accum ',accumVal,' arrayS ',arrayStart,' arrayE ',arrayEnd

      !--- ranksToComm 
      call cgp_array_write_f('ranksToComm',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,ranksToComm,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- commsMemPosInLoc
      call cgp_array_write_f('commsMemPosInLoc',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,commsMemPosInLoc,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      
      !--- commsMemSize
      call cgp_array_write_f('commsMemSize',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,commsMemSize,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- commsMemPosInNgb
      call cgp_array_write_f('commsMemPosInNgb',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,commsMemPosInNgb,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !---------------------------------------------------------------------------------------------------
      !--- numNodesToComm -> iArray = X
      arrayStart = mpi_rank+1
      arrayEnd   = mpi_rank+1
      array_size = mpi_size
      call cgp_array_write_f('numNodesToComm',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,numNodesToComm,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !---- READING THE VALUES FROM ALL THE PROCS TO KNOW THE TOTAL VECTOR PLACES REQUIRED AND STORE ASSOCIATED STUFF
      arrayStart = 1
      arrayEnd   = mpi_size
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,aux_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      array_size = accumVal

      arrayStart=1
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         arrayStart=arrayStart+aux_array(i)
      end do
      arrayEnd=arrayStart+aux_array(mpi_rank+1)-1
      !write(*,*) '[',mpi_rank,']nrwc aux: ',aux_array(:),' accum ',accumVal,' arrayS ',arrayStart,' arrayE ',arrayEnd

      !--- matrixCommScheme(:,1) iNodeL
      call cgp_array_write_f('matrixCommScheme_iNodeL',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,matrixCommScheme(:,1),cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- matrixCommScheme(:,2) iNodeGSrl
      call cgp_array_write_f('matrixCommScheme_iNodeGsrl',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,matrixCommScheme(:,2),cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- matrixCommScheme(:,3) ngbRank
      call cgp_array_write_f('matrixCommScheme_ngbRank',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,matrixCommScheme(:,3),cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f


      deallocate(aux_array)

   end subroutine save_CGNSmesh_parallel_data

   subroutine load_CGNSmesh_parallel_data()
      implicit none
      integer(cgsize_t) :: array_size,arrayStart,arrayEnd
      integer,dimension(:),allocatable :: aux_array
      integer :: i,iArray, accumVal

      !---- READING ALL THE DATA IN NODE 'Parallel_data'---------------------------------------------
      call cg_goto_f(index_file,index_base,cg_err,zonename,0,'Parallel_data',0,'end')
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      arrayStart = mpi_rank+1
      arrayEnd   = mpi_rank+1

      !rankNodeStart
      iArray=1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankNodeStart,cg_err)

      !rankNodeEnd
      iArray=2!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankNodeEnd,cg_err)

      !TO BE MOD!
      !rankElemStart
      iArray=3!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankElemStart,cg_err)

      !TO BE MOD
      !rankElemEnd
      iArray=4!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankElemEnd,cg_err)

      numNodesRankPar = rankNodeEnd - rankNodeStart + 1
      numElemsInRank  = rankElemEnd - rankElemStart + 1

      !write(*,*) '#rank ', mpi_rank, ' nodesInRank ', numNodesRankPar, ' iNodeS ', rankNodeStart, ' iNodeE ', rankNodeEnd
      !write(*,*) '#rank ', mpi_rank, ' elemsInRank ', numElemsInRank, ' iElemS ', rankElemStart, ' iElemE ', rankElemEnd

      ! numRanksWithComms -> iArray = 5
      iArray=5!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,numRanksWithComms,cg_err)

      !*now reading the values from all the procs to know the total vector places required and load associated stuff
      arrayStart = 1
      arrayEnd   = mpi_size
      allocate( aux_array(mpi_size) )
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,aux_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      
      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      array_size = accumVal

      arrayStart=1
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         arrayStart=arrayStart+aux_array(i)
      end do
      arrayEnd=arrayStart+aux_array(mpi_rank+1)-1

      allocate(ranksToComm(numRanksWithComms))
      allocate(commsMemPosInLoc(numRanksWithComms))
      allocate(commsMemPosInNgb(numRanksWithComms))
      allocate(commsMemSize(numRanksWithComms))

      !--- ranksToComm 
      iArray=6!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,ranksToComm,cg_err)

      !--- commsMemPosInLoc
      iArray=7!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,commsMemPosInLoc,cg_err)
      
      !--- commsMemSize
      iArray=8!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,commsMemSize,cg_err)

      !--- commsMemPosInNgb
      iArray=9!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,commsMemPosInNgb,cg_err)
      
      !write(*,*) '[',mpi_rank,']commMPIL ',commsMemPosInLoc(:)
      !write(*,*) '[',mpi_rank,']commMPIN ',commsMemPosInNgb(:)
      !write(*,*) '[',mpi_rank,']commMS ',commsMemSize(:)

      !---------------------------------------------------------------------------------------------------
      !--- numNodesToComm
      iArray=10!iArray+1
      arrayStart = mpi_rank+1
      arrayEnd   = mpi_rank+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,numNodesToComm,cg_err)

      !*now reading the values from all the procs to know the total vector places required and load associated stuff
      arrayStart = 1
      arrayEnd   = mpi_size
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,aux_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      accumVal=0
      do i=1,mpi_size
         accumVal=accumVal+aux_array(i)
      end do

      array_size = accumVal

      arrayStart=1
      do i=1,(mpi_rank) !from rank 0 mpi_rank-1
         arrayStart=arrayStart+aux_array(i)
      end do
      arrayEnd=arrayStart+aux_array(mpi_rank+1)-1
      !write(*,*) '[',mpi_rank,']nrwc aux: ',aux_array(:),' accum ',accumVal,' arrayS ',arrayStart,' arrayE ',arrayEnd

      allocate(matrixCommScheme(numNodesToComm,3))

      !--- ranksToComm 
      iArray=11!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,matrixCommScheme(:,1),cg_err)

      !--- commsMemPosInLoc
      iArray=12!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,matrixCommScheme(:,2),cg_err)
      
      !--- commsMemSize
      iArray=13!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,matrixCommScheme(:,3),cg_err)

      deallocate(aux_array)

   end subroutine load_CGNSmesh_parallel_data

!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------

   subroutine create_CGNSmesh_vertexSolution()
      implicit none
      !---- create a vertex solution
      call cg_sol_write_f(index_file,index_base,index_zone,'VertexSol',CGNS_ENUMV(Vertex),index_sol_V,cg_err)
      iF (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine create_CGNSmesh_vertexSolution

   subroutine create_CGNSmesh_cellCenterSolution()
      implicit none
      !---- create a cell-center solution
      call cg_sol_write_f(index_file,index_base,index_zone,'CellCenterSol',CGNS_ENUMV(CellCenter),index_sol_CC,cg_err)
      iF (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine create_CGNSmesh_cellCenterSolution

   integer function getFieldId_V(fieldName) result(fieldId)
      implicit none
      character(len=*), intent(in) :: fieldName
      integer :: i
      fieldId = -1
      do i=1,num_fields_V
         if(fieldNamesV(i).eq.fieldName) then
            fieldId = i
            exit
         end if
      end do
   end function getFieldId_V

   integer function getFieldId_CC(fieldName) result(fieldId)
      implicit none
      character(len=*), intent(in) :: fieldName
      integer :: i
      fieldId = -1
      do i=1,num_fields_CC
         if(fieldNamesCC(i).eq.fieldName) then
            fieldId = i
            exit
         end if
      end do
   end function getFieldId_CC

!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
   !### INTEGER | VERTEX  ###
   !-------------------------
   subroutine add_integerField_CGNSmesh_vertexSolution(fieldName_in,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(out), optional :: fieldId_out
      integer(int_size) :: fieldId

      !---- create integer field in the nodes/vertex
      call cgp_field_write_f(index_file,index_base,index_zone,index_sol_V,CGNS_ENUMV(Integer),fieldName_in,fieldId,cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f

      fieldNamesV(fieldId) = trim(fieldName_in)
      num_fields_V = fieldId
      if(mpi_rank.eq.0) write(*,*) 'Added Integer Field (V) ',fieldName_in,' with id ',fieldId
      if(present(fieldId_out)) fieldId_out = fieldId

   end subroutine add_integerField_CGNSmesh_vertexSolution

   subroutine write_integerField_CGNSmesh_vertexSolution(fieldId,data_field)
      implicit none
      integer(int_size), intent(in) :: fieldId
      integer(int_size), intent(in) :: data_field(numNodesRankPar)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_V,fieldId,&
                                  rankNodeStart, rankNodeEnd, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_integerField_CGNSmesh_vertexSolution

   subroutine add_write_integerField_CGNSmesh_vertexSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(in) :: data_field(numNodesRankPar)
      integer(int_size), intent(out), optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_integerField_CGNSmesh_vertexSolution(fieldName_in,fieldId)
      call write_integerField_CGNSmesh_vertexSolution(fieldId,data_field)
   end subroutine add_write_integerField_CGNSmesh_vertexSolution

   subroutine load_integerField_CGNSmesh_vertexSolution(fieldId,data_field)
      integer(int_size), intent(in) :: fieldId
      integer(int_size), intent(out) :: data_field(numNodesRankPar)

      call cgp_field_read_data_f(index_file,index_base,index_zone,index_sol_V,fieldId,&
                                 rankNodeStart,rankNodeEnd,data_field,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine load_integerField_CGNSmesh_vertexSolution

   subroutine load_integerField_CGNSmesh_vertexSolution_fieldName(fieldName,data_field)
      character(len=*), intent(in) :: fieldName
      integer(int_size), intent(out) :: data_field(numNodesRankPar)
      integer(int_size) :: fieldId

      fieldId = getFieldId_V(fieldName)
      if(fieldId.ne.-1) then
         call load_integerField_CGNSmesh_vertexSolution(fieldId,data_field)
      end if
   end subroutine load_integerField_CGNSmesh_vertexSolution_fieldName

!------------------------------------------------------------------------------------------------------------------------------
   !### DOUBLE | VERTEX  ###
   !------------------------
   subroutine add_doubleField_CGNSmesh_vertexSolution(fieldName_in,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      !---- create integer field in the nodes/vertex
      call cgp_field_write_f(index_file,index_base,index_zone,index_sol_V,CGNS_ENUMV(RealDouble),fieldName_in,fieldId,cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f

      fieldNamesV(fieldId) = trim(fieldName_in)
      num_fields_V = fieldId
      if(mpi_rank.eq.0) write(*,*) 'Added Double Field (V) ',fieldName_in,' with id ',fieldId
      if(present(fieldId_out)) fieldId_out = fieldId

   end subroutine add_doubleField_CGNSmesh_vertexSolution

   subroutine write_doubleField_CGNSmesh_vertexSolution(fieldId,data_field)
      implicit none
      integer(int_size), intent(in) :: fieldId
      real(8), intent(in) :: data_field(numNodesRankPar)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_V,fieldId,&
                                  rankNodeStart, rankNodeEnd, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_doubleField_CGNSmesh_vertexSolution

   subroutine add_write_doubleField_CGNSmesh_vertexSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      real(8), intent(in) :: data_field(numNodesRankPar)
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_doubleField_CGNSmesh_vertexSolution(fieldName_in,fieldId)
      call write_doubleField_CGNSmesh_vertexSolution(fieldId,data_field)

      if(present(fieldId_out)) fieldId_out = fieldId
   end subroutine add_write_doubleField_CGNSmesh_vertexSolution
!------------------------------------------------------------------------------------------------------------------------------
   !### INTEGER | CELL-CENTER  ###
   !------------------------------
   subroutine add_integerField_CGNSmesh_cellCenterSolution(fieldName_in,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      !---- create integer field in the nodes/vertex
      call cgp_field_write_f(index_file,index_base,index_zone,index_sol_CC,CGNS_ENUMV(Integer),fieldName_in,fieldId,cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f

      fieldNamesCC(fieldId) = trim(fieldName_in)
      num_fields_CC = fieldId
      if(mpi_rank.eq.0) write(*,*) 'Added Integer Field (CC) ',fieldName_in,' with id ',fieldId
      if(present(fieldId_out)) fieldId_out = fieldId

   end subroutine add_integerField_CGNSmesh_cellCenterSolution

   subroutine write_integerField_CGNSmesh_cellCenterSolution(fieldId,data_field)
      implicit none
      integer(int_size), intent(in) :: fieldId
      integer(int_size), intent(in) :: data_field(numElemsInRank)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_CC,fieldId,&
                                  rankElemStart, rankElemEnd, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_integerField_CGNSmesh_cellCenterSolution

    subroutine add_write_integerField_CGNSmesh_cellCenterSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(in) :: data_field(numElemsInRank)
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_integerField_CGNSmesh_cellCenterSolution(fieldName_in,fieldId)
      call write_integerField_CGNSmesh_cellCenterSolution(fieldId,data_field)

      if(present(fieldId_out)) fieldId_out = fieldId
   end subroutine add_write_integerField_CGNSmesh_cellCenterSolution

   subroutine load_integerField_CGNSmesh_cellCenterSolution(fieldId,data_field)
      integer(int_size), intent(in) :: fieldId
      integer(int_size), intent(out) :: data_field(numElemsInRank)

      call cgp_field_read_data_f(index_file,index_base,index_zone,index_sol_CC,fieldId,&
                                 rankElemStart,rankElemEnd,data_field,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine load_integerField_CGNSmesh_cellCenterSolution

   subroutine load_integerField_CGNSmesh_cellCenterSolution_fieldName(fieldName,data_field)
      character(len=*), intent(in) :: fieldName
      integer(int_size), intent(out) :: data_field(numElemsInRank)
      integer(int_size) :: fieldId

      fieldId = getFieldId_CC(fieldName)
      if(fieldId.ne.-1) then
         call load_integerField_CGNSmesh_cellCenterSolution(fieldId,data_field)
      end if
   end subroutine load_integerField_CGNSmesh_cellCenterSolution_fieldName

!------------------------------------------------------------------------------------------------------------------------------
   !### DOUBLE | CELL-CENTER  ###
   !------------------------------
   subroutine add_doubleField_CGNSmesh_cellCenterSolution(fieldName_in,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      !---- create integer field in the nodes/vertex
      call cgp_field_write_f(index_file,index_base,index_zone,index_sol_CC,CGNS_ENUMV(RealDouble),fieldName_in,fieldId,cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f

      fieldNamesCC(fieldId) = trim(fieldName_in)
      num_fields_CC = fieldId
      if(mpi_rank.eq.0) write(*,*) 'Added Double Field (CC) ',fieldName_in,' with id ',fieldId
      if(present(fieldId_out)) fieldId_out = fieldId

   end subroutine add_doubleField_CGNSmesh_cellCenterSolution

   subroutine write_doubleField_CGNSmesh_cellCenterSolution(fieldId,data_field)
      implicit none
      integer(int_size), intent(in) :: fieldId
      real(8), intent(in) :: data_field(numElemsInRank)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_CC,fieldId,&
                                  rankElemStart, rankElemEnd, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_doubleField_CGNSmesh_cellCenterSolution

    subroutine add_write_doubleField_CGNSmesh_cellCenterSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      real(8), intent(in) :: data_field(numElemsInRank)
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_doubleField_CGNSmesh_cellCenterSolution(fieldName_in,fieldId)
      call write_doubleField_CGNSmesh_cellCenterSolution(fieldId,data_field)

      if(present(fieldId_out)) fieldId_out = fieldId
   end subroutine add_write_doubleField_CGNSmesh_cellCenterSolution

! ################################################################################################
! -------------------------------- LOAD CGNS MESH FILE -------------------------------------------
! ################################################################################################

subroutine open_CGNSmesh_par(file_path,file_name)
   implicit none     
   character(500), intent(in)  :: file_path,file_name
   integer(cgsize_t) :: isize(ndime)
   integer(int_size) :: i,iArray,numOfArrays,fieldId
   character(len=32) :: field_name
   integer(cgenum_t) :: bocotype, datatype
   integer(int_size) :: nsols

   cgns_pathfile_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'.cgns'

   if(mpi_rank.eq.0) write(*,*) '# Opening cgns file: ', cgns_pathfile_name

   call cgp_open_f(cgns_pathfile_name,CG_MODE_MODIFY,index_file,cg_err)
   if (cg_err .ne. CG_OK) call cgp_error_exit_f

   index_base = 1
   index_zone = 1

   call cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,cg_err)
   if (cg_err .ne. CG_OK) call cgp_error_exit_f

   if(mpi_rank.eq.0) write(*,*) '# Readed CGNS zone ',trim(zonename),' zoneId ',index_zone

   totalNumNodesPar = isize(1) 
   totalNumElements = isize(2) 

   if(mpi_rank.eq.0) write(*,*) '  - totalNumNodesPar ',totalNumNodesPar,' totalNumElements ',totalNumElements

   call load_CGNSmesh_parallel_data()

   call load_CGNSmesh_coordinates()

   call load_CGNSmesh_connecArray()

   call cg_nsols_f(index_file,index_base,index_zone,nsols,cg_err)
   
   !checking that there are presents the two (2) FlowSolution_t nodes
   if(nsols.ne.2) call cgp_error_exit_f()

   index_sol_V = 1
   index_sol_CC = 2

   call cg_nfields_f(index_file,index_base,index_zone, index_sol_V, num_fields_V , cg_err)
   call cg_nfields_f(index_file,index_base,index_zone, index_sol_CC, num_fields_CC, cg_err)

   if(mpi_rank.eq.0) write(*,*) '# Loading the following number of Fields (V) ',num_fields_V
   fieldNamesV(:) = ""
   do i=1,num_fields_V
      call cg_field_info_f(index_file,index_base,index_zone,index_sol_V,i,datatype,field_name,cg_err)
      if(mpi_rank.eq.0) write(*,*) '  - field V ',i, ' name ',field_name, ' datatype ',datatype
      fieldNamesV(i) = trim(field_name)
   end do

   if(mpi_rank.eq.0) write(*,*) '# Loading the following number of Fields (CC) ',num_fields_CC
   fieldNamesCC(:) = ""
   do i=1,num_fields_CC
      call cg_field_info_f(index_file,index_base,index_zone,index_sol_CC,i,datatype,field_name,cg_err)
      if(mpi_rank.eq.0) write(*,*) '  - field CC ',i, ' name ',field_name, ' datatype ',datatype
      fieldNamesCC(i) = trim(field_name)
   end do

   allocate(globalIdSrl(numNodesRankPar))
   allocate(globalIdPar(numNodesRankPar))
   allocate(elemGid(numElemsInRank))

   call load_integerField_CGNSmesh_vertexSolution_fieldName('gidSrl',globalIdSrl)
   call load_integerField_CGNSmesh_vertexSolution_fieldName('gidPar',globalIdPar)

   call load_integerField_CGNSmesh_cellCenterSolution_fieldName('elemGid',elemGid)

   !i=10
   !write(*,*) '[',mpi_rank,']gidSrl(i) ',globalIdSrl(i),' gidPar(i) ',globalIdPar(i),' elemGid ',elemGid(i)

end subroutine open_CGNSmesh_par

subroutine close_CGNSmesh_par()
   implicit none

   call cgp_close_f(index_file,cg_err)

end subroutine close_CGNSmesh_par

end module mod_cgns_mesh
