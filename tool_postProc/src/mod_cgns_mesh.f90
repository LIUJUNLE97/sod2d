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

! haig de copiar les variables originals en int_size 4 a int_cgns 
integer(int_cgns) :: totalNumElements_int8
integer(int_cgns) :: numElemsInRank_int8, rankElemStart_int8, rankElemEnd_int8
integer(int_cgns) :: rankNodeStart_int8, rankNodeEnd_int8
integer(int_cgns), allocatable :: connecCGNS_int8(:)
!---------------------------------------------------------------------------------
real(8), allocatable :: coordParCGNS(:,:)

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

      !--------------------------------------------------------------------------
      !copying integers in required cgns format
      totalNumElements_int8 = totalNumElements
      numElemsInRank_int8   = numElemsInRank
      rankElemStart_int8    = rankElemStart
      rankElemEnd_int8      = rankElemEnd
      rankNodeStart_int8    = rankNodeStart
      rankNodeEnd_int8      = rankNodeEnd

      allocate(connecCGNS_int8(numElemsInRank_int8*nnode) )
      connecCGNS_int8(:) = connecCGNS(:)
      !-------------------------------------------------------------------------

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

#if 0
      call save_CGNSmesh_mpiRankVector()

      call save_CGNSmesh_connecArray()

      if(mpi_size.gt.1) then
         call save_CGNSmesh_parallel_data()
      end if
#endif

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
      isize(2) = totalNumElements_int8
      isize(3) = 0
   
      call cg_zone_write_f(index_file, index_base, zonename_in, isize, CGNS_ENUMV(Unstructured), zoneId, cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      if(mpi_rank.eq.0) write(*,*) ' # Created CGNS zone ',trim(zonename_in),' zoneId ',zoneId
   end subroutine


   subroutine add_CGNSmesh_nodesCoordinates()
      implicit none
      !integer :: iNodeL,iNodeGsrl

      !allocate(coord_x(numNodesRankPar))
      !allocate(coord_y(numNodesRankPar))
      !allocate(coord_z(numNodesRankPar))

      !!- create the coordinate data for this process
      !do iNodeL=1,numNodesRankPar
      !    iNodeGSrl=globalIdSrl(iNodeL)
      !    !write (*,*) 'iL ', iNodeL, ' iG ', iNodeG,' [', coord(iNodeG,1),']',' [', coord(iNodeG,2),']'
      !    coord_x(iNodeL) = coordGMSH(iNodeGSrl,1)
      !    coord_y(iNodeL) = coordGMSH(iNodeGSrl,2)
      !    coord_z(iNodeL) = coordGMSH(iNodeGSrl,3)
      !end do

      !-- create data nodes for coordinates
      call cgp_coord_write_f(index_file,index_base,index_zone,RealSingle,'CoordinateX',Cx,cg_err)
      call cgp_coord_write_f(index_file,index_base,index_zone,RealSingle,'CoordinateY',Cy,cg_err)
      call cgp_coord_write_f(index_file,index_base,index_zone,RealSingle,'CoordinateZ',Cz,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !- write the coordinate data in parallel
      call cgp_coord_write_data_f(index_file,index_base,index_zone,Cx,rankNodeStart_int8,rankNodeEnd_int8,coordPar(:,1),cg_err)
      call cgp_coord_write_data_f(index_file,index_base,index_zone,Cy,rankNodeStart_int8,rankNodeEnd_int8,coordPar(:,2),cg_err)
      call cgp_coord_write_data_f(index_file,index_base,index_zone,Cz,rankNodeStart_int8,rankNodeEnd_int8,coordPar(:,3),cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine add_CGNSmesh_nodesCoordinates

   subroutine add_CGNSmesh_elements()
      implicit none
      integer(cgsize_t) :: firstElem

      !---- create data node for elements
      firstElem=1
      call cgp_section_write_f(index_file, index_base, index_zone, 'Hex', CGNS_ENUMV(HEXA_64),&
                  firstElem, totalNumElements_int8, 0, index_elem, cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !---- write the element connectivity in parallel 
      call cgp_elements_write_data_f(index_file,index_base,index_zone,index_elem,&
                  rankElemStart_int8,rankElemEnd_int8,connecCGNS_int8,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine add_CGNSmesh_elements

   !------------------------------------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------------------------------------
#if 0
   subroutine save_CGNSmesh_mpiRankVector()
      integer,dimension(numElemsInRank_int8) :: mpiRankVec

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

      array_size = totalNumElements_int8*nnode

      call cgp_array_write_f('connecPar',CGNS_ENUMV(Integer),1,array_size,index_array,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f


      allocate(aux_array(numElemsInRank_int8*nnode))

      i=1
      do iElemL=1,numElemsInRank_int8
         do j=1,nnode
            iNodeL = connecPar(iElemL,j)

            aux_array(i) = iNodeL
            i=i+1
         end do
      end do

      arrayStart = (rankElemStart_int8-1)*nnode + 1
      arrayEnd = (rankElemEnd_int8)*nnode
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

      allocate(coordParCGNS(numNodesRankPar,3))

      call cgp_coord_read_data_f(index_file,index_base,index_zone,Cx,rankNodeStart_int8,rankNodeEnd_int8,coordParCGNS(:,1),cg_err)
      call cgp_coord_read_data_f(index_file,index_base,index_zone,Cy,rankNodeStart_int8,rankNodeEnd_int8,coordParCGNS(:,2),cg_err)
      call cgp_coord_read_data_f(index_file,index_base,index_zone,Cz,rankNodeStart_int8,rankNodeEnd_int8,coordParCGNS(:,3),cg_err)
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

      allocate( aux_array(numElemsInRank_int8*nnode) )
      allocate( connecPar(numElemsInRank_int8,nnode) )

      iArray = 1
      arrayStart = (rankElemStart_int8-1)*nnode + 1
      arrayEnd = (rankElemEnd_int8)*nnode
      
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,aux_array,cg_err)

      i=1
      do iElemL=1,numElemsInRank_int8
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

      !write rankNodeStart_int8,rankNodeEnd_int8
      !--------------------------------------------------------------------------------------------------
      array_size = mpi_size
      arrayStart = mpi_rank+1
      arrayEnd   = mpi_rank+1

      !--- rankNodeStart_int8 -> iArray = 1
      call cgp_array_write_f('rankNodeStart_int8',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankNodeStart_int8,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- rankNodeEnd_int8  -> iArray = 2
      call cgp_array_write_f('rankNodeEnd_int8',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankNodeEnd_int8,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- rankElemStart_int8 -> iArray = 3
      call cgp_array_write_f('rankElemStart_int8',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankElemStart_int8,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f

      !--- rankElemEnd_int8 -> iArray = 4
      call cgp_array_write_f('rankElemEnd_int8',CGNS_ENUMV(Integer),1,array_size,iArray,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
      call cgp_array_write_data_f(iArray,arrayStart,arrayEnd,rankElemEnd_int8,cg_err)
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

      !rankNodeStart_int8
      iArray=1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankNodeStart_int8,cg_err)

      !rankNodeEnd_int8
      iArray=2!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankNodeEnd_int8,cg_err)

      !TO BE MOD!
      !rankElemStart_int8
      iArray=3!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankElemStart_int8,cg_err)

      !TO BE MOD
      !rankElemEnd_int8
      iArray=4!iArray+1
      call cgp_array_read_data_f(iArray,arrayStart,arrayEnd,rankElemEnd_int8,cg_err)

      numNodesRankPar = rankNodeEnd_int8 - rankNodeStart_int8 + 1
      numElemsInRank_int8  = rankElemEnd_int8 - rankElemStart_int8 + 1

      !write(*,*) '#rank ', mpi_rank, ' nodesInRank ', numNodesRankPar, ' iNodeS ', rankNodeStart_int8, ' iNodeE ', rankNodeEnd_int8
      !write(*,*) '#rank ', mpi_rank, ' elemsInRank ', numElemsInRank_int8, ' iElemS ', rankElemStart_int8, ' iElemE ', rankElemEnd_int8

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
   totalNumElements_int8 = isize(2) 

   if(mpi_rank.eq.0) write(*,*) '  - totalNumNodesPar ',totalNumNodesPar,' totalNumElements_int8 ',totalNumElements_int8
   
   if(mpi_size.gt.1) then
      call load_CGNSmesh_parallel_data()
   end if

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
   allocate(elemGid(numElemsInRank_int8))

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

#endif
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
                                  rankNodeStart_int8, rankNodeEnd_int8, data_field, cg_err)
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
                                 rankNodeStart_int8,rankNodeEnd_int8,data_field,cg_err)
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
   !### FLOAT | VERTEX  ###
   !------------------------
   subroutine add_floatField_CGNSmesh_vertexSolution(fieldName_in,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      !---- create integer field in the nodes/vertex
      call cgp_field_write_f(index_file,index_base,index_zone,index_sol_V,CGNS_ENUMV(RealSingle),fieldName_in,fieldId,cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f

      fieldNamesV(fieldId) = trim(fieldName_in)
      num_fields_V = fieldId
      if(mpi_rank.eq.0) write(*,*) 'Added float Field (V) ',fieldName_in,' with id ',fieldId
      if(present(fieldId_out)) fieldId_out = fieldId

   end subroutine add_floatField_CGNSmesh_vertexSolution

   subroutine write_floatField_CGNSmesh_vertexSolution(fieldId,data_field)
      implicit none
      integer(int_size), intent(in) :: fieldId
      real(4), intent(in) :: data_field(numNodesRankPar)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_V,fieldId,&
                                  rankNodeStart_int8, rankNodeEnd_int8, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_floatField_CGNSmesh_vertexSolution

   subroutine add_write_floatField_CGNSmesh_vertexSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      real(4), intent(in) :: data_field(numNodesRankPar)
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_floatField_CGNSmesh_vertexSolution(fieldName_in,fieldId)
      call write_floatField_CGNSmesh_vertexSolution(fieldId,data_field)

      if(present(fieldId_out)) fieldId_out = fieldId
   end subroutine add_write_floatField_CGNSmesh_vertexSolution
!------------------------------------------------------------------------------------------------------------------------------

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
                                  rankNodeStart_int8, rankNodeEnd_int8, data_field, cg_err)
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
      integer(int_size), intent(in) :: data_field(numElemsInRank_int8)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_CC,fieldId,&
                                  rankElemStart_int8, rankElemEnd_int8, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_integerField_CGNSmesh_cellCenterSolution

    subroutine add_write_integerField_CGNSmesh_cellCenterSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      integer(int_size), intent(in) :: data_field(numElemsInRank_int8)
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_integerField_CGNSmesh_cellCenterSolution(fieldName_in,fieldId)
      call write_integerField_CGNSmesh_cellCenterSolution(fieldId,data_field)

      if(present(fieldId_out)) fieldId_out = fieldId
   end subroutine add_write_integerField_CGNSmesh_cellCenterSolution

   subroutine load_integerField_CGNSmesh_cellCenterSolution(fieldId,data_field)
      integer(int_size), intent(in) :: fieldId
      integer(int_size), intent(out) :: data_field(numElemsInRank_int8)

      call cgp_field_read_data_f(index_file,index_base,index_zone,index_sol_CC,fieldId,&
                                 rankElemStart_int8,rankElemEnd_int8,data_field,cg_err)
      if (cg_err .ne. CG_OK) call cgp_error_exit_f
   end subroutine load_integerField_CGNSmesh_cellCenterSolution

   subroutine load_integerField_CGNSmesh_cellCenterSolution_fieldName(fieldName,data_field)
      character(len=*), intent(in) :: fieldName
      integer(int_size), intent(out) :: data_field(numElemsInRank_int8)
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
      real(8), intent(in) :: data_field(numElemsInRank_int8)

      call cgp_field_write_data_f(index_file,index_base,index_zone,index_sol_CC,fieldId,&
                                  rankElemStart_int8, rankElemEnd_int8, data_field, cg_err)
      if (cg_err .NE. CG_OK) call cgp_error_exit_f
   end subroutine write_doubleField_CGNSmesh_cellCenterSolution

    subroutine add_write_doubleField_CGNSmesh_cellCenterSolution(fieldName_in,data_field,fieldId_out)
      implicit none
      character(len=*), intent(in) :: fieldName_in
      real(8), intent(in) :: data_field(numElemsInRank_int8)
      integer(int_size), intent(out),optional :: fieldId_out
      integer(int_size) :: fieldId

      call add_doubleField_CGNSmesh_cellCenterSolution(fieldName_in,fieldId)
      call write_doubleField_CGNSmesh_cellCenterSolution(fieldId,data_field)

      if(present(fieldId_out)) fieldId_out = fieldId
   end subroutine add_write_doubleField_CGNSmesh_cellCenterSolution



end module mod_cgns_mesh
