module mod_mpi_mesh
    use mod_mpi
    use mod_utils
    use iso_c_binding
    implicit none
!-----------------------------------   
#define _CHECK_ 0
#define int_cgns 8 
#define int_size 4
!-----------------------------------   
   interface
      subroutine gempa_do_partition(numElemsInRank,numRanksToPart,x,y,z,part) bind(c)
        import c_int, c_double
        integer(kind=c_int), intent(in), value :: numElemsInRank
        integer(kind=c_int), intent(in), value :: numRanksToPart
        real(kind=c_double), intent(in) :: x(numElemsInRank)
        real(kind=c_double), intent(in) :: y(numElemsInRank)
        real(kind=c_double), intent(in) :: z(numElemsInRank)
        integer(kind=c_int), intent(out) :: part(numElemsInRank)
      end subroutine gempa_do_partition
   end interface
!-----------------------------------   
   ! Dimensions -----------------------------
   integer(int_size), parameter :: ndime=3
   !-----------------------------------------
   ! Element characteristics ----------------
   integer(int_size), parameter :: nnode=64
   integer(int_size), parameter :: porder=3
   integer(int_size), parameter :: npbou=16
   integer(int_size), parameter :: ngaus=64
   !-----------------------------------------
   integer(int_size),parameter :: tV=4,tE=3,tF=2,tI=1
   integer(int_size),parameter :: tn2ijk(nnode) = [tV,tV,tE,tE,tV,tV,tE,tE,tE,tE,tF,tF,tE,tE,tF,tF,&
                                   tV,tV,tE,tE,tV,tV,tE,tE,tE,tE,tF,tF,tE,tE,tF,tF,&
                                   tE,tE,tF,tF,tE,tE,tF,tF,tF,tF,tI,tI,tF,tF,tI,tI,&
                                   tE,tE,tF,tF,tE,tE,tF,tF,tF,tF,tI,tI,tF,tF,tI,tI]  

   integer(int_size),parameter :: gmsh2ijk(nnode) = [1,4,11,12,2,3,15,16,9,20,33,34,10,19,36,35,&
               5,8,27,28,6,7,29,30,25,32,53,56,26,31,54,55,&
               13,23,41,44,17,21,45,46,37,50,57,60,38,49,58,59,&
               14,24,42,43,18,22,48,47,40,51,61,64,39,52,62,63]

   integer(int_size),parameter :: cgns2ijk(nnode)= [1,4,16,15,2,3,11,12,9,14,33,36,10,13,34,35,&
                 5,8,32,31,6,7,27,28,25,30,53,56,26,29,54,55,&
                 17,21,50,49,19,23,41,42,37,46,57,60,38,45,58,59,&
                 18,22,51,52,20,24,44,43,40,47,61,64,39,48,62,63]
 
!  according to cgns documentation....-----------------------------------------
!  integer(4),parameter :: cgns2ijk(nnode)= [1,4,16,15,2,3,11,12,9,14,33,36,10,13,34,35,&
!                5,8,32,31,6,7,27,28,25,30,53,56,26,29,54,55,&
!                17,23,50,49,19,21,41,42,37,46,57,60,38,45,58,59,&
!                18,24,51,52,20,22,44,43,40,47,61,64,39,48,62,63]

   integer(int_size),parameter :: dummy2ijk(nnode)= [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,&
                 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,&
                 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,&
                 49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64]


! ################################################################################################
! ----------------- VARS for alya/gmsh file reading ----------------------------------------------
! ################################################################################################
!------  -------------------------
integer(int_size), allocatable :: connecGMSH(:,:), boundGMSH(:,:)
real(8), allocatable :: coordGMSH(:,:)

! ################################################################################################
! ----------------- VARS for new Par mesh FORMAT -------------------------------------------------
! ################################################################################################

integer(int_cgns) :: numNodesRankPar, totalNumNodesPar, totalNumNodesSrl
integer(int_cgns) :: totalNumElements
integer(int_cgns) :: numElemsInRank, rankElemStart, rankElemEnd
integer(int_cgns) :: rankNodeStart, rankNodeEnd

integer(int_size), allocatable :: globalIdSrl(:), globalIdPar(:), elemGid(:)

real(8), allocatable :: coord_x(:), coord_y(:), coord_z(:)

integer(int_cgns), allocatable :: connecCGNS(:)
integer(int_size), allocatable :: connecPar(:,:)

! ################################################################################################
! ------------------------ VARS for MPI COMMS ----------------------------------------------------
! ################################################################################################

integer(int_size), allocatable :: matrixCommScheme(:,:),ranksToComm(:)
integer(int_size), allocatable :: commsMemPosInLoc(:),commsMemSize(:),commsMemPosInNgb(:)
integer(int_size) :: numNodesToComm,numRanksWithComms

! ################################################################################################
! --------------------------------- END VARS  ----------------------------------------------------
! ################################################################################################

contains

! ################################################################################################
! --------------------------- FOR READING ALYA MESH ----------------------------------------------
! ################################################################################################
   subroutine read_alya_mesh_files(file_path,file_name)
      implicit none     
      character(500), intent(in)  :: file_path, file_name
      integer(4) :: nelem, npoin, nboun, nbcodes

      call read_dims_file(file_path,file_name,npoin,nelem,nboun)
      allocate(connecGMSH(nelem,nnode))
      allocate(boundGMSH(nboun,npbou))
      allocate(coordGMSH(npoin,ndime))

      !if (nboun .ne. 0) then
      !   allocate(boundGMSH(nboun,npbou))
      !   allocate(bou_codes(nboun,2))
      !   call read_fixbou(file_path,file_name,nboun,nbcodes,bou_codes)
      !end if

      call read_geo_dat_file(file_path,file_name,npoin,nelem,nboun)

      totalNumElements = nelem
      totalNumNodesSrl = npoin

   end subroutine read_alya_mesh_files

   subroutine read_dims_file(file_path,file_name,npoin,nelem,nboun)
      implicit none     
      character(500), intent(in)  :: file_path, file_name
      integer(4)    , intent(out) :: npoin, nelem, nboun
      character(500)              :: file_type, line  

      write(file_type,*) ".dims.dat"

      if(mpi_rank.eq.0)write(*,*) "--| READING DIMS FILE..."
      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      read(99,*) line, npoin
      if(mpi_rank.eq.0)write(*,*) "--| NODES ON MESH : ",npoin
      read(99,*) line, nelem
      if(mpi_rank.eq.0)write(*,*) "--| ELEMENTS ON MESH : ",nelem
      read(99,*) line, nboun
      if(mpi_rank.eq.0)write(*,*) "--| BOUNDARY ELEMENTS ON MESH : ",nboun
      if(mpi_rank.eq.0)write(*,*) "--| END OF DIMS FILE!"
      close(99)
    
   end subroutine read_dims_file

   subroutine read_geo_dat_file(file_path,file_name,npoin,nelem,nboun)
      implicit none
      character(500), intent(in)  :: file_path, file_name
      integer(4)    , intent(in)  :: npoin, nelem, nboun
      integer(4)                  :: iline, int1, inode, idime, aux(nnode+1), bou_aux(npbou+1)
      character(2000)             :: file_type, line

      write(file_type,*) ".geo.dat"

      if(mpi_rank.eq.0) write(*,*) "--| READING GEO.DAT FILE..."
      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      !
      ! Nodes/element section, blank for now
      !
      read(99,*) ! Section header
      do iline = 1,nelem
         read(99,*)
      end do
      read(99,*) ! Section ender
      !
      ! Connectivity table section
      !
      if(mpi_rank.eq.0)write(*,*) "--| READING ELEMENT TABLE..."
      read(99,*) ! Section header
      do iline = 1,nelem
         read(99,'(a)') line
         read(line,*) (aux(inode), inode=1,nnode+1)
         connecGMSH(iline,1:nnode) = aux(2:nnode+1)
         !write(*,*) 'connecGMSH(',iline,'): ', connecGMSH(iline,:)
      end do
      read(99,*) ! Section ender
      !
      ! Nodal coordinates section
      !
      if(mpi_rank.eq.0)write(*,*) "--| READING COORDINATES..."
      read(99,*) ! Section header
      do iline = 1,npoin
         if (ndime == 2) then
            read(99,*) int1, coordGMSH(iline,1), coordGMSH(iline,2)
         else if (ndime == 3) then
            read(99,*) int1, coordGMSH(iline,1), coordGMSH(iline,2), coordGMSH(iline,3)
         end if
      end do
      read(99,*) ! Section ender
      !
      ! Boundary nodes section
      !
      if(mpi_rank.eq.0)write(*,*) "--| READING BOUNDARIES..."
      read(99,*) line! Section header
      do iline = 1,nboun
         read(99,'(a)') line
         read(line,*) (bou_aux(inode), inode=1,npbou+1)
         boundGMSH(iline,1:npbou) = bou_aux(2:npbou+1)
      end do
      close(99)
      if(mpi_rank.eq.0)write(*,*) "--| END OF GEO.DAT FILE!"
    
   end subroutine read_geo_dat_file

! ################################################################################################
! ---------------------------- CGNS MESH PARTITIONING --------------------------------------------
! ################################################################################################

   subroutine do_mesh_partitioning()
      implicit none

      call do_element_partitioning_serial()

      call do_element_partitioning_gempa()

      call do_node_partitioning_and_connectivity()

   end subroutine do_mesh_partitioning

   subroutine do_element_partitioning_serial()
      implicit none

      !---- number of elements and range this process will write
      numElemsInRank = (totalNumElements + mpi_size - 1) / mpi_size
      rankElemStart = numElemsInRank * mpi_rank + 1
      rankElemEnd   = numElemsInRank * (mpi_rank + 1)
      if (rankElemEnd .gt. totalNumElements) rankElemEnd = totalNumElements
      numElemsInRank = rankElemEnd - rankElemStart + 1

      !write(*,*) '#rank ',mpi_rank,' elemsInRank ',numElemsInRank,' iElemS ',rankElemStart,' iElemE ',rankElemEnd
      !write(*,*) '#rank ',mpi_rank,' totalNumElements ',totalNumElements,' totalNumNodesSrl ',totalNumNodesSrl

   end subroutine do_element_partitioning_serial

   subroutine do_element_partitioning_gempa()
      implicit none
      real(8), dimension(numElemsInRank) :: x,y,z
      integer, dimension(numElemsInRank,2) :: elemPart
      integer, parameter :: nodesToAvg(8) = [1,2,5,6,17,18,21,22]
      integer :: i,j,k,m,iElem,iNodeG,numCoords,iRank,iPos
      integer ::numElemsInRank_par,numElemsInRank_srl,numElems2get
      real(8) :: x_a,y_a,z_a
      integer,dimension(0:mpi_size-1) :: vecNumElemsRank
      integer,dimension(0:mpi_size-1,0:mpi_size-1) :: matNumElemsRank
      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size,target_displacement

      character(128) :: file_name, aux_string_rank

      numElemsInRank_srl = numElemsInRank

      i=1
      do iElem=rankElemStart,rankElemEnd
         elemPart(i,1) = iElem
         x_a=0.
         y_a=0.
         z_a=0.
         do j=1,8
            m = gmsh2ijk(nodesToAvg(j))
            iNodeG = connecGMSH(iElem,m) !de moment poso el primer, despres ja fare avg
            x_a = x_a + coordGMSH(iNodeG,1)
            y_a = y_a + coordGMSH(iNodeG,2)
            z_a = z_a + coordGMSH(iNodeG,3)
         end do

         x(i) = x_a/8.d0
         y(i) = y_a/8.d0
         z(i) = z_a/8.d0
         i=i+1
      end do

      numCoords=numElemsInRank_srl

      !----------------------------------------------------------------------------------
      !@TODO: en el cas de malles periodiques, aqui haurem de fer un trucu del almendrucu
      !       per associar els elements que siguin periodics
      !       un cop fet el particionament amb els elements associats, desfer aquesta
      !       'assosiacio' ficticia i posarlos on toca
      !----------------------------------------------------------------------------------


      !---- CALLING GEMPA in PARALLEL-----------------------
      call gempa_do_partition(numCoords,mpi_size,x,y,z,elemPart(:,2))
      !---------------------------------------------------------------------

      ! @now we have to info all the ranks about which are the elements that they will 'store'
      ! since we have done the partitioning in paralel...

      vecNumElemsRank(:)=0
      do i=1,numElemsInRank_srl
         iRank = elemPart(i,2)-1 !elemPart came with rank=1:mpi_size
         vecNumElemsRank(iRank) = vecNumElemsRank(iRank) + 1
      end do

#if _CHECK_
      write(aux_string_rank,'(I1)') mpi_rank
      file_name = 'elemPartition_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numElemsInRank_srl
         write(1,'(*(G0.7,:,","))') x(i),y(i),z(i),elemPart(i,1),elemPart(i,2)
      end do
      close(1)
#endif
      !write(*,*) 'vecNumElemsRank[',mpi_rank,'] ',vecNumElemsRank(:)

      !una finestra on comunicar el nou mpi_rank de cada node
      !es una primera fora de treballar, es pot mirar d'optimitzar per no requerir un vector tamany totalNumElments

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*mpi_size

      call MPI_Win_create(vecNumElemsRank,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iRank=0,mpi_size-1
         call MPI_Get(matNumElemsRank(:,iRank),mpi_size,MPI_INTEGER,iRank,target_displacement,&
                       mpi_size,MPI_INTEGER,window_id,mpi_err)
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      numElemsInRank_par=0

      do iRank=0,mpi_size-1
            numElemsInRank_par = numElemsInRank_par + matNumElemsRank(mpi_rank,iRank)
      end do

      !write(*,*) 'numElemsInRankP(,',mpi_rank,')->',numElemsInRank_par,' srl ',numElemsInRank_srl

      call quicksort_matrix_int(elemPart,2)

      i=1
      do iRank=0,mpi_size-1
         j=vecNumElemsRank(iRank)
         if(j.ne.0) then
            k=i+j-1
            call quicksort_matrix_int(elemPart,1,i,k)
            i=k+1
         endif
      end do

      allocate(elemGid(numElemsInRank_par))

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*numElemsInRank_srl

      call MPI_Win_create(elemPart(:,1),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      iPos=1
      do iRank=0,mpi_size-1
         numElems2get=matNumElemsRank(mpi_rank,iRank)
         

         if(numElems2get.ne.0) then

            target_displacement=0
            do i=0,mpi_rank-1
               target_displacement=target_displacement+matNumElemsRank(i,iRank)
            end do

            !if(mpi_rank.eq.2) write(*,*) 'iPos',iPos,' numE ',numElems2get,' td ',target_displacement

            call MPI_Get(elemGid(iPos),numElems2get,MPI_INTEGER,iRank,target_displacement,&
                       numElems2get,MPI_INTEGER,window_id,mpi_err)

            iPos = iPos + numElems2get
         end if

      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
!---------------------------------------------------------------------------------
      vecNumElemsRank(:)=0
      do i=0,mpi_size-1
         do j=0,mpi_size-1
            vecNumElemsRank(i)=vecNumElemsRank(i)+matNumElemsRank(i,j)
         end do
      end do

      !write(*,*) 'vecNumElemsRank[',mpi_rank,'] ',vecNumElemsRank(:)

      numElemsInRank = numElemsInRank_par

      rankElemStart=1
      do iRank=0,mpi_rank-1
         rankElemStart=rankElemStart+vecNumElemsRank(iRank)
      end do
      rankElemEnd = rankElemStart + vecNumElemsRank(iRank) - 1

      !write(*,*) 'vecNumElemsRank[',mpi_rank,'] ',vecNumElemsRank(:),' eStart ',rankElemStart,' eEnd ', rankElemEnd


!------------------------------------------------------------------------------

      !i=1    
      !do iElem=rankElemStart,rankElemEnd
      !   write(*,*) '[',mpi_rank,']iElem[',iElem,'] elemPart[',elemPart(i,1),' ',elemPart(i,2)
      !   i=i+1
      !end do

      !do i=1,numElemsInRank_par
      !   write(*,*) '[',mpi_rank,']i[',i,'] iElemG:',elemGid(i)
      !end do

#if _CHECK_
      file_name = 'elemGid_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numElemsInRank_par
         write(1,'(*(G0.7,:,","))') i,elemGid(i)
         !write(*,*) '[',mpi_rank,']i[',i,'] iElemG:',elemGid(i)
      end do
      close(1)

#endif
!------------------------------------------------------------------------------

   end subroutine do_element_partitioning_gempa

   subroutine do_node_partitioning_and_connectivity()
      implicit none
      integer, dimension(0:mpi_size-1) :: iNodeStartPar, iNodeEndPar
      integer, allocatable :: boundaryNodes(:)
      integer, allocatable :: vecSharedBN_full(:)

      !allocate(iNodeStartPar(0:mpi_size-1))
      !allocate(iNodeEndPar(0:mpi_size-1))

      call get_rankPartitionBoundaryNodes(boundaryNodes)

      call define_parallelNodePartitioning(iNodeStartPar,iNodeEndPar)

      call define_mpi_boundaries_inPar(boundaryNodes,vecSharedBN_full)

      call reorder_nodes_in_proc(iNodeStartPar)

      !aixo ho podria ficar en un altre modul... crec que quedaria mes 'endreçat'
      call generate_mpi_comm_scheme(vecSharedBN_full)

      deallocate(boundaryNodes)
      deallocate(vecSharedBN_full)

   end subroutine do_node_partitioning_and_connectivity

   subroutine get_rankPartitionBoundaryNodes(boundaryNodes)
      implicit none
      integer, allocatable, intent(out) :: boundaryNodes(:)

      integer, allocatable :: nodeType(:), nodeOwned(:)

      integer :: i,j,k,ind,nt,iElemL,iElemG,iNodeG,iRank
      integer :: numBNodesRankPar, numINodesRankPar, nodeCnt

      character(128) :: file_name, aux_string_rank

      allocate(nodeOwned(totalNumNodesSrl))
      allocate(nodeType(totalNumNodesSrl))

      nodeOwned=0
      nodeType=0

      do iElemL=1,numElemsInRank
         !iElemG = (iElemL-1) + rankElemStart
         iElemG = elemGid(iElemL)
         do k = 0,porder
            do i = 0,porder
               do j = 0,porder
                  ind = ((porder+1)**2)*k+(porder+1)*i+j+1
                  iNodeG = connecGMSH(iElemG,gmsh2ijk(ind))
                  nt = tn2ijk(ind)

                  nodeOwned(iNodeG) = nodeOwned(iNodeG) + 1

                 if(nodeType(iNodeG) .eq. 0) then
                    nodeType(iNodeG) = nt
                 else
                   if(nodeType(iNodeG) .ne. nt) write(*,*) 'fuck node', iNodeG,' nt ', nt, ' nt(ig) ', nodeType(iNodeG)
                 end if
               end do
            end do
         end do
      end do

      numNodesRankPar=0
      numBNodesRankPar=0
      numINodesRankPar=0
      do iNodeG=1,size(nodeOwned)
         nodeCnt = nodeOwned(iNodeG)
         if(nodeCnt .ne. 0) then
            numNodesRankPar = numNodesRankPar+1
            select case (nodeType(iNodeG))
               case (tV)
                  if(nodeCnt.ne.8) then 
                     numBNodesRankPar=numBNodesRankPar+1
                  else
                     numINodesRankPar=numINodesRankPar+1
                  end if
               case (tE)
                  if(nodeCnt.ne.4) then 
                     numBNodesRankPar=numBNodesRankPar+1
                  else
                     numINodesRankPar=numINodesRankPar+1
                  end if
               case (tF)
                  if(nodeCnt.ne.2) then 
                     numBNodesRankPar=numBNodesRankPar+1
                  else
                     numINodesRankPar=numINodesRankPar+1
                  end if
               case (tI)
                  numINodesRankPar=numINodesRankPar+1
               case default
                  write(*,*) "FUCKING ERROR!!!"
            end select
         end if
      end do

      !write(*,*) '#rank[',mpi_rank,']  nodesInRank ',numNodesRankPar, " bN ", numBNodesRankPar, " iN ",numINodesRankPar
     
      allocate(boundaryNodes(numBNodesRankPar))

      i=1
      do iNodeG=1,size(nodeOwned)
         nodeCnt = nodeOwned(iNodeG)
         if(nodeCnt .ne. 0) then
            select case (nodeType(iNodeG))
               case (tV)
                  if(nodeCnt.ne.8) then 
                     boundaryNodes(i) = iNodeG
                     i=i+1
                  end if
               case (tE)
                  if(nodeCnt.ne.4) then 
                     boundaryNodes(i) = iNodeG
                     i=i+1
                  end if
               case (tF)
                  if(nodeCnt.ne.2) then 
                     boundaryNodes(i) = iNodeG
                     i=i+1
                  end if
            end select
         end if
      end do

#if _CHECK_
      write(aux_string_rank,'(I1)') mpi_rank
      file_name = 'boundaryNodes_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numBNodesRankPar
         iNodeG=boundaryNodes(i)
         write(1,*) coordGMSH(iNodeG,1),',',coordGMSH(iNodeG,2),',',coordGMSH(iNodeG,3),',',iNodeG
      end do
      close(1)
#endif
   end subroutine get_rankPartitionBoundaryNodes

   subroutine get_serialNodePartitioning(numNodesRankSrl,iNodeStartSrl,iNodeEndSrl)
      integer, intent(out) :: numNodesRankSrl
      integer, intent(out) :: iNodeStartSrl(0:), iNodeEndSrl(0:)
      integer :: iRank

      numNodesRankSrl = (totalNumNodesSrl + mpi_size - 1) / mpi_size
      do iRank=0,mpi_size-1
          iNodeStartSrl(iRank) = numNodesRankSrl*iRank + 1
          iNodeEndSrl(iRank)   = numNodesRankSrl*(iRank + 1)
      end do
      if(iNodeEndSrl(mpi_size-1) .gt. totalNumNodesSrl) iNodeEndSrl(mpi_size-1) = totalNumNodesSrl
      numNodesRankSrl = iNodeEndSrl(mpi_rank) - iNodeStartSrl(mpi_rank)+1

      !do iRank=0,mpi_size-1
      !    write(*,*) ' ##rank ', iRank , ' iNS ', iNodeStartSrl(iRank), ' iNE ', iNodeEndSrl(iRank)
      !end do

   end subroutine get_serialNodePartitioning

   subroutine define_parallelNodePartitioning(iNodeStartPar,iNodeEndPar)
      integer,dimension(0:mpi_size-1),intent(out) :: iNodeStartPar, iNodeEndPar
      integer, allocatable :: vectorNumNodesRankPar(:)

      integer :: window_id
      integer :: iRank, auxCnt

      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      allocate(vectorNumNodesRankPar(0:mpi_size-1))
      vectorNumNodesRankPar=0

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*1

      call MPI_Win_create(numNodesRankPar, window_buffer_size, mpi_integer_size,&
                         MPI_INFO_NULL, MPI_COMM_WORLD, window_id, mpi_err)
      call MPI_Win_fence(0, window_id, mpi_err)

      target_displacement=0
      do iRank=0,mpi_size-1
         call MPI_Get(vectorNumNodesRankPar(iRank),1,MPI_INTEGER,iRank,target_displacement,&
                     1,MPI_INTEGER,window_id,mpi_err)
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      !write(*,*) 'rank[',mpi_rank,'] -> ', vectorNumNodesRankPar(:), ' nNrp ', numNodesRankPar

      auxCnt = 1
      do iRank=0,mpi_size-1
         iNodeStartPar(iRank) = auxCnt
         iNodeEndPar(iRank)   = iNodeStartPar(iRank) + vectorNumNodesRankPar(iRank) - 1
         auxCnt = iNodeEndPar(iRank) + 1
      end do

      totalNumNodesPar = iNodeEndPar(mpi_size-1)

      rankNodeStart = iNodeStartPar(mpi_rank)
      rankNodeEnd   = iNodeEndPar(mpi_rank)
      !write(*,*) '#rank ', mpi_rank, ' nodesInRank ', numNodesRankPar, ' iNodeS ', rankNodeStart, ' iNodeE ', rankNodeEnd
      !write(*,*) 'rank[',mpi_rank,'] -> start ',iNodeStartPar(:),' end ', iNodeEndPar(:), 'tNNP ', totalNumNodesPar

   end subroutine define_parallelNodePartitioning


   subroutine define_mpi_boundaries_inPar(boundaryNodes,vecSharedBN_full)
      implicit none
      integer, intent(in) :: boundaryNodes(:)
      integer, intent(out), allocatable :: vecSharedBN_full(:)

      integer :: numNodesRankSrl, numBndNodes, numBndNodesRank

      integer :: window_id,  origin_cnt
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      integer :: i,j,k,iNodeL,iNodeGSrl,iRank
      integer, allocatable :: vectorBN(:), matrixBN(:,:),vecSharedBN_part(:),auxVecDebug(:)
      integer, dimension(0:mpi_size-1) :: vecAuxCnt,iNodeStartSrl, iNodeEndSrl
      integer :: numRanksCnt,auxCnt,auxPos

      character(128) :: file_name, aux_string_rank

      ! getting a serial partitioning of the nodes without overlaping
      ! just to do the boundary calc in parallel 
      call get_serialNodePartitioning(numNodesRankSrl,iNodeStartSrl,iNodeEndSrl)

      !write(*,*) 'totalNumNodesSrl ', totalNumNodesSrl, ' numNodesRankSrl ', numNodesRankSrl,' mpi_size ', mpi_size

      allocate(vectorBN(totalNumNodesSrl))
      allocate(matrixBN(numNodesRankSrl,0:mpi_size-1))

      vectorBN = 0
      matrixBN = 0


      do i=1,size(boundaryNodes)
         iNodeGSrl = boundaryNodes(i)
         vectorBN(iNodeGSrl) = 1
      end do
      
      !---------------------------------------------------------------
      ! TO CHECK
#if _CHECK_
      write(aux_string_rank,'(I1)') mpi_rank
      file_name = 'boundaryNodes_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do iNodeGSrl=1,size(vectorBN)
           write(1,*) iNodeGSrl,',',vectorBN(iNodeGSrl)
      end do
      close(1)
#endif
      !---------------------------------------------------------------

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*totalNumNodesSrl
 
      call MPI_Win_create(vectorBN, window_buffer_size, mpi_integer_size, MPI_INFO_NULL, MPI_COMM_WORLD, window_id, mpi_err)
      call MPI_Win_fence(0, window_id, mpi_err)

      target_displacement = iNodeStartSrl(mpi_rank)-1
      !write(*,*) 'rank ', mpi_rank, ' targetdisp ', target_displacement
      do iRank=0,mpi_size-1
         call MPI_Get(matrixBN(:,iRank),numNodesRankSrl, MPI_INTEGER, irank, target_displacement,&
         numNodesRankSrl, MPI_INTEGER, window_id, mpi_err)
      end do
    
      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------

      !---------------------------------------------------------------
      ! TO CHECK
#if _CHECK_
      write(aux_string_rank,'(I1)') mpi_rank
      file_name = 'matrixBN_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numNodesRankSrl
           iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
           write(1,*) i,',',iNodeGSrl,',',matrixBN(i,0),',',matrixBN(i,1)
      end do
      close(1)
#endif
      !---------------------------------------------------------------
      !nou metode!
      !generating vector type
      !iNodeGsrl_i,numRanks(n),rank_1,rank_2,...,rank_n,iNodeGsrl_j,numRanks(n),rank_1,rank_2,...,rank_n,iNodeGsrl_k....
      auxCnt=0
      do i=1,numNodesRankSrl
         iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
         numRanksCnt=0
         do iRank=0,mpi_size-1
            if(matrixBN(i,iRank) .eq. 1) then
               numRanksCnt=numRanksCnt+1
            end if
         end do
         if(numRanksCnt.ge.2) then
            auxCnt=auxCnt+(2+numRanksCnt)
         end if
      end do

      !write(*,*) '[',mpi_rank,']auxCnt ',auxCnt
      
      !--------------------------------------------------------------------------------------
      !lets share how many memory position each rank needs!
      vecAuxCnt(mpi_rank) = auxCnt 

      window_buffer_size = mpi_integer_size*1
      call MPI_Win_create(vecAuxCnt(mpi_rank),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)
   
      do iRank=0,mpi_size-1
         target_displacement = 0
         if(iRank .ne. mpi_rank) then
            call MPI_Get(vecAuxCnt(iRank),1,MPI_INTEGER,iRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)
         else
         end if
      end do
   
      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      !if(mpi_rank.eq.0) write(*,*) 'vecAuxCnt ', vecAuxCnt(:)

      allocate(vecSharedBN_part(auxCnt))
      vecSharedBN_part = -1

      auxCnt=0
      do i=1,numNodesRankSrl
         iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
         numRanksCnt=0
         do iRank=0,mpi_size-1
            if(matrixBN(i,iRank) .eq. 1) then
               numRanksCnt=numRanksCnt+1
            end if
         end do
         if(numRanksCnt.ge.2) then           
            auxCnt=auxCnt+1
            vecSharedBN_part(auxCnt) = iNodeGSrl
            auxCnt=auxCnt+1
            vecSharedBN_part(auxCnt) = numRanksCnt
            do iRank=0,mpi_size-1
               if(matrixBN(i,iRank) .eq. 1) then
                  auxCnt=auxCnt+1
                  vecSharedBN_part(auxCnt) = iRank
               end if
            end do
         end if
      end do

      !if(mpi_rank.eq.0) write(*,*) vecSharedBN_part(:)

      auxCnt=0
      do iRank=0,mpi_size-1
         auxCnt=auxCnt+vecAuxCnt(iRank)
      end do
      
      allocate(vecSharedBN_full(auxCnt)) !dummy allocation for the moment to avoid crash

      !--------------------------------------------------------------------------------------
      !lets share between all the ranks the shared nodes/vertex!

      window_buffer_size = mpi_integer_size*vecAuxCnt(mpi_rank)
      call MPI_Win_create(vecSharedBN_part(1),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)
   
      auxCnt=1
      do iRank=0,mpi_size-1
         target_displacement = 0

         call MPI_Get(vecSharedBN_full(auxCnt),vecAuxCnt(iRank),MPI_INTEGER,iRank,target_displacement,&
                     vecAuxCnt(iRank),MPI_INTEGER,window_id,mpi_err)
         auxCnt=auxCnt+vecAuxCnt(iRank)
      end do
   
      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      !if(mpi_rank.eq.0) write(*,*) vecSharedBN_full(:)

#if _CHECK_
      file_name = 'vecSharedBN_full_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)

      i=0
      do while(i<size(vecSharedBN_full))
         i=i+1
         iNodeGsrl  = vecSharedBN_full(i)
         i=i+1
         numRanksCnt =vecSharedBN_full(i)
         allocate(auxVecDebug(numRanksCnt))
         do j=1,numRanksCnt
            i=i+1
            auxVecDebug(j)=vecSharedBN_full(i)
         end do

         write(1,'(*(G0.7,:,","))')iNodeGSrl,coordGMSH(iNodeGSrl,1),coordGMSH(iNodeGSrl,2),coordGMSH(iNodeGSrl,3),&
               auxVecDebug(:)

         deallocate(auxVecDebug)
      end do
      close(1)
#endif

#if 0
!JUST TO CHECK THE SHIT OF MPI_ACCUMULATE!!!
      i=0!mpi_rank
      j=mpi_rank
      write(*,*) '[',mpi_rank,']test PRE i',i
      window_buffer_size = mpi_integer_size*1
      call MPI_Win_create(i,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)
   
      do iRank=0,mpi_size-1
         target_displacement = 0
         call MPI_Accumulate(j,1,MPI_INTEGER,iRank,target_displacement,1,MPI_INTEGER,MPI_SUM,window_id,mpi_err)
      end do
   
      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)

      write(*,*) '[',mpi_rank,']test POST i',i

#endif

#if 0
!AIXO era el metode antic! i no permetia tenir mes de dos procs per node! vigila!!!!
      do i=1,numNodesRankSrl
         iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
         j=0
         do iRank=0,mpi_size-1
            if(j>2) then
               !####################################################################
               !OJO! QUE NO HO HAVIA PREVIST!! PERO AIXO SI POT PASSAR!!!! SITUACIO
               !  |       |       |
               !  |   0   |   1   |
               !  |-------·--------   <-- this vertex/node will share 4 PROCS!!!
               !  |       |       |
               !  |   3   |   4   |
               !
               write(*,*) 'FUCK ERROR in define_mpi_boundaries_inPar!!!!!! j CANNOT BE MORE THAN 2'
               write(*,*) 'This program is going to exit.'
               call exit(0)
            end if
            if(matrixBN(i,iRank) .eq. 1) then
               j=j+1
               matrixSharedBN_full(iNodeGSrl,j) = iRank
            end if
         end do
      end do

      !---------------------------------------------------------------
      ! TO CHECK
#if _CHECK_
      write(aux_string_rank,'(I1)') mpi_rank
      file_name = 'matrixBN_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numNodesRankSrl
           iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
           write(1,*) i,',',iNodeGSrl,',',matrixBN(i,0),',',matrixBN(i,1)
      end do
      close(1)

      file_name = 'matrixSharedBN_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numNodesRankSrl
         iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
         if(matrixSharedBN_full(iNodeGSrl,1).ne.-1 .and. matrixSharedBN_full(iNodeGSrl,2).ne.-1) then
            write(1,*) coordGMSH(iNodeGSrl,1),',',coordGMSH(iNodeGSrl,2),',',coordGMSH(iNodeGSrl,3),',',&
            iNodeGSrl,',',i,',',matrixSharedBN_full(iNodeGSrl,1),',',matrixSharedBN_full(iNodeGSrl,2)
         end if
      end do
      close(1)
#endif
      !---------------------------------------------------------------

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*numNodesRankSrl
 
      call MPI_Win_create(matrixSharedBN_full(iNodeStartSrl(mpi_rank),1), window_buffer_size,&
                                         mpi_integer_size, MPI_INFO_NULL, MPI_COMM_WORLD,& 
                                         window_id, mpi_err)
      call MPI_Win_fence(0, window_id, mpi_err)
    
      target_displacement = 0 

      do iRank=0,mpi_size-1
         if(mpi_rank .ne. iRank) then
            origin_cnt = iNodeEndSrl(iRank) - iNodeStartSrl(iRank) + 1

            !write(*,*) 'irank ', iRank, ' origin_cnt ', origin_cnt,' tar_disp ', target_displacement
            call MPI_Get(matrixSharedBN_full(iNodeStartSrl(iRank),1),origin_cnt, MPI_INTEGER,&
                         irank,target_displacement,origin_cnt, MPI_INTEGER, window_id, mpi_err)
         end if
      end do
    
      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)

      call MPI_Win_create(matrixSharedBN_full(iNodeStartSrl(mpi_rank),2), window_buffer_size,&
                                         mpi_integer_size, MPI_INFO_NULL, MPI_COMM_WORLD,& 
                                         window_id, mpi_err)
      call MPI_Win_fence(0, window_id, mpi_err)
    
      do iRank=0,mpi_size-1
         if(mpi_rank .ne. iRank) then
            origin_cnt = iNodeEndSrl(iRank) - iNodeStartSrl(iRank) + 1
            call MPI_Get(matrixSharedBN_full(iNodeStartSrl(iRank),2),origin_cnt, MPI_INTEGER,&
                         irank,target_displacement,origin_cnt, MPI_INTEGER, window_id, mpi_err)
         end if
      end do
    
      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)

      !---------------------------------------------------------------
      ! TO CHECK
#if _CHECK_
      file_name = 'matrixSharedBN_Glob_rank'// trim(aux_string_rank)//'.csv'
      do iNodeG=1,totalNumNodesSrl
         if(matrixSharedBN_full(iNodeG,1).ne.-1 .and. matrixSharedBN_full(iNodeG,2).ne.-1) then
            write(1,*) coordGMSH(iNodeG,1),',',coordGMSH(iNodeG,2),',',coordGMSH(iNodeG,3),',',&
            iNodeG,',',matrixSharedBN_full(iNodeG,1),',',matrixSharedBN_full(iNodeG,2)
         end if
      end do
      close(1)
#endif
      !---------------------------------------------------------------
#endif
   end subroutine define_mpi_boundaries_inPar

   subroutine reorder_nodes_in_proc(iNodeStartPar)
      integer,dimension(0:mpi_size-1),intent(in) :: iNodeStartPar
      integer :: iPos,m,indConn,indexIJK,indexCGNS,indexGMSH,indexNew
      integer :: iNodeL,iNodeGpar,iNodeGsrl,iElemL,iElemG

      integer,dimension(nnode) :: auxNodeNewOrderInElem,auxNewOrderIndex,auxCGNSorder
      integer,dimension(totalNumNodesSrl) :: isNodeAdded

      allocate(globalIdSrl(numNodesRankPar))
      allocate(globalIdPar(numNodesRankPar))

      allocate( connecCGNS(numElemsInRank*nnode) )
      allocate( connecPar(numElemsInRank,nnode) )

      isNodeAdded=-1
      iPos = 1

      indConn = -1

      call generate_new_nodeOrder_and_connectivity(cgns2ijk,auxNewOrderIndex,auxCGNSorder)
      !call generate_new_nodeOrder_and_connectivity(dummy2ijk,auxNewOrderIndex,auxCGNSorder)

      do iElemL=1,numElemsInRank
         !iElemG = (iElemL-1) + rankElemStart
         iElemG = elemGid(iElemL)
         auxNodeNewOrderInElem(:)=0

         do indexIJK=1,nnode
            indexGMSH = gmsh2ijk(indexIJK)
            iNodeGsrl = connecGMSH(iElemG,indexGMSH)

            indexNew = auxNewOrderIndex(indexIJK)
            
            auxNodeNewOrderInElem(indexNew) = iNodeGsrl
         end do

         do m=1,nnode
            iNodeGsrl = auxNodeNewOrderInElem(m)

            if(isNodeAdded(iNodeGsrl) < 0) then !node not added put it in the list
               iNodeL = iPos
               isNodeAdded(iNodeGSrl)  = iNodeL

               iNodeGPar = iNodeL + iNodeStartPar(mpi_rank) - 1

               globalIdSrl(iNodeL) = iNodeGsrl
               globalIdPar(iNodeL) = iNodeGPar

               iPos=iPos+1

            else 
               iNodeL = isNodeAdded(iNodeGSrl)
               iNodeGPar = globalIdPar(iNodeL)
            endif

            connecPar(iElemL,m) = iNodeL

            indexCGNS = auxCGNSorder(m)
            indConn = (iElemL-1)*nnode + indexCGNS
            
            connecCGNS(indConn) = iNodeGPar

         end do

         !write(*,*) '[',mpi_rank,']iElemG ',iElemG,' connecPar ',connecPar(iElemL,:)

      end do
   end subroutine reorder_nodes_in_proc

   subroutine generate_new_nodeOrder_and_connectivity(newOrderIJK,auxNewOrderIndex,auxCGNSorder)
      integer,intent(in)  :: newOrderIJK(:)
      integer,intent(out) :: auxNewOrderIndex(:),auxCGNSorder(:)
      integer :: i,j,k,indexIJK,indexNew,indexCGNS

      do k = 0,porder
         do i = 0,porder
            do j = 0,porder
               indexIJK = ((porder+1)**2)*k+(porder+1)*i+j+1

               indexNew = newOrderIJK(indexIJK)
               indexCGNS = cgns2ijk(indexIJK) !posicio requerida en el connec de cgns

               auxNewOrderIndex(indexIJK) = indexNew
               auxCGNSorder(indexNew) = indexCGNS

               !write(*,*) 'test->indexIJK ', indexIJK, ' iNew ', indexNew,' aux ', auxCGNSorder(indexNew)
            end do
         end do
      end do
   end subroutine generate_new_nodeOrder_and_connectivity

   subroutine generate_mpi_comm_scheme(vecSharedBN_full)
      !generate a matrix with the comm schemes for shared nodes between procs
      integer, intent(in)  :: vecSharedBN_full(:)
      integer, allocatable :: auxVecRanks(:)
      integer, dimension(0:mpi_size-1) :: commSchemeNumNodes
      integer,dimension(mpi_size*2) :: commSchemeStartEndNodes
      
      logical :: imIn
      integer :: i,j,k,iNodeL,iNodeGsrl,iRank,numRanksCnt
      integer :: window_id
      
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      character(128) :: file_name, aux_string_rank

      i=0
      numNodesToComm=0
      do while(i<size(vecSharedBN_full))
         i=i+1
         iNodeGsrl  = vecSharedBN_full(i)
         i=i+1
         numRanksCnt =vecSharedBN_full(i)

         imIn=.false.
         do j=1,numRanksCnt
            i=i+1
            if(vecSharedBN_full(i).eq.mpi_rank) imIn = .true.
         end do
         if(imIn) then
            numNodesToComm=numNodesToComm+(numRanksCnt-1)
         end if
      end do

      !write(*,*) '[',mpi_rank,']numNodesToComm ', numNodesToComm
      allocate(matrixCommScheme(numNodesToComm,3))
      matrixCommScheme(:,:)=-1

      i=0
      k=0
      do while(i<size(vecSharedBN_full))
         i=i+1
         iNodeGsrl  = vecSharedBN_full(i)
         i=i+1
         numRanksCnt =vecSharedBN_full(i)

         imIn=.false.
         allocate(auxVecRanks(numRanksCnt))
         do j=1,numRanksCnt
            i=i+1
            auxVecRanks(j)=vecSharedBN_full(i)
            if(auxVecRanks(j).eq.mpi_rank) imIn = .true.
         end do
         if(imIn) then
            iNodeL = gidSrl_to_lid(iNodeGSrl)
            do j=1,numRanksCnt
               if(auxVecRanks(j).ne.mpi_rank) then
                  k=k+1
                  matrixCommScheme(k,1) = iNodeL
                  matrixCommScheme(k,2) = iNodeGSrl
                  matrixCommScheme(k,3) = auxVecRanks(j)
                  !write(*,*) '[',mpi_rank,'] k ',k,',',matrixCommScheme(k,:)
               end if
            end do
         end if
         deallocate(auxVecRanks)
      end do

      !first, we sort the matrix depending on the rank
      call quicksort_matrix_int(matrixCommScheme,3)

      !second, determine how many nodes are shared with each other ranks
      commSchemeNumNodes(:)=0
      do i=1,numNodesToComm
         iRank=matrixCommScheme(i,3)
         commSchemeNumNodes(iRank)=commSchemeNumNodes(iRank)+1
      end do

      !third, sort the matrix depending on the global id for each rank
      !taking advantage of the loop, i'll generate the auxiliar vector commSchemeStartEndNodes
      !this vector will allow the other ranks know where (memory positions) they have to share the data with the shared data vector of this rank
      commSchemeStartEndNodes(:)=-1
      numRanksWithComms=0
      i=1
      do iRank=0,mpi_size-1
         j=commSchemeNumNodes(iRank)
         if(j.ne.0) then
            k=i+j-1
            !write(*,*) '[',mpi_rank,']sort rank',irank,' from ',i,' to ', k
            call quicksort_matrix_int(matrixCommScheme,2,i,k)
            !start position
            commSchemeStartEndNodes(2*iRank+1)=i
            !end position
            commSchemeStartEndNodes(2*iRank+2)=k
            i=k+1
            numRanksWithComms=numRanksWithComms+1
         endif
      end do

      !now I generate the vector ranksToComm, storing the ranks who with this rank will comm
      allocate(ranksToComm(numRanksWithComms))
      i=1
      do iRank=0,mpi_size-1
         j=commSchemeNumNodes(iRank)
         if(j.ne.0) then
         ranksToComm(i)=iRank
         i=i+1
         end if
      end do

      !now generate the vectors memPosIn.... to know the memory positions to where put/receive the info of ngb ranks
      allocate(commsMemPosInLoc(numRanksWithComms))
      allocate(commsMemPosInNgb(numRanksWithComms))
      allocate(commsMemSize(numRanksWithComms))
      !allocate(commsMemPosInNgbE(numRanksWithComms))
      do i=1,numRanksWithComms
         iRank=ranksToComm(i)
         commsMemPosInLoc(i)=commSchemeStartEndNodes(2*iRank+1)
         commsMemSize(i)=commSchemeStartEndNodes(2*iRank+2)-commsMemPosInLoc(i)+1
      end do

      !--------------------------------------------------------------------------------------

      window_buffer_size = mpi_integer_size*(mpi_size*2)
      call MPI_Win_create(commSchemeStartEndNodes,window_buffer_size,mpi_integer_size,&
                         MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)
   
      j=0
      do i=1,numRanksWithComms
         iRank=ranksToComm(i)
         j=j+1

         target_displacement = 2*mpi_rank
         call MPI_Get(commsMemPosInNgb(j),1,MPI_INTEGER,iRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)

      end do
   
      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

 
      !write(*,*) '[',mpi_rank,']csNumNodes->',commSchemeNumNodes(:)
      !write(*,*) '[',mpi_rank,']csStartEnd->',commSchemeStartEndNodes(:)
      !write(*,*) '[',mpi_rank,']numRanksWC->',numRanksWithComms
      !write(*,*) '[',mpi_rank,']ranksToComm->',ranksToComm
      !write(*,*) '[',mpi_rank,']commsMemPosInLoc->',commsMemPosInLoc
      !write(*,*) '[',mpi_rank,']commsMemPosInNgb->',commsMemPosInNgb
      !write(*,*) '[',mpi_rank,']commsMemSize->',commsMemSize

      !ara podria obtenir els ranks locals de les altres cpus... pero cal?

#if _CHECK_
      write(aux_string_rank,'(I1)') mpi_rank
      file_name = 'matrixCommScheme_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
         do i=1,numNodesToComm
             write(1,'(*(I5,:,","))') matrixCommScheme(i,:)
         end do
      close(1)
#endif

   end subroutine generate_mpi_comm_scheme
!----------------------------------------------------------------------------------------------------------

! ################################################################################################
! ----------------------------------- AUXILIAR FUNCS ---------------------------------------------
! ################################################################################################

   integer function gidSrl_to_lid(gidSrl) result(lid)
      implicit none
      integer, intent(in) :: gidSrl
      integer :: i

      lid = -1
      do i=1,numNodesRankPar
         if(gidSrl .eq. globalIdSrl(i)) then
            lid = i
            exit
         endif
      end do

   end function gidSrl_to_lid

   integer function gidPar_to_lid(gidPar) result(lid)
      implicit none
      integer, intent(in) :: gidPar
      integer :: i

      lid = -1
      do i=1,numNodesRankPar
         if(gidPar .eq. globalIdPar(i)) then
            lid = i
            exit
         endif
      end do

   end function gidPar_to_lid


end module mod_mpi_mesh
