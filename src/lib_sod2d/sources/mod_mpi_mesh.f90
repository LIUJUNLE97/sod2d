module mod_mpi_mesh
   use mod_constants
   use mod_mpi
   use mod_utils
   use mod_gmsh_indices !potser en el futur pot volar!
   use iso_c_binding
   implicit none
!-----------------------------------

! ################################################################################################
! ----------------- VARS for new Par mesh FORMAT -------------------------------------------------
! ################################################################################################

integer(4) :: numNodesRankPar,numElemsRankPar,totalNumElements
integer(4) :: rankNodeStart,rankNodeEnd,rankElemStart,rankElemEnd
integer(8) :: totalNumNodesPar, totalNumNodesSrl

integer(4),allocatable :: elemGid(:)
integer(8),allocatable :: globalIdSrl(:),globalIdPar(:)

real(rp), allocatable :: coordPar(:,:)

integer(4),allocatable :: connecVTK(:)
integer(4),allocatable :: connecParOrig(:,:),connecParWork(:,:)

integer(4),allocatable :: workingNodesPar(:)
integer(4) :: numWorkingNodesRankPar

integer(4),allocatable :: masSlaRankPar(:,:)
integer(4) :: nPerRankPar
logical :: isMeshPeriodic

integer(4) :: numBoundCodes, numBoundsRankPar, totalNumBoundsSrl
integer(4) :: ndofRankPar, numBoundaryNodesRankPar
integer(4), allocatable :: boundPar(:,:), boundParOrig(:,:),bouCodesPar(:), ldofPar(:), lbnodesPar(:), bouCodesNodesPar(:)
real(rp), allocatable :: boundNormalPar(:,:)
logical :: isMeshBoundaries

!For WallModels
integer(4) :: numBoundsWMRankPar
integer(4), allocatable :: listBoundsWallModel(:)

! ################################################################################################
! ------------------------ VARS for MPI COMMS ----------------------------------------------------
! ################################################################################################

integer(4),allocatable :: nodesToComm(:),ranksToComm(:)
integer(4),allocatable :: commsMemPosInLoc(:),commsMemSize(:),commsMemPosInNgb(:)
integer(4) :: numNodesToComm,numRanksWithComms

integer(4),allocatable :: bnd_nodesToComm(:),bnd_ranksToComm(:)
integer(4),allocatable :: bnd_commsMemPosInLoc(:),bnd_commsMemSize(:),bnd_commsMemPosInNgb(:)
integer(4) :: bnd_numNodesToComm,bnd_numRanksWithComms

! ################################################################################################
! --------------------------------- END VARS  ----------------------------------------------------
! ################################################################################################

character(*), parameter :: fmt_csv_msh = '(1x,*(g0,","))'

contains

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

!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------

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

   subroutine get_parallelNode_partitioning(iNodeStartPar,iNodeEndPar)
      integer,dimension(0:mpi_size-1),intent(out) :: iNodeStartPar, iNodeEndPar
      integer, allocatable :: vectorNumNodesRankPar(:)

      integer :: window_id
      integer :: iRank, auxCnt

      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      allocate(vectorNumNodesRankPar(0:mpi_size-1))

      vectorNumNodesRankPar(:)=0

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*1

      call MPI_Win_create(numNodesRankPar, window_buffer_size, mpi_integer_size,&
                         MPI_INFO_NULL, app_comm, window_id, mpi_err)
      call MPI_Win_fence(0, window_id, mpi_err)

      target_displacement=0
      do iRank=0,mpi_size-1
         call MPI_Get(vectorNumNodesRankPar(iRank),1,mpi_datatype_int,iRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
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

   end subroutine get_parallelNode_partitioning

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

      ! getting a serial partitioning of the nodes without overlaping
      ! just to do the boundary calc in parallel
      call get_serialNodePartitioning(numNodesRankSrl,iNodeStartSrl,iNodeEndSrl)

      !write(*,*) 'totalNumNodesSrl ', totalNumNodesSrl, ' numNodesRankSrl ', numNodesRankSrl,' mpi_size ', mpi_size

      allocate(vectorBN(totalNumNodesSrl))
      allocate(matrixBN(numNodesRankSrl,0:mpi_size-1))

      !$acc kernels
      vectorBN(:) = 0
      matrixBN(:,:) = 0
      !$acc end kernels

      !!!$acc parallel loop
      do i=1,size(boundaryNodes)
         iNodeGSrl = boundaryNodes(i)
         vectorBN(iNodeGSrl) = 1
      end do
      !!!$acc end parallel loop

      ! Create the window
      !--------------------------------------------------------------------------------------
      window_buffer_size = mpi_integer_size*totalNumNodesSrl

      call MPI_Win_create(vectorBN, window_buffer_size, mpi_integer_size, MPI_INFO_NULL, app_comm, window_id, mpi_err)
      call MPI_Win_fence(0, window_id, mpi_err)

      target_displacement = iNodeStartSrl(mpi_rank)-1
      !write(*,*) 'rank ', mpi_rank, ' targetdisp ', target_displacement
      do iRank=0,mpi_size-1
         call MPI_Get(matrixBN(:,iRank),numNodesRankSrl,mpi_datatype_int,iRank,target_displacement,numNodesRankSrl,mpi_datatype_int,window_id,mpi_err)
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------

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

      call MPI_Win_create(vecAuxCnt(mpi_rank),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      do iRank=0,mpi_size-1
         target_displacement = 0
         if(iRank .ne. mpi_rank) then
            call MPI_Get(vecAuxCnt(iRank),1,mpi_datatype_int,iRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
         else
         end if
      end do

      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      !if(mpi_rank.eq.0) write(*,*) 'vecAuxCnt ', vecAuxCnt(:)

      allocate(vecSharedBN_part(auxCnt))
      !$acc kernels
      vecSharedBN_part(:) = -1
      !$acc end kernels

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
            !if(numRanksCnt.eq.4) write(*,*) 'numRanksCnt ',numRanksCnt
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

      allocate(vecSharedBN_full(auxCnt))
      if(auxCnt.gt.0) then
         !--------------------------------------------------------------------------------------
         !lets share between all the ranks the shared nodes/vertex!

         window_buffer_size = mpi_integer_size*vecAuxCnt(mpi_rank)

         call MPI_Win_create(vecSharedBN_part,window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
         !call MPI_Win_create(vecSharedBN_part(1),window_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
         call MPI_Win_fence(0,window_id,mpi_err)

         auxCnt=1
         do iRank=0,mpi_size-1
            target_displacement = 0

            call MPI_Get(vecSharedBN_full(auxCnt),vecAuxCnt(iRank),mpi_datatype_int,iRank,target_displacement,vecAuxCnt(iRank),mpi_datatype_int,window_id,mpi_err)
            auxCnt=auxCnt+vecAuxCnt(iRank)
         end do

         !! Wait for the MPI_Get issued to complete before going any further
         call MPI_Win_fence(0,window_id,mpi_err)
         call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------
      end if
      !if(mpi_rank.eq.0) write(*,*) vecSharedBN_full(:)

   end subroutine define_mpi_boundaries_inPar

   subroutine generate_mpi_comm_scheme_i4(vecSharedBN_full)
      !generate a matrix with the comm schemes for shared nodes between procs
      integer(4), intent(in)  :: vecSharedBN_full(:)
      integer(4), allocatable :: auxVecRanks(:),matrixCommScheme(:,:)
      integer(4), dimension(0:mpi_size-1) :: commSchemeNumNodes
      integer(4), dimension(mpi_size*2) :: commSchemeStartEndNodes

      logical :: imIn
      integer(4) :: i,j,k,iNodeL,iRank,numRanksCnt
      integer(4) :: iNodeGsrl
      integer(4) :: window_id

      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      if(mpi_rank.eq.0) write(*,*) ' # Generating MPI Comm scheme...'

      i=0
      numNodesToComm=0
      do while(i<size(vecSharedBN_full))
         i=i+1
         iNodeGsrl  = vecSharedBN_full(i)
         i=i+1
         numRanksCnt = int(vecSharedBN_full(i),4)

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
         numRanksCnt = int(vecSharedBN_full(i),4)

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
                  !if(numRanksCnt.ge.3) write(*,*) '[',mpi_rank,'] k ',k,',',matrixCommScheme(k,:)
               end if
            end do
         end if
         deallocate(auxVecRanks)
      end do

      !first, we sort the matrix depending on the rank
      call quicksort_matrix_int4(matrixCommScheme,3)

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
            call quicksort_matrix_int4(matrixCommScheme,2,i,k)
            !start position
            commSchemeStartEndNodes(2*iRank+1)=i
            !end position
            commSchemeStartEndNodes(2*iRank+2)=k
            i=k+1
            numRanksWithComms=numRanksWithComms+1
         endif
      end do

      allocate(nodesToComm(numNodesToComm))

      nodesToComm(:) = matrixCommScheme(:,1)

      deallocate(matrixCommScheme)

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
                         MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      j=0
      do i=1,numRanksWithComms
         iRank=ranksToComm(i)
         j=j+1

         target_displacement = 2*mpi_rank
         call MPI_Get(commsMemPosInNgb(j),1,mpi_datatype_int,iRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)

      end do

      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      !write(*,*) '[',mpi_rank,']csNumNodes->',commSchemeNumNodes(:)
      !write(*,*) '[',mpi_rank,']csStartEnd->',commSchemeStartEndNodes(:)
      !write(*,*) '[',mpi_rank,']numRanksWC->',numRanksWithComms
      !write(*,*) '[',mpi_rank,']ranksToComm->',ranksToComm(:)
      !write(*,*) '[',mpi_rank,']commsMemPosInLoc->',commsMemPosInLoc(:)
      !write(*,*) '[',mpi_rank,']commsMemPosInNgb->',commsMemPosInNgb(:)
      !write(*,*) '[',mpi_rank,']commsMemSize->',commsMemSize(:)

      !ara podria obtenir els ranks locals de les altres cpus... pero cal?

   end subroutine generate_mpi_comm_scheme_i4
!----------------------------------------------------------------------------------------------------------

! ################################################################################################
! ----------------------------------- PLOTTING FUNCS ---------------------------------------------
! ################################################################################################

   subroutine print_csv_file_dfield(fileName,dfield)
      implicit none
      character(len=*), intent(in) :: fileName
      real(8), intent(in) :: dfield(:)
      integer :: iNodeL
      character(128) :: full_file_name, aux_string_rank

      write(aux_string_rank,'(I0)') mpi_rank
      full_file_name = trim(fileName) // trim(aux_string_rank)//'.csv'
      open(1, file=full_file_name)

      write(1,*) 'X,Y,Z,d_field'
      do iNodeL=1,numNodesRankPar
         write(1,fmt_csv_msh) coordPar(iNodeL,1),coordPar(iNodeL,2),coordPar(iNodeL,3),dfield(iNodeL)
      end do

      close(1)

   end subroutine print_csv_file_dfield

   subroutine print_csv_file_ffield(fileName,ffield)
      implicit none
      character(len=*), intent(in) :: fileName
      real(4), intent(in) :: ffield(:)
      integer :: iNodeL
      character(128) :: full_file_name, aux_string_rank

      write(aux_string_rank,'(I0)') mpi_rank
      full_file_name = trim(fileName) // trim(aux_string_rank)//'.csv'
      open(1, file=full_file_name)

      write(1,*) 'X,Y,Z,ffield'
      do iNodeL=1,numNodesRankPar
         write(1,fmt_csv_msh) coordPar(iNodeL,1),coordPar(iNodeL,2),coordPar(iNodeL,3),ffield(iNodeL)
      end do

      close(1)


   end subroutine print_csv_file_ffield

end module mod_mpi_mesh
