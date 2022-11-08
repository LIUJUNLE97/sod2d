module mod_mesh_boundaries
   use mod_constants
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms

contains

   subroutine splitBoundary_inPar()
      integer(4), allocatable    :: aux1(:),aux_ndof(:),aux_nbnodes(:)
      integer(4) :: ii,iNodeL,iNodeGSrl,iNodeGSrl_bound,iBound,iBoundL,ipbou
      integer(4) :: auxBoundCnt,aux_sum_NBRP
      integer(4) :: aux_boundList(totalNumBoundsSrl)

      if(mpi_rank.eq.0) write(*,*) ' # Splitting boundary nodes from DoFs in parallel...'

      !first, generate boundPar(:,:) & bouCodesPar(:,:)

      !1. how many boundaries my rank have?
      numBoundsRankPar = 0
      !$acc kernels
      aux_boundList(:) = 0
      !$acc end kernels

      !TODO: if possible this loop must be speedup using GPUS 
      loopBound: do iBound = 1,totalNumBoundsSrl
         auxBoundCnt = 0
         loopIp: do ipbou = 1,npbou
            iNodeGSrl_bound = boundGMSH(iBound,ipbou)
            loopG: do iNodeL = 1,numNodesRankPar
               iNodeGSrl = globalIdSrl(iNodeL)
               if(iNodeGSrl .eq. iNodeGSrl_bound) then
                  auxBoundCnt=auxBoundCnt+1
                  exit loopG
               end if
            end do loopG
            if(auxBoundCnt .eq. npbou) then
               !write(*,*) '[',mpi_rank,']iBound',iBound,'auxBoundCnt',auxBoundCnt
               aux_boundList(iBound) = 1
               numBoundsRankPar = numBoundsRankPar + 1
               exit loopIp
            end if
         end do loopIp
      end do loopBound
      !-----------------------------------------------------------------

      write(*,*) '[',mpi_rank,']numBoundsRankPar',numBoundsRankPar,'totalNumBoundsSrl',totalNumBoundsSrl
      
      call MPI_Allreduce(numBoundsRankPar,aux_sum_NBRP,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,mpi_err)

      if(aux_sum_NBRP.ne.totalNumBoundsSrl) then
         write(*,*) 'ERROR IN splitBoundary_inPar()->aux_sum_NBRP',aux_sum_NBRP,' not equal to totalNumBoundsSrl',totalNumBoundsSrl
         call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
      end if

      allocate(boundPar(numBoundsRankPar,npbou))
      allocate(bouCodesPar(numBoundsRankPar))
      
      ii=0
      !TODO: if possible this loop must be speedup using GPUS 
      do iBound = 1,totalNumBoundsSrl
         if(aux_boundList(iBound).eq.1) then
            ii=ii+1
            do ipbou = 1,npbou
               iNodeGSrl_bound = boundGMSH(iBound,ipbou)
               iNodeL = gidSrl_to_lid(iNodeGSrl_bound)
               boundPar(ii,ipbou) = iNodeL
            end do
            bouCodesPar(ii) = bou_codesGMSH(iBound,2)
         !write(*,*) '[',mpi_rank,']boundPar(',ii,')',boundPar(ii,:)
         end if
      end do
      !-----------------------------------------------------------------

      !------------------------------------------------------------------------
      allocate(aux1(numNodesRankPar))

      !ojo! crec que aqui ho podem fer diferent per tenir tots els nodes boundary reals!
!NEW METHOD
#if 1
      allocate(aux_ndof(numNodesRankPar))
      allocate(aux_nbnodes(numNodesRankPar))
      !$acc kernels
      aux1(:) = 0
      aux_ndof(:) = 0
      aux_nbnodes(:) = 0
      !$acc end kernels

      !possar tots els boundaries locals a 1
      !$acc parallel loop gang 
      do iBound = 1,numBoundsRankPar
         !$acc loop vector
         do ipbou = 1,npbou
            aux1(boundPar(iBound,ipbou)) = 1
         end do
      end do
      !$acc end parallel loop

      !fer mpi halo atomic
      call mpi_halo_atomic_update_int(aux1)
      !ara se nodes que son boundaries a altres cpus

      ! Determine how many nodes are boundary nodes
      numBoundaryNodesRankPar=0
      ndofRankPar = 0
      do iNodeL = 1,numNodesRankPar
         if (aux1(iNodeL) .ne. 0) then
            numBoundaryNodesRankPar = numBoundaryNodesRankPar+1
            aux_nbnodes(numBoundaryNodesRankPar) = iNodeL
         else
            ndofRankPar = ndofRankPar+1
            aux_ndof(ndofRankPar) = iNodeL
         end if
      end do

      if(ndofRankPar .ne. (numNodesRankPar - numBoundaryNodesRankPar)) then
         write(*,*) 'ERROR IN splitBoundary_inPar()-> ndofRankPar',ndofRankPar,' not equal to (nNRP-nBNRP) ',(numNodesRankPar - numBoundaryNodesRankPar)
         call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
      end if

      !ndofRankPar = numNodesRankPar - numBoundaryNodesRankPar
      write(*,*) '[',mpi_rank,'] ndof',ndofRankPar,'nbnodes',numBoundaryNodesRankPar

      !-------------------------------------------------------------------------------------
      ! Split aux1 into the 2 lists
      allocate(ldofPar(ndofRankPar))
      allocate(lbnodesPar(numBoundaryNodesRankPar))

      !$acc kernels copyout(lbnodesPar,ldofPar,aux2)
      ldofPar(:) = 0
      lbnodesPar(:) = 0
      !$acc end kernels

      !$acc parallel loop copyin(aux_nbnodes) copyout(lbnodesPar)
      do iNodeL = 1,numBoundaryNodesRankPar
         lbnodesPar(iNodeL) = aux_nbnodes(iNodeL)! iNodeL
      end do
      !$acc end parallel loop

      !$acc parallel loop copyin(aux_ndof) copyout(ldofPar)
      do iNodeL = 1,ndofRankPar
         ldofPar(iNodeL) = aux_ndof(iNodeL)! iNodeL
      end do
      !$acc end parallel loop

      deallocate(aux_nbnodes)
      deallocate(aux_ndof)

#else
TODO: DELETE IF EVERYTHING IS OK
!OLD METHOD
      ! Fill aux1 with all nodes in order
      !$acc parallel loop
      do iNodeL = 1,numNodesRankPar
         aux1(iNodeL) = iNodeL
      end do
      !$acc end parallel loop

      ! If node is on boundary, zero corresponding aux1 entry
      !$acc parallel loop gang 
      do iBound = 1,numBoundsRankPar
         !$acc loop vector
         do ipbou = 1,npbou
            aux1(boundPar(iBound,ipbou)) = 0
         end do
      end do
      !$acc end parallel loop

      
      ! Determine how many nodes are boundary nodes
      !
      numBoundaryNodesRankPar=0
      ndofRankPar = 0
      do iNodeL = 1,numNodesRankPar
         if (aux1(iNodeL) .eq. 0) then
            numBoundaryNodesRankPar = numBoundaryNodesRankPar+1
         end if
      end do

      !this%nbnodes = this%ndof    ! Nodes on boundaries
      !this%ndof = this%npoin-this%ndof ! Free nodes
      ndofRankPar = numNodesRankPar - numBoundaryNodesRankPar
      write(*,*) '[',mpi_rank,'] ndof',ndofRankPar,'nbnodes',numBoundaryNodesRankPar

      !-------------------------------------------------------------------------------------
      ! Split aux1 into the 2 lists
      allocate(ldofPar(ndofRankPar))
      allocate(lbnodesPar(numBoundaryNodesRankPar))

      idof = 0    ! Counter for free nodes
      ibnodes = 0 ! Counter for boundary nodes
      !$acc parallel loop reduction(+:idof,ibnodes)
      do iNodeL = 1,numNodesRankPar
         if (aux1(iNodeL) .eq. 0) then
            ibnodes = ibnodes+1
            lbnodesPar(ibnodes) = iNodeL
         else
            idof = idof+1
            ldofPar(idof) = aux1(iNodeL)
         end if
      end do
      !$acc end parallel loop
#endif
      deallocate(aux1)

   end subroutine splitBoundary_inPar

   subroutine generate_boundary_mpi_comm_scheme()
      !generate a matrix with the comm schemes for the boundary nodes shared nodes between procs
      integer,allocatable :: auxVecBndNodes(:)
      integer, dimension(0:mpi_size-1) :: bnd_commSchemeNumNodes
      integer,dimension(mpi_size*2) :: bnd_commSchemeStartEndNodes
      
      logical :: imIn
      integer :: i,j,k,iNodeL,bnd_iNodeL,iNodeGsrl,iRank,numRanksCnt
      integer :: window_id
      
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      character(128) :: file_name, aux_string_rank

      if(mpi_rank.eq.0) write(*,*) ' # Generating MPI Boundary Comm scheme...'

      allocate(auxVecBndNodes(numNodesToComm))
      auxVecBndNodes(:) = 0
      bnd_numNodesToComm = 0

      do k=1,numNodesToComm
         iNodeL = matrixCommScheme(k,1) 
         bndLoop: do i=1,numBoundaryNodesRankPar
            bnd_iNodeL = lbnodesPar(i)
            if(iNodeL .eq. bnd_iNodeL) then
               bnd_numNodesToComm = bnd_numNodesToComm + 1
               auxVecBndNodes(k) = 1
               exit bndLoop
            end if
         end do bndLoop
      end do
   
      allocate(bnd_matrixCommScheme(bnd_numNodesToComm,3))

      !write(*,*) '[',mpi_rank,'] b_nN2C',bnd_numNodesToComm,'nN2C',numNodesToComm

      i=0
      do k=1,numNodesToComm
         if(auxVecBndNodes(k).eq.1) then
            i=i+1
            bnd_matrixCommScheme(i,:) = matrixCommScheme(k,:)
         end if
      end do

      deallocate(auxVecBndNodes)

      !aqui ja tinc la matriu nomes amb els nodes boundary, si faig lo altre, hauria de ser pimpam

      !!!first, we sort the matrix depending on the rank
      !!call quicksort_matrix_int(matrixCommScheme,3)

      ! determine how many nodes are shared with each other ranks
      bnd_commSchemeNumNodes(:)=0
      do i=1,bnd_numNodesToComm
         iRank=bnd_matrixCommScheme(i,3)
         bnd_commSchemeNumNodes(iRank)=bnd_commSchemeNumNodes(iRank)+1
      end do

      !third, sort the matrix depending on the global id for each rank
      !taking advantage of the loop, i'll generate the auxiliar vector commSchemeStartEndNodes
      !this vector will allow the other ranks know where (memory positions) they have to share the data with the shared data vector of this rank
      bnd_commSchemeStartEndNodes(:)=-1
      bnd_numRanksWithComms=0
      i=1
      do iRank=0,mpi_size-1
         j=bnd_commSchemeNumNodes(iRank)
         if(j.ne.0) then
            k=i+j-1
            !write(*,*) '[',mpi_rank,']sort rank',irank,' from ',i,' to ', k
            !call quicksort_matrix_int(matrixCommScheme,2,i,k)
            !start position
            bnd_commSchemeStartEndNodes(2*iRank+1)=i
            !end position
            bnd_commSchemeStartEndNodes(2*iRank+2)=k
            i=k+1
            bnd_numRanksWithComms=bnd_numRanksWithComms+1
         endif
      end do

      !write(*,*) '[',mpi_rank,'] b_nRwC',bnd_numRanksWithComms,'nRwC',numRanksWithComms

      !now I generate the vector ranksToComm, storing the ranks who with this rank will comm
      allocate(bnd_ranksToComm(bnd_numRanksWithComms))
      i=1
      do iRank=0,mpi_size-1
         j=bnd_commSchemeNumNodes(iRank)
         if(j.ne.0) then
         bnd_ranksToComm(i)=iRank
         i=i+1
         end if
      end do

      !now generate the vectors memPosIn.... to know the memory positions to where put/receive the info of ngb ranks
      allocate(bnd_commsMemPosInLoc(bnd_numRanksWithComms))
      allocate(bnd_commsMemPosInNgb(bnd_numRanksWithComms))
      allocate(bnd_commsMemSize(bnd_numRanksWithComms))
      !allocate(commsMemPosInNgbE(numRanksWithComms))
      do i=1,bnd_numRanksWithComms
         iRank=bnd_ranksToComm(i)
         bnd_commsMemPosInLoc(i)=bnd_commSchemeStartEndNodes(2*iRank+1)
         bnd_commsMemSize(i)=bnd_commSchemeStartEndNodes(2*iRank+2)-bnd_commsMemPosInLoc(i)+1
      end do

      !--------------------------------------------------------------------------------------

      window_buffer_size = mpi_integer_size*(mpi_size*2)
      call MPI_Win_create(bnd_commSchemeStartEndNodes,window_buffer_size,mpi_integer_size,&
                         MPI_INFO_NULL,MPI_COMM_WORLD,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)
   
      j=0
      do i=1,bnd_numRanksWithComms
         iRank=bnd_ranksToComm(i)
         j=j+1

         target_displacement = 2*mpi_rank
         call MPI_Get(bnd_commsMemPosInNgb(j),1,MPI_INTEGER,iRank,target_displacement,1,MPI_INTEGER,window_id,mpi_err)

      end do
   
      !! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      !write(*,*) '[',mpi_rank,']bnd_csNumNodes->',bnd_commSchemeNumNodes(:)
      !write(*,*) '[',mpi_rank,']bnd_csStartEnd->',bnd_commSchemeStartEndNodes(:)
      !write(*,*) '[',mpi_rank,']bnd_numRanksWC->',bnd_numRanksWithComms
      !write(*,*) '[',mpi_rank,']bnd_ranksToComm->',bnd_ranksToComm(:)
      !write(*,*) '[',mpi_rank,']bnd_commsMemPosInLoc->',bnd_commsMemPosInLoc(:)
      !write(*,*) '[',mpi_rank,']bnd_commsMemPosInNgb->',bnd_commsMemPosInNgb(:)
      !write(*,*) '[',mpi_rank,']bnd_commsMemSize->',bnd_commsMemSize(:)

   end subroutine generate_boundary_mpi_comm_scheme


end module mod_mesh_boundaries