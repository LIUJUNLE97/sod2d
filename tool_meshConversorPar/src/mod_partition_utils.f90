module mod_partition_utils
   use mod_constants
   use mod_mpi
   use mod_utils
   implicit none

contains

   subroutine do_element_partitioning_serial(numElems2Par,iElemStart,iElemEnd,iElemsInRank,mpiRank,numRanks2Par)
      implicit none
      integer(4), intent(in) :: numElems2Par,mpiRank,numRanks2Par
      integer(4), intent(out) :: iElemStart,iElemEnd,iElemsInRank
      integer(4), dimension(0:numRanks2Par-1) :: vecElemsInMpiRank
      integer(4) :: iMpiRank

      if(numElems2Par.lt.numRanks2Par) then
         write(*,*) 'ERROR! The total number of elements to par:',numElems2Par,'is bigger than the number of cpus used to do the partition',numRanks2Par,&
                     'this is CRAZY BRO! The tool is going to crash! Use a reasonable number of CPUs (smaller than num of elements)'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      if(mpiRank.ge.numRanks2Par) then
         write(*,*) 'ERROR! The mpiRank',mpiRank,'cannot be bigger than numRanks2Par',numRanks2Par,&
                     'this is CRAZY BRO! The tool is going to crash! Bye!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      call distribution_algorithm(numElems2Par,numRanks2Par,vecElemsInMpiRank)

      iElemsInRank = vecElemsInMpiRank(mpiRank)
      iElemStart=1
      do iMpiRank=0,(mpiRank-1) !find the iElemStart
         iElemStart = iElemStart + vecElemsInMpiRank(iMpiRank)
      end do
      iElemEnd = iElemStart + (iElemsInRank - 1)

      !write(*,*) '#rank ',mpiRank,' iElemsInRank ',iElemsInRank,' iElemS ',iElemStart,' iElemE ',iElemEnd

   end subroutine do_element_partitioning_serial

   subroutine distribute_ranks2Part_in_mpiRank(numMshRanks2Part,iMshRankStart,iMshRankEnd,numRanksMpiRank,maxNumRanksMpiRank,ranksMpiRank,mapRanksToMpiRank)
      implicit none
      integer, intent(in) :: numMshRanks2Part
      integer, intent(out) :: iMshRankStart,iMshRankEnd,numRanksMpiRank,maxNumRanksMpiRank
      integer, allocatable, intent(inout) :: ranksMpiRank(:),mapRanksToMpiRank(:)
      integer, dimension(0:mpi_size-1) :: vecRanksMpiRank
      integer :: iMpiRank,iMshRankCnt,ii

      maxNumRanksMpiRank=1
      if(numMshRanks2Part.ge.mpi_size) then
         !---- number of elements and range this process will write
         call distribution_algorithm(numMshRanks2Part,mpi_size,vecRanksMpiRank)

         numRanksMpiRank = vecRanksMpiRank(mpi_rank)
         iMshRankStart=0
         do iMpiRank=0,(mpi_rank-1) !find the iMshRankStart
            iMshRankStart = iMshRankStart + vecRanksMpiRank(iMpiRank)
         end do
         iMshRankEnd = iMshRankStart + (numRanksMpiRank - 1)

         allocate(ranksMpiRank(numRanksMpiRank))
         do iMpiRank=1,numRanksMpiRank
            ranksMpiRank(iMpiRank)= iMshRankStart + (iMpiRank-1)
         end do

         allocate(mapRanksToMpiRank(numMshRanks2Part))
         iMshRankCnt=0
         do iMpiRank=0,(mpi_size-1)
            do ii=1,vecRanksMpiRank(iMpiRank)
               iMshRankCnt=iMshRankCnt+1
               mapRanksToMpiRank(iMshRankCnt)=iMpiRank
            end do
            maxNumRanksMpiRank = max(maxNumRanksMpiRank,vecRanksMpiRank(iMpiRank))
         end do

      else
         !possible first approach...
         !to use ONLY the number of CPUs=numMshRanks2Part and the others do nothing...
         !lets try

         if(mpi_rank.eq.0) then
            write(*,*) 'Tool not optimized for this case numMshRanks2Part<mpi_size -> ',numMshRanks2Part,'<',mpi_size
            write(*,*) 'Maybe it is better to run it with numMshRanks2Part=mpi_size... just suggesting'
         end if

         if(mpi_rank.lt.numMshRanks2Part) then
            numRanksMpiRank = 1
            iMshRankStart   = mpi_rank
            iMshRankEnd     = mpi_rank

            allocate(ranksMpiRank(numRanksMpiRank))
            ranksMpiRank(1) = mpi_rank
         else
            numRanksMpiRank = 0
            iMshRankStart   = 0
            iMshRankEnd     = 0

            allocate(ranksMpiRank(numRanksMpiRank))
         endif

         allocate(mapRanksToMpiRank(numMshRanks2Part))
         do iMpiRank=0,(numMshRanks2Part-1)
            mapRanksToMpiRank(iMpiRank+1)=iMpiRank
         end do

      end if
      !write(*,*) '#rank',mpi_rank,'numRanksMpiRank',numRanksMpiRank,'ranksInP',ranksMpiRank(:),'mapRanksToMpiRank',mapRanksToMpiRank(:)

   end subroutine distribute_ranks2Part_in_mpiRank

   subroutine define_parallelNodePartitioning(numMshRanks2Part,numMshRanksInMpiRank,numNodesMshRank,mshRanksInMpiRank,mapMshRankToMpiRank,mshRankNodeStart_i8,mshRankNodeEnd_i8,iNodeStartPar_i8,numNodesParTotal_i8)
      implicit none
      integer(4),intent(in) :: numMshRanks2Part
      integer(4),intent(in) :: numMshRanksInMpiRank,numNodesMshRank(numMshRanksInMpiRank),mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part)
      integer(8),intent(out) :: mshRankNodeStart_i8(numMshRanksInMpiRank),mshRankNodeEnd_i8(numMshRanksInMpiRank)
      integer(8),dimension(0:numMshRanks2Part-1),intent(out) :: iNodeStartPar_i8
      integer(8),intent(out) :: numNodesParTotal_i8
      integer(8),dimension(0:numMshRanks2Part-1) :: iNodeEndPar_i8
      integer(4),dimension(0:numMshRanks2Part-1) :: vecNumNodesMshRank
      integer(4) :: window_id
      integer(4) :: iMshRank,mshRank,mpiRank
      integer(8) :: auxNodeCnt_i8

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numNodesMshRank,vecNumNodesMshRank)
      !write(*,*) 'rank[',mpi_rank,']vNNR2s', vectorNumNodesRankPar2send(:),'vNNR2r',vectorNumNodesRankPar2rcv(:),'nNrp',numNodesMshRank(:)

      auxNodeCnt_i8 = 1
      do iMshRank=0,numMshRanks2Part-1
         iNodeStartPar_i8(iMshRank) = auxNodeCnt_i8
         iNodeEndPar_i8(iMshRank)   = iNodeStartPar_i8(iMshRank) + vecNumNodesMshRank(iMshRank) - 1
         auxNodeCnt_i8 = iNodeEndPar_i8(iMshRank) + 1
      end do

      numNodesParTotal_i8 = iNodeEndPar_i8(numMshRanks2Part-1)

      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         mshRankNodeStart_i8(iMshRank) = iNodeStartPar_i8(mshRank)
         mshRankNodeEnd_i8(iMshRank)   = iNodeEndPar_i8(mshRank)
         !write(*,*) 'mshRank',mshRank,'[',mpi_rank,']nodesInRank',numNodesMshRank(iMshRank),'iNodeS',mshRankNodeStart_i8(iMshRank),'iNodeE',mshRankNodeEnd_i8(iMshRank),'start',iNodeStartPar_i8(:),'end',iNodeEndPar_i8(:),'tNNP',numNodesParTotal_i8
         !write(*,*) 'mshRank',mshRank,'[',mpi_rank,']start',iNodeStartPar_i8(:),'end',iNodeEndPar_i8(:),'tNNP',numNodesParTotal_i8
      end do
   end subroutine define_parallelNodePartitioning

   subroutine get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,int2comm,vectorOut)
      implicit none
      integer(4),intent(in) :: numMshRanks2Part
      integer(4),intent(in) :: numMshRanksInMpiRank,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part),int2comm(numMshRanksInMpiRank)
      integer(4),dimension(0:numMshRanks2Part-1),intent(out) :: vectorOut

      integer(4),dimension(0:numMshRanks2Part-1) :: auxVector2send
      integer(4) :: window_id
      integer(4) :: iMshRank,mshRank,mpiRank,auxNodeCnt

      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      auxVector2send(:)=0
      vectorOut(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         auxVector2send(mshRank) = int2comm(iMshRank)
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(auxVector2send(0),win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0

      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vectorOut(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------

   end subroutine get_vector_with_mshRank_values_for_numMshRanks2Part


end module mod_partition_utils