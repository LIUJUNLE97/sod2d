module mod_mpi_mesh
   use mod_constants
   use mod_mpi
   use mod_utils
   use iso_c_binding
   implicit none
!-----------------------------------   
#define _CHECK_ 0
#define int_size 4
!-----------------------------------   
#ifndef GEMPAINTERFACE
   interface
      subroutine gempa_do_partition(numElemsInRank,numRanksToPart,x,y,z,weights,part) bind(c)
        import c_int, c_double
        integer(kind=c_int), intent(in), value :: numElemsInRank
        integer(kind=c_int), intent(in), value :: numRanksToPart
        real(kind=c_double), intent(in) :: x(numElemsInRank)
        real(kind=c_double), intent(in) :: y(numElemsInRank)
        real(kind=c_double), intent(in) :: z(numElemsInRank)
        integer(kind=c_int), intent(in) :: weights(numElemsInRank)
        integer(kind=c_int), intent(out) :: part(numElemsInRank)
      end subroutine gempa_do_partition
   end interface
#endif
!-----------------------------------   
   character(*), parameter :: fmt_csv = '(1x,*(g0,","))'
   ! Dimensions -----------------------------
   !integer(int_size), parameter :: ndime=3
   !-----------------------------------------
   ! Element characteristics ----------------
   !integer(int_size), parameter :: nnode=64
   !integer(int_size), parameter :: porder=3
   !integer(int_size), parameter :: npbou=16
   !integer(int_size), parameter :: ngaus=64
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
integer(int_size), allocatable :: connecGMSH(:,:), boundGMSH(:,:), bou_codesGMSH(:,:)
real(rp), allocatable :: coordGMSH(:,:)

! ################################################################################################
! ----------------- VARS for new Par mesh FORMAT -------------------------------------------------
! ################################################################################################

integer(int_size) :: numNodesRankPar, totalNumNodesPar, totalNumNodesSrl
integer(int_size) :: totalNumElements
integer(int_size) :: numElemsInRank, rankElemStart, rankElemEnd
integer(int_size) :: rankNodeStart, rankNodeEnd

integer(int_size), allocatable :: globalIdSrl(:), globalIdPar(:), elemGid(:)

real(rp), allocatable :: coordPar(:,:)

integer(int_size), allocatable :: connecCGNS(:)
integer(int_size), allocatable :: connecParOrig(:,:),connecParWork(:,:)

integer(int_size), allocatable :: workingNodesPar(:)
integer(int_size) :: numWorkingNodesRankPar

integer(int_size), allocatable :: masSlaSrl(:,:),masSlaRankPar(:,:)
integer(int_size) :: nPerSrl,nPerRankPar
logical :: isMeshPeriodic

integer(int_size) :: numBoundCodes, numBoundsRankPar, totalNumBoundsSrl
integer(int_size) :: ndofRankPar, numBoundaryNodesRankPar
integer(int_size), allocatable :: boundPar(:,:), bouCodesPar(:), ldofPar(:), lbnodesPar(:)
real(rp), allocatable :: boundNormalPar(:,:)
logical :: isMeshBoundaries

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

! ################################################################################################
! --------------------------- FOR READING ALYA MESH ----------------------------------------------
! ################################################################################################
   subroutine read_alya_mesh_files(file_path,file_name,isPeriodic)
      implicit none     
      character(500), intent(in)  :: file_path, file_name
      logical, optional :: isPeriodic
      integer(4) :: nelem,npoin,nboun,nper

      isMeshPeriodic = .false.
      if(present(isPeriodic)) then
         isMeshPeriodic = isPeriodic
      end if

      call read_dims_file(file_path,file_name,npoin,nelem,nboun)

      allocate(connecGMSH(nelem,nnode))
      allocate(coordGMSH(npoin,ndime))

      isMeshBoundaries = .false.
      if (nboun .ne. 0) then
         allocate(boundGMSH(nboun,npbou))
         allocate(bou_codesGMSH(nboun,2))

         isMeshBoundaries = .true.
         call read_fixbou_file(file_path,file_name,nboun,numBoundCodes,bou_codesGMSH)
      end if

      call read_geo_dat_file(file_path,file_name,npoin,nelem,nboun)
 
      if(isMeshPeriodic) then 
         call read_periodic_file(file_path,file_name,nper,masSlaSrl)

         !allocate some shit...
         !maybe we can do it in the func periodic_ops of mod_period, it would make sense...
         !allocate(connec_orig(this%nelem,nnode))
         !if (this%nboun .ne. 0) then
         !   allocate(bound_orig(this%nboun,npbou))
         !end if

      end if

      totalNumElements = nelem
      totalNumNodesSrl = npoin
      totalNumBoundsSrl = nboun
      nPerSrl = nper

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
      character(500), intent(in) :: file_path, file_name
      integer(4)    , intent(in) :: npoin, nelem, nboun
      integer(4)                 :: iline, int1, inode, idime, aux(nnode+1), bou_aux(npbou+1)
      character(2000)            :: file_type, line
      real(8)                    :: x,y,z

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
            if(rp == 4) then
               read(99,*) int1, x, y, z
               coordGMSH(iline,1)=real(x,rp)
               coordGMSH(iline,2)=real(y,rp) 
               coordGMSH(iline,3)=real(z,rp)
            else
               read(99,*) int1, coordGMSH(iline,1), coordGMSH(iline,2), coordGMSH(iline,3)
            end if
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

   subroutine read_periodic_file(file_path,file_name,nper,masSla)
      implicit none
      character(500), intent(in) :: file_path, file_name
      integer(4), intent(out) :: nper
      integer(4), intent(out), allocatable :: masSla(:,:)
      integer(4)                :: i,nlines,io_val
      character(500)            :: file_type, line

      if(mpi_rank.eq.0) write(*,*) "--| READING PERIODICITY FILE..."

      write(file_type,*) ".per"
      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")

      nlines = 0
      do
         read(99,*,iostat=io_val)
         if(io_val.ne.0) exit
         nlines = nlines + 1
      end do
      nper = nlines

      allocate(masSla(nper,2))

      rewind 99
      do i=1,nper
         read(99,*) masSla(i,1), masSla(i,2)
      end do

      if(mpi_rank.eq.0) write(*,*) "--| TOTAL PERIODIC NODES : ",nper

   end subroutine read_periodic_file

   subroutine read_fixbou_file(file_path,file_name,nboun,nbcodes,bou_codes)
      implicit none

      character(500), intent(in)  :: file_path, file_name
      integer(4), intent(in)      :: nboun
      integer(4), intent(inout)   :: nbcodes, bou_codes(nboun,2)
      integer(4)                  :: iboun, ii
      character(500)              :: file_type, line
      
      write(file_type,*) ".fix.bou"
      
      if(mpi_rank.eq.0) write(*,*) "--| READING FIXBOU FILE..."
      open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")

      read(99,*) ! Header
      do iboun = 1,nboun
         read(99,'(a)') line
         read(line,*) (bou_codes(iboun,ii), ii=1,2)
      end do
      nbcodes = maxval(bou_codes(:,2))
      if(mpi_rank.eq.0) write(*,*) "--| TOTAL BOUNDARY CODES : ",nbcodes

   end subroutine read_fixbou_file

! ################################################################################################
! ---------------------------- CGNS MESH PARTITIONING --------------------------------------------
! ################################################################################################

   subroutine do_mesh_partitioning()
      implicit none

      call do_element_partitioning_gempa()

      call do_node_partitioning_and_connectivity()

      call create_nodesCoordinates()

      if(isMeshPeriodic) then
         call create_masSla_parallel()
      end if

      call create_working_lists() !pot anar aqui no? Si, pero oju que potser em canvien altres merdes!

      if(totalNumBoundsSrl.ne.0) then
         call splitBoundary_inPar()
      end if

   end subroutine do_mesh_partitioning

   subroutine do_element_partitioning_serial(numElems2Par,iElemStart,iElemEnd,iElemsInRank)
      implicit none
      integer, intent(in) :: numElems2Par
      integer, intent(out) :: iElemStart,iElemEnd,iElemsInRank

      !---- number of elements and range this process will write
      iElemsInRank = (numElems2Par + mpi_size - 1) / mpi_size
      iElemStart = iElemsInRank * mpi_rank + 1
      iElemEnd   = iElemsInRank * (mpi_rank + 1)
      if (iElemEnd .gt. numElems2Par) iElemEnd = numElems2Par
      iElemsInRank = iElemEnd - iElemStart + 1

      !write(*,*) '#rank ',mpi_rank,' elemsInRank ',numElemsInRank,' iElemS ',rankElemStart,' iElemE ',rankElemEnd
      !write(*,*) '#rank ',mpi_rank,' totalNumElements ',totalNumElements,' totalNumNodesSrl ',totalNumNodesSrl

   end subroutine do_element_partitioning_serial

   recursive subroutine find_linked_elems(iRank,linkedElems,iElemG,unfoldedElems,numAddElemsInRank)
      implicit none
      integer, intent(in) :: iRank,linkedElems(:,:),iElemG
      integer, intent(inout) :: unfoldedElems(:),numAddElemsInRank
      integer :: numLinkedElems
      integer :: i,iEL1,iEL2
         
      numLinkedElems = size(linkedElems(:,1))

      do i=1,numLinkedElems
         iEL1 = linkedElems(i,1)
         if(iElemG.eq.iEL1) then   !#checking if the elem is in the linkedElems list
            iEL2 = linkedElems(i,2)
            !write(*,*) 'iEL1 ',iEL1, ' linked to iEL2 ', iEL2
            
            if(unfoldedElems(iEL2).le.0) then
               unfoldedElems(iEL2) = iRank
               numAddElemsInRank=numAddElemsInRank+1
               call find_linked_elems(iRank,linkedElems,iEL2,unfoldedElems,numAddElemsInRank)
            end if        
         end if
      end do

      numLinkedElems = size(linkedElems(:,1))

   end subroutine find_linked_elems

   subroutine do_element_partitioning_gempa()
      implicit none
      integer, parameter :: nodesToAvg(8) = [1,2,5,6,17,18,21,22]
      
      real(8), allocatable, dimension(:) :: x,y,z
      integer, allocatable :: listElems2Par(:),weightElems2Par(:),unfoldedElems(:)
      integer, allocatable :: linkedElems(:,:),elemPart(:,:),elemPartAux(:,:)

      integer :: i,j,k,m,iElem,iElemG,iEL1,iEL2,iEL2rep,iEL3,iNodeG,iRank,iPos
      integer :: numElemsInRank_par,numElems2get
      integer :: iElemStart,iElemEnd,iElemsInRank,numElems2Par,numLinkedElems,numAdditionalElemsInRank
      real(8) :: x_a,y_a,z_a
      integer,dimension(0:mpi_size-1) :: vecNumElemsRank
      integer,dimension(0:mpi_size-1,0:mpi_size-1) :: matNumElemsRank
      integer :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: window_buffer_size,target_displacement

      character(128) :: file_name, aux_string_rank

      if(mpi_rank.eq.0) write(*,*) ' # Doing element partitioning GEMPA...'

      !1. obtain the list element 2 part using gempa
      if(isMeshPeriodic) then
         call get_listElems2Par_Periodic(listElems2Par,weightElems2Par,linkedElems)
      else
         call get_listElems2Par(listElems2Par,weightElems2Par)
      end if
      
      numElems2Par = size(listElems2Par)
      !write(*,*) 'numElems2Par ',numElems2Par
      
      call do_element_partitioning_serial(numElems2Par,iElemStart,iElemEnd,iElemsInRank)

      allocate(elemPart(iElemsInRank,3))
      allocate(x(iElemsInRank))
      allocate(y(iElemsInRank))
      allocate(z(iElemsInRank))

      !numElemsInRank_srl = numElemsInRank

      i=1
      do iElem=iElemStart,iElemEnd
         iElemG=listElems2Par(iElem)
         elemPart(i,1) = iElemG 
         elemPart(i,3) = weightElems2Par(iElem)
         !write(*,*) '[',mpi_rank,'] iElem ',iElem,' iElemG ',iElemG,' weight ', elemPart(i,3)
         x_a=0.
         y_a=0.
         z_a=0.
         do j=1,8
            m = gmsh2ijk(nodesToAvg(j))
            iNodeG = connecGMSH(iElemG,m) !de moment poso el primer, despres ja fare avg
            x_a = x_a + coordGMSH(iNodeG,1)
            y_a = y_a + coordGMSH(iNodeG,2)
            z_a = z_a + coordGMSH(iNodeG,3)
         end do

         x(i) = x_a/8.d0
         y(i) = y_a/8.d0
         z(i) = z_a/8.d0
         i=i+1
      end do

      !numCoords=iElemsInRank !numElemsInRank_srl

      !----------------------------------------------------------------------------------
      !@TODO: en el cas de malles periodiques, aqui haurem de fer un trucu del almendrucu
      !       per associar els elements que siguin periodics
      !       un cop fet el particionament amb els elements associats, desfer aquesta
      !       'assosiacio' ficticia i posarlos on toca
      !----------------------------------------------------------------------------------

#ifndef GEMPAINTERFACE
      !---- CALLING GEMPA in PARALLEL-----------------------
      call gempa_do_partition(iElemsInRank,mpi_size,x,y,z,elemPart(:,3),elemPart(:,2))
      !---------------------------------------------------------------------
#endif

      ! if parallel, now we have to put the rank to the 'slave' elements who were not included in the partitioning process
      ! since they were linked to a 'master' element

      if(isMeshPeriodic) then
         numAdditionalElemsInRank=0
         allocate(unfoldedElems(totalNumElements))
         unfoldedElems(:)=-1
#if 1
         do iElem=1,iElemsInRank
            iElemG = elemPart(iElem,1)
            iRank = elemPart(iElem,2)
            call find_linked_elems(iRank,linkedElems,iElemG,unfoldedElems,numAdditionalElemsInRank)
         end do
         !do iElemG=1,size(unfoldedElems)
         !   iRank = unfoldedElems(iElemG)
         !   if(iRank.ge.0) then
         !      write(*,*) '[',mpi_rank,'] unfolded Elem ',iElemG,' rank ', iRank
         !   end if
         !end do
         !write(*,*) 'unfolded ',unfoldedElems(:)
#else
         numLinkedElems = size(linkedElems(:,1))
         !write(*,*) 'numLinkedElems ',numLinkedElems
         do iElem=1,iElemsInRank
            iElemG = elemPart(iElem,1)
            iRank = elemPart(iElem,2)
            do j=1,numLinkedElems 
               iEL1 = linkedElems(j,1)
               if(iElemG.eq.iEL1) then   !#checking if the elem is in the linkedElems list
                  iEL2 = linkedElems(j,2)
                  write(*,*) 'iEL1 ',iEL1, ' linked to iEL2 ', iEL2
                  
                  if(unfoldedElems(iEL2).le.0) then
                     unfoldedElems(iEL2) = iRank
                     numAdditionalElemsInRank=numAdditionalElemsInRank+1
                  end if

                  do k=1,numLinkedElems
                     iEL2rep = linkedElems(k,1)
                     if(iEL2.eq.iEL2rep) then 
                        iEL3 = linkedElems(k,2)
                        write(*,*) 'iEL2 ',iEL2, ' linked to iEL3 ', iEL3
                        if(unfoldedElems(iEL3).le.0) then
                           unfoldedElems(iEL3) = iRank
                           numAdditionalElemsInRank=numAdditionalElemsInRank+1
                        end if
                     end if
                  end do
               end if
               !if(iElemG.eq.iEL2) write(*,*) 'iEL2 ',iEL2, ' linked to iEL1 ', iEL1
            end do
         end do
#endif

         !now we have the 'unfolded elems' and the rank that they belong!
         !we need to add it to the elemPart matrix(:,2)!!!
         !write(*,*) 'numAdditional ',numAdditionalElemsInRank
         allocate(elemPartAux(iElemsInRank,2))
         !$acc kernels
         elemPartAux(:,1)=elemPart(:,1)
         elemPartAux(:,2)=elemPart(:,2)
         !$acc end kernels

         deallocate(elemPart)
         iElemsInRank = iElemsInRank + numAdditionalElemsInRank
         allocate(elemPart(iElemsInRank,2))

         i=1
         !!!$acc kernels
         do j=1,size(elemPartAux(:,1))
            elemPart(i,:)=elemPartAux(i,:)
            !!!$acc atomic update
            i=i+1
            !!!$acc end atomic
         end do
         !!!$acc end kernels

         do iElemG=1,size(unfoldedElems)
            iRank = unfoldedElems(iElemG)
            if(iRank.ge.0) then
               elemPart(i,1) = iElemG
               elemPart(i,2) = iRank
               !write(*,*) '[',mpi_rank,'] unf.Elem ',iElemG,' rank ', iRank,' i ',i
               i=i+1
            end if
         end do

         deallocate(unfoldedElems)
         deallocate(linkedElems)
         deallocate(elemPartAux)
      end if

      ! @now we have to info all the ranks about which are the elements that they will 'store'
      ! since we have done the partitioning in paralel...

      vecNumElemsRank(:)=0
      !!!$acc parallel loop
      do i=1,iElemsInRank
         iRank = elemPart(i,2)-1 !elemPart came with rank=1:mpi_size
         !!!$acc atomic update
         vecNumElemsRank(iRank) = vecNumElemsRank(iRank) + 1
         !!!$acc end atomic
      end do
      !!!$acc end parallel loop

#if _CHECK_
      write(aux_string_rank,'(I0)') mpi_rank
      file_name = 'elemPartition_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      write(1,*) 'X,Y,Z,iElemG,rank'
      do i=1,iElemsInRank ! numElemsInRank_srl
         iElemG=elemPart(i,1)
         x_a=0.
         y_a=0.
         z_a=0.
         do j=1,8
            m = gmsh2ijk(nodesToAvg(j))
            iNodeG = connecGMSH(iElemG,m) !de moment poso el primer, despres ja fare avg
            x_a = x_a + coordGMSH(iNodeG,1)
            y_a = y_a + coordGMSH(iNodeG,2)
            z_a = z_a + coordGMSH(iNodeG,3)
         end do
         x_a = x_a/8.d0
         y_a = y_a/8.d0
         z_a = z_a/8.d0

         !write(1,fmt_csv) x(i),y(i),z(i),elemPart(i,1),elemPart(i,2)
         write(1,fmt_csv) x_a,y_a,z_a,iElemG,elemPart(i,2)
      end do
      close(1)
#endif
      !write(*,*) 'vecNumElemsRank[',mpi_rank,'] ',vecNumElemsRank(:)

      !una finestra on comunicar el nou mpi_rank de cada node
      !es una primera fora de treballar, es pot mirar d'optimitzar per no requerir un vector tamany totalNumElements

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

      !write(*,*) 'numElemsInRankP(,',mpi_rank,')->',numElemsInRank_par!,' srl ',numElemsInRank_srl

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
      window_buffer_size = mpi_integer_size*iElemsInRank ! numElemsInRank_srl

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
      deallocate(listElems2Par)
      deallocate(elemPart)
      deallocate(x)
      deallocate(y)
      deallocate(z)
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
         write(1,fmt_csv) i,',',elemGid(i)
      end do
      close(1)

#endif
!------------------------------------------------------------------------------

   end subroutine do_element_partitioning_gempa

   subroutine get_listElems2Par(listElems2Par,weightElems2Par)
      implicit none
      integer, allocatable, intent(out) :: listElems2Par(:),weightElems2Par(:)
      integer iElem

      allocate(listElems2Par(totalNumElements))
      allocate(weightElems2Par(totalNumElements))
      !$acc parallel loop
      do iElem=1,totalNumElements
         listElems2Par(iElem) = iElem
         weightElems2Par(iElem) = 1
      end do
      !$acc end parallel loop

   end subroutine get_listElems2Par

   subroutine get_listElems2Par_Periodic(listElems2Par,weightElems2Par,linkedElems)
      implicit none
      integer, allocatable, intent(out) :: listElems2Par(:),weightElems2Par(:)
      integer, allocatable, intent(out) :: linkedElems(:,:)
      integer, parameter :: nodesToPer(6) = [11,27,36,40,41,42]
      integer, allocatable :: linkerNodesAll(:),listAllElems(:),listElems2NotPar(:)
      integer :: numLinkerNodes,numLinkedElems,numElems2Par,numElems2NotPar
      integer :: iElem,iElemG,i,j,m,iNodeGSrl,weightCnt
      integer :: iNode1,iNode2,iElem1,iElem2

      allocate(linkerNodesAll(totalNumNodesSrl))

      !$acc kernels
      linkerNodesAll(:) = -1
      !$acc end kernels

      numLinkerNodes=0
      do iElem=1,totalNumElements
         do j=1,6
            m = gmsh2ijk(nodesToPer(j))
            iNodeGSrl = connecGMSH(iElem,m) !de moment poso el primer, despres ja fare avg

            if(linkerNodesAll(iNodeGSrl).lt.0) then
               linkerNodesAll(iNodeGSrl) = iElem
               numLinkerNodes=numLinkerNodes+1
            else !the linker node is inner, not boundary
               linkerNodesAll(iNodeGSrl) = -2
               numLinkerNodes=numLinkerNodes-1
            end if
         end do
      end do

      numLinkedElems = numLinkerNodes/2
      !nPer = size(masSla(:,1))
      !write(*,*) 'nPer ',nPer

      allocate(linkedElems(numLinkedElems,2))
      allocate(listAllElems(totalNumElements))
      !$acc kernels
      listAllElems(:)=1
      !$acc end kernels
      numElems2Par=totalNumElements
      numElems2NotPar=0

      j=1
      do i=1,nPerSrl
         iNode1 = masSlaSrl(i,1)
         iNode2 = masSlaSrl(i,2)
         iElem1 = linkerNodesAll(iNode1)
         iElem2 = linkerNodesAll(iNode2)
         if(iElem1.ge.0) then
            linkedElems(j,1) = iElem1
            linkedElems(j,2) = iElem2
            j=j+1

            if(listAllElems(iElem2).ge.0) then
               listAllElems(iElem2) = -1
               numElems2Par = numElems2Par-1
               numElems2NotPar = numElems2NotPar+1
            end if
            !if(mpi_rank.eq.0) write(*,*) 'iElem1 ',iElem1,' iElem2 ',iElem2
         end if
         !if(linkerNodesAll(iNodeGSrl).gt.0) then
         !   linkerNodes(j)=linkerNodesAll(iNodeGSrl)
      end do

      !write(*,*) 'numElems2Par: ',numElems2Par,' numElem2NotPar ',numElems2NotPar
      allocate(listElems2Par(numElems2Par))
      allocate(weightElems2Par(numElems2Par))
      allocate(listElems2NotPar(numElems2NotPar))

      !$acc kernels
      listElems2Par(:)=-1
      listElems2NotPar(:)=-1
      !$acc end kernels
      i=1
      j=1
      do iElem=1,totalNumElements
         if(listAllElems(iElem).ge.0) then
            listElems2Par(i) = iElem
            weightElems2Par(i) = 1
            i=i+1
         else
            listElems2NotPar(j) = iElem
            j=j+1
         end if
      end do

      do iElem=1,numElems2Par
         iElemG = listElems2Par(iElem)
         weightCnt=1
         do i=1,numLinkedElems
            if(linkedElems(i,1).eq.iElemG) weightCnt=weightCnt+1
         end do
         weightElems2Par(iElem) = weightCnt
         !write(*,*) 'iElemG ',iElemG,' weight ',weightCnt
      end do

      !give more weight to the elements that account for more nodes
      !---------------------------------------------------------------------

      !---------------------------------------------------------------------
      !write(*,*) 'numElems2Par: ',numElems2Par,' i ',i,' j ',j
      !if(mpi_rank.eq.0) write(*,*) 'listElems2Par ',listElems2Par
      !if(mpi_rank.eq.0) write(*,*) 'listElemsNot2Par ',listElems2NotPar

      deallocate(linkerNodesAll)
      deallocate(listAllElems)

      !write(*,*) 'numLinkerNodes ',numLinkerNodes

   end subroutine get_listElems2Par_Periodic

   subroutine do_node_partitioning_and_connectivity()
      implicit none
      integer, dimension(0:mpi_size-1) :: iNodeStartPar, iNodeEndPar
      integer, allocatable :: boundaryNodes(:)
      integer, allocatable :: vecSharedBN_full(:)

      if(mpi_rank.eq.0) write(*,*) ' # Doing node partitioning...'

      call get_rankPartitionBoundaryNodes(boundaryNodes)

      call define_parallelNodePartitioning(iNodeStartPar,iNodeEndPar)

      call define_mpi_boundaries_inPar(boundaryNodes,vecSharedBN_full)

      deallocate(boundaryNodes)

      call reorder_nodes_in_proc(iNodeStartPar)

      !aixo ho podria ficar en un altre modul... crec que quedaria mes 'endreçat'
      call generate_mpi_comm_scheme(vecSharedBN_full)

      deallocate(vecSharedBN_full)

   end subroutine do_node_partitioning_and_connectivity

   subroutine create_nodesCoordinates()
      implicit none
      integer :: iNodeL,iNodeGsrl

      if(mpi_rank.eq.0) write(*,*) ' # Creating nodes coordinates...'
      !------------------------------------------------------
      allocate(coordPar(numNodesRankPar,ndime)) !only works for ndime=3

      !- create the coordinate data for this process
      !$acc parallel loop
      do iNodeL=1,numNodesRankPar
          iNodeGSrl=globalIdSrl(iNodeL)
          !write (*,*) 'iL ', iNodeL, ' iG ', iNodeG,' [', coord(iNodeG,1),']',' [', coord(iNodeG,2),']'
          coordPar(iNodeL,1) = coordGMSH(iNodeGSrl,1)
          coordPar(iNodeL,2) = coordGMSH(iNodeGSrl,2)
          coordPar(iNodeL,3) = coordGMSH(iNodeGSrl,3)
      end do
      !$acc end parallel loop

   end subroutine create_nodesCoordinates

   subroutine create_masSla_parallel()
      implicit none
      integer :: iNodeL,iNodeL1,iNodeL2,iNodeG,iPer,iPerPar

      if(mpi_rank.eq.0) write(*,*) ' # Creting master-slave rels in parallel...'

      nPerRankPar=0
      !$acc parallel loop reduction(+:nPerRankPar)
      do iNodeL = 1,numNodesRankPar
         iNodeG = globalIdSrl(iNodeL)
         !$acc looop vector
         do iPer = 1,nPerSrl
            if (iNodeG .eq. masSlaSrl(iPer,2)) then
               nPerRankPar=nPerRankPar+1
            end if
         end do
      end do
      !$acc end parallel loop

      !TRY TO ACC THIS PART OF THE CODE: MAIN ISSUE WITH gidSrl_to_lid func in device!!
      allocate(masSlaRankPar(nPerRankPar,2))
      iPerPar=0
      do iNodeL = 1,numNodesRankPar
         iNodeG = globalIdSrl(iNodeL)
         do iPer = 1,nPerSrl
            if (iNodeG .eq. masSlaSrl(iPer,2)) then
               iPerPar=iPerPar+1
               !masSlaRankPar(iPerPar,:) = masSlaSrl(iPer,:)
               iNodeL1 = gidSrl_to_lid(masSlaSrl(iPer,1))
               iNodeL2 = gidSrl_to_lid(masSlaSrl(iPer,2))
               masSlaRankPar(iPerPar,1) = iNodeL1
               masSlaRankPar(iPerPar,2) = iNodeL2
               !write(*,*) 'iPerPar',iPerPar,'masSla',masSlaRankPar(iPerPar,:)
            end if
         end do
      end do

      write(*,*) '[',mpi_rank,']nPerSrl',nPerSrl,'nPerRankPar',nPerRankPar,'iPerPar',iPerPar

   end subroutine create_masSla_parallel

   subroutine create_working_lists()
      implicit none
      integer :: iNodeL,iElem,iPer,iAux
      integer :: iNodeL_Per,iNodeL_Per_Pair
      integer, allocatable :: aux_workingNodesPar(:)

      if(mpi_rank.eq.0) write(*,*) ' # Creating working lists...'

      allocate( connecParWork(numElemsInRank,nnode) )
      !$acc kernels
      connecParWork(:,:) = connecParOrig(:,:)
      !$acc end kernels

      if(isMeshPeriodic) then !do all stuff in case mesh is periodic

         !----------------------------------------------------------------
         !-------------  CONNEC   ----------------------------------------

        !allocate( connecParOrig(numElemsInRank,nnode) )
        !connecParOrig(:,:) = connecPar(:,:)

         !$acc kernels
         do iElem = 1,numElemsInRank
            do iAux = 1,nnode
               iNodeL = connecParWork(iElem,iAux)
               !iNodeG = globalIdSrl(iNodeL)
               do iPer = 1,nPerRankPar
                  iNodeL_Per = masSlaRankPar(iPer,2)
                  if (iNodeL .eq. iNodeL_Per) then
                     iNodeL_Per_Pair = masSlaRankPar(iPer,1)
                     !iNodeL_Per_Pair = gidSrl_to_lid(iNodeG_Per_Pair)
                     connecParWork(iElem,iAux) = iNodeL_Per_Pair
                  end if
               end do
            end do
         end do
         !$acc end kernels

         !if(present(bound) .and. present(bound_orig)) then
         !   !!!$acc kernels
         !   bound_orig(:,:) = bound(:,:)
         !   !!!$acc end kernels
         !end if

         !if(present(bound)) then
         !   !!!$acc parallel loop gang vector_length(vecLength)
         !   do iboun = 1,nboun
         !      !!!$acc loop vector
         !      do ipbou = 1,npbou
         !         !!!$acc loop seq
         !         do iper = 1,nper
         !            if (bound(iboun,ipbou) .eq. masSla(iper,2)) then
         !               bound(iboun,ipbou) = masSla(iper,1)
         !            end if
         !         end do
         !      end do
         !   end do
         !   !!!$acc end parallel loop
         !end if

         numWorkingNodesRankPar = numNodesRankPar - nPerRankPar

         allocate(aux_workingNodesPar(numNodesRankPar))
         allocate(workingNodesPar(numWorkingNodesRankPar))

         !$acc parallel loop
         do iNodeL = 1,numNodesRankPar
            aux_workingNodesPar(iNodeL) = iNodeL
         end do
         !$acc end parallel loop

         do iPer = 1,nPerRankPar
            iNodeL_Per = masSlaRankPar(iPer,2)
            do iNodeL = 1,numNodesRankPar
               !iNodeG = globalIdSrl(iNodeL)
               if (iNodeL_Per .eq. iNodeL) then
                  aux_workingNodesPar(iNodeL) = 0
                  exit
               end if
            end do
         end do

         !TODO: GPU-PARALLELIZE THIS LOOP
         iAux = 0
         do iNodeL = 1,numNodesRankPar
            if (aux_workingNodesPar(iNodeL) .ne. 0) then
               iAux = iAux+1
               workingNodesPar(iAux) = aux_workingNodesPar(iNodeL)
               !write(*,*)'wNP(',iAux,')',workingNodesPar(iAux)
            end if
         end do

         deallocate(aux_workingNodesPar)
      else !non-periodic meshes

         numWorkingNodesRankPar = numNodesRankPar
         allocate(workingNodesPar(numWorkingNodesRankPar))

         !$acc parallel loop
         do iNodeL = 1,numWorkingNodesRankPar
            workingNodesPar(iNodeL) = iNodeL
         end do
         !$acc end parallel loop

      end if

   end subroutine create_working_lists

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
      write(aux_string_rank,'(I0)') mpi_rank
      file_name = 'boundaryNodes_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numBNodesRankPar
         iNodeG=boundaryNodes(i)
         write(1,fmt_csv) coordGMSH(iNodeG,1),coordGMSH(iNodeG,2),coordGMSH(iNodeG,3),iNodeG
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

      vectorNumNodesRankPar(:)=0

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

      !---------------------------------------------------------------
      ! TO CHECK
#if _CHECK_
      write(aux_string_rank,'(I0)') mpi_rank
      file_name = 'boundaryNodes_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do iNodeGSrl=1,size(vectorBN)
           write(1,fmt_csv) iNodeGSrl,vectorBN(iNodeGSrl)
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
      file_name = 'matrixBN_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
      do i=1,numNodesRankSrl
           iNodeGSrl = iNodeStartSrl(mpi_rank)-1 + i
           write(1,fmt_csv) i,iNodeGSrl,matrixBN(i,0),matrixBN(i,1)
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

         write(1,fmt_csv)iNodeGSrl,coordGMSH(iNodeGSrl,1),coordGMSH(iNodeGSrl,2),coordGMSH(iNodeGSrl,3),&
               auxVecDebug(:)

         !write(1,'(*(G0.7,:,","))')iNodeGSrl,coordGMSH(iNodeGSrl,1),coordGMSH(iNodeGSrl,2),coordGMSH(iNodeGSrl,3),&
         !      auxVecDebug(:)

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
      write(aux_string_rank,'(I0)') mpi_rank
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
      allocate( connecParOrig(numElemsInRank,nnode) )

      isNodeAdded=-1
      iPos = 1

      indConn = -1

      !@TODO SUPER VERY FUCKING IMPORTANT
      !--------------------------------------------------------------------------------------------------------
      ! now this only works using gmsh2ijk, because the old part of the code is using the atoijk set in 
      ! the subroutine set_hex64_list in mod_elem.f90 file
      ! atoijk is the same than in gmsh2ijk
      ! first do a first migration of the code and once it work, think how to allow different node ordering
      ! BUT FIX IT!!!!!!

      call generate_new_nodeOrder_and_connectivity(gmsh2ijk,auxNewOrderIndex,auxCGNSorder)
      !call generate_new_nodeOrder_and_connectivity(cgns2ijk,auxNewOrderIndex,auxCGNSorder)
      !call generate_new_nodeOrder_and_connectivity(dummy2ijk,auxNewOrderIndex,auxCGNSorder)

      !----------------------------------------------------------------------------------------------------------

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

            connecParOrig(iElemL,m) = iNodeL

            indexCGNS = auxCGNSorder(m)
            indConn = (iElemL-1)*nnode + indexCGNS
            
            connecCGNS(indConn) = iNodeGPar

         end do

         !write(*,*) '[',mpi_rank,']iElemG ',iElemG,' connecParOrig ',connecParOrig(iElemL,:)

      end do
   end subroutine reorder_nodes_in_proc

   subroutine generate_new_nodeOrder_and_connectivity(newOrderIJK,auxNewOrderIndex,auxCGNSorder)
      integer,intent(in)  :: newOrderIJK(:)
      integer,intent(out) :: auxNewOrderIndex(:),auxCGNSorder(:)
      integer :: i,j,k,indexIJK,indexNew,indexCGNS

      !!!$acc kernels
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
      !!!$acc end kernels
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
                  !if(numRanksCnt.ge.3) write(*,*) '[',mpi_rank,'] k ',k,',',matrixCommScheme(k,:)
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
      !write(*,*) '[',mpi_rank,']ranksToComm->',ranksToComm(:)
      !write(*,*) '[',mpi_rank,']commsMemPosInLoc->',commsMemPosInLoc(:)
      !write(*,*) '[',mpi_rank,']commsMemPosInNgb->',commsMemPosInNgb(:)
      !write(*,*) '[',mpi_rank,']commsMemSize->',commsMemSize(:)

      !ara podria obtenir els ranks locals de les altres cpus... pero cal?

#if _CHECK_
      write(aux_string_rank,'(I0)') mpi_rank
      file_name = 'matrixCommScheme_rank'// trim(aux_string_rank)//'.csv'
      open(1, file=file_name)
         do i=1,numNodesToComm
             write(1,fmt_csv) matrixCommScheme(i,:)
         end do
      close(1)
#endif

   end subroutine generate_mpi_comm_scheme
!----------------------------------------------------------------------------------------------------------

   subroutine splitBoundary_inPar()
      integer(4), allocatable    :: aux1(:)
      integer(4) :: ii,iNodeL,iNodeGSrl,iNodeGSrl_bound,iBound,iBoundL,ipbou,idof,ibnodes
      integer(4) :: auxBoundCnt,aux_sum_NBRP
      integer(4) :: aux_boundList(totalNumBoundsSrl)

      if(mpi_rank.eq.0) write(*,*) ' # Splitting boundary nodes from DoFs in parallel...'

      !first, generate boundPar(:,:) & bouCodesPar(:,:)

      !1. how many boundaries my rank have?
      numBoundsRankPar = 0
      aux_boundList(:) = 0

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

      write(*,*) '[',mpi_rank,']numBoundsRankPar',numBoundsRankPar,'totalNumBoundsSrl',totalNumBoundsSrl
      
      call MPI_Allreduce(numBoundsRankPar,aux_sum_NBRP,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,mpi_err)

      if(aux_sum_NBRP.ne.totalNumBoundsSrl) then
         write(*,*) 'ERROR IN splitBoundary_inPar()->aux_sum_NBRP',aux_sum_NBRP,' not equal to totalNumBoundsSrl',totalNumBoundsSrl
         call MPI_Abort(MPI_COMM_WORLD, -1, mpi_err)
      end if

      allocate(boundPar(numBoundsRankPar,npbou))
      allocate(bouCodesPar(numBoundsRankPar))
      
      ii=0
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

      !------------------------------------------------------------------------
      ! Fill aux1 with all nodes in order
      allocate(aux1(numNodesRankPar))
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

      deallocate(aux1)

   end subroutine splitBoundary_inPar



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
         !write(1,fmt_csv) coord_x(iNodeL),coord_y(iNodeL),coord_z(iNodeL),dfield(iNodeL)
         write(1,fmt_csv) coordPar(iNodeL,1),coordPar(iNodeL,2),coordPar(iNodeL,3),dfield(iNodeL)
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
         write(1,fmt_csv) coordPar(iNodeL,1),coordPar(iNodeL,2),coordPar(iNodeL,3),ffield(iNodeL)
      end do

      close(1)


   end subroutine print_csv_file_ffield

end module mod_mpi_mesh
