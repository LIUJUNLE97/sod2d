module mod_meshConversorTool
   use mod_constants
   use mod_mpi
   use mod_utils
   use mod_ijk_indices
   use elem_hex
   use quadrature_rules
   use mod_hdf5
   use iso_c_binding
   use mod_partition_utils
   use mod_mesh_quality
   implicit none

!-----------------------------------
#define _CHECK_ 0
!-----------------------------------
#ifndef GEMPAINTERFACE
   interface
      subroutine gempa_do_partition(numElemsMpiRankToPart,numMshRanksToPart,x,y,z,weights,part) bind(c)
        import c_int, c_double
        integer(kind=c_int), intent(in), value :: numElemsMpiRankToPart
        integer(kind=c_int), intent(in), value :: numMshRanksToPart
        real(kind=c_double), intent(in) :: x(numElemsMpiRankToPart)
        real(kind=c_double), intent(in) :: y(numElemsMpiRankToPart)
        real(kind=c_double), intent(in) :: z(numElemsMpiRankToPart)
        integer(kind=c_int), intent(in) :: weights(numElemsMpiRankToPart)
        integer(kind=c_int), intent(out) :: part(numElemsMpiRankToPart)
      end subroutine gempa_do_partition
   end interface
#endif
!-----------------------------------

! ################################################################################################
integer, parameter :: fileId_geo = 99, fileId_bou = 98, fileId_dims = 97
! ################################################################################################
!     FOR DEBUG
character(128) :: file_name_deb,aux_string_mpirank_deb,aux_string_mshrank_deb
character(*), parameter :: fmt_csv_deb = '(1x,*(g0,","))'
! ################################################################################################

contains

   subroutine read_gmsh_h5_file_and_do_partitioning_in_parallel(gmsh_filePath,gmsh_fileName,mesh_h5_filePath,mesh_h5_fileName,numMshRanks2Part,isLinealOutput,evalMeshQuality)
      implicit none
      character(len=*),intent(in)  :: gmsh_filePath,gmsh_fileName,mesh_h5_filePath,mesh_h5_fileName
      integer(4),intent(in) :: numMshRanks2Part
      logical,intent(in) :: isLinealOutput, evalMeshQuality
      logical :: isPeriodic=.false.,isBoundaries=.false.,isMapFaces=.false.

      character(1024) :: meshFile_h5_full_name

      integer(4) :: numElemsMpiRank,numNodesMpiRank,numBoundsMpiRank,numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4), allocatable :: listElemsMpiRank(:)
      integer(4), allocatable :: linkedPerElemsSrl(:,:),listPerElemsSrl(:)
      integer(4), allocatable :: linkedPerToMapElemsAndFacesSrl(:,:)
      integer(8), allocatable :: listNodesMpiRank_i8(:),connecMpiRank_i8(:,:),masSlaNodesSrl_i8(:,:)

      real(8), allocatable :: coordNodesMpiRank(:,:)
      integer(4),allocatable :: numBoundaryNodesMshRank(:),numMpiBoundaryNodes(:),numMshBoundaryNodes(:),numInnerNodes(:),numBoundElemsMshRank(:)
      integer(4),allocatable :: numDoFMshRank(:),numElemsVTKMshRank(:),sizeConnecVTKMshRank(:)

      type(jagged_vector_int4) :: listBoundElemsMshRank_jv,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,boundaryNodes_jv,dofNodes_jv
      type(jagged_vector_int8) :: listNodesMshRank_i8_jv,mshBoundaryNodes_i8_jv,mpiBoundaryNodes_i8_jv
      type(jagged_matrix_int4) :: listElemsBoundFacesMshRank_jm
      type(jagged_matrix_int8) :: connecMshRank_i8_jm,connecOrigBoundFacesMshRank_i8_jm,globalIdSrlOrdered_i8_jm,masSlaRankPar_i8_jm,perMapLinkedNodesRankPar_i8_jm
      type(jagged_matrix_int4) :: connecOrigBoundFacesMshRank_jm,connecBoundFacesMshRank_jm
      type(jagged_matrix_real8) :: coordMshRank_jm
      integer(8),dimension(0:numMshRanks2Part-1) :: iNodeStartPar_i8
      integer(4),dimension(0:numMshRanks2Part-1) :: vecNumWorkingNodes,vecNumNodesToCommMshRank,vecNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecBndNumMshRanksWithComms
      integer(4),dimension(0:numMshRanks2Part-1) :: vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank,vecNumPerMapLinkedNodesMshRank
      real(8),allocatable :: Ngp_l(:,:),dNgp(:,:,:),wgp(:)
      integer(4) :: iMshRank,mshRank
      integer(hid_t) :: gmsh_h5_fileId, sod2dmsh_h5_fileId
      real(8),dimension(10) :: start_time,end_time,elapsed_time_r

      integer(4) :: mporder,mnnode,mngaus,mnpbou,mnnodeVTK,numVTKElemsPerMshElem,connecChunkSize
      integer(4),dimension(:),allocatable :: gmsh2ijk,vtk2ijk,gmsh2ij,vtk2ij,a2ijk,a2ij
      integer(4),dimension(:),allocatable :: auxNewOrderIndex,auxVTKorder
      integer(4),dimension(:,:,:),allocatable :: ijk2a,ijk2gmsh,ijk2vtk

      ! ################################################################################################
      ! ----------------- VARS for new Par mesh FORMAT -------------------------------------------------
      ! ################################################################################################

      integer(4) :: numElemsGmsh,numBoundFacesGmsh,numPerFacesGmsh,numPerLinkedNodesGmsh,numMapFacesGmsh
      integer(8) :: numNodesGmsh_i8,numNodesParTotal_i8
      integer(4), allocatable :: numElemsMshRank(:),mshRankElemStart(:),mshRankElemEnd(:)
      integer(8), allocatable :: mshRankNodeStart_i8(:), mshRankNodeEnd_i8(:)
      integer(4), allocatable :: numBoundFacesMshRank(:),numNodesMshRank(:)
      integer(4) :: maxBoundCode

      type(jagged_vector_int4) :: elemGidMshRank_jv
      type(jagged_vector_int8) :: globalIdSrl_i8_jv,globalIdPar_i8_jv

      type(jagged_vector_int4) :: connecVTK_jv
      type(jagged_matrix_int4) :: connecParOrig_jm,connecParWork_jm

      type(jagged_matrix_real8) :: coordPar_jm,coordVTK_jm

      type(jagged_vector_int4) :: workingNodesPar_jv
      integer(4),allocatable :: numWorkingNodesMshRank(:)

      type(jagged_matrix_int4) :: masSlaRankPar_jm
      integer(4),allocatable :: numPerNodesMshRank(:)

      type(jagged_matrix_int4) :: perMapLinkedNodesRankPar_jm
      integer(4),allocatable :: numPerMapLinkedNodesMshRank(:)

      integer(4) :: mshRankInMpiRankStart,mshRankInMpiRankEnd,numMshRanksInMpiRank,maxNumMshRanks
      integer(4), allocatable :: mshRanksInMpiRank(:),mapMshRankToMpiRank(:)

      integer(4) :: mapFaceDir
      real(8) :: mapFaceGapCoord

! ################################################################################################
! ------------------------ VARS for MPI COMMS ----------------------------------------------------
! ################################################################################################

      type(jagged_vector_int4) :: nodesToComm_jv,commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv
      integer(4),allocatable :: numNodesToCommMshRank(:),numMshRanksWithComms(:)

      type(jagged_vector_int4) :: bnd_nodesToComm_jv,bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv
      integer(4),allocatable :: bnd_numNodesToCommMshRank(:),bnd_numMshRanksWithComms(:)

! ################################################################################################
! ------------------------ VARS for Mesh Quality ----------------------------------------------------
! ################################################################################################
      type(jagged_matrix_real8) :: quality_jm
      real(8) :: quality_elem,auxAvg,elemVertOrig(3)
      real(8) :: quality_vec(2)
      integer(4) :: iElem,iNode,iElemVTK,iElemGid
      integer(4) :: numTangledElemsMpiRank,numTangledElemsMshRank,numTangledElemsTotal
      real(8) :: maxQuality(2), minQuality(2), avgQuality(2)
      real(8) :: rankMaxQuality(2), rankMinQuality(2), rankAvgQuality(2)
      character(len=12) :: aux_numRanks,aux_mpiRank
      character(len=512) :: meshquality_filename,meshqualitySum_filename

! ################################################################################################
! ################################################################################################

      start_time(1) = MPI_Wtime()

      call init_hdf5_interface()

#if _CHECK_
       if(mpi_rank.eq.0) write(*,*) "Be Careful! _CHECK_ flag is ON! Debug files (*.csv) will be created!"
#endif

      call open_gmsh_h5_file(gmsh_filePath,gmsh_fileName,gmsh_h5_fileId)

      call read_dims_gmsh_h5_file(gmsh_h5_fileId,mporder,numElemsGmsh,numNodesGmsh_i8,numBoundFacesGmsh,numPerFacesGmsh,numPerLinkedNodesGmsh,numMapFacesGmsh,isPeriodic,isBoundaries,isMapFaces)

      call initGmshIJKTables(mporder)
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)
      call set_allocate_hexahedronHO_ijk_indices(mporder,gmsh2ijk,vtk2ijk)
      call set_allocate_quadrilateralHO_ij_indices(mporder,gmsh2ij,vtk2ij)

      allocate(a2ijk(mnnode))
      allocate(a2ij(mnpbou))
      allocate(ijk2a(0:porder,0:porder,0:porder))

      allocate(Ngp_l(mngaus,mnnode))
      allocate(dNgp(ndime,mnnode,mngaus))
      allocate(wgp(mngaus))

      !for the moment we set that the a2ijk is the gmsh2ijk
      a2ijk(:) = gmsh2ijk(:)
      a2ij(:)  = gmsh2ij(:)
      !ijk2a(:,:,:) = ijk2gmsh(:,:,:)

      start_time(2) = MPI_Wtime()
      call read_elems_nodes_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,gmsh_h5_fileId,isPeriodic,isMapFaces,numElemsGmsh,numNodesGmsh_i8,numPerFacesGmsh,numPerLinkedNodesGmsh,numMapFacesGmsh,&
                     numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank_i8,connecMpiRank_i8,coordNodesMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8,&
                     linkedPerToMapElemsAndFacesSrl,mapFaceDir,mapFaceGapCoord)
      end_time(2) = MPI_Wtime()
      elapsed_time_r(2) = end_time(2) - start_time(2)

      call distribute_ranks2Part_in_mpiRank(numMshRanks2Part,mshRankInMpiRankStart,mshRankInMpiRankEnd,numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank,mapMshRankToMpiRank)

      start_time(3) = MPI_Wtime()
      call do_element_partitioning_gempa_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,gmsh_h5_fileId,isPeriodic,isMapFaces,numElemsGmsh,numNodesGmsh_i8,numBoundFacesGmsh,&
                     numMshRanks2Part,numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank,mapMshRankToMpiRank,numElemsMpiRank,listElemsMpiRank,&
                     numNodesMpiRank,listNodesMpiRank_i8,connecMpiRank_i8,coordNodesMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerFacesGmsh,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8,&
                     numMapFacesGmsh,linkedPerToMapElemsAndFacesSrl,&
                     numElemsMshRank,numBoundElemsMshRank,mshRankElemStart,mshRankElemEnd,numNodesMshRank,numBoundFacesMshRank,numPerNodesMshRank,numPerMapLinkedNodesMshRank,maxBoundCode,&
                     elemGidMshRank_jv,listNodesMshRank_i8_jv,connecMshRank_i8_jm,listBoundElemsMshRank_jv,listElemsBoundFacesMshRank_jm,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,&
                     connecOrigBoundFacesMshRank_i8_jm,coordMshRank_jm,masSlaRankPar_i8_jm,perMapLinkedNodesRankPar_i8_jm)
      end_time(3) = MPI_Wtime()
      elapsed_time_r(3) = end_time(3) - start_time(3)

      call close_hdf5_file(gmsh_h5_fileId)

      !-----------------------------------------------------------------------
      allocate(boundaryNodes_jv%vector(numMshRanksInMpiRank))
      allocate(dofNodes_jv%vector(numMshRanksInMpiRank))
      allocate(mpiBoundaryNodes_i8_jv%vector(numMshRanksInMpiRank))
      allocate(mshBoundaryNodes_i8_jv%vector(numMshRanksInMpiRank))

      allocate(numMpiBoundaryNodes(numMshRanksInMpiRank))
      allocate(numMshBoundaryNodes(numMshRanksInMpiRank))
      allocate(numInnerNodes(numMshRanksInMpiRank))

      allocate(numBoundaryNodesMshRank(numMshRanksInMpiRank))
      allocate(numDoFMshRank(numMshRanksInMpiRank))

      allocate(mshRankNodeStart_i8(numMshRanksInMpiRank))
      allocate(mshRankNodeEnd_i8(numMshRanksInMpiRank))

      numMpiBoundaryNodes(:) = 0
      numMshBoundaryNodes(:) = 0
      numInnerNodes(:) = 0
      numBoundaryNodesMshRank(:) = 0
      numDoFMshRank(:) = 0
      mshRankNodeStart_i8(:) = 0
      mshRankNodeEnd_i8(:) = 0

      !-----------------------------------------------------------------------

      start_time(4) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Reading Msh Partitions Boundary Nodes"
      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)

         call get_rankPartitionBoundaryNodes_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,&
                  mshRank,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),&
                  listNodesMshRank_i8_jv%vector(iMshRank)%elems,connecMshRank_i8_jm%matrix(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems,&
                  listElemsBoundFacesMshRank_jm%matrix(iMshRank)%elems,listBoundFacesMshRank_jv%vector(iMshRank)%elems,&
                  connecOrigBoundFacesMshRank_i8_jm%matrix(iMshRank)%elems,&
                  numPerNodesMshRank(iMshRank),masSlaRankPar_i8_jm%matrix(iMshRank)%elems,&
                  numMshBoundaryNodes(iMshRank),numMpiBoundaryNodes(iMshRank),numInnerNodes(iMshRank),&
                  mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems,mshBoundaryNodes_i8_jv%vector(iMshRank)%elems)
      end do
      end_time(4) = MPI_Wtime()
      elapsed_time_r(4) = end_time(4) - start_time(4)

      call define_parallelNodePartitioning(numMshRanks2Part,numMshRanksInMpiRank,numNodesMshRank,mshRanksInMpiRank,mapMshRankToMpiRank,mshRankNodeStart_i8,mshRankNodeEnd_i8,iNodeStartPar_i8,numNodesParTotal_i8)

      !--------------------------------------------------------------------------------------
      allocate(auxNewOrderIndex(mnnode))
      allocate(auxVTKorder(mnnode))

      call generate_new_nodeOrder_and_connectivity(mporder,mnnode,gmsh2ijk,vtk2ijk,a2ijk,auxNewOrderIndex,auxVTKorder)
      call evalShapeFunctions_Ngp_l(Ngp_l,dNgp,wgp,mporder,mnnode,mngaus,a2ijk)

      allocate(globalIdSrl_i8_jv%vector(numMshRanksInMpiRank))
      allocate(globalIdSrlOrdered_i8_jm%matrix(numMshRanksInMpiRank))
      allocate(globalIdPar_i8_jv%vector(numMshRanksInMpiRank))
      allocate(connecVTK_jv%vector(numMshRanksInMpiRank))
      allocate(connecParOrig_jm%matrix(numMshRanksInMpiRank))
      allocate(connecParWork_jm%matrix(numMshRanksInMpiRank))
      allocate(coordPar_jm%matrix(numMshRanksInMpiRank))
      allocate(coordVTK_jm%matrix(numMshRanksInMpiRank))
      allocate(numWorkingNodesMshRank(numMshRanksInMpiRank))
      allocate(workingNodesPar_jv%vector(numMshRanksInMpiRank))
      allocate(connecBoundFacesMshRank_jm%matrix(numMshRanksInMpiRank))

      allocate(connecOrigBoundFacesMshRank_jm%matrix(numMshRanksInMpiRank))
      allocate(masSlaRankPar_jm%matrix(numMshRanksInMpiRank))
      allocate(perMapLinkedNodesRankPar_jm%matrix(numMshRanksInMpiRank))

      allocate(numElemsVTKMshRank(numMshRanksInMpiRank))
      allocate(sizeConnecVTKMshRank(numMshRanksInMpiRank))

      start_time(5) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Reordering nodes, interpolating coordinates and creating working lists..."

      if(isLinealOutput) then
         mnnodeVTK = 8
         numVTKElemsPerMshElem = (mporder**3)
      else
         mnnodeVTK = mnnode      
         numVTKElemsPerMshElem = 1
      end if

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)

         allocate(globalIdSrl_i8_jv%vector(iMshRank)%elems(numNodesMshRank(iMshRank)))
         allocate(globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(numNodesMshRank(iMshRank),2))
         allocate(globalIdPar_i8_jv%vector(iMshRank)%elems(numNodesMshRank(iMshRank)))

         numElemsVTKMshRank(iMshRank)   = numElemsMshRank(iMshRank)*numVTKElemsPerMshElem
         sizeConnecVTKMshRank(iMshRank) = numElemsVTKMshRank(iMshRank)*mnnodeVTK

         allocate(connecVTK_jv%vector(iMshRank)%elems(sizeConnecVTKMshRank(iMshRank)))
         allocate(connecParOrig_jm%matrix(iMshRank)%elems(numElemsMshRank(iMshRank),mnnode))
         allocate(connecParWork_jm%matrix(iMshRank)%elems(numElemsMshRank(iMshRank),mnnode))

         allocate(coordPar_jm%matrix(iMshRank)%elems(numNodesMshRank(iMshRank),3))
         allocate(coordVTK_jm%matrix(iMshRank)%elems(numNodesMshRank(iMshRank),3))

         allocate(connecBoundFacesMshRank_jm%matrix(iMshRank)%elems(numBoundFacesMshRank(iMshRank),mnpbou))
         allocate(connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems(numBoundFacesMshRank(iMshRank),mnpbou))
         allocate(masSlaRankPar_jm%matrix(iMshRank)%elems(numPerNodesMshRank(iMshRank),2))
         allocate(perMapLinkedNodesRankPar_jm%matrix(iMshRank)%elems(numPerMapLinkedNodesMshRank(iMshRank),2))

         !$acc kernels
         globalIdSrl_i8_jv%vector(iMshRank)%elems(:) = -1
         globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(:,:) = -1
         globalIdPar_i8_jv%vector(iMshRank)%elems(:) = -1
         connecVTK_jv%vector(iMshRank)%elems(:) = 0
         connecParOrig_jm%matrix(iMshRank)%elems(:,:) = 0
         connecParWork_jm%matrix(iMshRank)%elems(:,:) = 0
         coordPar_jm%matrix(iMshRank)%elems(:,:) = 0.0_rp
         coordVTK_jm%matrix(iMshRank)%elems(:,:) = 0.0_rp
         connecBoundFacesMshRank_jm%matrix(iMshRank)%elems(:,:) = 0

         connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems(:,:) = 0
         masSlaRankPar_jm%matrix(iMshRank)%elems(:,:) = 0
         perMapLinkedNodesRankPar_jm%matrix(iMshRank)%elems(:,:) = 0
         !$acc end kernels

         !reordering nodes
         call reorder_nodes_in_mshRank(mporder,mnnode,mnpbou,gmsh2ijk,vtk2ijk,auxNewOrderIndex,auxVTKorder,&
                        mshRank,numMshRanks2Part,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),&
                        numPerNodesMshRank(iMshRank),numPerMapLinkedNodesMshRank(iMshRank),sizeConnecVTKMshRank(iMshRank),isLinealOutput,&
                        elemGidMshRank_jv%vector(iMshRank)%elems,listNodesMshRank_i8_jv%vector(iMshRank)%elems,connecMshRank_i8_jm%matrix(iMshRank)%elems,iNodeStartPar_i8,&
                        connecOrigBoundFacesMshRank_i8_jm%matrix(iMshRank)%elems,masSlaRankPar_i8_jm%matrix(iMshRank)%elems,perMapLinkedNodesRankPar_i8_jm%matrix(iMshRank)%elems,&
                        globalIdSrl_i8_jv%vector(iMshRank)%elems,globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems,globalIdPar_i8_jv%vector(iMshRank)%elems,&
                        connecVTK_jv%vector(iMshRank)%elems,connecParOrig_jm%matrix(iMshRank)%elems,&
                        connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,masSlaRankPar_jm%matrix(iMshRank)%elems,perMapLinkedNodesRankPar_jm%matrix(iMshRank)%elems)

         !------------------------------------------------------------------------------
         deallocate(connecMshRank_i8_jm%matrix(iMshRank)%elems)
         deallocate(connecOrigBoundFacesMshRank_i8_jm%matrix(iMshRank)%elems)
         deallocate(masSlaRankPar_i8_jm%matrix(iMshRank)%elems)
         deallocate(perMapLinkedNodesRankPar_i8_jm%matrix(iMshRank)%elems)
         !------------------------------------------------------------------------------

         call set_nodesCoordinates(mnnode,mnpbou,mngaus,mshRank,isLinealOutput,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),globalIdSrl_i8_jv%vector(iMshRank)%elems,&
            listNodesMshRank_i8_jv%vector(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems,Ngp_l,connecParOrig_jm%matrix(iMshRank)%elems,&
            coordPar_jm%matrix(iMshRank)%elems,coordVTK_jm%matrix(iMshRank)%elems)

         !------------------------------------------------------------------------------
         deallocate(listNodesMshRank_i8_jv%vector(iMshRank)%elems)
         !------------------------------------------------------------------------------

         call create_working_lists_parallel(mnnode,mnpbou,isPeriodic,mshRank,numElemsMshRank(iMshRank),numNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),&
               connecParOrig_jm%matrix(iMshRank)%elems,numPerNodesMshRank(iMshRank),masSlaRankPar_jm%matrix(iMshRank)%elems,&
               connecParWork_jm%matrix(iMshRank)%elems,connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,connecBoundFacesMshRank_jm%matrix(iMshRank)%elems,&
               numWorkingNodesMshRank(iMshRank),workingNodesPar_jv%vector(iMshRank)%elems)

      end do

      deallocate(connecMshRank_i8_jm%matrix)
      deallocate(connecOrigBoundFacesMshRank_i8_jm%matrix)
      deallocate(masSlaRankPar_i8_jm%matrix)
      deallocate(perMapLinkedNodesRankPar_i8_jm%matrix)
      deallocate(listNodesMshRank_i8_jv%vector)

      end_time(5) = MPI_Wtime()
      elapsed_time_r(5) = end_time(5) - start_time(5)

      start_time(6) = MPI_Wtime()
      call generate_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numMpiBoundaryNodes,&
                                    mpiBoundaryNodes_i8_jv,globalIdSrl_i8_jv,globalIdSrlOrdered_i8_jm,&
                                    nodesToComm_jv,commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv,numNodesToCommMshRank,numMshRanksWithComms)

      end_time(6) = MPI_Wtime()
      elapsed_time_r(6) = end_time(6) - start_time(6)

      start_time(7) = MPI_Wtime()

      call generate_dof_and_boundary_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,&
               numNodesMshRank,numMshBoundaryNodes,mshBoundaryNodes_i8_jv,numMpiBoundaryNodes,mpiBoundaryNodes_i8_jv,globalIdSrl_i8_jv,globalIdSrlOrdered_i8_jm,&
               numBoundaryNodesMshRank,boundaryNodes_jv,numDoFMshRank,dofNodes_jv,&
               bnd_nodesToComm_jv,bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv,bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms)

      do iMshRank=1,numMshRanksInMpiRank
         deallocate(globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems)
      end do
      deallocate(globalIdSrlOrdered_i8_jm%matrix)

      end_time(7) = mpi_wtime()
      elapsed_time_r(7) = end_time(7) - start_time(7)

      !-------------------------------------------------------------------------------------------

      deallocate(Ngp_l)

      call MPI_Barrier(app_comm,mpi_err)
      !----------------------------------------------------------------------------------------------
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numWorkingNodesMshRank,vecNumWorkingNodes)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numNodesToCommMshRank,vecNumNodesToCommMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numMshRanksWithComms,vecNumMshRanksWithComms)

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numBoundFacesMshRank,vecNumBoundFacesMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numDoFMshRank,vecNumDoFMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numBoundaryNodesMshRank,vecNumBoundaryNodesMshRank)

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,bnd_numNodesToCommMshRank,vecBndNumNodesToCommMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,bnd_numMshRanksWithComms,vecBndNumMshRanksWithComms)

      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numPerNodesMshRank,vecNumPerNodesMshRank)
      call get_vector_with_mshRank_values_for_numMshRanks2Part(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numPerMapLinkedNodesMshRank,vecNumPerMapLinkedNodesMshRank)

      !write(*,*) 'l1[',mpi_rank,']vnwn',vecNumWorkingNodes(:),'vnntcmr',vecNumNodesToCommMshRank(:),'vnmrwc',vecNumMshRanksWithComms(:)
      !write(*,*) 'l2[',mpi_rank,']vnbf',vecNumBoundFacesMshRank(:),'vndof',vecNumDoFMshRank(:),'vnbn',vecNumBoundaryNodesMshRank(:)

      !----------------------------------------------------------------------------------------------
      allocate(quality_jm%matrix(numMshRanksInMpiRank))
      if(evalMeshQuality) then
         if(numMshRanksInMpiRank .gt. 0) then
            ! Open a file for outputting realted information
            write(aux_mpiRank,'(I0)') mpi_rank
            write(aux_numRanks,'(I0)') numMshRanks2Part
            !-----------------------------------------------------------------------------------------------
            meshquality_filename = "meshQuality_"//trim(adjustl(mesh_h5_fileName))//'-'//trim(aux_numRanks)//'_rank'//trim(aux_mpiRank)//'.dat'
            open(unit=555, file=meshquality_filename, status="unknown", action="write", form="formatted")
            write(555,*) "--| Evaluating mesh quality"
            write(555,*) "----| List of tangled elements (GMSH global numeration)"
            call flush(555)
         end if
         !-----------------------------------------------------------------------------------------------
         if(mpi_rank.eq.0) then
            meshqualitySum_filename = "meshQualitySummary_"//trim(adjustl(mesh_h5_fileName))//'-'//trim(aux_numRanks)//'.dat'
            open(unit=554, file=meshqualitySum_filename, status="unknown", action="write", form="formatted")
            write(554,*) "--| Evaluating mesh quality (summary)"
            call flush(554)
         end if
         !-----------------------------------------------------------------------------------------------
         call MPI_Barrier(app_comm,mpi_err)

         ! Initialize vars to get high-level info about mesh quality
         maxQuality(:) = 0.0d0          ! Max absolute
         minQuality(:) = 0.0d0          ! Min absolute
         avgQuality(:) = 0.0d0          ! Avg across all ranks
         numTangledElemsMpiRank=0

         call eval_MeshQuality(numMshRanksInMpiRank,numElemsVTKMshRank,numElemsMshRank,numNodesMshRank,mnnode,mngaus, coordPar_jm,connecParOrig_jm,dNgp,wgp,mshRanksInMpiRank,elemGidMshRank_jv,numVTKElemsPerMshElem,quality_jm,rankMaxQuality, rankMinQuality, rankAvgQuality)

         call MPI_Barrier(app_comm,mpi_err)
         write(555,*) "----| End of list of tangled elements"
         write(555,*) "--| Mesh quality evaluated"
         write(555,*) "--| Mesh statistics in rank"
         write(555,*) "----| Num.Tangled Elems:",numTangledElemsMpiRank
         write(555,*) "----| Quality:          Anisotropic                Isotropic"
         write(555,*) "----| Max quality:",rankMaxQuality(:)
         write(555,*) "----| Min quality:",rankMinQuality(:)
         write(555,*) "----| Avg quality:",rankAvgQuality(:)
         close(555)
         ! Compute max, min and avg across all ranks
         call MPI_Reduce(numTangledElemsMpiRank,numTangledElemsTotal,1,mpi_datatype_int4,MPI_SUM,0,app_comm,mpi_err)
         call MPI_Reduce(rankMaxQuality(1),maxQuality(1),1,mpi_datatype_real8,MPI_MAX,0,app_comm,mpi_err)
         call MPI_Reduce(rankMinQuality(1),minQuality(1),1,mpi_datatype_real8,MPI_MIN,0,app_comm,mpi_err)
         call MPI_Reduce(rankAvgQuality(1),avgQuality(1),1,mpi_datatype_real8,MPI_SUM,0,app_comm,mpi_err)
         call MPI_Reduce(rankMaxQuality(2),maxQuality(2),1,mpi_datatype_real8,MPI_MAX,0,app_comm,mpi_err)
         call MPI_Reduce(rankMinQuality(2),minQuality(2),1,mpi_datatype_real8,MPI_MIN,0,app_comm,mpi_err)
         call MPI_Reduce(rankAvgQuality(2),avgQuality(2),1,mpi_datatype_real8,MPI_SUM,0,app_comm,mpi_err)
         avgQuality = avgQuality / mpi_size
         ! Write high-level data to file
         if(mpi_rank.eq.0) then
            write(554,*) "--| Mesh statistics"
            write(554,*) "----| Num.Tangled Elems:",numTangledElemsTotal
            write(554,*) "----| Quality:          Anisotropic                Isotropic"
            write(554,*) "----| Max quality:", maxQuality(:)
            write(554,*) "----| Min quality:", minQuality(:)
            write(554,*) "----| Avg quality:", avgQuality(:)
            close(554)
         end if
      else
         do iMshRank=1,numMshRanksInMpiRank
            allocate(quality_jm%matrix(iMshRank)%elems(numElemsVTKMshRank(iMshRank),2))
         end do
      end if
      !----------------------------------------------------------------------------------------------
      start_time(8) = MPI_Wtime()
      if(mpi_rank.eq.0) write(*,*) "--| Saving Partitioned HDF5 mesh file:",mesh_h5_fileName

      connecChunkSize = 10000000

      call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName,numMshRanks2Part,meshFile_h5_full_name)
      call MPI_Barrier(app_comm,mpi_err)

      call create_hdf5_file(meshFile_h5_full_name,sod2dmsh_h5_fileId)

      call create_hdf5_groups_datasets_in_meshFile_from_tool(mnnode,mnpbou,sod2dmsh_h5_fileId,isPeriodic,isBoundaries,isMapFaces,isLinealOutput,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8,&
               vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank,vecNumPerMapLinkedNodesMshRank)

      call create_groups_datasets_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,sod2dmsh_h5_fileId,isLinealOutput,evalMeshQuality,numMshRanks2Part,numElemsGmsh,numNodesParTotal_i8,mnnodeVTK,numVTKElemsPerMshElem)

      call MPI_Barrier(app_comm,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         call write_mshRank_data_in_hdf5_meshFile_from_tool(mporder,mnnode,mnpbou,sod2dmsh_h5_fileId,mshRank,numMshRanks2Part,isPeriodic,isBoundaries,isMapFaces,isLinealOutput,numElemsGmsh,numBoundFacesGmsh,&
            numElemsMshRank(iMshRank),mshRankElemStart(iMshRank),mshRankElemEnd(iMshRank),mshRankNodeStart_i8(iMshRank),&
            mshRankNodeEnd_i8(iMshRank),numNodesMshRank(iMshRank),numWorkingNodesMshRank(iMshRank),numBoundFacesMshRank(iMshRank),numBoundaryNodesMshRank(iMshRank),numDoFMshRank(iMshRank),maxBoundCode,&
            numElemsVTKMshRank(iMshRank),sizeConnecVTKMshRank(iMshRank),mnnodeVTK,numVTKElemsPerMshElem,mapFaceDir,mapFaceGapCoord,&
            a2ijk,a2ij,gmsh2ijk,gmsh2ij,vtk2ijk,vtk2ij,&
            elemGidMshRank_jv%vector(iMshRank)%elems,globalIdSrl_i8_jv%vector(iMshRank)%elems,globalIdPar_i8_jv%vector(iMshRank)%elems,&
            connecParOrig_jm%matrix(iMshRank)%elems,connecParWork_jm%matrix(iMshRank)%elems,coordPar_jm%matrix(iMshRank)%elems,workingNodesPar_jv%vector(iMshRank)%elems,&
            boundaryNodes_jv%vector(iMshRank)%elems,dofNodes_jv%vector(iMshRank)%elems,boundFacesCodesMshRank_jv%vector(iMshRank)%elems,connecOrigBoundFacesMshRank_jm%matrix(iMshRank)%elems,connecBoundFacesMshRank_jm%matrix(iMshRank)%elems,&
            numPerNodesMshRank(iMshRank),masSlaRankPar_jm%matrix(iMshRank)%elems,numPerMapLinkedNodesMshRank(iMshRank),perMapLinkedNodesRankPar_jm%matrix(iMshRank)%elems,&
            numNodesToCommMshRank(iMshRank),numMshRanksWithComms(iMshRank),nodesToComm_jv%vector(iMshRank)%elems,commsMemPosInLoc_jv%vector(iMshRank)%elems,&
            commsMemSize_jv%vector(iMshRank)%elems,commsMemPosInNgb_jv%vector(iMshRank)%elems,ranksToComm_jv%vector(iMshRank)%elems,&
            bnd_numNodesToCommMshRank(iMshRank),bnd_numMshRanksWithComms(iMshRank),bnd_nodesToComm_jv%vector(iMshRank)%elems,bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems,&
            bnd_commsMemSize_jv%vector(iMshRank)%elems,bnd_commsMemPosInNgb_jv%vector(iMshRank)%elems,bnd_ranksToComm_jv%vector(iMshRank)%elems,&
            vecNumWorkingNodes,vecNumMshRanksWithComms,vecNumNodesToCommMshRank,vecBndNumMshRanksWithComms,vecBndNumNodesToCommMshRank,vecNumBoundFacesMshRank,vecNumDoFMshRank,vecNumBoundaryNodesMshRank,vecNumPerNodesMshRank,vecNumPerMapLinkedNodesMshRank)

        call write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,sod2dmsh_h5_fileId,evalMeshQuality,mshRank,numMshRanks2Part,numElemsMshRank(iMshRank),&
            numElemsVTKMshRank(iMshRank),sizeConnecVTKMshRank(iMshRank),mnnodeVTK,numVTKElemsPerMshElem,&
            mshRankElemStart(iMshRank),mshRankElemEnd(iMshRank),mshRankNodeStart_i8(iMshRank),mshRankNodeEnd_i8(iMshRank),numNodesMshRank(iMshRank),&
            coordVTK_jm%matrix(iMshRank)%elems,connecVTK_jv%vector(iMshRank)%elems,quality_jm%matrix(iMshRank)%elems,connecChunkSize)
      end do

      do iMshRank=(numMshRanksInMpiRank+1),maxNumMshRanks
        call dummy_write_mshRank_data_in_hdf5_meshFile_from_tool(sod2dmsh_h5_fileId,numMshRanks2Part,isPeriodic,isBoundaries,isMapFaces,isLinealOutput)

        call dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(sod2dmsh_h5_fileId,evalMeshQuality,numMshRanks2Part,connecChunkSize)
      end do

      call MPI_Barrier(app_comm,mpi_err)
      !--------------------------------------------------------------------------------------------------------------------------------


      call MPI_Barrier(app_comm,mpi_err)
      !--------------------------------------------------------------------------------------------------------------------------------

      call close_hdf5_file(sod2dmsh_h5_fileId)

      deallocate(a2ijk)
      deallocate(a2ij)
      deallocate(gmsh2ijk)
      deallocate(gmsh2ij)
      deallocate(vtk2ijk)
      deallocate(vtk2ij)

      end_time(8) = mpi_wtime()
      elapsed_time_r(8) = end_time(8) - start_time(8)

      end_time(1) = MPI_Wtime()
      elapsed_time_r(1) = end_time(1) - start_time(1)
      if(mpi_rank.eq.0) then
         write(*,*) "--| Everything Done!"
         write(*,*) "------------------------------------------------"
         write(*,*) 'Total Mesh generation time',elapsed_time_r(1)
         write(*,*) ' -ReadElemsNodesAndPer',elapsed_time_r(2)
         write(*,*) ' -ElementPartition',elapsed_time_r(3)
         write(*,*) ' -RankBoundaries',elapsed_time_r(4)
         write(*,*) ' -ReOrder,Coords,Working',elapsed_time_r(5)
         write(*,*) ' -MpiCommScheme',elapsed_time_r(6)
         write(*,*) ' -DoF & BoundCommScheme',elapsed_time_r(7)
         write(*,*) ' -SavingHDF5',elapsed_time_r(8)
      end if

      call end_hdf5_interface()

   end subroutine read_gmsh_h5_file_and_do_partitioning_in_parallel
!---------------------------------------------------------------------------------------------------------------------------
   subroutine open_gmsh_h5_file(gmsh_h5_filePath,gmsh_h5_fileName,gmsh_h5_fileId)
      implicit none
      character(len=*), intent(in) :: gmsh_h5_filePath,gmsh_h5_fileName
      integer(hid_t), intent(out) :: gmsh_h5_fileId
      character(512) :: gmsh_h5_fullFileName
      integer(hid_t) :: file_id,plist_id
      integer :: h5err

      gmsh_h5_fullFileName = trim(adjustl(gmsh_h5_filePath))//trim(adjustl(gmsh_h5_fileName))//'.h5'

      ! Setup file access property list with parallel I/O access.
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,h5err)
      call h5pset_fapl_mpio_f(plist_id,app_comm,MPI_INFO_NULL,h5err)

      if(mpi_rank.eq.0) write(*,*) '# Opening gmsh h5 mesh file:', gmsh_h5_fullFileName

      call h5fopen_f(gmsh_h5_fullFileName, H5F_ACC_RDWR_F,gmsh_h5_fileId,h5err,access_prp=plist_id)
      if(h5err .ne. 0) then
         write(*,*) 'FATAL ERROR! Cannot open h5 meshfile ',trim(adjustl(gmsh_h5_fullFileName))
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      call h5pclose_f(plist_id, h5err)

   end subroutine open_gmsh_h5_file

   subroutine read_dims_gmsh_h5_file(gmsh_h5_fileId,mporder,numElems,numNodes_i8,numBoundFaces,numPerFaces,numPerLinks,numMapFaces,isPeriodic,isBoundaries,isMapFaces)
      implicit none
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      integer(4),intent(out) :: mporder,numElems,numBoundFaces,numPerFaces,numPerLinks,numMapFaces
      integer(8),intent(out) :: numNodes_i8
      logical, intent(out) :: isPeriodic,isBoundaries,isMapFaces

      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hid_t) :: dtype
      integer(4) :: h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(8),allocatable :: aux_array_i8(:)

      !write(*,*) 'Loading parallel data hdf5...'

      dtype = h5_datatype_int8
      ms_dims(1) = 1
      ms_offset(1) = 0
      allocate(aux_array_i8(1))

      dsetname = '/dims/order'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      mporder=int(aux_array_i8(1),4)

      dsetname = '/dims/numNodes'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      numNodes_i8=aux_array_i8(1)

      dsetname = '/dims/numElements'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      numElems=int(aux_array_i8(1),4)

      dsetname = '/dims/numBoundaryFaces'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      numBoundFaces=int(aux_array_i8(1),4)

      dsetname = '/dims/numPeriodicFaces'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      numPerFaces=int(aux_array_i8(1),4)

      dsetname = '/dims/numPeriodicLinks'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      numPerLinks=int(aux_array_i8(1),4)

      dsetname = '/dims/numMappedFaces'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      numMapFaces=int(aux_array_i8(1),4)

      if(numBoundFaces.ne.0) isBoundaries = .true.
      if(numPerFaces.ne.0) isPeriodic = .true.
      if(numMapFaces.ne.0) isMapFaces = .true.

      if(mpi_rank.eq.0) write(*,*) 'mporder',mporder,'numNodes',numNodes_i8,'numElems',numElems,'numBoundFaces',numBoundFaces,'numPerFaces',numPerFaces,'numPerLinks',numPerLinks,'numMapFaces',numMapFaces

      deallocate(aux_array_i8)

   end subroutine read_dims_gmsh_h5_file


   subroutine read_elem_connec_and_nodes_coords_from_gmsh_h5_file(mnnode,gmsh_h5_fileId,numElemsSrl,numNodesSrl_i8,numElemsInRank,listElemsRank,numNodesInRank,connecRank_i8,listNodesRank_i8,coordNodesRank)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      integer(4),intent(in) :: numElemsSrl,numElemsInRank,listElemsRank(numElemsInRank)
      integer(8),intent(in) :: numNodesSrl_i8
      integer(4),intent(out) :: numNodesInRank
      integer(8),intent(out) :: connecRank_i8(numElemsInRank,mnnode)
      integer(8),allocatable,intent(inout) :: listNodesRank_i8(:)
      real(8),allocatable,intent(inout) :: coordNodesRank(:,:)
      integer(4) :: iElem,iNode,nodeCnt,iAux,jAux,iDim
      integer(8) :: iNodeG_i8,prevNodeG_i8
      !-----------------------------------------------------------
      character(128) :: dsetname
      integer(4),parameter :: ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      !------------------------------------------------------------
      integer(8),allocatable :: rawNodeListRank_i8(:)
      real(8), dimension(5) :: start_time,end_time,elapsed_time,elapsed_time_m

      if(mpi_rank.eq.0) write(*,*) ' -numElems2read',numElemsInRank
      start_time(1) = MPI_Wtime()

      dsetname = '/connec'
      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      !haig de jugar amb el ratio fs_dims i ms_dims

      !call read_elem_connec_by_coords_selection(dset_id,fspace_id,plist_id,numElemsInRank,listElemsRank,connecRank_i8)
      call read_elem_connec_by_chunks(mnnode,dset_id,fspace_id,plist_id,numElemsSrl,numElemsInRank,listElemsRank,connecRank_i8)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      !----------------------------------------------------------------------------------------------------------
      start_time(2) = MPI_Wtime()
      !----------------------------------------------------------------------------------------------------------
      ! REORDERING NODE LIST

      allocate(rawNodeListRank_i8(numElemsInRank*mnnode))
      !$acc kernels
      rawNodeListRank_i8(:) = 0
      !$acc end kernels

      ! add all the iNodeGs of this mpirank
      nodeCnt=0
      do iAux = 1,numElemsInRank
         do jAux = 1,mnnode
            iNodeG_i8 = connecRank_i8(iAux,jAux)
            nodeCnt=nodeCnt+1
            rawNodeListRank_i8(nodeCnt) = iNodeG_i8
         end do
      end do

      !sorting the nodes id
      call quicksort_array_int8(rawNodeListRank_i8)

      prevNodeG_i8=0
      nodeCnt=0
      do iAux = 1,(numElemsInRank*mnnode)
         if((rawNodeListRank_i8(iAux)).ne.(prevNodeG_i8)) then
            nodeCnt=nodeCnt+1
            prevNodeG_i8 = rawNodeListRank_i8(iAux)
         end if
      end do
      numNodesInRank = nodeCnt
      !write(*,*) 'newMethod.[',mpi_rank,'] numNodesInRank ',numNodesInRank

      allocate(listNodesRank_i8(numNodesInRank))
      allocate(coordNodesRank(numNodesInRank,ndime))

      !$acc kernels
      coordNodesRank(:,:) = 0
      listNodesRank_i8(:) = 0
      !$acc end kernels

      prevNodeG_i8=0
      nodeCnt=0
      do iAux = 1,(numElemsInRank*mnnode)
         if((rawNodeListRank_i8(iAux)).ne.(prevNodeG_i8)) then
            nodeCnt=nodeCnt+1
            listNodesRank_i8(nodeCnt) = rawNodeListRank_i8(iAux)
            prevNodeG_i8 = rawNodeListRank_i8(iAux)
         end if
      end do
      !if(mpi_rank.eq.0) write(*,*) 'listNodesRank',listNodesRank_i8(:)

      deallocate(rawNodeListRank_i8)

      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      start_time(3) = MPI_Wtime()
      ! Nodal coordinates section
      !------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0)write(*,*) "--| Reading coordinates..."

      dsetname = '/coords'

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      !call read_nodes_coords_by_coords_selection(dset_id,fspace_id,plist_id,numNodesInRank,listNodesRank_i8,coordNodesRank)
      call read_nodes_coords_by_chunks(dset_id,fspace_id,plist_id,numNodesSrl_i8,numNodesInRank,listNodesRank_i8,coordNodesRank)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(3) = MPI_Wtime()
      elapsed_time(3) = end_time(3) - start_time(3)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),3,mpi_datatype_real8,MPI_MAX,0,app_comm,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '   1.Read Elems',elapsed_time_m(1)
         write(*,*) '   2.Sort Nodes',elapsed_time_m(2)
         write(*,*) '   3.Read Coord',elapsed_time_m(3)
      end if

   end subroutine read_elem_connec_and_nodes_coords_from_gmsh_h5_file
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_elem_connec_by_chunks(mnnode,dset_id,fspace_id,plist_id,numElemsSrl,numElemsInRank,listElemsRank,connecRank_i8)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer(4),intent(in) :: numElemsSrl,numElemsInRank,listElemsRank(numElemsInRank)
      integer(8),intent(out) :: connecRank_i8(numElemsInRank,mnnode)
      integer(4) :: nextElemToRead,elemCnt,iElem,iElemG
      integer(4) :: maxRows2read,numChunks,iChunk,chunk2sing
      !-----------------------------------------------------------
      integer(4),parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: mspace_id
      integer(4) :: h5err
      integer(hid_t) :: dtype,test_dtype
      integer(4), allocatable :: vecChunks(:)
      integer(8), allocatable :: auxConnec_i8(:,:)
      !----------------------------------------------------------------------------------------------------------
      !TESTING ZONE READ ELEMENTS BY CHUNKS!
      !maxRows2read = 20000
      maxRows2read = 100000
      chunk2sing = 50

      numChunks = ceiling(real(numElemsSrl)/real(maxRows2read))
      if(mpi_rank.eq.0) write(*,*) ' #Reading elems by chunks | numChunks',numChunks,'maxRows2read',maxRows2read


      ms_dims(1) = int(mnnode,hsize_t)
      ms_dims(2) = 1
      ms_offset(1) = 0
      ms_offset(2) = 0
      dtype = h5_datatype_int8

      allocate(vecChunks(numChunks))
      call distribution_algorithm(numElemsSrl,numChunks,vecChunks)

      elemCnt=0
      iElemG=0
      nextElemToRead=listElemsRank(elemCnt+1)
      ms_offset(2) = 0
      do iChunk=1,numChunks
         ms_dims(1) = int(mnnode,hsize_t)
         ms_dims(2) = int(vecChunks(iChunk),hsize_t)
         allocate(auxConnec_i8(ms_dims(1),ms_dims(2)))
         if((mpi_rank.eq.0).and.(mod(iChunk,chunk2sing).eq.0)) write(*,*) '  -chunk(elem)',iChunk,'/',numChunks
         !write(*,*) 'iChunk',iChunk,'ms_dims',ms_dims,'ms_offset',ms_offset(2),'+',ms_offset(2) + ms_dims(2)!,'fsdims',fs_dims

         call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         call h5dread_f(dset_id,dtype,auxConnec_i8,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5sclose_f(mspace_id,h5err)

         do iElem=1,ms_dims(2)
            iElemG=iElemG+1
            if(iElemG.eq.nextElemToRead) then
               elemCnt=elemCnt+1
               connecRank_i8(elemCnt,1:mnnode) = auxConnec_i8(:,iElem)
               !write(*,*) 'connecRank_i8[',mpi_rank,'](',elemCnt,')',connecRank(elemCnt,:)
               if(elemCnt.lt.numElemsInRank) nextElemToRead=listElemsRank(elemCnt+1)
            end if
         end do

         deallocate(auxConnec_i8)
         ms_offset(2) = ms_offset(2) + ms_dims(2)
      end do
      !write(*,*) 'elemCnt',elemCnt,'numElemsRank',numElemsInRank

      deallocate(vecChunks)
      !----------------------------------------------------------------------------------------------------------

   end subroutine read_elem_connec_by_chunks
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_elem_connec_by_coords_selection(mnnode,dset_id,fspace_id,plist_id,numElemsInRank,listElemsRank,connecRank_i8)
      implicit none
      integer(4),intent(in) :: mnnode
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer(4),intent(in) :: numElemsInRank,listElemsRank(numElemsInRank)
      integer(8),intent(out) :: connecRank_i8(numElemsInRank,mnnode)
      !-----------------------------------------------------------------------
      integer(8), allocatable :: auxConnec_i8(:,:)
      integer(4) :: iAux,iElem,iNode,elemToRead
      !-----------------------------------------------------------------------
      integer(4),parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hid_t) :: mspace_id
      integer(4) :: h5err
      integer(hid_t) :: dtype
      integer(hsize_t) :: ms_numElems
      integer(hsize_t),allocatable :: ms_coords(:,:)
      !-----------------------------------------------------------------------

      dtype = h5_datatype_int8

      ms_dims(1) = int(mnnode,hsize_t)
      ms_dims(2) = int(numElemsInRank,hsize_t)
      ms_numElems = mnnode*numElemsInRank
      !write(*,*) 'ms_numElems(',mpi_rank,')',ms_numElems
      allocate(ms_coords(ms_rank,mnnode*numElemsInRank))
      allocate(auxConnec_i8(mnnode,numElemsInRank))

      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.

      !default
      ms_coords(1,:) = 1
      ms_coords(2,:) = 1

      iAux=0
      do iElem=1,numElemsInRank
         elemToRead=listElemsRank(iElem)
         do iNode=1,mnnode
            iAux=iAux+1
            ms_coords(1,iAux) = iNode
            ms_coords(2,iAux) = elemToRead
         end do
      end do

      call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err)
      call h5dread_f(dset_id,dtype,auxConnec_i8,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      do iElem=1,numElemsInRank
         connecRank_i8(iElem,:) = auxConnec_i8(:,iElem)
      end do

      deallocate(auxConnec_i8)
      deallocate(ms_coords)

      call h5sclose_f(mspace_id,h5err)

   end subroutine read_elem_connec_by_coords_selection
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_nodes_coords_by_chunks(dset_id,fspace_id,plist_id,numNodesSrl_i8,numNodesInRank,listNodesRank_i8,coordNodesRank)
      implicit none
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer(8),intent(in) :: numNodesSrl_i8,listNodesRank_i8(numNodesInRank)
      integer(4),intent(in) :: numNodesInRank
      real(8),intent(inout) :: coordNodesRank(:,:)
      !-----------------------------------------------------------------------
      real(8),allocatable :: auxCoords(:,:)
      integer(8) :: nextNodeToRead,nodeCnt,iNode,iNodeG
      !-----------------------------------------------------------------------
      integer,parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: mspace_id
      integer(4) :: h5err
      integer(hid_t) :: dtype
      integer(8), allocatable :: vecChunks(:),auxConnec(:,:)
      !----------------------------------------------------------------------------------------------------------
      integer(8) :: maxRows2read,numChunks,iChunk,chunk2sing

      !maxRows2read = 1000000
      maxRows2read = 100000
      chunk2sing = 50

      numChunks = ceiling(real(numNodesSrl_i8)/real(maxRows2read))
      if(mpi_rank.eq.0) write(*,*) ' #Reading nodescoords by chunks | numChunks',numChunks,'maxRows2read',maxRows2read

      ms_dims(1) = int(ndime,hsize_t)
      ms_dims(2) = 1
      ms_offset(1) = 0
      ms_offset(2) = 0
      dtype = h5_datatype_real8

      allocate(vecChunks(numChunks))
      call distribution_algorithm_int8(numNodesSrl_i8,numChunks,vecChunks)

      nodeCnt=0
      iNodeG=0
      nextNodeToRead=listNodesRank_i8(nodeCnt+1)
      ms_offset(2) = 0
      do iChunk=1,numChunks
         ms_dims(1) = int(ndime,hsize_t)
         ms_dims(2) = int(vecChunks(iChunk),hsize_t)
         allocate(auxCoords(ms_dims(1),ms_dims(2)))
         if((mpi_rank.eq.0).and.(mod(iChunk,chunk2sing).eq.0)) write(*,*) '  -chunk(crds)',iChunk,'/',numChunks
         !write(*,*) 'iChunk',iChunk,'ms_dims',ms_dims,'ms_offset',ms_offset(2),'+',ms_offset(2) + ms_dims(2),'fsdims',fs_dims

         call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err)
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         call h5dread_f(dset_id,dtype,auxCoords,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5sclose_f(mspace_id,h5err)

         do iNode=1,ms_dims(2)
            iNodeG=iNodeG+1
            if(iNodeG.eq.nextNodeToRead) then
               nodeCnt=nodeCnt+1
               coordNodesRank(nodeCnt,1:ndime) = auxCoords(:,iNode)

               if(nodeCnt.lt.numNodesInRank) nextNodeToRead=listNodesRank_i8(nodeCnt+1)
            end if
         end do

         deallocate(auxCoords)
         ms_offset(2) = ms_offset(2) + ms_dims(2)
      end do

      deallocate(vecChunks)
      !----------------------------------------------------------------------------------------------------------

   end subroutine read_nodes_coords_by_chunks
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine read_nodes_coords_by_coords_selection(dset_id,fspace_id,plist_id,numNodesInRank,listNodesRank_i8,coordNodesRank)
      implicit none
      integer(hid_t),intent(in) :: dset_id,fspace_id,plist_id
      integer(4),intent(in) :: numNodesInRank
      integer(8),intent(in) :: listNodesRank_i8(numNodesInRank)
      real(rp),intent(inout) :: coordNodesRank(:,:)
      !-----------------------------------------------------------------------
      real(8),allocatable :: auxCoords(:,:)
      integer(4) :: iAux,iElem,iNode,iDim
      integer(8) :: nodeToRead_i8
      !-----------------------------------------------------------------------
      integer,parameter :: ms_rank=2
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hid_t) :: mspace_id
      integer :: h5err
      integer(hid_t) :: dtype
      integer(hsize_t) :: ms_numElems
      integer(hsize_t),allocatable :: ms_coords(:,:)
      !-----------------------------------------------------------------------

      dtype = h5_datatype_real8

      ms_dims(1) = int(ndime,hsize_t)
      ms_dims(2) = int(numNodesInRank,hsize_t)

      ms_numElems = ms_dims(1)*ms_dims(2)
      allocate(ms_coords(ms_rank,ms_dims(1)*ms_dims(2)))
      allocate(auxCoords(ms_dims(1),ms_dims(2)))

      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.

      !default
      ms_coords(1,:) = 1
      ms_coords(2,:) = 1

      iAux=0
      do iNode=1,numNodesInRank
         nodeToRead_i8=listNodesRank_i8(iNode)
         do iDim=1,ndime
            iAux=iAux+1
            ms_coords(1,iAux) = iDim
            ms_coords(2,iAux) = nodeToRead_i8
         end do
      end do

      call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err)
      call h5dread_f(dset_id,dtype,auxCoords,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      do iNode=1,numNodesInRank
         coordNodesRank(iNode,:) = auxCoords(:,iNode)
      end do

      deallocate(auxCoords)
      deallocate(ms_coords)

      call h5sclose_f(mspace_id,h5err)

   end subroutine read_nodes_coords_by_coords_selection
!-------------------------------------------------------------------------------------------------------------------------------------------

   subroutine read_periodic_faces_and_links_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,&
                  numElemsInRank,listElemsInRank,connecInRank_i8,numBoundElemsMpiRank,listBoundElemsMpiRank,numLinkedPerElemsSrl,linkedPerElemsSrl,&
                  numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8,&
                  numPerFacesInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8)
      implicit none
      integer(4), intent(in) :: mporder,mnnode,mnpbou
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      integer(4),intent(in) :: numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,listElemsInRank(numElemsInRank)
      integer(4),intent(in) :: numBoundElemsMpiRank,listBoundElemsMpiRank(numBoundElemsMpiRank)
      integer(8),intent(in) :: connecInRank_i8(numElemsInRank,mnnode)
      integer(4),intent(out) :: numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl,numPerFacesInRank
      integer(4),allocatable,intent(inout) :: linkedPerElemsSrl(:,:),listPerElemsSrl(:)
      integer(8),allocatable,intent(inout) :: masSlaNodesSrl_i8(:,:)
      integer(4),allocatable,intent(inout) :: listElemsPerFacesInRank(:,:),listPerFacesInRank(:)
      integer(8),allocatable,intent(inout) :: connecPerFacesInRank_i8(:,:)
      integer(4) :: h5err

      character(128) :: dsetname2read

      real(8), dimension(2) :: start_time,end_time,elapsed_time,elapsed_time_m

      if(mpi_rank.eq.0) write(*,*) "--| Reading periodic boundaries..."

      start_time(1) = MPI_Wtime()
      dsetname2read = '/periodicFaces'
      call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,.true.,gmsh_h5_fileId,dsetname2read,&
                           numPerFacesSrl,numElemsInRank,listElemsInRank,numBoundElemsMpiRank,listBoundElemsMpiRank,connecInRank_i8,&
                           numPerFacesInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8)
      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      if(mpi_rank.eq.0) write(*,*) "--| Reading periodic links..."

      start_time(2) = MPI_Wtime()
      call read_periodic_links_from_gmsh_h5_file_in_parallel(mporder,mnpbou,gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,numPerFacesInRank,&
                     listElemsInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8,numLinkedPerElemsSrl,linkedPerElemsSrl,&
                     numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8)
      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),2,mpi_datatype_real8,MPI_MAX,0,app_comm,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '   1.Time ReadPerBounds',elapsed_time(1)
         write(*,*) '   2.Time ReadPerLinks',elapsed_time(2)
      end if
      !-----------------------------------------------

   end subroutine read_periodic_faces_and_links_from_gmsh_h5_file_in_parallel

   subroutine read_mapped_faces_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,gmsh_h5_fileId,numMapFacesSrl,numElemsInRank,listElemsInRank,&
                     numNodesInRank,numPerFacesInRank,connecInRank_i8,listNodesInRank_i8,coordNodesInRank,numBoundElemsMpiRank,listBoundElemsMpiRank,&
                     listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8,linkedPerToMapElemsAndFacesSrl,mapFaceDir,mapFaceGapCoord)
      implicit none
      integer(4), intent(in) :: mporder,mnnode,mnpbou
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      integer(4),intent(in) :: numMapFacesSrl,numElemsInRank,numNodesInRank,numPerFacesInRank,listElemsInRank(numElemsInRank)
      integer(8),intent(in) :: connecInRank_i8(numElemsInRank,mnnode),listNodesInRank_i8(numNodesInRank)
      real(8),intent(in) :: coordNodesInRank(numNodesInRank,ndime)
      integer(4),intent(in) :: numBoundElemsMpiRank,listBoundElemsMpiRank(numBoundElemsMpiRank)
      integer(4),intent(in) :: listElemsPerFacesInRank(:,:),listPerFacesInRank(:)
      integer(8),intent(in) :: connecPerFacesInRank_i8(:,:)
      integer(4),intent(out) :: linkedPerToMapElemsAndFacesSrl(numMapFacesSrl,4)
      integer(4),intent(out) :: mapFaceDir
      real(8),intent(out)  :: mapFaceGapCoord
      integer(4) :: numMapFacesInRank,h5err
      real(8) :: mapFaceGap
      integer(4),allocatable :: listElemsMapFacesInRank(:,:),listMapFacesInRank(:)
      integer(8),allocatable :: connecMapFacesInRank_i8(:,:)


      character(128) :: dsetname2read

      real(8), dimension(1) :: start_time,end_time,elapsed_time,elapsed_time_m

      if(mpi_rank.eq.0) write(*,*) "--| Reading mapped boundaries..."

      start_time(1) = MPI_Wtime()

      dsetname2read = '/mappedFaces'
      call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,.false.,gmsh_h5_fileId,dsetname2read,&
                           numMapFacesSrl,numElemsInRank,listElemsInRank,numBoundElemsMpiRank,listBoundElemsMpiRank,connecInRank_i8,&
                           numMapFacesInRank,listElemsMapFacesInRank,listMapFacesInRank,connecMapFacesInRank_i8)

      call read_mapped_dir_and_gap_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,mapFaceDir,mapFaceGap)

      call link_mapfaces_and_perfaces_in_parallel(mporder,mnpbou,mapFaceDir,mapFaceGap,&
               numMapFacesSrl,numElemsInRank,numMapFacesInRank,numNodesInRank,numPerFacesInRank,listElemsInRank,&
               listElemsMapFacesInRank,listMapFacesInRank,connecMapFacesInRank_i8,&
               listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8,&
               listNodesInRank_i8,coordNodesInRank,linkedPerToMapElemsAndFacesSrl,mapFaceGapCoord)

      !changing mapFaceDir depeneding on the mapFaceGap (pos or neg)
      if(mapFaceGap<0.0_rp) mapFaceDir = -1*mapFaceDir

      !--------------------------------------------------------------------------------------------------------------------
      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),1,mpi_datatype_real8,MPI_MAX,0,app_comm,mpi_err)

      if(mpi_rank.eq.0) then
         write(*,*) '   1.Time Read Map Faces',elapsed_time(1)
      end if
      !-----------------------------------------------

   end subroutine read_mapped_faces_from_gmsh_h5_file_in_parallel

   subroutine read_mapped_dir_and_gap_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,mapFaceDir,mapFaceGap)
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      integer(4),intent(out) :: mapFaceDir
      real(8),intent(out) :: mapFaceGap

      character(128) :: dsetname
      integer(hsize_t), dimension(1) :: ms_dims
      integer(hid_t) :: dtype
      integer(4) :: h5err
      integer(HSSIZE_T), dimension(1) :: ms_offset
      integer(8),allocatable :: aux_array_i8(:)
      real(8),allocatable :: aux_array_r8(:)

      dtype = h5_datatype_int8
      ms_dims(1) = 1
      ms_offset(1) = 0
      allocate(aux_array_i8(1))

      dsetname = '/mappedInlet/dir'
      call read_dataspace_1d_int8_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_i8)
      mapFaceDir=int(aux_array_i8(1),4)

      deallocate(aux_array_i8)
      allocate(aux_array_r8(1))
      dtype = h5_datatype_real8

      dsetname = '/mappedInlet/gap'
      call read_dataspace_1d_fp64_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_dims,ms_offset,aux_array_r8)
      mapFaceGap=aux_array_r8(1)

      deallocate(aux_array_r8)

      if(mpi_rank.eq.0) write(*,*) ' ## mapFaceDir',mapFaceDir,'mapFaceGap',mapFaceGap

   end subroutine read_mapped_dir_and_gap_from_gmsh_h5_file_in_parallel

   subroutine read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,checkForDuplicatedLinkedFaces,gmsh_h5_fileId,dsetname,&
                     numFacesSrl,numElemsInRank,listElemsInRank,numBoundElemsInRank,listBoundElemsInRank,&
                     connecInRank_i8,numFacesInRank,listElemsFacesInRank,listFacesInRank,connecFacesInRank_i8)
      implicit none
      integer(4), intent(in) :: mporder,mnnode,mnpbou
      logical, intent(in) :: checkForDuplicatedLinkedFaces
      integer(hid_t), intent(in) :: gmsh_h5_fileId
      character(*), intent(in) :: dsetname
      integer(4), intent(in) :: numFacesSrl,numElemsInRank,listElemsInRank(numElemsInRank)
      integer(4), intent(in) :: numBoundElemsInRank,listBoundElemsInRank(numBoundElemsInRank)
      integer(8), intent(in) :: connecInRank_i8(numElemsInRank,mnnode)
      integer(4), intent(out) :: numFacesInRank
      integer(4), allocatable, intent(inout) :: listElemsFacesInRank(:,:),listFacesInRank(:)
      integer(8), allocatable, intent(inout) :: connecFacesInRank_i8(:,:)

      integer(4) :: iFace,iChunk,iFaceG,numFacesToRead,iElem,iVert,ind_gmsh,iElemL,iElemG,iAux,jAux,nodeCnt
      integer(4) :: faces2readInChunk,maxFaces2read,face2sing=50000
      integer(8) :: iNodeG_inFace,iNodeG_inElem,iFaceNodes_i8(mnpbou)
      integer(8),allocatable :: auxFacesInRank_i8(:,:)
      logical :: vertexFound

      integer(4),allocatable :: arrayFacesInRank(:),arrayFacesLinkedInRank(:),arrayFacesLinkedTotal(:),matrixFacesLinkedTotal(:,:)
      integer(4) :: iMpiRank,factor_comms_mpi
      integer(4),dimension(0:mpi_size-1) :: vecNumFacesLinkedInRank,counts_recv_mpi,displacements_mpi!,vecNumPerElemsInRank
      integer(4) :: numFacesLinkedInRank,numFacesLinkedTotal,counts_send_mpi

      integer(4),parameter :: numQuadVert=4,numHexaVert=8
      integer(8) :: f_iNodeG_i8(numQuadVert),e_iNodeG_i8(numHexaVert) !for the vertex of squares of the face and the element
      integer(4) :: gmshQuadVertInd(numQuadVert),gmshHexVertInd(numHexaVert)

      integer(4),parameter :: ms_rank=2,ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dtype

      real(8), dimension(3) :: start_time,end_time,elapsed_time,elapsed_time_m

      !for the moment I will use this auxiliar array, but if we see that is allocating too much memory
      !maybe we can try another strategy, without allocating an array of numBoundFaceSrl
      !but doing a double loop
      !1st loop: determine the bounds that are in the msh rank
      !2nd loop: store only the nodes in the rank
      !OR ANOTHER POSSIBLE STRATEGY:
      ! Determine a maxNumBoundsToRead and then read by chunks
      !

      if(mpi_rank.eq.0) write(*,*) ' # Reading faces & matching elements for dsetname',dsetname

      allocate(listElemsFacesInRank(numElemsInRank,maxBoundsPerElem))
      listElemsFacesInRank(:,:) = 0

      !---------------------------------------------------

      !do the reading by chunks-------
      maxFaces2read =1000

      dtype = h5_datatype_int8
      ms_dims(1) = int(mnpbou,hsize_t)
      ms_dims(2) = 0
      ms_offset(1) = 0
      ms_offset(2) = 0
      !do iChunk=1,numChunks

      faces2readInChunk=numFacesSrl
      allocate(auxFacesInRank_i8(mnpbou+2,faces2readInChunk))
      auxFacesInRank_i8(:,:) = 0

      ms_dims(2) = int(faces2readInChunk,hsize_t)

      start_time(1) = MPI_Wtime()

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      call h5dread_f(dset_id,dtype,auxFacesInRank_i8(1:mnpbou,1:faces2readInChunk),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)

      call get_gmshQuadHOVertIndex(mporder,gmshQuadVertInd)
      call get_gmshHexHOVertIndex(mporder,gmshHexVertInd)

      start_time(2) = MPI_Wtime()

      numFacesLinkedInRank=0

      !read the chunk
      iChunk=0
      do iFace=1,faces2readInChunk
         if((mpi_rank.eq.0).and.(mod(iFace,face2sing).eq.0)) write(*,*) '  - iFace',iFace,'/',faces2readInChunk

         iFaceG=iFace !REVISAR AIXO QUAN IMPLEMENTI ELS CHUNKS A FULL
      !do iFace = 1,numFacesSrl
         !ms_offset(2) = iFace-1
         !call read_dataspace_int4_hyperslab_parallel(gmsh_h5_fileId,dsetname,ms_rank,ms_dims,ms_offset,iFaceNodes)
         !-----------------------------------------------------
         iFaceNodes_i8(:)=auxFacesInRank_i8(1:mnpbou,iFace)
         !fill the corners of the face to check
         do iVert=1,numQuadVert
            ind_gmsh = gmshQuadVertInd(iVert) !gmsh2ij_vertices(iVert)
            f_iNodeG_i8(iVert) = iFaceNodes_i8(ind_gmsh)
         end do
         !write(*,*) '[',mpi_rank,']iFace',iFace,' -> f_iNodeG ',f_iNodeG_i8(:)

         !Search to which element this face belongs
         elemLoop : do iElem=1,numBoundElemsInRank!numElemsInRank
            nodeCnt=0
            iElemL=listBoundElemsInRank(iElem)
            iElemG=listElemsInRank(iElemL)
            !fill the corners of the element
            do iVert=1,numHexaVert
               ind_gmsh = gmshHexVertInd(iVert)
               e_iNodeG_i8(iVert) = connecInRank_i8(iElemL,ind_gmsh)
            end do

            fLoop: do iVert=1,numQuadVert
               vertexFound=.false.
               iNodeG_inFace = f_iNodeG_i8(iVert)
               eLoop: do jAux=1,numHexaVert
                  iNodeG_inElem = e_iNodeG_i8(jAux)
                  if(iNodeG_inFace .eq. iNodeG_inElem) then
                     nodeCnt=nodeCnt+1
                     vertexFound=.true.
                     exit eLoop
                  end if
               end do eLoop
               if(.not.(vertexFound)) exit fLoop
            end do fLoop
            if(nodeCnt.ge.4) then
               numFacesLinkedInRank=numFacesLinkedInRank+1
               !write(*,*) '[',mpi_rank,']iFace',iFace,' -> elem ',iElemG
               auxFacesInRank_i8(mnpbou+1,iFace)=iFaceG
               auxFacesInRank_i8(mnpbou+2,iFace)=iElemL
               exit elemLoop !only adding to first elem (in case inner face)
#if 0
               addLoop : do iAux = 1,maxBoundsPerElem
                  if(listElemsFacesInRank(iElem,iAux).eq.0) then
                     !write(*,*) '[',mpi_rank,'] adding iFace',iFace,'to elem',iElemG
                     listElemsFacesInRank(iElem,iAux) = iFace

                     !exit addLoop
                     exit elemLoop
                  end if
               end do addLoop
#endif
            end if

         end do elemLoop
      end do
      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      if(checkForDuplicatedLinkedFaces) then 
         if(mpi_rank.eq.0) write(*,*) 'Checking for duplicated linked faces (for internal faces)'
         !only to be used when reading initially the periodic faces!!! when mesh is partitioned this must be FALSE
         !write(*,*) '[',mpi_rank,']numFacesLinkedInRank',numFacesLinkedInRank
         !---- --------------------------------------------------------------------------------------

         !WHEN WE NEED TO CHECK THAT IF INTERNAL FACES DUE TO INITIAL PARTITIONING, THIS ONE APPEARS TWICE, LINKED TO TWO ELEMENTS
         !If so, force assigning it to element with smaller ID
         call MPI_Allgather(numFacesLinkedInRank,1,mpi_datatype_int4,vecNumFacesLinkedInRank(0),1,mpi_datatype_int4,app_comm,mpi_err)

         numFacesLinkedTotal = 0
         do iMpiRank=0,mpi_size-1
            numFacesLinkedTotal = numFacesLinkedTotal + vecNumFacesLinkedInRank(iMpiRank)
         end do
         !write(*,*) '[',mpi_rank,']numFacesLinkedTotal',numFacesLinkedTotal

         !AQUI POSAR IF, si numFacesLinkedTotal es com el serial, ja esta, si es mes, s'ha de 'separar el gra de la palla'
         if(numFacesLinkedTotal.eq.numFacesSrl) then !no dupicated faces, we can go on
            numFacesInRank=numFacesLinkedInRank
         else
            if(mpi_rank.eq.0) write(*,*) 'Detected INNER Faces! numFacesLinkedTotal',numFacesLinkedTotal,'numFacesSrl',numFacesSrl,'-> Doing smart delete algorithm!'

            factor_comms_mpi = 2
            counts_send_mpi = numFacesLinkedInRank*factor_comms_mpi
            counts_recv_mpi(0) = vecNumFacesLinkedInRank(0)*factor_comms_mpi
            displacements_mpi(0) = 0
            do iMpiRank=1,mpi_size-1
               counts_recv_mpi(iMpiRank)   = vecNumFacesLinkedInRank(iMpiRank)*factor_comms_mpi
               displacements_mpi(iMpiRank) = displacements_mpi(iMpiRank-1) + counts_recv_mpi(iMpiRank-1)
               !write(*,*) 'counts(',iMpiRank,')',counts_recv_mpi(iMpiRank),'disp',displacements_mpi(iMpiRank)
            end do

            allocate( arrayFacesLinkedInRank(numFacesLinkedInRank*factor_comms_mpi))
            allocate( arrayFacesLinkedTotal(numFacesLinkedTotal*factor_comms_mpi))

            arrayFacesLinkedInRank(:) = 0
            arrayFacesLinkedTotal(:) = 0

            jAux=0
            do iFace=1,faces2readInChunk
               iFaceG = auxFacesInRank_i8(mnpbou+1,iFace)
               if(iFaceG.ne.0) then
                  jAux=jAux+1
                  arrayFacesLinkedInRank(jAux) = iFaceG
                  jAux=jAux+1
                  iElem = auxFacesInRank_i8(mnpbou+2,iFace)
                  arrayFacesLinkedInRank(jAux) = iElem
               end if
            end do

            call MPI_Allgatherv(arrayFacesLinkedInRank,counts_send_mpi,mpi_datatype_int4,arrayFacesLinkedTotal,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)

            deallocate(arrayFacesLinkedInRank)
            allocate(matrixFacesLinkedTotal(numFacesLinkedTotal,3))

            iMpiRank=0
            jAux=0
            do iAux=1,numFacesLinkedTotal
               !if(mpi_rank.eq.0) write(*,*) 'jAux',jAux,'impirank',iMpiRank
               jAux=jAux+1
               matrixFacesLinkedTotal(iAux,1) = arrayFacesLinkedTotal(jAux)
               jAux=jAux+1                
               matrixFacesLinkedTotal(iAux,2) = arrayFacesLinkedTotal(jAux)
               matrixFacesLinkedTotal(iAux,3) = iMpiRank

               !if(mpi_rank.eq.2) write(*,*) 'iAux',iAux,'iMpiRank',matrixFacesLinkedTotal(iAux,:)
               if(jAux.ge.counts_recv_mpi(iMpiRank)+displacements_mpi(iMpiRank)) iMpiRank = iMpiRank + 1
            end do

            deallocate(arrayFacesLinkedTotal)

            call quicksort_matrix_int4(matrixFacesLinkedTotal,1)

            !do iAux=1,numFacesLinkedTotal
            !   if(mpi_rank.eq.2) write(*,*) 'iAux',iAux,'iMpiRank',matrixFacesLinkedTotal(iAux,:)
            !end do

            !counting the number of non-duplicated faces (also the duplicated ones)
            allocate(arrayFacesLinkedInRank(numFacesLinkedInRank)) !the max number of faces in the rank

            jAux=0
            numFacesInRank=0

            !fer el loop de 3/4 fins numFacesLinkedTotal
            iFace=1
            do while (iFace<numFacesLinkedTotal)
               if(matrixFacesLinkedTotal(iFace,1).eq.matrixFacesLinkedTotal(iFace+1,1)) then
                  if(matrixFacesLinkedTotal(iFace,2)<matrixFacesLinkedTotal(iFace+1,2)) then
                     if(mpi_rank.eq.matrixFacesLinkedTotal(iFace,3)) then
                        numFacesInRank=numFacesInRank+1
                        arrayFacesLinkedInRank(numFacesInRank) = matrixFacesLinkedTotal(iFace,1)
                     !write(*,*) '[',mpi_rank,']addsFace',arrayFacesLinkedInRank(numFacesInRank),'iFace(a)',iFace
                     end if
                  else
                     if(mpi_rank.eq.matrixFacesLinkedTotal(iFace+1,3)) then
                        numFacesInRank=numFacesInRank+1
                        arrayFacesLinkedInRank(numFacesInRank) = matrixFacesLinkedTotal(iFace+1,1)
                     !write(*,*) '[',mpi_rank,']addsFace',arrayFacesLinkedInRank(numFacesInRank),'iFace(b)',iFace
                     end if
                  end if
                  jAux=jAux+1
                  iFace=iFace+2 !jump to the following, face cannot be shared by more than two elems
               else
                  if(mpi_rank.eq.matrixFacesLinkedTotal(iFace,3)) then
                     numFacesInRank=numFacesInRank+1
                     arrayFacesLinkedInRank(numFacesInRank) = matrixFacesLinkedTotal(iFace,1)
                     !write(*,*) '[',mpi_rank,']addsFace',arrayFacesLinkedInRank(numFacesInRank),'iFace(c)',iFace
                  end if
                  iFace=iFace+1
               end if
               !if(mpi_rank.eq.0) write(*,*) 'iAux',iAux,'iMpiRank',matrixFacesLinkedTotal(iAux,:)
            end do

            !if the last face is not duplicated, must be added to its processor
            iFace=numFacesLinkedTotal
            if(matrixFacesLinkedTotal(iFace,1).ne.matrixFacesLinkedTotal(iFace-1,1)) then
               if(mpi_rank.eq.matrixFacesLinkedTotal(iFace,3)) then
                  numFacesInRank=numFacesInRank+1
                  arrayFacesLinkedInRank(numFacesInRank) = matrixFacesLinkedTotal(iFace,1)
                  !write(*,*) '[',mpi_rank,']addsFace',arrayFacesLinkedInRank(numFacesInRank),'iFace(end)',iFace
               end if
            end if

            write(*,*) '[',mpi_rank,'] # total duplicated:',jAux

            allocate(arrayFacesInRank(numFacesInRank))

            arrayFacesInRank(:) = arrayFacesLinkedInRank(1:numFacesInRank)
            call quicksort_array_int4(arrayFacesInRank)

            deallocate(arrayFacesLinkedInRank)

            jAux=0
            do iFace=1,faces2readInChunk
               iFaceG = auxFacesInRank_i8(mnpbou+1,iFace)
               if(iFaceG.ne.0) then
                  iAux = binarySearch_int_i4(arrayFacesInRank,iFaceG)
                  if(iAux.eq.0) then
                     !write(*,*) '[',mpi_rank,']deleteFace',iFace
                     auxFacesInRank_i8(:,iFace) = 0
                     jAux=jAux+1
                     !write(*,*) 'face',iFaceG,'now does not exists in rank',mpi_rank
                  end if
               end if
            end do
            write(*,*) '[',mpi_rank,'] # deleted(in rank):',jAux,'-> numFacesInRank',numFacesInRank,'numFacesLinkedInRank',numFacesLinkedInRank

            deallocate(arrayFacesInRank)

         end if

         !-------------------------------------------------------------------------------------------
         call MPI_Barrier(app_comm,h5err)
         call MPI_Allgather(numFacesInRank,1,mpi_datatype_int4,vecNumFacesLinkedInRank(0),1,mpi_datatype_int4,app_comm,mpi_err)

         numFacesLinkedTotal = 0
         do iMpiRank=0,mpi_size-1
            numFacesLinkedTotal = numFacesLinkedTotal + vecNumFacesLinkedInRank(iMpiRank)
         end do
         !write(*,*) '[',mpi_rank,']numFacesLinkedTotal',numFacesLinkedTotal

         if(numFacesLinkedTotal.ne.numFacesSrl) then
            write(*,*) 'ERROR! Crashing in read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel! EXIT!'
            write(*,*) 'numFacesLinkedTotal',numFacesLinkedTotal,'NOT EQUAL TO numFacesSrl',numFacesSrl
            call MPI_Abort(app_comm,-1,mpi_err)
         end if
      else
         numFacesInRank=numFacesLinkedInRank
      end if
!-------------------------------------------------------------------------------------------

      allocate(listFacesInRank(numFacesInRank))
      allocate(connecFacesInRank_i8(numFacesInRank,mnpbou))

      start_time(3) = MPI_Wtime()

      iAux=0
      do iFace=1,faces2readInChunk
         iFaceG = auxFacesInRank_i8(mnpbou+1,iFace)
         if(iFaceG.ne.0) then
            iAux=iAux+1
            listFacesInRank(iAux)        = iFaceG
            connecFacesInRank_i8(iAux,:) = auxFacesInRank_i8(1:mnpbou,iFace)
            iElem                        = auxFacesInRank_i8(mnpbou+2,iFace)
            addLoop : do jAux = 1,maxBoundsPerElem
               if(listElemsFacesInRank(iElem,jAux).eq.0) then
                  !write(*,*) '[',mpi_rank,'] adding iFace',iFace,'to elem',iElem,'iAux',iAux
                  listElemsFacesInRank(iElem,jAux) = iFaceG
                  exit addLoop
               end if
            end do addLoop

         end if
      end do

      end_time(3) = MPI_Wtime()
      elapsed_time(3) = end_time(3) - start_time(3)

      deallocate(auxFacesInRank_i8)

      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),3,mpi_datatype_real8,MPI_MAX,0,app_comm,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '     1.Bounds(readHDF5)',elapsed_time_m(1)
         write(*,*) '     2.Bounds(linkElemAndFace)',elapsed_time_m(2)
         write(*,*) '     3.Bounds(genListAndConn)',elapsed_time_m(3)
      end if

   end subroutine read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel

   subroutine read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname,numBoundsInRank,listBoundFacesInRank,boundFacesCodesInRank,maxBoundCodeInRank)
      implicit none
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      character(*),intent(in) :: dsetname
      integer(4),intent(in) :: numBoundsInRank,listBoundFacesInRank(numBoundsInRank)
      integer(4),intent(out) :: maxBoundCodeInRank
      integer,allocatable,intent(inout) :: boundFacesCodesInRank(:)

      integer(8),allocatable :: boundFacesCodesInRank_i8(:)
      integer(4) :: iBound

      !-----------------------------------------------------------
      integer(4),parameter :: ms_rank=1,ds_rank=1
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      integer(hsize_t) :: ms_numElems
      integer(hid_t) :: dtype
      integer(hsize_t),allocatable :: ms_coords(:,:)
      !------------------------------------------------------------

      allocate(boundFacesCodesInRank(numBoundsInRank))
      allocate(boundFacesCodesInRank_i8(numBoundsInRank))
      boundFacesCodesInRank(:) = 0
      boundFacesCodesInRank_i8(:) = 0
      !---------------------------------------------------------------
      dtype = h5_datatype_int8

      if(numBoundsInRank.ne.0) then
         ms_dims(1) = int(numBoundsInRank,hsize_t)
      else
         ms_dims(1) = 1 !dummy acces case numBoundsInRank=0
      end if
      ms_numElems = ms_dims(1)
      allocate(ms_coords(ms_rank,ms_dims(1)))

      !default
      ms_coords(1,:) = 1

      do iBound=1,numBoundsInRank
         ms_coords(1,iBound) = listBoundFacesInRank(iBound)
      end do

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err)
      call h5dread_f(dset_id,dtype,boundFacesCodesInRank_i8,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      deallocate(ms_coords)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      !---------------------------------------------------------------

      boundFacesCodesInRank(:) = int(boundFacesCodesInRank_i8(:),4)
      deallocate(boundFacesCodesInRank_i8)

      maxBoundCodeInRank = 0
      do iBound = 1,numBoundsInRank
         maxBoundCodeInRank = max(maxBoundCodeInRank,boundFacesCodesInRank(iBound))
      end do
      !write(*,*) '(',mpi_rank,')maxBC',maxBoundCodeInRank,'numBiR',numBoundsInRank

   end subroutine read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel

   subroutine read_periodic_links_from_gmsh_h5_file_in_parallel(mporder,mnpbou,gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,numPerFacesInRank,&
               listElemsInRank,listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8,numLinkedPerElemsSrl,linkedPerElemsAll,&
               numPerElemsSrl,listPerElemsAll,numMasSlaNodesSrl,masSlaNodesAll_i8)
      implicit none
      integer(4),intent(in) :: mporder,mnpbou
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      integer(4),intent(in) :: numPerFacesSrl,numPerLinkedNodesSrl,numElemsInRank,numPerFacesInRank
      integer(4),intent(in) :: listElemsInRank(numElemsInRank),listElemsPerFacesInRank(numElemsInRank,maxBoundsPerElem),listPerFacesInRank(numPerFacesInRank)
      integer(8),intent(in) :: connecPerFacesInRank_i8(numPerFacesInRank,mnpbou)
      integer(4),intent(out) :: numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4),allocatable,intent(inout) :: linkedPerElemsAll(:,:),listPerElemsAll(:)
      integer(8),allocatable,intent(inout) :: masSlaNodesAll_i8(:,:)

      integer(4), allocatable ::listPerElemsInRank(:)

      character(128) :: dsetname
      integer(4),allocatable :: perFaceElemInRank(:),perFaceElemAll(:),listPerFacesAll(:)
      integer(8),allocatable :: refNodesPerFacesInRank_i8(:),refNodesPerFacesAll_i8(:)
      logical,allocatable :: masterPerFacesInRank(:),masterPerFacesAll(:)
      integer(4),dimension(0:mpi_size-1) :: vecNumPerFacesInRank,vecNumPerElemsInRank

      integer(4),parameter :: numQuadVert=4
      integer(4):: gmshQuadInnerVertInd(numQuadVert)

      integer(4) :: iLink,iElem,iPerFace,iFaceG,iBound,ind_gmsh,iElemG,iPos,iPosFace
      integer(4) :: iFaceMaster,iFaceSlave,iFaceGmstr,iFaceGslv,elemGidPerFaceM,elemGidPerFaceS
      integer(4) :: iAux,jAux,iVert,jVert
      integer(4) :: vertNodeCnt,linkedPerElemsCnt,numPerElemsInRank
      integer(8) :: iNodeG,iNodeGtoReplace,iNodeGM,iNodeGS
      logical :: vertNodeFound
      integer(8),allocatable :: matPerLinkNodes_i8(:,:),matPerLinkNodesT_i8(:,:)

      integer(4),parameter :: ms_rank=2,ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dtype

      integer(4) :: memPos,memSize,window_id1,window_id2,mpiRankTrgt,iMpiRank
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size1,win_buffer_size2,trgt_disp

      real(8), dimension(8) :: start_time,end_time,elapsed_time,elapsed_time_m
      !----------------------------------------------------------------------------------------------------------
      !  0. Allocate vectors

      allocate(refNodesPerFacesInRank_i8(numPerFacesInRank*numQuadVert))
      allocate(perFaceElemInRank(numPerFacesInRank))
      allocate(masterPerFacesInRank(numPerFacesInRank))
      allocate(refNodesPerFacesAll_i8(numPerFacesSrl*numQuadVert))
      allocate(perFaceElemAll(numPerFacesSrl))
      allocate(masterPerFacesAll(numPerFacesSrl))
      allocate(listPerFacesAll(numPerFacesSrl))

      !$acc kernels
      refNodesPerFacesInRank_i8(:) = -1
      refNodesPerFacesAll_i8(:) = -1
      perFaceElemInRank(:) = -1
      perFaceElemAll(:) = -1
      masterPerFacesInRank(:) = .true.
      masterPerFacesAll(:) = .true.
      !$acc end kernels

      !----------------------------------------------------------------------------------------------------------

      call MPI_Allgather(numPerFacesInRank,1,mpi_datatype_int,vecNumPerFacesInRank,1,mpi_datatype_int,app_comm,mpi_err)
      !write(*,*) '2.vecNumPerFaces',vecNumPerFacesInRank(:)
      !write(*,*) '[',mpi_rank,']',listPerFacesInRank(:)
      !if(mpi_rank.eq.0) write(*,*) '1.refNodesPerFacesAll',refNodesPerFacesAll_i8(:)
      !if(mpi_rank.eq.0) write(*,*) '1.masterPerFacesAll',masterPerFacesAll(:)

      !----------------------------------------------------------------------------------------------------------
      !  2. Second read the periodic links from h5 file
      if(mpi_rank.eq.0) write(*,*) "  |-> Reading periodic links from h5 file..."
      start_time(2) = MPI_Wtime()
      !IMPORTANT: for the moment I will this auxiliar matrix with the full size numPerLinkedNodesSrl for performance reasons
      !if in a future this can be a memory bound problem (which should not be, since the numPerLinkedNodes should not be very big compared with the inner nodes)
      !rethink the method, maybe do the reading line by line of the hdf5 for all the faces or do the method by chunks

      allocate(matPerLinkNodes_i8(numPerLinkedNodesSrl,2))
      allocate(matPerLinkNodesT_i8(2,numPerLinkedNodesSrl))

      dsetname = 'periodicLinks'

      dtype = h5_datatype_int8
      ms_dims(1) = 2
      ms_dims(2) = int(numPerLinkedNodesSrl,hsize_t)
      ms_offset(1) = 0
      ms_offset(2) = 0

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      call h5dread_f(dset_id,dtype,matPerLinkNodesT_i8(1:2,1:numPerLinkedNodesSrl),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      do iLink=1,numPerLinkedNodesSrl
         matPerLinkNodes_i8(iLink,:) = matPerLinkNodesT_i8(:,iLink)
      end do

      deallocate(matPerLinkNodesT_i8)

      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

      call get_gmshQuadHOInnerVertIndex(mporder,gmshQuadInnerVertInd)

      !----------------------------------------------------------------------------------------------------------
      !  3. Fill the refNodesPerFacesAll_i8 matrix
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling refNodesPerFacesAll..."
      start_time(3) = MPI_Wtime()
      call quicksort_matrix_int8(matPerLinkNodes_i8,1)

      do iPerFace=1,numPerFacesInRank
         perFaceLoop : do iVert=1,numQuadVert
            ind_gmsh = gmshQuadInnerVertInd(iVert)!gmsh2ij_vertInnerNodes(iVert)
            iNodeG = connecPerFacesInRank_i8(iPerFace,ind_gmsh)
            iAux = (iPerFace-1)*numQuadVert + iVert
            iPos = binarySearch_int_i8(matPerLinkNodes_i8(:,1),iNodeG)
            !write(*,*) 'ipf',iPerFace,'iAux',iAux,'iNodeG',iNodeG,'iPos',iPos
            if(iPos.ne.0) then !found
               refNodesPerFacesInRank_i8(iAux) = iNodeG
            else
               masterPerFacesInRank(iPerFace) = .false.
               exit perFaceLoop
            end if
         end do perFaceLoop
      end do

      !write(*,*) '[',mpi_rank,']list',listPerFacesInRank(:),'masterPerFace',masterPerFacesInRank(:)
      !then, the one not founds in master, must be in the 'slave column'

      call quicksort_matrix_int8(matPerLinkNodes_i8,2)

      do iPerFace=1,numPerFacesInRank
         if(.not.(masterPerFacesInRank(iPerFace))) then
            do iVert=1,numQuadVert
               ind_gmsh = gmshQuadInnerVertInd(iVert)!gmsh2ij_vertInnerNodes(iVert)
               iNodeG = connecPerFacesInRank_i8(iPerFace,ind_gmsh)
               iAux = (iPerFace-1)*numQuadVert + iVert
               iPos = binarySearch_int_i8(matPerLinkNodes_i8(:,2),iNodeG)
               refNodesPerFacesInRank_i8(iAux) = matPerLinkNodes_i8(iPos,1)
            end do
         end if
      end do
      end_time(3) = MPI_Wtime()
      elapsed_time(3) = end_time(3) - start_time(3)
      if(mpi_rank.eq.0) write(*,*) "  |-> Looking for and setting the periodic master nodes..."
      start_time(4) = MPI_Wtime()

      !do iPerFace=1,numPerFacesInRank
      !   write(*,*) '[',mpi_rank,']perFace(',listPerFacesInRank(iPerFace),'own',ownedPerFacesInRank(iPerFace,:),'oth',otherPerFacesInRank(iPerFace,:),'M/S',masterPerFacesInRank(iPerFace)
      !end do

      do iLink=1,numPerLinkedNodesSrl
         iNodeG = matPerLinkNodes_i8(iLink,1)
         call find_masterNode_recursive(iNodeG,numPerLinkedNodesSrl,matPerLinkNodes_i8,iNodeGtoReplace)

         if(iNodeGtoReplace.ne.0) then
            matPerLinkNodes_i8(iLink,1) = iNodeGtoReplace
         end if
      end do

      numMasSlaNodesSrl = numPerLinkedNodesSrl
      do iLink=1,numPerLinkedNodesSrl-1
         if((matPerLinkNodes_i8(iLink,1).eq.matPerLinkNodes_i8(iLink+1,1)).and.(matPerLinkNodes_i8(iLink,2).eq.matPerLinkNodes_i8(iLink+1,2))) then
            numMasSlaNodesSrl = numMasSlaNodesSrl-1
         end if
      end do
      !write(*,*) 'numMasSlaNodesSrl',numMasSlaNodesSrl,'numPLN',numPerLinkedNodesSrl

      allocate(masSlaNodesAll_i8(numMasSlaNodesSrl,2))

      iAux  = 1
      iLink = 1
      masSlaNodesAll_i8(iAux,:) = matPerLinkNodes_i8(iLink,:)
      !write(*,*) 'iAux',iAux,'msns',masSlaNodesAll(iAux,:)
      do iLink=2,numPerLinkedNodesSrl
         if(.not.((masSlaNodesAll_i8(iAux,1).eq.matPerLinkNodes_i8(iLink,1)).and.(masSlaNodesAll_i8(iAux,2).eq.matPerLinkNodes_i8(iLink,2)))) then
            iAux=iAux+1
            masSlaNodesAll_i8(iAux,:) = matPerLinkNodes_i8(iLink,:)
            !write(*,*) 'iAux',iAux,'msns',masSlaNodesAll(iAux,:)
         end if
      end do

      deallocate(matPerLinkNodes_i8)

      !-----------------------------------------------------------------------------------------------------------
      ! X. fill the vector perFaceElemInRank with the local elems
      end_time(4) = MPI_Wtime()
      elapsed_time(4) = end_time(4) - start_time(4)
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling vector perFaceElemInRank with local Elems..."
      start_time(5) = MPI_Wtime()

      numPerElemsInRank=0
      do iElem=1,numElemsInRank
         do iBound=1,maxBoundsPerElem
            iElemG = listElemsInRank(iElem)
            iFaceG = listElemsPerFacesInRank(iElem,iBound)
            iPosFace = binarySearch_int_i4(listPerFacesInRank,iFaceG)
            if(iFaceG.ne.0) then
               perFaceElemInRank(iPosFace) = iElemG
            end if
         end do
         if(listElemsPerFacesInRank(iElem,1).ne.0) then
            numPerElemsInRank=numPerElemsInRank+1
         end if
      end do

      call MPI_Allgather(numPerElemsInRank,1,mpi_datatype_int,vecNumPerElemsInRank,1,mpi_datatype_int,app_comm,mpi_err)

      numPerElemsSrl=0
      do iMpiRank=0,mpi_size-1
         numPerElemsSrl=numPerElemsSrl+vecNumPerElemsInRank(iMpiRank)
      end do

#if 0
      write(*,*) '[',mpi_rank,']nPFiR',numPerFacesInRank,'nPFS',numPerFacesSrl
      write(*,*) '[',mpi_rank,']nPEiR',numPerElemsInRank,'nPES',numPerElemsSrl
      call MPI_Barrier(app_comm,h5err)
      if(mpi_rank.eq.0) then
         write(*,*) 'list(:)',listPerFacesInRank
         write(*,*) 'perFaceElemInRank(:)',perFaceElemInRank
      end if
      call MPI_Barrier(app_comm,h5err)
      if(mpi_rank.eq.3) then
         write(*,*) 'list(:)',listPerFacesInRank
         write(*,*) 'perFaceElemInRank(:)',perFaceElemInRank
      end if
      call MPI_Barrier(app_comm,h5err)
      call MPI_Abort(app_comm,-1,mpi_err)
#endif

      end_time(5) = MPI_Wtime()
      elapsed_time(5) = end_time(5) - start_time(5)
      !----------------------------------------------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling vectors listPerElems (Rank & All)..."
      start_time(6) = MPI_Wtime()


      allocate(listPerElemsInRank(numPerElemsInRank))
      allocate(listPerElemsAll(numPerElemsSrl))
      iAux=0
      do iElem=1,numElemsInRank
         if(listElemsPerFacesInRank(iElem,1).ne.0) then
            iElemG = listElemsInRank(iElem)
            iAux=iAux+1
            listPerElemsInRank(iAux) = iElemG
         end if
      end do
      !write(*,*) '[',mpi_rank,']lPEiR',listPerElemsInRank(:)
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerElemsInRank
      call MPI_Win_create(listPerElemsInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerElemsInRank(iMpiRank)
         call MPI_Get(listPerElemsAll(memPos),memSize,mpi_datatype_int,iMpiRank,trgt_disp,memSize,mpi_datatype_int,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------

      deallocate(listPerElemsInRank)
      call quicksort_array_int4(listPerElemsAll) !i think is not necessary... but just in case
      !write(*,*) '[',mpi_rank,']lPESrl',listPerElemsAll(:)

      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_int8_size*(numPerFacesInRank*numQuadVert)
      call MPI_Win_create(refNodesPerFacesInRank_i8,win_buffer_size1,mpi_int8_size,MPI_INFO_NULL,app_comm,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)*numQuadVert
         call MPI_Get(refNodesPerFacesAll_i8(memPos),memSize,mpi_datatype_int8,iMpiRank,trgt_disp,memSize,mpi_datatype_int8,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesInRank
      call MPI_Win_create(perFaceElemInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)
         call MPI_Get(perFaceElemAll(memPos),memSize,mpi_datatype_int,iMpiRank,trgt_disp,memSize,mpi_datatype_int,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'perFaceElemAll',perFaceElemAll
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesInRank
      call MPI_Win_create(masterPerFacesInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)
         call MPI_Get(masterPerFacesAll(memPos),memSize,mpi_datatype_int,iMpiRank,trgt_disp,memSize,mpi_datatype_int,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'masterPerFacesAll',masterPerFacesAll

      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------------------------------------------------
      win_buffer_size1 = mpi_integer_size*numPerFacesInRank
      call MPI_Win_create(listPerFacesInRank,win_buffer_size1,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id1,mpi_err)
      call MPI_Win_fence(0,window_id1,mpi_err)

      memPos = 1
      trgt_disp = 0
      do iMpiRank=0,mpi_size-1
         memSize   = vecNumPerFacesInRank(iMpiRank)
         call MPI_Get(listPerFacesAll(memPos),memSize,mpi_datatype_int,iMpiRank,trgt_disp,memSize,mpi_datatype_int,window_id1,mpi_err)
         memPos = memPos + memSize
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id1,mpi_err)
      call MPI_Win_free(window_id1,mpi_err)
      !---------------------------------------------------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'listPerFacesAll',listPerFacesAll

      !-------------------------------------------------
      !--------------------------------------------------------------------------------
      !--------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) '2.refNodesPerFacesAll_i8',refNodesPerFacesAll_i8(:)
      !if(mpi_rank.eq.0) write(*,*) '2.masterPerFacesAll',masterPerFacesAll(:)
      !if(mpi_rank.eq.0) write(*,*) '2.perFaceElemAll',perFaceElemAll(:)

      end_time(6) = MPI_Wtime()
      elapsed_time(6) = end_time(6) - start_time(6)
      !----------------------------------------------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) "  |-> Finding the matches: the linked periodic faces!..."
      start_time(7) = MPI_Wtime()
      !------------------------------------------------------------------------------------------
      ! finding the 'matches', ie the linked periodic faces!
      numLinkedPerElemsSrl = numPerFacesSrl/2
      allocate(linkedPerElemsAll(numLinkedPerElemsSrl,2))

      linkedPerElemsCnt=0
      do iFaceMaster=1,numPerFacesSrl
         if(masterPerFacesAll(iFaceMaster)) then
            id2loop: do iFaceSlave=1,numPerFacesSrl
               if(.not.(masterPerFacesAll(iFaceSlave))) then
                  vertNodeCnt=0
                  iVertLoop: do iVert=1,numQuadVert
                     iAux = (iFaceMaster-1)*numQuadVert + iVert
                     vertNodeFound = .false.
                     iNodeGM = refNodesPerFacesAll_i8(iAux)
                     jVertLoop: do jVert=1,numQuadVert
                        jAux = (iFaceSlave-1)*numQuadVert + jVert
                        iNodeGS = refNodesPerFacesAll_i8(jAux)
                        if(iNodeGM.eq.iNodeGS) then
                           vertNodeFound = .true.
                           vertNodeCnt = vertNodeCnt+1
                           exit jVertLoop
                        end if
                     end do jVertLoop
                     if(.not.(vertNodeFound)) exit iVertLoop
                  end do iVertLoop
                  if(vertNodeCnt.eq.4) then
                     !if(mpi_rank.eq.0) write(*,*) '(M)pFId',perFaceIdM,'(S)pFId',perFaceIdS
                     linkedPerElemsCnt=linkedPerElemsCnt+1
                     iFaceGmstr = listPerFacesAll(iFaceMaster)
                     iFaceGslv  = listPerFacesAll(iFaceSlave)
                     !if(mpi_rank.eq.0) write(*,*) 'matching iFaceGmstr',iFaceGmstr,'iFaceGslv',iFaceGslv

                     linkedPerElemsAll(linkedPerElemsCnt,1) = perFaceElemAll(iFaceMaster)!elemGidPerFaceM
                     linkedPerElemsAll(linkedPerElemsCnt,2) = perFaceElemAll(iFaceSlave)!elemGidPerFaceS

                     exit id2loop
                  end if
               end if
            end do id2loop
         end if
      end do

      if(linkedPerElemsCnt.ne.numLinkedPerElemsSrl) then
         write(*,*) 'ERROR! Crashing in read_periodic_links_from_gmsh_h5_file_in_parallel! EXIT!'
         write(*,*) 'linkedPerElemsCnt',linkedPerElemsCnt,'NOT EQUAL TO numLinkedPerElemsSrl',numLinkedPerElemsSrl
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      call quicksort_matrix_int4(linkedPerElemsAll,1)

#if 0
      if(mpi_rank.eq.0) then
         do iPos=1,numLinkedPerElemsSrl
            write(*,*) '0.lPES',linkedPerElemsAll(iPos,:)
         end do
      end if
      call MPI_Barrier(app_comm, mpi_err)
      if(mpi_rank.eq.1) then
         do iPos=1,numLinkedPerElemsSrl
            write(*,*) '1.lPES',linkedPerElemsAll(iPos,:)
         end do
      end if
#endif

      !----------------------------------------------------------------------------------------------------------
      end_time(7) = MPI_Wtime()
      elapsed_time(7) = end_time(7) - start_time(7)
      if(mpi_rank.eq.0) write(*,*) "  |-> Done periodic links stuff!!..."

      deallocate(refNodesPerFacesAll_i8)
      deallocate(listPerFacesAll)
      deallocate(perFaceElemAll)
      deallocate(masterPerFacesAll)

      elapsed_time(1) = end_time(1) - start_time(1)
      elapsed_time(2) = end_time(2) - start_time(2)
      elapsed_time(3) = end_time(3) - start_time(3)
      elapsed_time(4) = end_time(4) - start_time(4)
      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),7,mpi_datatype_real8,MPI_MAX,0,app_comm,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) '     1.PerLinks(s1)',elapsed_time_m(1)
         write(*,*) '     2.PerLinks(readHDF5)',elapsed_time_m(2)
         write(*,*) '     3.refNodesPerFacesAll',elapsed_time_m(3)
         write(*,*) '     4.Doing PerMaster',elapsed_time_m(4)
         write(*,*) '     5.Vecs.perFaceElemAll',elapsed_time_m(5)
         write(*,*) '     6.Vecs.listPerElems',elapsed_time_m(6)
         write(*,*) '     7.FindingMatches',elapsed_time_m(7)
      end if

   end subroutine read_periodic_links_from_gmsh_h5_file_in_parallel

   subroutine link_mapfaces_and_perfaces_in_parallel(mporder,mnpbou,mapFaceDir,mapFaceGap,&
               numMapFacesSrl,numElemsInRank,numMapFacesInRank,numNodesInRank,numPerFacesInRank,listElemsInRank,&
               listElemsMapFacesInRank,listMapFacesInRank,connecMapFacesInRank_i8,&
               listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8,&
               listNodesInRank_i8,coordNodesInRank,linkedPerToMapElemsAndFacesSrl,mapFaceGapCoord)
      implicit none
      integer(4),intent(in) :: mporder,mnpbou
      integer(4),intent(in) :: mapFaceDir,numMapFacesSrl,numElemsInRank,numMapFacesInRank,numNodesInRank,numPerFacesInRank
      integer(4),intent(in) :: listElemsInRank(numElemsInRank),listElemsMapFacesInRank(numElemsInRank,maxBoundsPerElem),listMapFacesInRank(numMapFacesInRank)
      integer(4),intent(in) :: listElemsPerFacesInRank(numElemsInRank,maxBoundsPerElem),listPerFacesInRank(numPerFacesInRank)
      integer(8),intent(in) :: connecMapFacesInRank_i8(numMapFacesInRank,mnpbou),listNodesInRank_i8(numNodesInRank),connecPerFacesInRank_i8(numPerFacesInRank,mnpbou)
      real(8),intent(in)   :: mapFaceGap,coordNodesInRank(numNodesInRank,ndime)
      integer(4),intent(out) :: linkedPerToMapElemsAndFacesSrl(numMapFacesSrl,4)
      real(8),intent(out)  :: mapFaceGapCoord

      real(8) :: avgCoordInRank,avgCoordAll,targetCoord,minTargetCoord,maxTargetCoord
      real(8) :: minCentroidDist,centroidDistAbs,centroidDist(ndime),faceCentroid(ndime)

      integer(4) :: numPerFacesToMapInRank,numPerFacesToMapSrl,numMapElemsInRank
      real(8), allocatable :: coordsInMapFaceInRank(:),coordsInMapFaceSrl(:)

      integer(4),allocatable :: mapFaceElemInRank(:),mapFaceElemSrl(:),listMapFacesSrl(:)
      integer(4),allocatable :: listPerFacesToMapFacesInRank(:,:),listPerElemsToMapElemsInRank(:,:)
      integer(4),allocatable :: auxListPerFacesInRank(:),auxListMapFacesInRank(:),auxListPerElemsInRank(:),auxListMapElemsInRank(:)
      integer(4),allocatable :: auxListPerFacesSrl(:),auxListMapFacesSrl(:),auxListPerElemsSrl(:),auxListMapElemsSrl(:)

      integer(4),dimension(0:mpi_size-1) :: vecNumMapFacesInRank,vecNumPerFacesToMapInRank
      real(8),dimension(0:mpi_size-1) :: vecAvgCoordInRank

      integer(4) :: counts_send_mpi,factor_comms_mpi,iMpiRank
      integer(4),dimension(0:mpi_size-1) :: counts_recv_mpi,displacements_mpi

      integer(4),parameter :: numQuadVert=4
      integer(4):: gmshQuadVertInd(numQuadVert)

      integer(4) :: m,iFace,iElem,iBound,iNode,iRank,idime,minDist_mapFaceG,minDist_mapElemG
      integer(4) :: iLink,iMapFace,iPerFace,iFaceG,ind_gmsh,iElemG,iPosFace,iPerFaceG,iPerElemG
      integer(4) :: iAux,jAux,iVert,jVert
      integer(4) :: vertNodeCnt,linkedMapElemsCnt
      integer(8) :: iNodeG
      integer(8),allocatable :: matMapLinkNodes_i8(:,:),matMapLinkNodesT_i8(:,:)

      real(8), dimension(8) :: start_time,end_time,elapsed_time,elapsed_time_m
      !----------------------------------------------------------------------------------------------------------
      !  0. Allocate vectors

      allocate(mapFaceElemInRank(numMapFacesInRank))
      allocate(mapFaceElemSrl(numMapFacesSrl))
      allocate(listMapFacesSrl(numMapFacesSrl))
      allocate(coordsInMapFaceInRank(numMapFacesInRank*ndime))
      allocate(coordsInMapFaceSrl(numMapFacesSrl*ndime))

      !$acc kernels
      mapFaceElemInRank(:) = -1
      mapFaceElemSrl(:) = -1
      listMapFacesSrl(:) = -1
      coordsInMapFaceInRank(:) = 0.0d0
      coordsInMapFaceSrl(:) = 0.0d0
      !$acc end kernels

      !----------------------------------------------------------------------------------------------------------
      call get_gmshQuadHOVertIndex(mporder,gmshQuadVertInd)

      avgCoordInRank = 0.0d0
      do iMapFace=1,numMapFacesInRank
         faceCentroid(:) = 0.0d0
         do iVert=1,numQuadVert
            ind_gmsh = gmshQuadVertInd(iVert)
            iNodeG = connecMapFacesInRank_i8(iMapFace,ind_gmsh)
            iNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)

            do idime=1,ndime
               faceCentroid(idime) =  faceCentroid(idime) + coordNodesInRank(iNode,idime)
            end do

         end do

         do idime=1,ndime
            faceCentroid(idime) = faceCentroid(idime)/real(numQuadVert,8)
            coordsInMapFaceInRank((iMapFace-1)*ndime+idime) = faceCentroid(idime)
         end do

         avgCoordInRank = (avgCoordInRank*(iMapFace-1) + faceCentroid(mapFaceDir))/real(iMapFace,8)
         !write(*,*) 'iMF',iMapFace,'faceC',faceCentroid(:)
      end do

      call MPI_Allgather(numMapFacesInRank,1,mpi_datatype_int,vecNumMapFacesInRank,1,mpi_datatype_int4,app_comm,mpi_err)
      call MPI_Allgather(avgCoordInRank,1,mpi_datatype_real8,vecAvgCoordInRank,1,mpi_datatype_real8,app_comm,mpi_err)

      avgCoordAll = 0.0d0
      do iRank=0,(mpi_size-1)
         avgCoordAll = avgCoordAll + vecNumMapFacesInRank(iRank)*vecAvgCoordInRank(iRank)
      end do
      avgCoordAll = avgCoordAll/real(numMapFacesSrl,8)

      targetCoord = avgCoordAll - mapFaceGap
      minTargetCoord = targetCoord - abs(mapFaceGap)*0.01
      maxTargetCoord = targetCoord + abs(mapFaceGap)*0.01
      !write(*,*) '[',mpi_rank,']avgCoordAll',avgCoordAll,'numNodesInMapFace',numNodesInMapFaceRawAll
      !write(*,*) '[',mpi_rank,']targetCoord',targetCoord
      mapFaceGapCoord =  0.5*(targetCoord + avgCoordAll) !targetCoord !to comment with oriol
      !-----------------------------------------------------------------------------------------------------------

      !-----------------------------------------------------------------------------------------------------------
      ! X. fill the vector mapFaceElemInRank with the local elems
      if(mpi_rank.eq.0) write(*,*) "  |-> Filling vector mapFaceElemInRank with local Elems..."

      numMapElemsInRank=0
      do iElem=1,numElemsInRank
         do iBound=1,maxBoundsPerElem
            iElemG = listElemsInRank(iElem)
            iFaceG = listElemsMapFacesInRank(iElem,iBound)
            iPosFace = binarySearch_int_i4(listMapFacesInRank,iFaceG)
            if(iFaceG.ne.0) then
               mapFaceElemInRank(iPosFace) = iElemG
            end if
         end do
         if(listElemsMapFacesInRank(iElem,1).ne.0) then
            numMapElemsInRank=numMapElemsInRank+1
         end if
      end do

      if(numMapElemsInRank.ne.numMapFacesInRank) then
         write(*,*) 'ERROR! Crashing in link_mapfaces_and_perfaces_in_parallel! EXIT!'
         write(*,*) 'numMapElemsInRank',numMapElemsInRank,'NOT EQUAL TO numMapFacesInRank',numMapFacesInRank
         call MPI_Abort(app_comm,-1,mpi_err)
      end if
      !-----------------------------------------------------------------------------------------------------------

      factor_comms_mpi = 1
      counts_send_mpi = numMapFacesInRank*factor_comms_mpi
      counts_recv_mpi(0) = vecNumMapFacesInRank(0)*factor_comms_mpi
      displacements_mpi(0) = 0
      do iMpiRank=1,mpi_size-1
         counts_recv_mpi(iMpiRank)   = vecNumMapFacesInRank(iMpiRank)*factor_comms_mpi
         displacements_mpi(iMpiRank) = displacements_mpi(iMpiRank-1) + counts_recv_mpi(iMpiRank-1)
      end do

      call MPI_Allgatherv(listMapFacesInRank,counts_send_mpi,mpi_datatype_int4,listMapFacesSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)
      call MPI_Allgatherv(mapFaceElemInRank,counts_send_mpi,mpi_datatype_int4,mapFaceElemSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)

      factor_comms_mpi = 3
      counts_send_mpi = numMapFacesInRank*factor_comms_mpi
      counts_recv_mpi(0) = vecNumMapFacesInRank(0)*factor_comms_mpi
      displacements_mpi(0) = 0
      do iMpiRank=1,mpi_size-1
         counts_recv_mpi(iMpiRank)   = vecNumMapFacesInRank(iMpiRank)*factor_comms_mpi
         displacements_mpi(iMpiRank) = displacements_mpi(iMpiRank-1) + counts_recv_mpi(iMpiRank-1)
      end do
      call MPI_Allgatherv(coordsInMapFaceInRank,counts_send_mpi,mpi_datatype_real8,coordsInMapFaceSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_real8,app_comm,mpi_err)

      !------------------------------------------------------------------------------------------------------------

      allocate(listPerFacesToMapFacesInRank(numPerFacesInRank,2))
      allocate(listPerElemsToMapElemsInRank(numPerFacesInRank,2))
      listPerFacesToMapFacesInRank(:,:) = 0
      listPerElemsToMapElemsInRank(:,:) = 0 

      numPerFacesToMapInRank = 0

      do iElem=1,numElemsInRank
         do iBound=1,maxBoundsPerElem
            iPerElemG = listElemsInRank(iElem)
            iPerFaceG = listElemsPerFacesInRank(iElem,iBound)
            if(iPerFaceG.ne.0) then
               iPerFace = binarySearch_int_i4(listPerFacesInRank,iPerFaceG)
               faceCentroid(:) = 0.0d0
               do iVert=1,numQuadVert
                  ind_gmsh = gmshQuadVertInd(iVert)
                  iNodeG = connecPerFacesInRank_i8(iPerFace,ind_gmsh)

                  iNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
                  do idime=1,ndime
                     faceCentroid(idime) =  faceCentroid(idime) + coordNodesInRank(iNode,idime)
                  end do
               end do

               do idime=1,ndime
                  faceCentroid(idime) = faceCentroid(idime)/real(numQuadVert,8)
               end do
               if( faceCentroid(mapFaceDir).gt.minTargetCoord .and. faceCentroid(mapFaceDir).lt.maxTargetCoord ) then
                  numPerFacesToMapInRank = numPerFacesToMapInRank + 1

                  minCentroidDist = 1.e9
                  do iMapFace = 1,numMapFacesSrl
                     centroidDistAbs = 0.0_rp
                     do idime=1,ndime
                        centroidDist(idime) = coordsInMapFaceSrl((iMapFace-1)*ndime+idime) - faceCentroid(idime)
                        centroidDistAbs = centroidDistAbs + centroidDist(idime)*centroidDist(idime)
                     end do
                     centroidDistAbs = sqrt(centroidDistAbs)
                     !write(*,*) ' iMP',iMapFace,'faceG',listMapFacesSrl(iMapFace),'centroidDistAbs',centroidDistAbs
                     !write(*,*) ' fCPF',faceCentroid(:)
                     !write(*,*) ' fCMP',coordsInMapFaceSrl((iMapFace-1)*ndime+1),coordsInMapFaceSrl((iMapFace-1)*ndime+2),coordsInMapFaceSrl((iMapFace-1)*ndime+3)

                     if(centroidDistAbs.lt.minCentroidDist) then
                        minCentroidDist = centroidDistAbs
                        minDist_mapFaceG = listMapFacesSrl(iMapFace)
                        minDist_mapElemG = mapFaceElemSrl(iMapFace)
                     end if

                  end do

                  listPerFacesToMapFacesInRank(iPerFace,1) = iPerFaceG
                  listPerFacesToMapFacesInRank(iPerFace,2) = minDist_mapFaceG

                  listPerElemsToMapElemsInRank(iPerFace,1) = iPerElemG
                  listPerElemsToMapElemsInRank(iPerFace,2) = minDist_mapElemG

                  !write(*,*) 'iPF',iPerFace,'iPerFaceG',iPerFaceG,'mapFaceG',minDist_mapFaceG,'iPerElemG',iPerElemG,'mapElemG',minDist_mapElemG,'faceC',faceCentroid(:),'minCD',minCentroidDist
                  !write(*,*) 'iPF',iPerFace,'iFaceG',iFaceG,'faceC',faceCentroid(:),'minCD',minCentroidDist,'mapFaceG',minDist_mapFaceG,'mapElemG',minDist_mapElemG
               end if
            end if
         end do
      end do

      call MPI_Allreduce(numPerFacesToMapInRank,numPerFacesToMapSrl,1,mpi_datatype_int4,MPI_SUM,app_comm,mpi_err)
      call MPI_Allgather(numPerFacesToMapInRank,1,mpi_datatype_int,vecNumPerFacesToMapInRank,1,mpi_datatype_int4,app_comm,mpi_err)
      !write(*,*) 'numPerFacesToMap',numPerFacesToMapInRank,'Total',numPerFacesToMapSrl

      if(numMapFacesSrl.ne.numPerFacesToMapSrl) then
         write(*,*) 'ERROR! Crashing in link_mapfaces_and_perfaces_in_parallel! EXIT!'
         write(*,*) 'numMapFacesSrl',numMapFacesSrl,'NOT EQUAL TO numPerFacesToMapSrl',numPerFacesToMapSrl
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      !----------------------------------------------------------------------------------------------------------------------------

      allocate(auxListPerFacesInRank(numPerFacesToMapInRank),auxListMapFacesInRank(numPerFacesToMapInRank),auxListPerElemsInRank(numPerFacesToMapInRank),auxListMapElemsInRank(numPerFacesToMapInRank))
      allocate(auxListPerFacesSrl(numMapFacesSrl),auxListMapFacesSrl(numMapFacesSrl),auxListPerElemsSrl(numMapFacesSrl),auxListMapElemsSrl(numMapFacesSrl))
      !write(*,*) '[',mpi_rank,']numPerFacesInRank',numPerFacesInRank
      !write(*,*) '[',mpi_rank,']numMapFacesInRank',numMapFacesInRank
      !write(*,*) '[',mpi_rank,']numPerFacesToMapInRank',numPerFacesToMapInRank

      iMapFace=0
      do iPerFace = 1,numPerFacesInRank
         if(listPerFacesToMapFacesInRank(iPerFace,1).ne.0) then
            iMapFace = iMapFace + 1
            !write(*,*) '[',mpi_rank,']iPerFace',iPerFace,'iMapFace',iMapFace

            auxListPerFacesInRank(iMapFace) = listPerFacesToMapFacesInRank(iPerFace,1)
            auxListMapFacesInRank(iMapFace) = listPerFacesToMapFacesInRank(iPerFace,2)
            auxListPerElemsInRank(iMapFace) = listPerElemsToMapElemsInRank(iPerFace,1)
            auxListMapElemsInRank(iMapFace) = listPerElemsToMapElemsInRank(iPerFace,2)            

         end if
      end do

      deallocate(listPerFacesToMapFacesInRank)
      deallocate(listPerElemsToMapElemsInRank)

      factor_comms_mpi = 1
      counts_send_mpi = numPerFacesToMapInRank*factor_comms_mpi
      counts_recv_mpi(0) = vecNumPerFacesToMapInRank(0)*factor_comms_mpi
      displacements_mpi(0) = 0
      do iMpiRank=1,mpi_size-1
         counts_recv_mpi(iMpiRank)   = vecNumPerFacesToMapInRank(iMpiRank)*factor_comms_mpi
         displacements_mpi(iMpiRank) = displacements_mpi(iMpiRank-1) + counts_recv_mpi(iMpiRank-1)
      end do

      call MPI_Allgatherv(auxListPerFacesInRank,counts_send_mpi,mpi_datatype_int4,auxListPerFacesSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)
      call MPI_Allgatherv(auxListMapFacesInRank,counts_send_mpi,mpi_datatype_int4,auxListMapFacesSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)
      call MPI_Allgatherv(auxListPerElemsInRank,counts_send_mpi,mpi_datatype_int4,auxListPerElemsSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)
      call MPI_Allgatherv(auxListMapElemsInRank,counts_send_mpi,mpi_datatype_int4,auxListMapElemsSrl,counts_recv_mpi,displacements_mpi,mpi_datatype_int4,app_comm,mpi_err)

      do iMapFace = 1,numMapFacesSrl
         linkedPerToMapElemsAndFacesSrl(iMapFace,1) = auxListPerElemsSrl(iMapFace)
         linkedPerToMapElemsAndFacesSrl(iMapFace,2) = auxListMapElemsSrl(iMapFace)
         linkedPerToMapElemsAndFacesSrl(iMapFace,3) = auxListPerFacesSrl(iMapFace)
         linkedPerToMapElemsAndFacesSrl(iMapFace,4) = auxListMapFacesSrl(iMapFace)
      end do

      deallocate(auxListPerFacesInRank,auxListMapFacesInRank,auxListPerElemsInRank,auxListMapElemsInRank)
      deallocate(auxListPerFacesSrl,auxListMapFacesSrl,auxListPerElemsSrl,auxListMapElemsSrl)

      call quicksort_matrix_int4(linkedPerToMapElemsAndFacesSrl,1)

   end subroutine link_mapfaces_and_perfaces_in_parallel

   subroutine unfold_linked_elems(numElems2PartInRank,numLinkedPerElemsSrl,numMapFacesSrl,numPerElemsSrl,linkedPerElemsSrl,linkedPerToMapElemsAndFacesSrl,listPerElemsSrl,elemPart)
      implicit none
      integer(4),intent(inout) :: numElems2PartInRank
      integer(4),intent(in) :: numLinkedPerElemsSrl,numMapFacesSrl,numPerElemsSrl
      integer(4),intent(in) :: linkedPerElemsSrl(numLinkedPerElemsSrl,2),linkedPerToMapElemsAndFacesSrl(numMapFacesSrl,4),listPerElemsSrl(numPerElemsSrl)
      integer,allocatable,intent(inout) :: elemPart(:,:)
      integer,allocatable :: elemPartAuxPer(:,:),elemPartAuxMap(:,:),elemPartAuxOrig(:,:)
      integer(4) :: numAdditionalPerElemsInRank,numAdditionalMapElemsInRank,numElems2PartInRankOrig
      integer(4) :: unfoldedPerElems(numPerElemsSrl),unfoldedMapElems(numMapFacesSrl)
      integer(4) :: iElem,iPosMapElem,iElemG,mshRank,iCnt,iStart,iEnd,i,j

      numAdditionalPerElemsInRank=0
      unfoldedPerElems(:)=-1
      unfoldedMapElems(:)=-1

      do iElem=1,numElems2PartInRank
         iElemG  = elemPart(iElem,1)
         mshRank = elemPart(iElem,2)

         call find_linked_elems_inPar(mshRank,iElemG,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,unfoldedPerElems,numAdditionalPerElemsInRank)
      end do

      !do iElemG=1,size(unfoldedPerElems)
      !   iRank = unfoldedPerElems(iElemG)
      !   if(iRank.ge.0) then
      !      write(*,*) '[',mpi_rank,'] unfolded Elem ',iElemG,' rank ', iRank
      !   end if
      !end do
      !write(*,*) 'unfolded[',mpi_rank,']',unfoldedPerElems(:)

      !now we have the 'unfolded elems' and the rank that they belong!
      !we need to add it to the elemPart matrix(:,2)!!!
      !write(*,*) 'numAdditionalPer[',mpi_rank,']',numAdditionalPerElemsInRank
      allocate(elemPartAuxPer(numAdditionalPerElemsInRank,2))
      iCnt=0
      do iElem=1,numPerElemsSrl
         mshRank = unfoldedPerElems(iElem)
         if(mshRank.ge.0) then
            iCnt=iCnt+1
            iElemG = listPerElemsSrl(iElem)
            elemPartAuxPer(iCnt,1) = iElemG
            elemPartAuxPer(iCnt,2) = mshRank
            !write(*,*) '[',mpi_rank,'] unf.Elem ',iElemG,' rank ', mshRank,' i ',i
         end if
      end do

      !-----------------------------------------------------------------------------------------
      !now do the mapped elems
      numAdditionalMapElemsInRank = 0
      do iElem=1,numElems2PartInRank
         iElemG  = elemPart(iElem,1)
         mshRank = elemPart(iElem,2)

         iPosMapElem = binarySearch_int_i4(linkedPerToMapElemsAndFacesSrl(:,1),iElemG)
         if(iPosMapElem.ne.0) then
            numAdditionalMapElemsInRank = numAdditionalMapElemsInRank + 1
            unfoldedMapElems(iPosMapElem)=mshRank
            !write(*,*) '[',mpi_rank,'] iPosMapElem',iPosMapElem,'iElemG',iElemG,'mshRank',mshRank
         end if
      end do

      do iElem=1,numAdditionalPerElemsInRank
         iElemG  = elemPartAuxPer(iElem,1)
         mshRank = elemPartAuxPer(iElem,2)

         iPosMapElem = binarySearch_int_i4(linkedPerToMapElemsAndFacesSrl(:,1),iElemG)
         if(iPosMapElem.ne.0) then
            numAdditionalMapElemsInRank = numAdditionalMapElemsInRank + 1
            unfoldedMapElems(iPosMapElem)=mshRank
            !write(*,*) '[',mpi_rank,'] iPosMapElem',iPosMapElem,'iElemG',iElemG,'mshRank',mshRank
         end if
      end do

      !write(*,*) 'numAdditionalMap[',mpi_rank,']',numAdditionalMapElemsInRank
      allocate(elemPartAuxMap(numAdditionalMapElemsInRank,2))
      iCnt=0
      do iElem=1,numMapFacesSrl
         mshRank = unfoldedMapElems(iElem)
         if(mshRank.ge.0) then
            iCnt=iCnt+1
            iElemG = linkedPerToMapElemsAndFacesSrl(iElem,2)
            elemPartAuxMap(iCnt,1) = iElemG
            elemPartAuxMap(iCnt,2) = mshRank
            !write(*,*) '[',mpi_rank,'] unfMapElem',iElemG,'rank',mshRank,'iCnt',iCnt
         end if
      end do

      !--------------------------------------------------------------------------------------------------

      numElems2PartInRankOrig = numElems2PartInRank
      allocate(elemPartAuxOrig(numElems2PartInRankOrig,2))
      !$acc kernels
      elemPartAuxOrig(:,1)=elemPart(:,1)
      elemPartAuxOrig(:,2)=elemPart(:,2)
      !$acc end kernels

      deallocate(elemPart)
      numElems2PartInRank = numElems2PartInRankOrig + numAdditionalPerElemsInRank + numAdditionalMapElemsInRank
      !write(*,*) '[',mpi_rank,']numElems2PartInRank',numElems2PartInRank
      allocate(elemPart(numElems2PartInRank,2))

      iStart = 1
      iEnd   = numElems2PartInRankOrig
      if(numElems2PartInRankOrig.ge.1) then
         !write(*,*) '[',mpi_rank,']iStart',iStart,'iEnd',iEnd
         elemPart(iStart:iEnd,:) = elemPartAuxOrig(1:numElems2PartInRankOrig,:)
      end if

      iStart = iEnd + 1
      iEnd   = iEnd + numAdditionalPerElemsInRank
      if(numAdditionalPerElemsInRank.ge.1) then
         !write(*,*) '[',mpi_rank,']iStart',iStart,'iEnd',iEnd
         elemPart(iStart:iEnd,:) = elemPartAuxPer(1:numAdditionalPerElemsInRank,:)
      end if

      iStart = iEnd + 1
      iEnd   = iEnd + numAdditionalMapElemsInRank
      if(numAdditionalMapElemsInRank.ge.1) then
         !write(*,*) '[',mpi_rank,']iStart',iStart,'iEnd',iEnd
         elemPart(iStart:iEnd,:) = elemPartAuxMap(1:numAdditionalMapElemsInRank,:)
      end if

      deallocate(elemPartAuxOrig,elemPartAuxPer,elemPartAuxMap)

   end subroutine unfold_linked_elems

   recursive subroutine find_linked_elems_inPar(mshRank,iElemG,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,unfoldedElems,numAddElemsInRank)
      implicit none
      integer, intent(in) :: mshRank,iElemG,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl(numLinkedPerElemsSrl,2),listPerElemsSrl(numPerElemsSrl)
      integer, intent(inout) :: unfoldedElems(numPerElemsSrl),numAddElemsInRank
      integer :: iPosEL1ref,iPosEL1ini,iPosEL1end,iPosEL1,iPosEL2,iElemGLink!,iEL2
      logical :: keepLooping

      iPosEL1ref = binarySearch_int_i4(linkedPerElemsSrl(:,1),iElemG)

      if(iPosEL1ref.ne.0) then   !#checking if the elem is in the linkedElems list
         iPosEL1ini=iPosEL1ref
         iPosEL1end=iPosEL1ref
         keepLooping = .true.
         do while(keepLooping)
            if(iPosEL1ini.eq.1) then
               keepLooping = .false.
            else
               if(iElemG.eq.linkedPerElemsSrl(iPosEL1ini-1,1)) then
                  iPosEL1ini = iPosEL1ini - 1
                  keepLooping = .true.
               else
                  keepLooping = .false.
               end if
            end if
         end do

         keepLooping = .true.
         do while(keepLooping)
            if(iPosEL1end.eq.numLinkedPerElemsSrl) then
               keepLooping = .false.
            else
               if(iElemG.eq.linkedPerElemsSrl(iPosEL1end+1,1)) then
                  iPosEL1end = iPosEL1end + 1
                  keepLooping = .true.
               else
                  keepLooping = .false.
               end if
            end if
         end do
         !write(*,*) '[',mpi_rank,']iELini',iPosEL1ini,'iELend',iPosEL1end

         do iPosEL1=iPosEL1ini,iPosEL1end
            iElemGLink = linkedPerElemsSrl(iPosEL1,2)
            !write(*,*) 'iElemG',iElemG,'iPosEl1',iPosEL1,'El',linkedPerElemsSrl(iPosEL1,1),'+1',linkedPerElemsSrl(iPosEL1+1,1)
            !write(*,*) 'iElemG',iElemG,'linked to', iElemGLink

            iPosEL2 = binarySearch_int_i4(listPerElemsSrl,iElemGLink)
            if(unfoldedElems(iPosEL2).le.0) then
               unfoldedElems(iPosEL2) = mshRank
               numAddElemsInRank=numAddElemsInRank+1
               call find_linked_elems_inPar(mshRank,iElemGLink,numLinkedPerElemsSrl,numPerElemsSrl,linkedPerElemsSrl,listPerElemsSrl,unfoldedElems,numAddElemsInRank)
            end if
         end do

      end if

   end subroutine find_linked_elems_inPar

   recursive subroutine find_masterNode_recursive(iNodeG,numPerLinkedNodesSrl,matPerLinkNodes_i8,iNodeGtoReplace)
      implicit none
      integer(8),intent(in) :: iNodeG,matPerLinkNodes_i8(numPerLinkedNodesSrl,2) !this matPerLinksNodes must be sorted by slave col
      integer(4),intent(in) :: numPerLinkedNodesSrl
      integer(8),intent(out) :: iNodeGtoReplace
      integer(4) :: iPos,iPosNew
      integer(8) :: iNodeGnew

      iNodeGtoReplace = 0
      iPos = binarySearch_int_i8(matPerLinkNodes_i8(:,2),iNodeG)
      if(iPos.ne.0) then
         !if(mpi_rank.eq.0) write(*,*) 'iNG(M)',iNodeG,'is in slave',matPerLinkNodes(iPos,:)
         iNodeGnew = matPerLinkNodes_i8(iPos,1)
         iPosNew   = binarySearch_int_i8(matPerLinkNodes_i8(:,2),iNodeGnew)
         if(iPosNew.ne.0) then
            !if(mpi_rank.eq.0) write(*,*) 'iNGnew(M)',iNodeGnew,'is in slave too!',matPerLinkNodes(iPosNew,:)
            call find_masterNode_recursive(iNodeGnew,numPerLinkedNodesSrl,matPerLinkNodes_i8,iNodeGtoReplace)
         else
            iNodeGtoReplace = iNodeGnew
         end if
      end if

   end subroutine find_masterNode_recursive


   subroutine read_elems_nodes_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,gmsh_h5_fileId,isPeriodic,isMapFaces,numElemsSrl,numNodesSrl_i8,numPerFacesSrl,numPerLinkedNodesSrl,numMapFacesSrl,&
                     numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,listNodesMpiRank_i8,connecMpiRank_i8,coordNodesMpiRank,&
                     numLinkedPerElemsSrl,linkedPerElemsSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8,&
                     linkedPerToMapElemsAndFacesSrl,mapFaceDir,mapFaceGapCoord)
      implicit none
      integer(4),intent(in)     :: mporder,mnnode,mnpbou,gmsh2ijk(mnnode)
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      logical,intent(in)        :: isPeriodic,isMapFaces
      integer(4),intent(in)     :: numElemsSrl,numPerFacesSrl,numPerLinkedNodesSrl,numMapFacesSrl
      integer(8),intent(in)     :: numNodesSrl_i8
      integer(4),intent(out)    :: numElemsMpiRank,numNodesMpiRank,numLinkedPerElemsSrl,numPerElemsSrl,numMasSlaNodesSrl
      integer(4),intent(out)    :: mapFaceDir
      real(8),intent(out)       :: mapFaceGapCoord
      integer(8),allocatable,intent(inout) :: listNodesMpiRank_i8(:),connecMpiRank_i8(:,:),masSlaNodesSrl_i8(:,:)
      integer(4),allocatable,intent(inout) :: listElemsMpiRank(:),linkedPerElemsSrl(:,:),listPerElemsSrl(:),linkedPerToMapElemsAndFacesSrl(:,:)
      real(8), allocatable, intent(inout) :: coordNodesMpiRank(:,:)

      integer(4),allocatable :: listElemsPerFacesInRank(:,:),listPerFacesInRank(:),listBoundElemsMpiRank(:)
      integer(8),allocatable :: connecPerFacesInRank_i8(:,:)
      integer(4)             :: numPerFacesMpiRank,numBoundElemsMpiRank,iElemMpiRankStart,iElemMpiRankEnd
      integer(4)             :: auxCnt,iElem,h5err

      real(8),dimension(5) :: start_time,end_time,elapsed_time,elapsed_time_m

      !----------------------------------------------------------------------------------------------

      call do_element_partitioning_serial(numElemsSrl,iElemMpiRankStart,iElemMpiRankEnd,numElemsMpiRank,mpi_rank,mpi_size)

      !----------------------------------------------------------------------------------------------

      ! Connectivity table section
      !------------------------------------------------------------------------------------------
      allocate(connecMpiRank_i8(numElemsMpiRank,mnnode))
      allocate(listElemsMpiRank(numElemsMpiRank))

      !$acc kernels
      connecMpiRank_i8(:,:) = 0
      listElemsMpiRank(:) = 0
      !$acc end kernels

      auxCnt=0
      do iElem = iElemMpiRankStart,iElemMpiRankEnd
         auxCnt=auxCnt+1
         listElemsMpiRank(auxCnt) = iElem
      end do

      call MPI_Barrier(app_comm, mpi_err)
      start_time(1) = MPI_Wtime()
      if(mpi_rank.eq.0)write(*,*) "--| Reading element table and coordinates for Mpi Ranks"
      call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(mnnode,gmsh_h5_fileId,numElemsSrl,numNodesSrl_i8,numElemsMpiRank,listElemsMpiRank,&
                                                               numNodesMpiRank,connecMpiRank_i8,listNodesMpiRank_i8,coordNodesMpiRank)
      end_time(1) = MPI_Wtime()

      !------------------------------------------------------------------------------------------
      ! Periodic faces

      if(isPeriodic) then

         call get_listBoundaryElementsInRank_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,numElemsMpiRank,numNodesMpiRank,&
                              listNodesMpiRank_i8,connecMpiRank_i8,numBoundElemsMpiRank,listBoundElemsMpiRank)


         start_time(2) = MPI_Wtime()
         call read_periodic_faces_and_links_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,gmsh_h5_fileId,numPerFacesSrl,numPerLinkedNodesSrl,&
                  numElemsMpiRank,listElemsMpiRank,connecMpiRank_i8,numBoundElemsMpiRank,listBoundElemsMpiRank,numLinkedPerElemsSrl,linkedPerElemsSrl,&
                  numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8,numPerFacesMpiRank,listElemsPerFacesInRank,&
                  listPerFacesInRank,connecPerFacesInRank_i8)
         end_time(2) = MPI_Wtime()

         if(isMapFaces) then
            start_time(3) = MPI_Wtime()
            
            allocate(linkedPerToMapElemsAndFacesSrl(numMapFacesSrl,4))

            call read_mapped_faces_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,gmsh_h5_fileId,numMapFacesSrl,numElemsMpiRank,listElemsMpiRank,&
                     numNodesMpiRank,numPerFacesMpiRank,connecMpiRank_i8,listNodesMpiRank_i8,coordNodesMpiRank,numBoundElemsMpiRank,listBoundElemsMpiRank,&
                     listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8,linkedPerToMapElemsAndFacesSrl,mapFaceDir,mapFaceGapCoord)
            end_time(3) = MPI_Wtime()
         else
            mapFaceDir = 0
            mapFaceGapCoord = 0.0d0
            allocate(linkedPerToMapElemsAndFacesSrl(0,0))
         end if

         deallocate(listElemsPerFacesInRank,listPerFacesInRank,connecPerFacesInRank_i8)


      else
         numLinkedPerElemsSrl = 0
         numPerElemsSrl = 0
         numMasSlaNodesSrl = 0
         allocate(linkedPerElemsSrl(numLinkedPerElemsSrl,0))
         allocate(listPerElemsSrl(numPerElemsSrl))
         allocate(masSlaNodesSrl_i8(numMasSlaNodesSrl,0))

         allocate(linkedPerToMapElemsAndFacesSrl(0,0))
      end if

      elapsed_time(1) = end_time(1) - start_time(1)
      elapsed_time(2) = end_time(2) - start_time(2)
      elapsed_time(3) = end_time(3) - start_time(3)
      call MPI_Reduce(elapsed_time(1),elapsed_time_m(1),3,mpi_datatype_real8,MPI_MAX,0,app_comm,h5err)

      if(mpi_rank.eq.0) then
         write(*,*) ' 1.Time ReadElemNode',elapsed_time_m(1)
         write(*,*) ' 2.Time ReadPeriodic',elapsed_time_m(2)
         write(*,*) ' 3.Time ReadMapFaces',elapsed_time_m(3)
      end if

   end subroutine read_elems_nodes_gmsh_h5_file_in_parallel


   subroutine do_element_partitioning_gempa_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,gmsh_h5_fileId,isPeriodic,isMapFaces,numElemsSrl,numNodesSrl_i8,numBoundFacesSrl,&
                  numMshRanks2Part,numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank,mapMshRankToMpiRank,numElemsMpiRank,listElemsMpiRank,&
                  numNodesMpiRank,listNodesMpiRank_i8,connecMpiRank_i8,coordNodesMpiRank,&
                  numLinkedPerElemsSrl,linkedPerElemsSrl,numPerFacesSrl,numPerElemsSrl,listPerElemsSrl,numMasSlaNodesSrl,masSlaNodesSrl_i8,&
                  numMapFacesSrl,linkedPerToMapElemsAndFacesSrl,&
                  numElemsMshRank,numBoundElemsMshRank,mshRankElemStart,mshRankElemEnd,numNodesMshRank,numBoundFacesMshRank,numPerNodesMshRank,numPerMapLinkedNodesMshRank,maxBoundCode,&
                  elemGid_jv,listNodesMshRank_i8_jv,connecMshRank_i8_jm,listBoundElemsMshRank_jv,listElemsBoundFacesMshRank_jm,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,&
                  connecBoundFacesMshRank_i8_jm,coordMshRank_jm,masSlaRankPar_i8_jm,perMapLinkedNodesRankPar_i8_jm)
      implicit none
      integer(4),intent(in) :: mporder,mnnode,mnpbou,gmsh2ijk(mnnode)
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      logical,intent(in) :: isPeriodic,isMapFaces
      integer(4),intent(in) :: numElemsSrl,numBoundFacesSrl
      integer(8),intent(in) :: numNodesSrl_i8
      integer(4),intent(in) :: numElemsMpiRank,numNodesMpiRank,numMshRanks2Part,numLinkedPerElemsSrl,numPerFacesSrl,numPerElemsSrl
      integer(4),intent(inout) :: numMasSlaNodesSrl,numMapFacesSrl
      integer(4),intent(in) :: numMshRanksInMpiRank,maxNumMshRanks,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part)
      integer(4),intent(inout),allocatable :: listElemsMpiRank(:)
      integer(4),intent(in) :: linkedPerElemsSrl(numLinkedPerElemsSrl,2),listPerElemsSrl(numPerElemsSrl)
      integer(8),intent(inout),allocatable :: listNodesMpiRank_i8(:),connecMpiRank_i8(:,:),masSlaNodesSrl_i8(:,:)
      integer(4),intent(inout),allocatable :: linkedPerToMapElemsAndFacesSrl(:,:)
      real(8),intent(inout),allocatable :: coordNodesMpiRank(:,:)
      integer(4),intent(out) :: maxBoundCode
      integer(4),intent(out),allocatable :: numElemsMshRank(:),numBoundElemsMshRank(:),mshRankElemStart(:),mshRankElemEnd(:),numNodesMshRank(:),numBoundFacesMshRank(:),numPerNodesMshRank(:),numPerMapLinkedNodesMshRank(:)
      type(jagged_vector_int4),intent(out) :: elemGid_jv,listBoundFacesMshRank_jv,boundFacesCodesMshRank_jv,listBoundElemsMshRank_jv
      type(jagged_vector_int8),intent(out) :: listNodesMshRank_i8_jv
      type(jagged_matrix_int8),intent(out) :: connecMshRank_i8_jm,connecBoundFacesMshRank_i8_jm,masSlaRankPar_i8_jm,perMapLinkedNodesRankPar_i8_jm
      type(jagged_matrix_int4),intent(out) :: listElemsBoundFacesMshRank_jm
      type(jagged_matrix_real8),intent(out) :: coordMshRank_jm

      integer(4),allocatable :: dummyListElems(:),dummyListBoundFaces(:),dummyBoundFacesCodes(:),dummyListElemsBounds(:,:),dummyListBoundElems(:)
      integer(8),allocatable :: dummyListNodes(:),dummyConnec(:,:),dummyConnecBoundFaces(:,:)
      real(8),allocatable :: dummyCoordNodes(:,:)
      integer(4) :: dummyNumNodes,dummyNumBoundFaces,dummyNumBoundElems

      real(8), allocatable, dimension(:) :: x,y,z
      integer(4),allocatable :: listElems2Part(:),weightElems2Part(:),elemPart(:,:)
      integer(8),allocatable :: connecPerFacesSrl_i8(:,:),connecMapFacesSrl_i8(:,:)

      integer(4),parameter :: numHexaVert=8
      integer(4) :: gmshHexVertInd(numHexaVert)

      integer(4) :: ii,jj,kk,m,iVert,iElem,iElem2Part,iNode,iElemG,iBound,iPos
      integer(8) :: iNodeG
      integer(4) :: mshRank,mpiRank,iMshRank,iMpiRank
      integer(4) :: numElems2Part,numValues2get,aux_numElemsMshRank,aux_sum_NBP,aux_sum_NB_mpiRank
      integer(4) :: maxBoundCodeMpiRank,maxBoundCodeMshRank
      character(128) :: dsetname2read

      real(8) :: x_a,y_a,z_a
      integer(4),dimension(0:numMshRanks2Part-1) :: numElems2mshRankInMpiRank,numElemsInMshRank
      integer(4),dimension(0:numMshRanks2Part-1,0:mpi_size-1) :: matNumElems2mshRankInMpiRank

      integer(4) :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size,target_displacement

      if(mpi_rank.eq.0) write(*,*) '--| Doing element partitioning GEMPA in parallel...'

      !1. obtain the list element 2 part using gempa
      if(isPeriodic) then
         call get_listElems2Par_periodic_in_parallel(numElemsMpiRank,numLinkedPerElemsSrl,numMapFacesSrl,numPerElemsSrl,listElemsMpiRank,linkedPerElemsSrl,linkedPerToMapElemsAndFacesSrl,listElems2Part,weightElems2Part,numElems2Part)
      else
         call get_listElems2Par_in_parallel(numElemsMpiRank,listElemsMpiRank,listElems2Part,weightElems2Part,numElems2Part)
      end if

      allocate(elemPart(numElems2Part,4))
      allocate(x(numElems2Part))
      allocate(y(numElems2Part))
      allocate(z(numElems2Part))

      call get_gmshHexHOVertIndex(mporder,gmshHexVertInd)

      do iElem2Part=1,numElems2Part
         elemPart(iElem2Part,1) = listElems2Part(iElem2Part)
         elemPart(iElem2Part,3) = weightElems2Part(iElem2Part)
         iElemG = listElems2Part(iElem2Part)
         iElem = binarySearch_int_i4(listElemsMpiRank,iElemG)
         elemPart(iElem2Part,4) = iElem

         x_a=0.0d0
         y_a=0.0d0
         z_a=0.0d0
         do iVert=1,numHexaVert
            m = gmshHexVertInd(iVert)
            iNodeG = connecMpiRank_i8(iElem,m)

            iPos = binarySearch_int_i8(listNodesMpiRank_i8,iNodeG)

            x_a = x_a + coordNodesMpiRank(iPos,1)
            y_a = y_a + coordNodesMpiRank(iPos,2)
            z_a = z_a + coordNodesMpiRank(iPos,3)
         end do

         x(iElem2Part) = x_a/8.0d0
         y(iElem2Part) = y_a/8.0d0
         z(iElem2Part) = z_a/8.0d0
         !write(*,*) '[',mpi_rank,']iElemG',elemPart(iElem2Part,1),'w',elemPart(iElem2Part,3),'x',x(iElem2Part),'y',y(iElem2Part),'z',z(iElem2Part)
         !write(*,*) '[',mpi_rank,']iElemG',elemPart(iElem2Part,1),'z',z(iElem2Part)
      end do

      deallocate(listElems2Part)
      deallocate(listElemsMpiRank)

#ifndef GEMPAINTERFACE
      !---- CALLING GEMPA in PARALLEL-----------------------
      call gempa_do_partition(numElems2Part,numMshRanks2Part,x,y,z,elemPart(:,3),elemPart(:,2))
      !---------------------------------------------------------------------
#else
      write(*,*) 'FATAL ERROR! Trying to call do_element_partitioning_gempa() for program compiled with flag GEMPAINTERFACE!'
      call MPI_Abort(app_comm,-1,mpi_err)
#endif

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'gempaPartitionInfo_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,weight,iElem'
      do ii=1,numElems2Part
         write(1,fmt_csv_deb) x(ii),y(ii),z(ii),elemPart(ii,1),elemPart(ii,2),elemPart(ii,3),elemPart(ii,4)
      end do
      close(1)
#endif

      deallocate(x)
      deallocate(y)
      deallocate(z)

      if(isPeriodic) then
         ! if parallel, now we have to put the rank to the 'slave' elements who were not included in the partitioning process
         ! since they were linked to a 'master' element
         call unfold_linked_elems(numElems2Part,numLinkedPerElemsSrl,numMapFacesSrl,numPerElemsSrl,linkedPerElemsSrl,linkedPerToMapElemsAndFacesSrl,listPerElemsSrl,elemPart)
      end if

      ! @now we have to info all the ranks about which are the elements that they will 'store'
      ! since we have done the partitioning in paralel...

      numElems2mshRankInMpiRank(:)=0
      !!!$acc parallel loop
      do iElem=1,numElems2Part
         mshRank = elemPart(iElem,2)-1 !elemPart came with rank=1:mpi_size
         !!!$acc atomic update
         numElems2mshRankInMpiRank(mshRank) = numElems2mshRankInMpiRank(mshRank) + 1
         !!!$acc end atomic
      end do
      !!!$acc end parallel loop

      !write(*,*) 'numElems2mshRankInMpiRank[',mpi_rank,'] ',numElems2mshRankInMpiRank(:)

#if _CHECK_
      if(.not.(isPeriodic)) then
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'elemPartition_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,iElem'
      do ii=1,numElems2Part
         iElemG   = elemPart(ii,1)
         mshRank  = elemPart(ii,2)-1
         !iElem    = elemPart(ii,4)

         x_a=0.0d0
         y_a=0.0d0
         z_a=0.0d0
         do jj=1,numHexaVert
            m = gmshHexVertInd(jj)
            iNodeG = connecMpiRank_i8(iElem,m)

            iPos = binarySearch_int_i8(listNodesMpiRank_i8,iNodeG)

            x_a = x_a + coordNodesMpiRank(iPos,1)
            y_a = y_a + coordNodesMpiRank(iPos,2)
            z_a = z_a + coordNodesMpiRank(iPos,3)
         end do
         x_a = x_a/8.d0
         y_a = y_a/8.d0
         z_a = z_a/8.d0

         write(1,fmt_csv_deb) x_a,y_a,z_a,iElemG,mshRank,iElem
      end do
      close(1)
      end if
#endif

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part
      call MPI_Win_create(numElems2mshRankInMpiRank,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMpiRank=0,mpi_size-1
         call MPI_Get(matNumElems2mshRankInMpiRank(:,iMpiRank),numMshRanks2Part,mpi_datatype_int,iMpiRank,target_displacement,&
                                               numMshRanks2Part,mpi_datatype_int,window_id,mpi_err)
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------

      call quicksort_matrix_int4(elemPart,2)

      ii=1
      do mshRank=0,numMshRanks2Part-1
         jj=numElems2mshRankInMpiRank(mshRank)
         if(jj.ne.0) then
            kk=ii+jj-1
            call quicksort_matrix_int4(elemPart,1,ii,kk)
            ii=kk+1
         endif
      end do

#if _CHECK_
      if(.not.(isPeriodic)) then
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      file_name_deb = 'elemPartitionOrdered_mpiRank'// trim(aux_string_mpirank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iElemG,rank,iElem'
      do ii=1,numElems2Part
         iElemG  = elemPart(ii,1)
         mshRank = elemPart(ii,2)-1
         iElem   = elemPart(ii,4)

         x_a=0.0d0
         y_a=0.0d0
         z_a=0.0d0
         do jj=1,numHexaVert
            m = gmshHexVertInd(jj)
            iNodeG = connecMpiRank_i8(iElem,m)

            iPos = binarySearch_int_i8(listNodesMpiRank_i8,iNodeG)

            x_a = x_a + coordNodesMpiRank(iPos,1)
            y_a = y_a + coordNodesMpiRank(iPos,2)
            z_a = z_a + coordNodesMpiRank(iPos,3)
         end do
         x_a = x_a/8.d0
         y_a = y_a/8.d0
         z_a = z_a/8.d0

         write(1,fmt_csv_deb) x_a,y_a,z_a,iElemG,mshRank,iElem
      end do
      close(1)
      end if
#endif

      !--------------------------------------------------------------------------------------
      ! Deallocating arrays
      deallocate(listNodesMpiRank_i8)
      deallocate(coordNodesMpiRank)
      deallocate(connecMpiRank_i8)

      !--------------------------------------------------------------------------------------

      numElemsInMshRank(:)=0

      do mpiRank=0,mpi_size-1
         do mshRank=0,numMshRanks2Part-1
            numElemsInMshRank(mshRank) = numElemsInMshRank(mshRank) + matNumElems2mshRankInMpiRank(mshRank,mpiRank)
         end do
      end do

      allocate(elemGid_jv%vector(numMshRanksInMpiRank))
      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         allocate(elemGid_jv%vector(iMshRank)%elems(numElemsInMshRank(mshRank)))
         !$acc kernels
         elemGid_jv%vector(iMshRank)%elems(:) = -1
         !$acc end kernels
      end do

      ! Create the window
      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numElems2Part

      call MPI_Win_create(elemPart(:,1),win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         iPos=1
         mshRank=mshRanksInMpiRank(iMshRank)

         do mpiRank=0,mpi_size-1
            numValues2get=matNumElems2mshRankInMpiRank(mshRank,mpiRank)

            if(numValues2get.ne.0) then
               target_displacement=0
               do jj=0,mshRank-1
                  target_displacement=target_displacement+matNumElems2mshRankInMpiRank(jj,mpiRank)
               end do

               call MPI_Get(elemGid_jv%vector(iMshRank)%elems(iPos),numValues2get,mpi_datatype_int,mpiRank,target_displacement,&
                          numValues2get,mpi_datatype_int,window_id,mpi_err)

               iPos = iPos + numValues2get
            end if

         end do
      end do

      !!! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         write(aux_string_mshrank_deb,'(I0)') mshRank
         file_name_deb = 'elemGid_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
         open(1,file=file_name_deb)
         do iElem=1,numElemsInMshRank(iMshRank)
            !write(1,fmt_csv_deb) iElem,elemGid(iElem,iMshRank)
            write(1,fmt_csv_deb) iElem,elemGid_jv%vector(iMshRank)%elems(iElem)
         end do
         close(1)
      end do
#endif
!---------------------------------------------------------------------------------
      deallocate(elemPart)

      allocate(numElemsMshRank(numMshRanksInMpiRank))
      allocate(mshRankElemStart(numMshRanksInMpiRank))
      allocate(mshRankElemEnd(numMshRanksInMpiRank))
      allocate(numBoundElemsMshRank(numMshRanksInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank= mshRanksInMpiRank(iMshRank)
         numElemsMshRank(iMshRank) = numElemsInMshRank(mshRank)
      end do
      !write(*,*) 'numElemsMshRank[',mpi_rank,']',numElemsMshRank(:)

      mshRankElemStart(:)=1
      do iMshRank=1,numMshRanksInMpiRank
         mshRank= mshRanksInMpiRank(iMshRank)
         do jj=0,mshRank-1
            mshRankElemStart(iMshRank)=mshRankElemStart(iMshRank)+numElemsInMshRank(jj)
         end do
         mshRankElemEnd(iMshRank) = mshRankElemStart(iMshRank) + numElemsInMshRank(jj)-1
      end do
      !write(*,*) 'numElemsMshRank[',mpi_rank,']',numElemsMshRank(:),'mshRanksInMpiRank',mshRanksInMpiRank(:),'eStart',mshRankElemStart(:),'eEnd ',mshRankElemEnd(:)

      numBoundElemsMshRank(:) = 0
!------------------------------------------------------------------------------

      allocate(numNodesMshRank(numMshRanksInMpiRank))
      allocate(numBoundFacesMshRank(numMshRanksInMpiRank))
      allocate(listNodesMshRank_i8_jv%vector(numMshRanksInMpiRank))
      allocate(connecMshRank_i8_jm%matrix(numMshRanksInMpiRank))!,boundMshRank_jm,bou_codesMshRank_jm
      allocate(coordMshRank_jm%matrix(numMshRanksInMpiRank))

      allocate(numPerNodesMshRank(numMshRanksInMpiRank))
      allocate(numPerMapLinkedNodesMshRank(numMshRanksInMpiRank))
      allocate(masSlaRankPar_i8_jm%matrix(numMshRanksInMpiRank))
      allocate(perMapLinkedNodesRankPar_i8_jm%matrix(numMshRanksInMpiRank))

      !if(numBoundFacesSrl.ne.0) then
      allocate(listBoundFacesMshRank_jv%vector(numMshRanksInMpiRank))
      allocate(boundFacesCodesMshRank_jv%vector(numMshRanksInMpiRank))
      allocate(listElemsBoundFacesMshRank_jm%matrix(numMshRanksInMpiRank))
      allocate(connecBoundFacesMshRank_i8_jm%matrix(numMshRanksInMpiRank))
      allocate(listBoundElemsMshRank_jv%vector(numMshRanksInMpiRank))
      !end if

      numNodesMshRank(:)  = 0
      numBoundFacesMshRank(:) = 0
      numPerNodesMshRank(:) = 0
      numPerMapLinkedNodesMshRank(:) = 0

      maxBoundCodeMshRank = 0
      maxBoundCodeMpiRank = 0
      maxBoundCode=0

      if(mpi_rank.eq.0) then
         write(*,*) '--| Reading element table and coordinates for Msh Ranks | numMshRanks(0)',numMshRanksInMpiRank,'max',maxNumMshRanks
         write(*,*) '  |-> numMshRanks(rank=0)',numMshRanksInMpiRank,'max',maxNumMshRanks
      end if

      do iMshRank=1,numMshRanksInMpiRank
         !write(*,*) 'REAL-rank[',mpi_rank,']doing',iMshRank,'max',maxNumMshRanks
         aux_numElemsMshRank=numElemsMshRank(iMshRank)
         mshRank= mshRanksInMpiRank(iMshRank)
         !write(*,*) '#Reading stuff for mshRank',mshRank,'in mpi_rank',mpi_rank,'(iMshRank',iMshRank,')','aux_numElemsMshRank',aux_numElemsMshRank

         allocate(connecMshRank_i8_jm%matrix(iMshRank)%elems(aux_numElemsMshRank,mnnode))

         call quicksort_array_int4(elemGid_jv%vector(iMshRank)%elems)

         call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(mnnode,gmsh_h5_fileId,numElemsSrl,numNodesSrl_i8,aux_numElemsMshRank,elemGid_jv%vector(iMshRank)%elems,numNodesMshRank(iMshRank),&
                                          connecMshRank_i8_jm%matrix(iMshRank)%elems,listNodesMshRank_i8_jv%vector(iMshRank)%elems,coordMshRank_jm%matrix(iMshRank)%elems)

         if(numBoundFacesSrl.ne.0) then

            call get_listBoundaryElementsInRank_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,aux_numElemsMshRank,numNodesMshRank(iMshRank),&
                              listNodesMshRank_i8_jv%vector(iMshRank)%elems,connecMshRank_i8_jm%matrix(iMshRank)%elems,&
                              numBoundElemsMshRank(iMshRank),listBoundElemsMshRank_jv%vector(iMshRank)%elems)

            dsetname2read = '/boundFaces'
            call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,.false.,gmsh_h5_fileId,dsetname2read,&
                     numBoundFacesSrl,aux_numElemsMshRank,elemGid_jv%vector(iMshRank)%elems,&
                     numBoundElemsMshRank(iMshRank),listBoundElemsMshRank_jv%vector(iMshRank)%elems,connecMshRank_i8_jm%matrix(iMshRank)%elems,&
                     numBoundFacesMshRank(iMshRank),listElemsBoundFacesMshRank_jm%matrix(iMshRank)%elems,&
                     listBoundFacesMshRank_jv%vector(iMshRank)%elems,connecBoundFacesMshRank_i8_jm%matrix(iMshRank)%elems)


            dsetname2read = '/boundFacesId'
            call read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,numBoundFacesMshRank(iMshRank),&
                     listBoundFacesMshRank_jv%vector(iMshRank)%elems,boundFacesCodesMshRank_jv%vector(iMshRank)%elems,maxBoundCodeMshRank)

            maxBoundCodeMpiRank = max(maxBoundCodeMpiRank,maxBoundCodeMshRank)
            !if(mpi_rank.eq.0) write(*,*) ' maxBoundCode',maxBoundCodeMpiRank
         else
            allocate(listBoundFacesMshRank_jv%vector(iMshRank)%elems(numBoundFacesMshRank(iMshRank)))
            allocate(boundFacesCodesMshRank_jv%vector(iMshRank)%elems(numBoundFacesMshRank(iMshRank)))
            allocate(connecBoundFacesMshRank_i8_jm%matrix(iMshRank)%elems(numBoundFacesMshRank(iMshRank),mnpbou))
            allocate(listElemsBoundFacesMshRank_jm%matrix(iMshRank)%elems(numElemsMshRank(iMshRank),maxBoundsPerElem))
            listElemsBoundFacesMshRank_jm%matrix(iMshRank)%elems(:,:) = 0
         end if

      end do

      aux_numElemsMshRank=1
      allocate(dummyListElems(1))
      dummyListElems(1)=1
      allocate(dummyConnec(1,mnnode))
      dummyNumBoundElems=1
      allocate(dummyListBoundElems(dummyNumBoundElems))

      do iMshRank=(numMshRanksInMpiRank+1),maxNumMshRanks
         !write(*,*) 'FAKE-rank[',mpi_rank,']doing',iMshRank,'max',maxNumMshRanks

         call read_elem_connec_and_nodes_coords_from_gmsh_h5_file(mnnode,gmsh_h5_fileId,numElemsSrl,numNodesSrl_i8,aux_numElemsMshRank,dummyListElems,dummyNumNodes,&
                                          dummyConnec,dummyListNodes,dummyCoordNodes)

         if(numBoundFacesSrl.ne.0) then
            dummyListElems(1) = 0 !to avoid possible fake matches!
            dummyConnec(1,:) = 0  !to avoid possible fake matches!
            dummyListBoundElems(1) = 0

            dsetname2read = '/boundFaces'
            call read_boundaries_for_elemsInRank_from_gmsh_h5_file_in_parallel(mporder,mnnode,mnpbou,.false.,gmsh_h5_fileId,dsetname2read,&
                     numBoundFacesSrl,aux_numElemsMshRank,dummyListElems,dummyNumBoundElems,dummyListBoundElems,dummyConnec,&
                     dummyNumBoundFaces,dummyListElemsBounds,dummyListBoundFaces,dummyConnecBoundFaces)

            dsetname2read = '/boundFacesId'
            call read_boundaries_boundCodes_from_gmsh_h5_file_in_parallel(gmsh_h5_fileId,dsetname2read,dummyNumBoundFaces,&
                     dummyListBoundFaces,dummyBoundFacesCodes,maxBoundCodeMshRank)

            deallocate(dummyListElemsBounds)
            deallocate(dummyListBoundFaces)
            deallocate(dummyConnecBoundFaces)
            deallocate(dummyBoundFacesCodes)
         end if

         deallocate(dummyListNodes)
         deallocate(dummyCoordNodes)

      end do
      deallocate(dummyListElems)
      deallocate(dummyConnec)
      deallocate(dummyListBoundElems)

      call MPI_Allreduce(maxBoundCodeMpiRank,maxBoundCode,1,MPI_INT,MPI_MAX,app_comm,mpi_err)

      !fer el check aqui de que s'han llegit correctament totes les boundaries
      if(numBoundFacesSrl.ne.0) then
         aux_sum_NB_mpiRank=0
         do iMshRank=1,numMshRanksInMpiRank
            aux_sum_NB_mpiRank = aux_sum_NB_mpiRank + numBoundFacesMshRank(iMshRank)
         end do

         call MPI_Allreduce(aux_sum_NB_mpiRank,aux_sum_NBP,1,MPI_INT,MPI_SUM,app_comm,mpi_err)

         if(aux_sum_NBP.ne.numBoundFacesSrl) then
            write(*,*) '[',mpi_rank,']Error in do_element_partitioning_gempa_in_parallel(..) after read_boundaries(..): aux_sum_NBP',aux_sum_NBP,' not equal to  numBoundFacesSrl',numBoundFacesSrl
            call MPI_Abort(app_comm, -1, mpi_err)
         end if
      end if

      !create here masSlaNodesMshRank

      if(isPeriodic) then
         do iMshRank=1,numMshRanksInMpiRank
            mshRank= mshRanksInMpiRank(iMshRank)
            call generate_masSlaRankPar(mshRank,numNodesMshRank(iMshRank),listNodesMshRank_i8_jv%vector(iMshRank)%elems,numMasSlaNodesSrl,masSlaNodesSrl_i8,&
                        numPerNodesMshRank(iMshRank),masSlaRankPar_i8_jm%matrix(iMshRank)%elems)
         end do

         !since now i have the masSlaRankPar_i8_jm with global nodes ids i can deallocate the serial one
         deallocate(masSlaNodesSrl_i8)
         numMasSlaNodesSrl=0

         if(isMapFaces) then
            allocate(connecPerFacesSrl_i8(numPerFacesSrl,mnpbou))
            allocate(connecMapFacesSrl_i8(numMapFacesSrl,mnpbou))

            call read_mapPerBounds_from_gmsh_h5_file_in_parallel(mporder,mnpbou,gmsh_h5_fileId,&
                     numPerFacesSrl,numMapFacesSrl,connecPerFacesSrl_i8,connecMapFacesSrl_i8)
            
            do iMshRank=1,numMshRanksInMpiRank
               mshRank= mshRanksInMpiRank(iMshRank)
               call generate_mappedNodeLinksRankPar(mporder,mnnode,mnpbou,mshRank,numPerFacesSrl,numMapFacesSrl,&
                     numElemsMshRank(iMshRank),elemGid_jv%vector(iMshRank)%elems,numNodesMshRank(iMshRank),listNodesMshRank_i8_jv%vector(iMshRank)%elems,&
                     coordMshRank_jm%matrix(iMshRank)%elems,connecPerFacesSrl_i8,connecMapFacesSrl_i8,linkedPerToMapElemsAndFacesSrl,&
                     numPerMapLinkedNodesMshRank(iMshRank),perMapLinkedNodesRankPar_i8_jm%matrix(iMshRank)%elems)
            end do

            deallocate(connecPerFacesSrl_i8,connecMapFacesSrl_i8)
            deallocate(linkedPerToMapElemsAndFacesSrl)
            numMapFacesSrl=0
         else
            do iMshRank=1,numMshRanksInMpiRank
               numPerMapLinkedNodesMshRank(iMshRank)=0
               allocate(perMapLinkedNodesRankPar_i8_jm%matrix(iMshRank)%elems(0,0))
            end do
         end if

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         write(aux_string_mshrank_deb,'(I0)') mshRank
         file_name_deb = 'masSla_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
         open(1,file=file_name_deb)
         write(1,*) 'X,Y,Z,iNodeG,MS'
         do ii=1,numPerNodesMshRank(iMshRank)
            iNodeG=masSlaRankPar_i8_jm%matrix(iMshRank)%elems(ii,1)
            iPos = binarySearch_int_i8(listNodesMshRank_i8_jv%vector(iMshRank)%elems,iNodeG)
            x_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,1)
            y_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,2)
            z_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,3)
            write(1,fmt_csv_deb) x_a,y_a,z_a,iNodeG,'1'
         end do
         do ii=1,numPerNodesMshRank(iMshRank)
            iNodeG=masSlaRankPar_i8_jm%matrix(iMshRank)%elems(ii,2)
            iPos = binarySearch_int_i8(listNodesMshRank_i8_jv%vector(iMshRank)%elems,iNodeG)
            x_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,1)
            y_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,2)
            z_a = coordMshRank_jm%matrix(iMshRank)%elems(iPos,3)
            write(1,fmt_csv_deb) x_a,y_a,z_a,iNodeG,'2'
         end do
         close(1)
      end do
#endif
      else
         do iMshRank=1,numMshRanksInMpiRank
            numPerNodesMshRank(iMshRank)=0
            allocate(masSlaRankPar_i8_jm%matrix(iMshRank)%elems(0,0))
            numPerMapLinkedNodesMshRank(iMshRank)=0
            allocate(perMapLinkedNodesRankPar_i8_jm%matrix(iMshRank)%elems(0,0))
         end do
      end if

   end subroutine do_element_partitioning_gempa_in_parallel
!------------------------------------------------------------------------------------------------

   subroutine generate_masSlaRankPar(mshRank,numNodesInRank,listNodesInRank_i8,numMasSlaNodesSrl,masSlaNodesSrl_i8,numMasSlaNodesInRank,masSlaNodesInRank_i8)
      implicit none
      integer(4),intent(in) :: mshRank,numNodesInRank,numMasSlaNodesSrl
      integer(8),intent(in) :: listNodesInRank_i8(numNodesInRank),masSlaNodesSrl_i8(numMasSlaNodesSrl,2)
      integer(4),intent(out) :: numMasSlaNodesInRank
      integer(8),allocatable,intent(out) :: masSlaNodesInRank_i8(:,:)
      integer(4) :: iLink,iPos,iNode,iAux
      integer(8) :: iNodeGM,iNodeGS
      integer(4) :: numPerNodesInRank

      numMasSlaNodesInRank = 0
      do iLink = 1,numMasSlaNodesSrl
         iNodeGS = masSlaNodesSrl_i8(iLink,2)

         iPos = binarySearch_int_i8(listNodesInRank_i8,iNodeGS)
         if(iPos.ne.0) numMasSlaNodesInRank = numMasSlaNodesInRank + 1! nodeIsPerInRank(iPosS) = .true.
      end do

      allocate(masSlaNodesInRank_i8(numMasSlaNodesInRank,2))
      iAux = 0
      do iLink = 1,numMasSlaNodesSrl
         iNodeGM = masSlaNodesSrl_i8(iLink,1)
         iNodeGS = masSlaNodesSrl_i8(iLink,2)

         !iPosM = binarySearch_int_i(listNodesInRank,iNodeGM)
         !if(iPosM.ne.0) nodeIsPerInRank(iPosM) = .true.

         iPos = binarySearch_int_i8(listNodesInRank_i8,iNodeGS)
         if(iPos.ne.0) then
            iAux=iAux+1
            masSlaNodesInRank_i8(iAux,1) = iNodeGM
            masSlaNodesInRank_i8(iAux,2) = iNodeGS
            !write(*,*) 'mshRank',mshRank,'MS',masSlaNodesInRank(iAux,:)
         end if
      end do

      !write(*,*) 'numMasSlaRank(',mshRank,')',numMasSlaNodesInRank,'numMasSlaSrl',numMasSlaNodesSrl

   end subroutine generate_masSlaRankPar

   subroutine read_mapPerBounds_from_gmsh_h5_file_in_parallel(mporder,mnpbou,gmsh_h5_fileId,numPerFacesSrl,numMapFacesSrl,&
         connecPerFacesSrl_i8,connecMapFacesSrl_i8)
      implicit none
      integer(4), intent(in) :: mporder,mnpbou,numPerFacesSrl,numMapFacesSrl
      integer(hid_t),intent(in) :: gmsh_h5_fileId
      integer(8),intent(out) :: connecPerFacesSrl_i8(mnpbou,numPerFacesSrl),connecMapFacesSrl_i8(mnpbou,numMapFacesSrl)

      character(256) :: dsetname
      integer(4),parameter :: ms_rank=2,ds_rank=2
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id
      integer(hsize_t),dimension(ds_rank) :: fs_dims,fs_maxdims
      integer(4) :: h5err
      integer(hsize_t), dimension(ms_rank) :: ms_dims
      integer(hssize_t), dimension(ms_rank) :: ms_offset
      integer(hid_t) :: dtype

      real(8), dimension(2) :: start_time,end_time,elapsed_time,elapsed_time_m

      start_time(1) = MPI_Wtime()

      dtype = h5_datatype_int8
      ms_dims(1) = int(mnpbou,hsize_t)
      ms_dims(2) = 0
      ms_offset(1) = 0
      ms_offset(2) = 0

      dsetname = '/periodicFaces'
      connecPerFacesSrl_i8(:,:) = 0
      ms_dims(2) = int(numPerFacesSrl,hsize_t)

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      call h5dread_f(dset_id,dtype,connecPerFacesSrl_i8(1:mnpbou,1:numPerFacesSrl),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(1) = MPI_Wtime()
      elapsed_time(1) = end_time(1) - start_time(1)
      start_time(2) = MPI_Wtime()

      dsetname = '/mappedFaces'
      connecMapFacesSrl_i8(:,:) = 0
      ms_dims(2) = int(numMapFacesSrl,hsize_t)

      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)

      call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      call h5dread_f(dset_id,dtype,connecMapFacesSrl_i8(1:mnpbou,1:numMapFacesSrl),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      end_time(2) = MPI_Wtime()
      elapsed_time(2) = end_time(2) - start_time(2)

   end subroutine read_mapPerBounds_from_gmsh_h5_file_in_parallel

   subroutine generate_mappedNodeLinksRankPar(mporder,mnnode,mnpbou,mshRank,numPerFacesSrl,numMapFacesSrl,numElemsInRank,elemGidInRank,numNodesInRank,&
         listNodesInRank_i8,coordNodesInRank,connecPerFacesSrl_i8,connecMapFacesSrl_i8,linkedPerToMapElemsAndFacesSrl,numLinkedMapNodesInRank,perMapLinkedNodes_i8)
      implicit none
      integer(4),intent(in) :: mporder,mnnode,mnpbou,mshRank,numPerFacesSrl,numMapFacesSrl
      integer(4),intent(in) :: numElemsInRank,numNodesInRank
      integer(4),intent(in) :: elemGidInRank(numElemsInRank)
      integer(8),intent(in) :: listNodesInRank_i8(numNodesInRank)
      integer(8),intent(in) :: connecPerFacesSrl_i8(mnpbou,numPerFacesSrl),connecMapFacesSrl_i8(mnpbou,numMapFacesSrl)
      real(8),intent(in)   :: coordNodesInRank(numNodesInRank,ndime)
      integer(4),intent(in) :: linkedPerToMapElemsAndFacesSrl(numMapFacesSrl,4)
      integer(4),intent(out) :: numLinkedMapNodesInRank
      integer(8),allocatable,intent(inout) :: perMapLinkedNodes_i8(:,:)

      integer(4) :: numLinkedMapFacesInRank,numRawLinkedMapNodesInRank
      integer(4) :: iElem,iElemPerG,iElemMapG,iFacePerG,iFaceMapG,elemPos
      integer(4) :: idime,ip,jp,iPosNodePer,iPosNodeMap,rawCnt
      integer(8) :: iNodePerG,iNodeMapG,minDist_iNodeMapG,iNodePerG2
      real(8) :: coordPerNode(ndime),coordMapNode(ndime),minNodeDist,nodeDistAbs,nodeDist(ndime)
      integer(8),allocatable :: rawPerMapLinkedNodes_i8(:,:)
      logical :: keepLooping
      !------------------------------------------------------------------------------------------------------------------

      numLinkedMapFacesInRank = 0
      do iElem = 1,numElemsInRank
         iElemPerG = elemGidInRank(iElem)
         elemPos = binarySearch_int_i4(linkedPerToMapElemsAndFacesSrl(:,1),iElemPerG)
         if(elemPos.ne.0) then
            numLinkedMapFacesInRank = numLinkedMapFacesInRank + 1
         end if
      end do
      !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']numLinkedMapFacesInRank',numLinkedMapFacesInRank

      numRawLinkedMapNodesInRank = numLinkedMapFacesInRank*mnpbou
      allocate(rawPerMapLinkedNodes_i8(numRawLinkedMapNodesInRank,2))

      rawCnt = 0
      do iElem = 1,numElemsInRank
         iElemPerG = elemGidInRank(iElem)
         elemPos = binarySearch_int_i4(linkedPerToMapElemsAndFacesSrl(:,1),iElemPerG)
         if(elemPos.ne.0) then
            iElemMapG = linkedPerToMapElemsAndFacesSrl(elemPos,2)
            iFacePerG = linkedPerToMapElemsAndFacesSrl(elemPos,3)
            iFaceMapG = linkedPerToMapElemsAndFacesSrl(elemPos,4)
            !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']iEpG',iElemPerG,'iEmG',iElemMapG,'iFpG',iFacePerG,'iFmG',iFaceMapG

            do ip=1,mnpbou
               iNodePerG = connecPerFacesSrl_i8(ip,iFacePerG)
               iPosNodePer = binarySearch_int_i8(listNodesInRank_i8,iNodePerG)
               coordPerNode(:) = coordNodesInRank(iPosNodePer,:) 
               !write(*,*) '  ip',ip,'coord',coordPerNode(:)
               
               minNodeDist = 1.0e9
               do jp=1,mnpbou
                  iNodeMapG = connecMapFacesSrl_i8(jp,iFaceMapG)
                  iPosNodeMap = binarySearch_int_i8(listNodesInRank_i8,iNodeMapG)
                  coordMapNode(:) = coordNodesInRank(iPosNodeMap,:) 

                  nodeDistAbs = 0.0d0
                  do idime=1,ndime
                     nodeDist(idime) = coordPerNode(idime) - coordMapNode(idime)
                     nodeDistAbs     = nodeDistAbs + nodeDist(idime)*nodeDist(idime)
                  end do
                  nodeDistAbs = sqrt(nodeDistAbs)

                  if(nodeDistAbs.lt.minNodeDist) then
                     minNodeDist = nodeDistAbs
                     minDist_iNodeMapG = iNodeMapG 
                  end if
               end do
               rawCnt = rawCnt + 1
               rawPerMapLinkedNodes_i8(rawCnt,1) = iNodePerG
               rawPerMapLinkedNodes_i8(rawCnt,2) = minDist_iNodeMapG
               !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']iEpG',iElemPerG,'iNPG',iNodePerG,'iNMG',minDist_iNodeMapG,'minDist',minNodeDist
            end do

         end if
      end do

      call quicksort_matrix_int8(rawPerMapLinkedNodes_i8,1)
      !if(mpi_rank.eq.0) write(*,*) '[',mpi_rank,']mshRank[',mshRank,']rawPerMapLinkedNodes_i8(:,:)',rawPerMapLinkedNodes_i8(:,1)
      !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']FLAG numRawLinked',numRawLinkedMapNodesInRank

      ip=1
      numLinkedMapNodesInRank=numRawLinkedMapNodesInRank
      do while (ip<numRawLinkedMapNodesInRank)
         iNodePerG = rawPerMapLinkedNodes_i8(ip,1)
         jp=ip+1
         keepLooping = .true.
         do while(keepLooping.and.(jp.le.numRawLinkedMapNodesInRank))
            iNodePerG2 = rawPerMapLinkedNodes_i8(jp,1)
            ip = jp
            !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']iNPG',iNodePerG,'iNPG2',iNodePerG2,'ip',ip,'jp',jp
            if(iNodePerG.eq.iNodePerG2) then
               !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']rep-iNPG',iNodePerG,'ip',ip,'jp',jp
               numLinkedMapNodesInRank = numLinkedMapNodesInRank - 1
               rawPerMapLinkedNodes_i8(jp,:) = 0
               jp=jp+1
               keepLooping = .true.
            else
               keepLooping = .false.
            end if
         end do
      end do
      !write(*,*) '[',mpi_rank,']mshRank[',mshRank,']numLinkedMapNodesInRank',numLinkedMapNodesInRank,'numRaw',numRawLinkedMapNodesInRank

      allocate(perMapLinkedNodes_i8(numLinkedMapNodesInRank,2))

      jp=0
      do ip=1,numRawLinkedMapNodesInRank
         if(rawPerMapLinkedNodes_i8(ip,1).ne.0) then
            jp=jp+1
            perMapLinkedNodes_i8(jp,:) = rawPerMapLinkedNodes_i8(ip,:)
         end if
      end do
      !if(mpi_rank.eq.0) write(*,*) '[',mpi_rank,']mshRank[',mshRank,']perMapLinkedNodes_i8(:,1)',perMapLinkedNodes_i8(:,1)

      deallocate(rawPerMapLinkedNodes_i8)

      !------------------------------------------------------------------------------------------------------------------

   end subroutine generate_mappedNodeLinksRankPar

   subroutine get_listElems2Par_in_parallel(numElemsMpiRank,listElemsMpiRank,listElems2Part,weightElems2Part,numElems2Part)
      implicit none
      integer(4), intent(in)  :: numElemsMpiRank,listElemsMpiRank(numElemsMpiRank)
      integer, allocatable, intent(out) :: listElems2Part(:),weightElems2Part(:)
      integer, intent(out) :: numElems2Part
      integer :: iElem, iElemG

      allocate(listElems2Part(numElemsMpiRank))
      allocate(weightElems2Part(numElemsMpiRank))
      !$acc parallel loop
      do iElem=1,numElemsMpiRank
         listElems2Part(iElem) = listElemsMpiRank(iElem)
         weightElems2Part(iElem) = 1
      end do
      !$acc end parallel loop

      numElems2Part = numElemsMpiRank

   end subroutine get_listElems2Par_in_parallel

   subroutine get_listElems2Par_periodic_in_parallel(numElemsMpiRank,numLinkedPerElemsSrl,numMapFacesSrl,numPerElemsSrl,listElemsMpiRank,linkedPerElemsSrl,linkedPerToMapElemsAndFacesSrl,listElems2Part,weightElems2Part,numElems2Part)
      implicit none
      integer(4),intent(in) :: numElemsMpiRank,numLinkedPerElemsSrl,numMapFacesSrl,numPerElemsSrl
      integer(4),intent(in) :: listElemsMpiRank(numElemsMpiRank),linkedPerElemsSrl(numLinkedPerElemsSrl,2),linkedPerToMapElemsAndFacesSrl(numMapFacesSrl,4)
      integer(4),allocatable,intent(out) :: listElems2Part(:),weightElems2Part(:)
      integer(4),intent(out) :: numElems2Part
      logical :: isElem2PartMpiRank(numElemsMpiRank)

      integer :: iElem,iElemG,iPos,iPerElem,iMapElem,auxCnt

      isElem2PartMpiRank(:) = .true.

      !Periodic elems
      do iElem=1,numLinkedPerElemsSrl
         iElemG = linkedPerElemsSrl(iElem,2)

         iPos = binarySearch_int_i4(listElemsMpiRank,iElemG)
         !write(*,*) '[',mpi_rank,']iElem',iElem,'iElemG',iElemG,'iPos',iPos
         if(iPos.ne.0) then !the element is in the slave rank, so not to part
            isElem2PartMpiRank(iPos) = .false.
         end if
      end do

      !Mapped elems
      do iElem=1,numMapFacesSrl
         iElemG = linkedPerToMapElemsAndFacesSrl(iElem,2)
         iPos = binarySearch_int_i4(listElemsMpiRank,iElemG)
         if(iPos.ne.0) then !the element is a mapped face (so will be 'slave'), so not to part
            isElem2PartMpiRank(iPos) = .false.
         end if
      end do

      numElems2Part=0
      do iElem=1,numElemsMpiRank
         if(isElem2PartMpiRank(iElem)) then
            numElems2Part=numElems2Part+1
         end if
      end do

      allocate(listElems2Part(numElems2Part))
      allocate(weightElems2Part(numElems2Part))

      auxCnt=0
      do iElem=1,numElemsMpiRank
         if(isElem2PartMpiRank(iElem)) then
            auxCnt=auxCnt+1
            listElems2Part(auxCnt)   = listElemsMpiRank(iElem)
            weightElems2Part(auxCnt) = 1
            do iPerElem = 1,numLinkedPerElemsSrl
               iElemG = linkedPerElemsSrl(iPerElem,1)
               if(listElems2Part(auxCnt).eq.iElemG) then
                  weightElems2Part(auxCnt) = weightElems2Part(auxCnt) + 1
               end if
            end do
            do iMapElem = 1,numMapFacesSrl
               iElemG = linkedPerToMapElemsAndFacesSrl(iMapElem,1)
               if(listElems2Part(auxCnt).eq.iElemG) then
                  weightElems2Part(auxCnt) = weightElems2Part(auxCnt) + 1
               end if
            end do
         end if
      end do

      !if(mpi_rank.eq.1) then
      !   write(*,*) 'listElems2Part',listElems2Part(:)
      !   write(*,*) 'weightElems2Part',weightElems2Part(:)
      !end if
      !write(*,*) '[',mpi_rank,']numElems2Part',numElems2Part,'auxCnt',auxCnt

   end subroutine get_listElems2Par_periodic_in_parallel

   subroutine get_rankPartitionBoundaryNodes_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,mshRank,numElemsInRank,numNodesInRank,numBoundsInRank,listNodesInRank_i8,connecInRank_i8,coordInRank,&
               listElemsBoundsInRank,listBoundFacesInRank,connecBoundFacesInRank_i8,numMasSlaNodesInRank,masSlaNodesInRank_i8,&
               numMshBoundNodesRankPar,numMpiBoundNodesRankPar,numInnerNodesRankPar,mpiBoundaryNodes_i8,mshBoundaryNodes_i8)
      implicit none
      integer(4), intent(in) :: mporder,mnnode,mnpbou
      integer(4), intent(in) :: gmsh2ijk(mnnode)
      integer(4), intent(in) :: mshRank,numElemsInRank,numNodesInRank,numBoundsInRank,numMasSlaNodesInRank
      integer(4), intent(in) :: listElemsBoundsInRank(numElemsInRank,maxBoundsPerElem),listBoundFacesInRank(numBoundsInRank)
      integer(8), intent(in) :: listNodesInRank_i8(numNodesInRank),connecInRank_i8(numElemsInRank,mnnode),connecBoundFacesInRank_i8(numBoundsInRank,mnpbou),masSlaNodesInRank_i8(numMasSlaNodesInRank,2)
      real(8),intent(in) :: coordInRank(numNodesInRank,3)

      integer(4), intent(out) :: numMshBoundNodesRankPar,numMpiBoundNodesRankPar,numInnerNodesRankPar
      integer(8), allocatable, intent(inout) :: mpiBoundaryNodes_i8(:),mshBoundaryNodes_i8(:)

      integer(4) :: nodeOwnedCnt(numNodesInRank)
      logical :: nodeInBoundary(numNodesInRank),nodeInMshBoundary(numNodesInRank),nodeInMpiBoundary(numNodesInRank)
      integer(4) :: i,ind,ind_f,iElemL,iElemG,iNode,indexNode,indexNode_inFace
      integer(8) :: iNodeG,iNodeG_inFace
      integer(4) :: auxCnt,mpiAuxCnt,mshAuxCnt
      integer(4) :: iBound,numBoundFacesToCheck,boundFacesToCheck(maxBoundsPerElem)
      integer(4),dimension(mnpbou) :: gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd
      integer(4):: gmshQuadInnerVertInd(4)
      !-------------------------------------------

      !$acc kernels
      nodeOwnedCnt(:)=0
      nodeInBoundary(:)=.false.
      nodeInMshBoundary(:)=.false.
      nodeInMpiBoundary(:)=.false.
      !$acc end kernels

      do iElemL=1,numElemsInRank
         do ind = 1,mnnode
            iNodeG = connecInRank_i8(iElemL,ind)
            indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
            nodeOwnedCnt(indexNode) = nodeOwnedCnt(indexNode) + 1
         end do
      end do

      do ind=1,numMasSlaNodesInRank
         iNodeG = masSlaNodesInRank_i8(ind,1)
         indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
         nodeOwnedCnt(indexNode) = 2

         iNodeG = masSlaNodesInRank_i8(ind,2)
         indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
         nodeOwnedCnt(indexNode) = 0
      end do

      !check if face is in boundary

      call get_gmshHexHOFacesIndex(mporder,mnpbou,gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd)
      call get_gmshQuadHOInnerVertIndex(mporder,gmshQuadInnerVertInd)

      do iElemL=1,numElemsInRank

         numBoundFacesToCheck=0
         do iBound=1,maxBoundsPerElem
            if(listElemsBoundsInRank(iElemL,iBound).ne.0) then
               numBoundFacesToCheck=numBoundFacesToCheck+1
               boundFacesToCheck(numBoundFacesToCheck) = listElemsBoundsInRank(iElemL,iBound)
            end if
         end do

         call find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,gmshHexFaceFrontInd, listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)
         call find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,gmshHexFaceLeftInd,  listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)
         call find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,gmshHexFaceTopInd,   listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)
         call find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,gmshHexFaceBackInd,  listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)
         call find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,gmshHexFaceRightInd, listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)
         call find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,gmshHexFaceBottomInd,listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundFacesInRank,connecBoundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)

      end do

      numMshBoundNodesRankPar=0
      numMpiBoundNodesRankPar=0
      numInnerNodesRankPar=0

      do iNode=1,numNodesInRank
         if(nodeInBoundary(iNode)) then !node is boundary
            if(nodeInMshBoundary(iNode)) numMshBoundNodesRankPar=numMshBoundNodesRankPar+1
            if(nodeInMpiBoundary(iNode)) numMpiBoundNodesRankPar=numMpiBoundNodesRankPar+1
         else !node is inner
            numInnerNodesRankPar=numInnerNodesRankPar+1
         end if
      end do

      !write(*,*) '1.#rank[',mpi_rank,'] numNodesInRank',numNodesInRank,'bMshN',numMshBoundNodesRankPar,'bMpiN',numMpiBoundNodesRankPar,'iN',numInnerNodesRankPar
      allocate(mpiBoundaryNodes_i8(numMpiBoundNodesRankPar))
      allocate(mshBoundaryNodes_i8(numMshBoundNodesRankPar))

      auxCnt=0
      mpiAuxCnt=0
      mshAuxCnt=0
      do iNode=1,numNodesInRank
         iNodeG=listNodesInRank_i8(iNode)
         if(nodeInBoundary(iNode)) then !node is boundary
            auxCnt=auxCnt+1
            if(nodeInMshBoundary(iNode)) then
               mshAuxCnt=mshAuxCnt+1
               mshBoundaryNodes_i8(mshAuxCnt) = iNodeG
            end if
            if(nodeInMpiBoundary(iNode)) then
               mpiAuxCnt=mpiAuxCnt+1
               mpiBoundaryNodes_i8(mpiAuxCnt) = iNodeG
            end if
            !write(*,*) '1.#rank[',mpi_rank,']mshRank(',mshRank,')iNodeG',iNodeG,'auxCnt',auxCnt
         end if
      end do

      !write(*,*) '#rank[',mpi_rank,']  nodesInRank ',numNodesMshRank, " bN ", numBNodesRankPar, " iN ",numINodesRankPar,'i',i

#if _CHECK_
      write(aux_string_mpirank_deb,'(I0)') mpi_rank
      write(aux_string_mshrank_deb,'(I0)') mshRank

      file_name_deb = 'boundaryMpiNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numMpiBoundNodesRankPar
         iNodeG=mpiBoundaryNodes_i8(iNode)
         indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)

         write(1,fmt_csv_deb) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

      file_name_deb = 'boundaryMshNodes_mpiRank'//trim(aux_string_mpirank_deb)//'_mshRank_'//trim(aux_string_mshrank_deb)//'.csv'
      open(1, file=file_name_deb)
      write(1,*) 'X,Y,Z,iNodeG'
      do iNode=1,numMshBoundNodesRankPar
         iNodeG=mshBoundaryNodes_i8(iNode)
         indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)

         write(1,fmt_csv_deb) coordInRank(indexNode,1),coordInRank(indexNode,2),coordInRank(indexNode,3),iNodeG
      end do
      close(1)

#endif
   end subroutine get_rankPartitionBoundaryNodes_in_parallel

   subroutine find_mpi_and_msh_boundaries(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,numBoundsInRank,iElemL,gmsh2ijk,face2ijk,listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,listBoundsInRank,boundFacesInRank_i8,numBoundFacesToCheck,boundFacesToCheck,nodeInBoundary,nodeInMshBoundary,nodeInMpiBoundary,gmshQuadInnerVertInd)
      implicit none
      integer(4),intent(in) :: mporder,mnnode,mnpbou
      integer(4),intent(in) :: numElemsInRank,numNodesInRank,numBoundsInRank
      integer(4),intent(in) :: iElemL,gmsh2ijk(mnnode),face2ijk(mnpbou),nodeOwnedCnt(numNodesInRank),listBoundsInRank(numBoundsInRank)
      integer(8),intent(in) :: listNodesInRank_i8(numNodesInRank),connecInRank_i8(numElemsInRank,mnnode),boundFacesInRank_i8(numBoundsInRank,mnpbou)
      integer(4),intent(in) :: numBoundFacesToCheck,boundFacesToCheck(maxBoundsPerElem)
      integer(4),intent(in) :: gmshQuadInnerVertInd(4)
      logical,intent(out) :: nodeInBoundary(numNodesInRank),nodeInMshBoundary(numNodesInRank),nodeInMpiBoundary(numNodesInRank)
      integer(4) :: ii,jj,ind,ind_f,iFace,indexNode,indexNode_inFace,ind_gmsh
      integer(8) :: iNodeG,iNodeG_inFace,iNodeG_ToCheck
      integer(4) :: iBound,iBoundG,iElemG,indexBound
      logical :: isMshBound

      !# 1.Whatever face ------------------------------------------------------
      ind = face2ijk(mporder+3)
      iNodeG_inFace = connecInRank_i8(iElemL,ind)
      indexNode_inFace = binarySearch_int_i8(listNodesInRank_i8,iNodeG_inFace)
      if(nodeOwnedCnt(indexNode_inFace).eq.1) then !this face is boundary

         !----------------------------------------------------------------
         isMshBound=.false.
         do iBound=1,numBoundFacesToCheck
            iBoundG = boundFacesToCheck(iBound)
            indexBound = binarySearch_int_i4(listBoundsInRank,iBoundG)
            checkLoop: do jj=1,4
               ind_gmsh = gmshQuadInnerVertInd(jj)
               iNodeG_ToCheck=boundFacesInRank_i8(indexBound,ind_gmsh)
               if(iNodeG_inFace.eq.iNodeG_ToCheck) then
                  isMshBound = .true.
                  exit checkLoop
               end if
               !write(*,*) '[',mpi_rank,']iElemL',iElemL,'check iBoundG',iBoundG,'jj',jj,'iNodeG2c',iNodeG_ToCheck,'iNodeGinF',iNodeG_inFace
            end do checkLoop
            !write(*,*) '[',mpi_rank,']iElemL',iElemL,'check iBoundG',iBoundG,'isMshBound',isMshBound
         end do
         !----------------------------------------------------------------
         !----------------------------------------------------------------
         do iFace=1,mnpbou
            ind_f = face2ijk(iFace)
            iNodeG = connecInRank_i8(iElemL,ind_f)
            indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
            nodeInBoundary(indexNode) = .true.
            if(isMshBound) then
               nodeInMshBoundary(indexNode) = .true.
            else
               nodeInMpiBoundary(indexNode) = .true.
            end if
         end do
         !----------------------------------------------------------------
      end if

   end subroutine find_mpi_and_msh_boundaries

   subroutine reorder_nodes_in_mshRank(mporder,mnnode,mnpbou,gmsh2ijk,vtk2ijk,auxNewOrderIndex,auxVTKorder,mshRank,numMshRanks2Part,&
               numElemsInRank,numNodesInRank,numBoundsInRank,numMasSlaNodesInRank,numPerMapinkedNodesInRank,sizeConnecVTK,linOutput,&
               elemGidInRank,listNodesInRank_i8,connecInRank_i8,iNodeStartPar_i8,boundFacesInRank_i8,masSlaNodesInRank_i8,perMapLinkedNodesInRank_i8,&
               globalIdSrlInRank_i8,globalIdSrlOrderedInRank_i8,globalIdParInRank_i8,connecVTKInRank,connecParOrigInRank,boundFacesInRank,&
               masSlaNodesInRank,perMapLinkedNodesInRank)
      implicit none
      integer(4),intent(in) :: mporder,mnnode,mnpbou
      integer(4),intent(in),dimension(mnnode) :: gmsh2ijk,vtk2ijk,auxNewOrderIndex,auxVTKorder
      integer(4),intent(in) :: mshRank,numMshRanks2Part,numElemsInRank,numNodesInRank,numBoundsInRank,numMasSlaNodesInRank,numPerMapinkedNodesInRank,sizeConnecVTK
      logical,intent(in)    :: linOutput
      integer(4),intent(in) :: elemGidInRank(numElemsInRank)
      integer(8),intent(in) :: listNodesInRank_i8(numNodesInRank),connecInRank_i8(numElemsInRank,mnnode)
      integer(8),dimension(0:numMshRanks2Part-1),intent(in) :: iNodeStartPar_i8
      integer(8),intent(in) :: boundFacesInRank_i8(numBoundsInRank,mnpbou),masSlaNodesInRank_i8(numMasSlaNodesInRank,2),perMapLinkedNodesInRank_i8(numPerMapinkedNodesInRank,2)

      integer(8),intent(out) :: globalIdSrlInRank_i8(numNodesInRank),globalIdParInRank_i8(numNodesInRank)
      integer(8),intent(out),dimension(numNodesInRank,2) :: globalIdSrlOrderedInRank_i8
      integer(4),intent(out) :: connecVTKInRank(sizeConnecVTK),connecParOrigInRank(numElemsInRank,mnnode)
      integer(4),intent(out) :: boundFacesInRank(numBoundsInRank,mnpbou),masSlaNodesInRank(numMasSlaNodesInRank,2),perMapLinkedNodesInRank(numPerMapinkedNodesInRank,2)

      integer(4) :: m,indConn,indexIJK,indexVTK,indexGMSH,indexNew,nodeIndexCnt,indPosListNodes
      integer(4) :: iNodeL,iElemL,iElemG,iNode,iBound,iLink
      integer(8) :: iNodeGsrl,iNodeGpar
      integer(4) :: ii,jj,kk,igp,jgp,kgp
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      integer(8),dimension(mnnode) :: auxNodeNewOrderInElem_i8
      integer(4),dimension(numNodesInRank) :: isNodeAdded
      integer(4) :: aux_ibound_GtoL(mnpbou), aux_iMasSla_GtoL(2)

      !$acc kernels
      isNodeAdded(:)=-1
      globalIdSrlOrderedInRank_i8(:,:)=-1
      !$acc end kernels
      nodeIndexCnt = 0
      indConn = -1

      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)
      !----------------------------------------------------------------------------------------------------------

      do iElemL=1,numElemsInRank
         iElemG = elemGidInRank(iElemL)
         auxNodeNewOrderInElem_i8(:)=0

         do indexIJK=1,mnnode
            indexGMSH = gmsh2ijk(indexIJK)
            iNodeGsrl = connecInRank_i8(iElemL,indexGMSH)

            indexNew = auxNewOrderIndex(indexIJK)

            auxNodeNewOrderInElem_i8(indexNew) = iNodeGsrl
         end do

         do m=1,mnnode
            iNodeGsrl = auxNodeNewOrderInElem_i8(m)
            indPosListNodes = binarySearch_int_i8(listNodesInRank_i8,iNodeGsrl)

            if(isNodeAdded(indPosListNodes) < 0) then !node not added put it in the list
               nodeIndexCnt=nodeIndexCnt+1
               iNodeL = nodeIndexCnt

               isNodeAdded(indPosListNodes)  = iNodeL

               iNodeGPar = int(iNodeL,8) + iNodeStartPar_i8(mshRank) - 1

               globalIdSrlInRank_i8(iNodeL) = iNodeGsrl
               globalIdParInRank_i8(iNodeL) = iNodeGPar

               globalIdSrlOrderedInRank_i8(iNodeL,1) = iNodeGsrl
               globalIdSrlOrderedInRank_i8(iNodeL,2) = iNodeL
            else
               iNodeL = isNodeAdded(indPosListNodes)
            endif

            connecParOrigInRank(iElemL,m) = iNodeL

            if(.not.(linOutput)) then
               indexVTK = auxVTKorder(m)
               indConn = (iElemL-1)*mnnode + indexVTK
               connecVTKInRank(indConn) = iNodeL
            end if
         end do

         if(linOutput) then
            indConn = (iElemL-1)*((mporder**3)*8)
            do ii=0,(mporder-1)
               do jj=0,(mporder-1)
                  do kk=0,(mporder-1)
                     
                     !1.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !2.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !3.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !4.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !5.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !6.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !7.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !8.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = auxNewOrderIndex(indexIJK)
                     iNodeL = connecParOrigInRank(iElemL,indexNew)
                     indConn = indConn+1
                     connecVTKInRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                  end do
               end do
            end do
         end if

         !write(*,*) '[',mpi_rank,']iElemG ',iElemG,' connecParOrig ',connecParOrig(iElemL,:)
      end do

      call quicksort_matrix_int8(globalIdSrlOrderedInRank_i8,1)

      !@TODO: REORDER THE BOUNDARIES ALSO!
      do iBound=1,numBoundsInRank
         aux_ibound_GtoL(:) = 0
         !write(*,*) '1.reorder ibound[',mpi_rank,']npbouG',boundFacesInRank(iBound,:)
         do iNode=1,mnpbou
            iNodeGsrl = boundFacesInRank_i8(iBound,iNode)
            indexNew = binarySearch_int_i8(globalIdSrlOrderedInRank_i8(:,1),iNodeGsrl)
            iNodeL = globalIdSrlOrderedInRank_i8(indexNew,2)
            !if(mpi_rank.eq.2)write(*,*) 'ibound[',mpi_rank,']',iBound,'iNodeGsrl',iNodeGsrl,'iNodeL',iNodeL
            aux_ibound_GtoL(iNode) = iNodeL
         end do

         boundFacesInRank(iBound,:) = aux_ibound_GtoL(:)
         !write(*,*) '2.reorder ibound[',mpi_rank,']npbouG',boundFacesInRank(iBound,:)
      end do

      !----------------------------------------------------------------
      do iLink=1,numMasSlaNodesInRank
         !master
         iNodeGsrl = masSlaNodesInRank_i8(iLink,1)
         indexNew = binarySearch_int_i8(globalIdSrlOrderedInRank_i8(:,1),iNodeGsrl)
         iNodeL = globalIdSrlOrderedInRank_i8(indexNew,2)
         aux_iMasSla_GtoL(1) = iNodeL
         !slave
         iNodeGsrl = masSlaNodesInRank_i8(iLink,2)
         indexNew = binarySearch_int_i8(globalIdSrlOrderedInRank_i8(:,1),iNodeGsrl)
         iNodeL = globalIdSrlOrderedInRank_i8(indexNew,2)
         aux_iMasSla_GtoL(2) = iNodeL

         masSlaNodesInRank(iLink,:) = aux_iMasSla_GtoL(:)
      end do

      !reordeno el vector per la columna slaves, ja que son unics
      call quicksort_matrix_int4(masSlaNodesInRank,2)
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      do iLink=1,numPerMapinkedNodesInRank
         !periodic face
         iNodeGsrl = perMapLinkedNodesInRank_i8(iLink,1)
         indexNew = binarySearch_int_i8(globalIdSrlOrderedInRank_i8(:,1),iNodeGsrl)
         iNodeL = globalIdSrlOrderedInRank_i8(indexNew,2)
         aux_iMasSla_GtoL(1) = iNodeL
         !map face
         iNodeGsrl = perMapLinkedNodesInRank_i8(iLink,2)
         indexNew = binarySearch_int_i8(globalIdSrlOrderedInRank_i8(:,1),iNodeGsrl)
         iNodeL = globalIdSrlOrderedInRank_i8(indexNew,2)
         aux_iMasSla_GtoL(2) = iNodeL

         perMapLinkedNodesInRank(iLink,:) = aux_iMasSla_GtoL(:)
      end do

      !----------------------------------------------------------------

   end subroutine reorder_nodes_in_mshRank

   subroutine generate_new_nodeOrder_and_connectivity(mporder,mnnode,gmsh2ijk,vtk2ijk,newOrderIJK,auxNewOrderIndex,auxVTKorder)
      integer(4),intent(in) :: mporder,mnnode
      integer(4),intent(in) :: gmsh2ijk(mnnode),vtk2ijk(mnnode),newOrderIJK(mnnode)
      integer(4),intent(out) :: auxNewOrderIndex(mnnode),auxVTKorder(mnnode)
      integer(4) :: i,j,k,indexIJK,indexNew,indexVTK,indexCGNS

      do k = 0,mporder
         do i = 0,mporder
            do j = 0,mporder
               !indexIJK = ((mporder+1)**2)*k+(mporder+1)*i+j+1
               indexIJK = get_indexIJK_sod2d(mporder,i,j,k)

               indexNew = newOrderIJK(indexIJK)
               indexVTK = vtk2ijk(indexIJK)

               auxNewOrderIndex(indexIJK) = indexNew
               auxVTKorder(indexNew) = indexVTK

               !write(*,*) 'test->indexIJK ', indexIJK, ' iNew ', indexNew,' aux ', auxCGNSorder(indexNew)
            end do
         end do
      end do

   end subroutine generate_new_nodeOrder_and_connectivity

   subroutine generate_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,numMpiBoundaryNodes,mpiBoundaryNodes_i8_jv,globalIdSrl_i8_jv,globalIdSrlOrdered_i8_jm,&
               nodesToComm_jv,commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv,numNodesToCommMshRank,numMshRanksWithComms)
   !generate a matrix with the comm schemes for shared nodes between procs
      implicit none
      integer(4),intent(in) :: numMshRanks2Part,numMpiBoundaryNodes(numMshRanksInMpiRank)
      integer(4),intent(in) :: numMshRanksInMpiRank,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part)
      type(jagged_vector_int8),intent(in) :: mpiBoundaryNodes_i8_jv,globalIdSrl_i8_jv
      type(jagged_matrix_int8),intent(in) :: globalIdSrlOrdered_i8_jm
      type(jagged_vector_int4),intent(inout) :: nodesToComm_jv,commsMemPosInLoc_jv,commsMemSize_jv,commsMemPosInNgb_jv,ranksToComm_jv
      integer(4),allocatable,intent(inout) :: numNodesToCommMshRank(:),numMshRanksWithComms(:)

      integer(4) :: iMshRank,mshRank,iMshRankTrgt,mpiRank,iAux,memPos,memSize,memDisp
      integer(4) :: mshRankOrig,mshRankTrgt,numBNOrig,numBNTrgt,mpiRankTrgt
      integer(4) :: iNode,iPos,iNodeL,iNodeLPos
      integer(8) :: iNodeGSrl

      integer(4), dimension(0:numMshRanks2Part-1) :: auxCommSchemeNumNodes
      integer(4),allocatable :: auxCommSchemeMemPos(:,:),auxCommSchemeMemPosAll(:)

      integer(4), dimension(0:numMshRanks2Part-1) :: vecNumMpiBN,vecMemDispMpiBN
      integer(8),allocatable :: arrayMpiBNInMpiRank_i8(:),auxMpiBNMshRank_i8(:)
      integer(4) :: numMpiBNInMpiRank,auxMemDispMpiBN

      integer(4) :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      !------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) "--| Generating MPI Comm scheme"

      !NEW METHOD
      !---------------------------------------------------------------------

      vecNumMpiBN(:) = 0
      vecMemDispMpiBN(:) = 0
      numMpiBNInMpiRank = 0

      auxMemDispMpiBN = 0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vecNumMpiBN(mshRank) = numMpiBoundaryNodes(iMshRank)
         numMpiBNInMpiRank = numMpiBNInMpiRank + numMpiBoundaryNodes(iMshRank)
         if(iMshRank.ne.1) auxMemDispMpiBN = auxMemDispMpiBN + numMpiBoundaryNodes(iMshRank-1)
         vecMemDispMpiBN(mshRank) = auxMemDispMpiBN
      end do

      !if(mpi_rank.eq.7) then
      !   write(*,*) 'numMpiBNInMpiRank',numMpiBNInMpiRank
      !   write(*,*) 'vecNumMpiBN(',mpi_rank,')',vecNumMpiBN(:)
      !   write(*,*) 'vecMemDispMpiBN(',mpi_rank,')',vecMemDispMpiBN(:)
      !endif

      allocate(arrayMpiBNInMpiRank_i8(numMpiBNInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         do iAux=1,numMpiBoundaryNodes(iMshRank)
            memPos = vecMemDispMpiBN(mshRank) + iAux
            arrayMpiBNInMpiRank_i8(memPos) = mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      !if(mpi_rank.eq.7) then
      !   write(*,*) 'arrayMpiBNInMpiRank_i8(',mpi_rank,')',arrayMpiBNInMpiRank_i8(:)
      !endif

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vecNumMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecNumMpiBN(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      call MPI_Win_create(vecMemDispMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecMemDispMpiBN(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.7) write(*,*) 'vecNumMpiBN',vecNumMpiBN(:)
      !if(mpi_rank.eq.7) write(*,*) 'vecMemDispMpiBN',vecMemDispMpiBN(:)

      allocate(numNodesToCommMshRank(numMshRanksInMpiRank))
      allocate(numMshRanksWithComms(numMshRanksInMpiRank))
      allocate(nodesToComm_jv%vector(numMshRanksInMpiRank))
      allocate(ranksToComm_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemSize_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemPosInLoc_jv%vector(numMshRanksInMpiRank))
      allocate(commsMemPosInNgb_jv%vector(numMshRanksInMpiRank))
      allocate(auxCommSchemeMemPos(numMshRanksInMpiRank,numMshRanks2Part))

      !--------------------------------------------------------------------------------------------
      call MPI_Barrier(app_comm,mpi_err)

      win_buffer_size = mpi_int8_size*numMpiBNInMpiRank
      call MPI_Win_create(arrayMpiBNInMpiRank_i8,win_buffer_size,mpi_int8_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)

      call MPI_Barrier(app_comm,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numBNOrig = numMpiBoundaryNodes(iMshRank)

         numNodesToCommMshRank(iMshRank)=0
         auxCommSchemeNumNodes(:)=0
         auxCommSchemeMemPos(iMshRank,:)=0

         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then

               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispMpiBN(mshRankTrgt)
               memSize = vecNumMpiBN(mshRankTrgt)
               allocate(auxMpiBNMshRank_i8(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxMpiBNMshRank_i8,memSize,mpi_datatype_int8,mpiRankTrgt,target_displacement,memSize,mpi_datatype_int8,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               !if(mpi_rank.eq.0 .and. mshRankTrgt.eq.1) write(*,*) 'memSize',memSize,'td',vecMemDispMpiBN(mshRankTrgt),'auxMpiBN(',mpi_rank,')',auxMpiBNMshRank
               !if(mshRankOrig.eq.7 .and. mshRankTrgt .eq. 10) then
               !   write(*,*) 'mshRankOrig 7 & mshRankTrgt 10 (7):'
               !   do iNode=1,numBNOrig
               !      write(*,*) iNode,mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)
               !   end do
               !   write(*,*) 'mshRankOrig 7 & mshRankTrgt 10 (10):'
               !   do iNode=1,vecNumMpiBN(mshRankTrgt)
               !      write(*,*) iNode,auxMpiBNMshRank(iNode)
               !   end do
               !end if

               do iNode=1,numBNOrig
                  iNodeGSrl = mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)

                  iPos = binarySearch_int_i8(auxMpiBNMshRank_i8,iNodeGSrl)
                  if(iPos.ne.0) then
                     numNodesToCommMshRank(iMshRank) = numNodesToCommMshRank(iMshRank) + 1
                     auxCommSchemeNumNodes(mshRankTrgt) = auxCommSchemeNumNodes(mshRankTrgt)+1

                     !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBNO',numBNOrig,'mshRankT',mshRankTrgt,'memPos',memPos,'numBNT',numBNTrgt
                     !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'mshRankT',mshRankTrgt,'iNGO',iNodeGOrig,'iNGT',iNodeGTrgt
                  end if
                  !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBNO',numBNOrig,'mshRankT',mshRankTrgt,'memPos',memPos,'numBNT',numBNTrgt
               end do

               deallocate(auxMpiBNMshRank_i8)
            end if

         end do

         ! setting numRanksWithsComms
         numMshRanksWithComms(iMshRank) = 0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               numMshRanksWithComms(iMshRank) = numMshRanksWithComms(iMshRank) + 1
            end if
         end do

         allocate(ranksToComm_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))
         allocate(commsMemSize_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))
         allocate(commsMemPosInLoc_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))
         allocate(commsMemPosInNgb_jv%vector(iMshRank)%elems(numMshRanksWithComms(iMshRank)))

         iAux=0
         memPos=1
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            auxCommSchemeMemPos(iMshRank,iMshRankTrgt)=memPos
            if(auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               iAux=iAux+1
               mshRankTrgt = iMshRankTrgt-1
               ranksToComm_jv%vector(iMshRank)%elems(iAux) = mshRankTrgt

               commsMemSize_jv%vector(iMshRank)%elems(iAux) = auxCommSchemeNumNodes(mshRankTrgt)
               commsMemPosInLoc_jv%vector(iMshRank)%elems(iAux) = memPos
               memPos=memPos+auxCommSchemeNumNodes(mshRankTrgt)
            end if
         end do

         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'csNumNodes->',auxCommSchemeNumNodes(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'numRanksWC->',numMshRanksWithComms(iMshRank)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'ranksToComm->',ranksToComm_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInLoc->',commsMemPosInLoc_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemSize->',commsMemSize_jv%vector(iMshRank)%elems(:)

         !write(*,*) 'rank[',mpi_rank,']mshRankOrig',mshRankOrig,'numNodesToCommMshRank',numNodesToCommMshRank(iMshRank),'numBNOrig',numBNOrig,'cSNN',auxCommSchemeNumNodes(:),&
         !          'numMshRanksWithComms',numMshRanksWithComms(iMshRank),'r2C',ranksToComm_jv%vector(iMshRank)%elems(:)

         allocate(nodesToComm_jv%vector(iMshRank)%elems(numNodesToCommMshRank(iMshRank)))
         iAux=0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispMpiBN(mshRankTrgt)
               memSize = vecNumMpiBN(mshRankTrgt)
               allocate(auxMpiBNMshRank_i8(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxMpiBNMshRank_i8,memSize,mpi_datatype_int8,mpiRankTrgt,target_displacement,memSize,mpi_datatype_int8,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numBNOrig
                  iNodeGSrl = mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)

                  iPos = binarySearch_int_i8(auxMpiBNMshRank_i8,iNodeGSrl)
                  if(iPos.ne.0) then
                     iAux=iAux+1

                     iNodeLPos = binarySearch_int_i8(globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
                     iNodeL = globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(iNodeLPos,2)
                     !write(*,*) 'iNodeL',iNodeL,'iNodeG',iNodeGSrl,'mshRank',mshRankTrgt
                     nodesToComm_jv%vector(iMshRank)%elems(iAux) = iNodeL
                  end if

               end do

               deallocate(auxMpiBNMshRank_i8)
            end if
         end do

      end do

      call MPI_Win_free(window_id,mpi_err)

      !---------------------------------------------------------------------

      allocate(auxCommSchemeMemPosAll(numMshRanks2Part**2))

      auxCommSchemeMemPosAll(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numMshRanks2Part
            iPos = mshRank*numMshRanks2Part+iMshRankTrgt
            auxCommSchemeMemPosAll(iPos) = auxCommSchemeMemPos(iMshRank,iMshRankTrgt)
         end do
      end do
      !write(*,*) '1.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      ! Create the window
      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*(numMshRanks2Part**2)

      call MPI_Win_create(auxCommSchemeMemPosAll,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mpiRank = mapMshRankToMpiRank(iMshRank)
         memSize = numMshRanks2Part
         memPos = (iMshRank-1)*numMshRanks2Part+1
         target_displacement = (iMshRank-1)*numMshRanks2Part

         call MPI_Get(auxCommSchemeMemPosAll(memPos),memSize,mpi_datatype_int,mpiRank,target_displacement,memSize,mpi_datatype_int,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0,window_id,mpi_err)
      call MPI_Win_free(window_id,mpi_err)
      !--------------------------------------------------------------------------------------
      !write(*,*) '2.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numMshRanksWithComms(iMshRank)
            mshRankTrgt=ranksToComm_jv%vector(iMshRank)%elems(iMshRankTrgt)
            iPos = mshRankTrgt*numMshRanks2Part+mshRankOrig+1
            commsMemPosInNgb_jv%vector(iMshRank)%elems(iMshRankTrgt) = auxCommSchemeMemPosAll(iPos)
            !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'iMRT',iMshRankTrgt,'mRT',mshRankTrgt,'iPos',iPos
         end do
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInNgb->',commsMemPosInNgb_jv%vector(iMshRank)%elems(:)
      end do

      deallocate(auxCommSchemeMemPos)
      deallocate(auxCommSchemeMemPosAll)

   end subroutine generate_mpi_comm_scheme_parallel

   subroutine generate_dof_and_boundary_mpi_comm_scheme_parallel(numMshRanks2Part,numMshRanksInMpiRank,mshRanksInMpiRank,mapMshRankToMpiRank,&
               numNodesMshRank,numMshBoundaryNodes,mshBoundaryNodes_i8_jv,numMpiBoundaryNodes,mpiBoundaryNodes_i8_jv,globalIdSrl_i8_jv,globalIdSrlOrdered_i8_jm,&
               numBoundaryNodesMshRank,boundaryNodes_jv,numDoFMshRank,dofNodes_jv,&
               bnd_nodesToComm_jv,bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv,bnd_numNodesToCommMshRank,bnd_numMshRanksWithComms)
      !generate a matrix with the comm schemes for shared nodes between procs
      integer,intent(in) :: numMshRanks2Part,numNodesMshRank(numMshRanksInMpiRank),numMpiBoundaryNodes(numMshRanksInMpiRank),numMshBoundaryNodes(numMshRanksInMpiRank)
      integer,intent(in) :: numMshRanksInMpiRank,mshRanksInMpiRank(numMshRanksInMpiRank),mapMshRankToMpiRank(numMshRanks2Part)
      type(jagged_vector_int8),intent(in) :: mshBoundaryNodes_i8_jv,mpiBoundaryNodes_i8_jv,globalIdSrl_i8_jv
      type(jagged_matrix_int8),intent(in) :: globalIdSrlOrdered_i8_jm
      integer,dimension(numMshRanksInMpiRank),intent(out) :: numBoundaryNodesMshRank,numDoFMshRank
      type(jagged_vector_int4),intent(inout) :: boundaryNodes_jv,dofNodes_jv
      type(jagged_vector_int4),intent(inout) :: bnd_nodesToComm_jv,bnd_commsMemPosInLoc_jv,bnd_commsMemSize_jv,bnd_commsMemPosInNgb_jv,bnd_ranksToComm_jv
      integer(4),allocatable,intent(inout) :: bnd_numNodesToCommMshRank(:),bnd_numMshRanksWithComms(:)

      integer(4) :: iMshRank,mshRank,iMshRankTrgt,mpiRank,memPos,memSize,memDisp
      integer(4) :: iAux,auxCntBoundaryNode,auxCnt1,auxCnt2
      integer(4) :: mshRankOrig,mshRankTrgt,mpiRankTrgt,numMpiBNOrig,numMpiBNTrgt,numMshBNOrig,numMshBNTrgt,numBndMpiNodeOrig,numBndMpiNodeTrgt
      integer(4) :: iNode,iNodeBnd,iNodeMpi,iPos,iNodeL,iNodeLPos,iPosMshBN
      integer(8) :: iNodeGSrl

      integer(4) :: numBndMpiBoundaryNodes(numMshRanksInMpiRank),numBndMpiBoundaryNodesAll(numMshRanks2Part)
      integer(4), dimension(0:numMshRanks2Part-1) :: bnd_auxCommSchemeNumNodes
      integer(4),allocatable :: bnd_auxCommSchemeMemPos(:,:),bnd_auxCommSchemeMemPosAll(:)
      logical,allocatable :: auxListBN(:)

      type(jagged_vector_int8) :: bndMpiBoundaryNodes_i8_jv
      integer(4), dimension(0:numMshRanks2Part-1) :: vecNumMshBN,vecMemDispMshBN,vecNumBndMpiBN,vecMemDispBndMpiBN
      integer(8),allocatable :: arrayMshBNInMpiRank_i8(:),auxMshBNMshRank_i8(:),arrayBndMpiBNInMpiRank_i8(:),auxBndMpiBNMshRank_i8(:)
      integer(4) :: numMshBNInMpiRank,numBndMpiBNInMpiRank,auxMemDispMshBN,auxMemDispBndMpiBN

      integer(4) :: window_id
      integer(KIND=MPI_ADDRESS_KIND) :: win_buffer_size
      integer(KIND=MPI_ADDRESS_KIND) :: target_displacement

      !------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) "--| Generating DoF & Boundary MPI Comm scheme"

      !NEW APPROACH
      !---------------------------------------------------------------------
      vecNumMshBN(:) = 0
      vecMemDispMshBN(:) = 0
      numMshBNInMpiRank = 0

      auxMemDispMshBN = 0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vecNumMshBN(mshRank) = numMshBoundaryNodes(iMshRank)
         numMshBNInMpiRank = numMshBNInMpiRank + numMshBoundaryNodes(iMshRank)
         !if(iMshRank.ne.1) vecMemDispMshBN(mshRank) = vecMemDispMshBN(mshRank) + numMshBoundaryNodes(iMshRank-1)
         if(iMshRank.ne.1) auxMemDispMshBN = auxMemDispMshBN + numMshBoundaryNodes(iMshRank-1)
         vecMemDispMshBN(mshRank) = auxMemDispMshBN
      end do
      !write(*,*) 'numMshBNInMpiRank',numMshBNInMpiRank
      !write(*,*) 'vecNumMshBN(',mpi_rank,')',vecNumMshBN(:)
      !write(*,*) 'vecMemDispMshBN(',mpi_rank,')',vecMemDispMshBN(:)
      allocate(arrayMshBNInMpiRank_i8(numMshBNInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         do iAux=1,numMshBoundaryNodes(iMshRank)
            memPos = vecMemDispMshBN(mshRank) + iAux
            arrayMshBNInMpiRank_i8(memPos) = mshBoundaryNodes_i8_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vecNumMshBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecNumMshBN(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      call MPI_Win_create(vecMemDispMshBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecMemDispMshBN(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------------------------------
      !lets determine ALL the boundary nodes in each rank
      !the boundary nodes in a rank will be:
      ! 1. the boundary nodes in the container mshBoundaryNodes_jv of each rank
      ! 2. plus the nodes that are in mpiBoundaryNode_jv in the rank and in the boundary nodes container mshBoundaryNodes_jv in the other ranks

      !--------------------------------------------------------------------------------------------
      call MPI_Barrier(app_comm,mpi_err)
      win_buffer_size = mpi_int8_size*numMshBNInMpiRank
      call MPI_Win_create(arrayMshBNInMpiRank_i8,win_buffer_size,mpi_int8_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Barrier(app_comm,mpi_err)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numMshBNOrig = numMshBoundaryNodes(iMshRank)
         numMpiBNOrig = numMpiBoundaryNodes(iMshRank)

         auxCntBoundaryNode=0

         allocate(auxListBN(numNodesMshRank(iMshRank)))
         auxListBN(:)=.false.
         !1.
         do iNode=1,numMshBNOrig
            iNodeGSrl = mshBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)
            iNodeLPos = binarySearch_int_i8(globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
            iNodeL    = globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(iNodeLPos,2)

            auxListBN(iNodeL) = .true.
            auxCntBoundaryNode = auxCntBoundaryNode + 1
         end do

         !2.
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispMshBN(mshRankTrgt)
               memSize = vecNumMshBN(mshRankTrgt)
               allocate(auxMshBNMshRank_i8(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxMshBNMshRank_i8,memSize,mpi_datatype_int8,mpiRankTrgt,target_displacement,memSize,mpi_datatype_int8,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numMpiBNOrig
                  iNodeGSrl = mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)
                  iPosMshBN = binarySearch_int_i8(mshBoundaryNodes_i8_jv%vector(iMshRank)%elems,iNodeGSrl)
                  if(iPosMshBN.eq.0) then !if is not already in a boundary face
                     iPos = binarySearch_int_i8(auxMshBNMshRank_i8,iNodeGSrl)
                     if(iPos.ne.0) then
                        iNodeLPos = binarySearch_int_i8(globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
                        iNodeL    = globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(iNodeLPos,2)
                        if(.not.(auxListBN(iNodeL))) then
                           auxListBN(iNodeL) = .true.
                           auxCntBoundaryNode = auxCntBoundaryNode + 1
                        end if
                     end if
                  end if
               end do

               deallocate(auxMshBNMshRank_i8)

            end if
         end do

         !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numMshBN',numMshBNOrig,'numMpiBN',numMpiBNOrig,'addCntBN',auxCntBoundaryNode
         numBoundaryNodesMshRank(iMshRank) = auxCntBoundaryNode
         numDoFMshRank(iMshRank) = numNodesMshRank(iMshRank) - numBoundaryNodesMshRank(iMshRank)
         !write(*,*) 'rank[',mpi_rank,']mshRankO',mshRankOrig,'numBN',numBoundaryNodesMshRank(iMshRank),'numDoF',numDoFMshRank(iMshRank),'numN',numNodesMshRank(iMshRank)

         allocate(boundaryNodes_jv%vector(iMshRank)%elems(numBoundaryNodesMshRank(iMshRank)))
         allocate(dofNodes_jv%vector(iMshRank)%elems(numDoFMshRank(iMshRank)))

         auxCnt1=0
         auxCnt2=0
         do iNodeL=1,numNodesMshRank(iMshRank)
            if(auxListBN(iNodeL)) then
               auxCnt1=auxCnt1+1
               boundaryNodes_jv%vector(iMshRank)%elems(auxCnt1) = iNodeL
            else
               auxCnt2=auxCnt2+1
               dofNodes_jv%vector(iMshRank)%elems(auxCnt2) = iNodeL
            end if
         end do

         deallocate(auxListBN)
      end do

      call MPI_Win_free(window_id,mpi_err)
      !---------------------------------------------------------------------------------------------------------
      ! END DOF SECTION!!!
      !---------------------------------------------------------------------------------------------------------

      deallocate(arrayMshBNInMpiRank_i8) !new line, be careful
      allocate(bndMpiBoundaryNodes_i8_jv%vector(numMshRanksInMpiRank))
      !generar el la jagged matrix amb els nou mpi boundary nodes, a partir d'aqui el proces es identic que abans
      !vaja, straightforward! :)

      do iMshRank=1,numMshRanksInMpiRank

         allocate(auxListBN(numBoundaryNodesMshRank(iMshRank)))
         auxListBN(:) = .false.
         auxCnt1=0
         do iNodeBnd=1,numBoundaryNodesMshRank(iMshRank)
            iNodeL = boundaryNodes_jv%vector(iMshRank)%elems(iNodeBnd)
            iNodeGSrl = globalIdSrl_i8_jv%vector(iMshRank)%elems(iNodeL)
            iPos = binarySearch_int_i8(mpiBoundaryNodes_i8_jv%vector(iMshRank)%elems,iNodeGSrl)
            if(iPos.ne.0) then
               auxCnt1=auxCnt1+1
               auxListBN(iNodeBnd) = .true.
            end if
         end do
         numBndMpiBoundaryNodes(iMshRank) = auxCnt1
         !write(*,*) 'mshRankO',mshRankOrig,'numBndMpi',numBndMpiBoundaryNodes(iMshRank)

         allocate(bndMpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(numBndMpiBoundaryNodes(iMshRank)))

         auxCnt1=0
         do iNodeBnd=1,numBoundaryNodesMshRank(iMshRank)
            if(auxListBN(iNodeBnd)) then
               iNodeL = boundaryNodes_jv%vector(iMshRank)%elems(iNodeBnd)
               iNodeGSrl = globalIdSrl_i8_jv%vector(iMshRank)%elems(iNodeL)
               auxCnt1=auxCnt1+1
               !if(mpi_rank.eq.3)write(*,*) 'rank[',mpi_rank,']aux',auxCnt1,'iNodeL',iNodeL,'iNodeGSrl',iNodeGSrl
               bndMpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(auxCnt1) = iNodeGSrl
            end if
         end do
         deallocate(auxListBN)

         call quicksort_array_int8(bndMpiBoundaryNodes_i8_jv%vector(iMshRank)%elems)
      end do

!--------------------------------------------------------------------------------------------

      vecNumBndMpiBN(:) = 0
      vecMemDispBndMpiBN(:) = 0
      numBndMpiBNInMpiRank = 0

      auxMemDispBndMpiBN=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         vecNumBndMpiBN(mshRank) = numBndMpiBoundaryNodes(iMshRank)
         numBndMpiBNInMpiRank = numBndMpiBNInMpiRank + numBndMpiBoundaryNodes(iMshRank)
         !if(iMshRank.ne.1) vecMemDispBndMpiBN(mshRank) = vecMemDispBndMpiBN(mshRank) + numBndMpiBoundaryNodes(iMshRank-1)
         if(iMshRank.ne.1) auxMemDispBndMpiBN = auxMemDispBndMpiBN + numBndMpiBoundaryNodes(iMshRank-1)
         vecMemDispBndMpiBN(mshRank) = auxMemDispBndMpiBN
      end do
      allocate(arrayBndMpiBNInMpiRank_i8(numBndMpiBNInMpiRank))

      do iMshRank=1,numMshRanksInMpiRank
         mshRank = mshRanksInMpiRank(iMshRank)
         do iAux=1,numBndMpiBoundaryNodes(iMshRank)
            memPos = vecMemDispBndMpiBN(mshRank) + iAux
            arrayBndMpiBNInMpiRank_i8(memPos) = bndMpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iAux)
         end do
      end do

      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*numMshRanks2Part

      call MPI_Win_create(vecNumBndMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecNumBndMpiBN(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      call MPI_Win_create(vecMemDispBndMpiBN,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mshRank = iMshRank-1
         mpiRank = mapMshRankToMpiRank(iMshRank)
         target_displacement = mshRank

         call MPI_Get(vecMemDispBndMpiBN(mshRank),1,mpi_datatype_int,mpiRank,target_displacement,1,mpi_datatype_int,window_id,mpi_err)
      end do

      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------------
      !if(mpi_rank.eq.0) write(*,*) 'vecNumBndMpiBN',vecNumBndMpiBN(:)
      !if(mpi_rank.eq.0) write(*,*) 'vecMemDispBndMpiBN',vecMemDispBndMpiBN(:)
      !---------------------------------------------------------------------------------------------------------
      allocate(bnd_numNodesToCommMshRank(numMshRanksInMpiRank))
      allocate(bnd_numMshRanksWithComms(numMshRanksInMpiRank))
      allocate(bnd_nodesToComm_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_ranksToComm_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_commsMemSize_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_commsMemPosInLoc_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_commsMemPosInNgb_jv%vector(numMshRanksInMpiRank))
      allocate(bnd_auxCommSchemeMemPos(numMshRanksInMpiRank,numMshRanks2Part))

      !---------------------------------------------------------------------------------------------------------------
      call MPI_Barrier(app_comm,mpi_err)
      win_buffer_size = mpi_int8_size*numBndMpiBNInMpiRank
      call MPI_Win_create(arrayBndMpiBNInMpiRank_i8,win_buffer_size,mpi_int8_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Barrier(app_comm,mpi_err)
      !---------------------------------------------------------------------------------------------------------------

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         numBndMpiNodeOrig = numBndMpiBoundaryNodes(iMshRank)

         bnd_numNodesToCommMshRank(iMshRank)=0
         bnd_auxCommSchemeNumNodes(:)=0
         bnd_auxCommSchemeMemPos(iMshRank,:)=0

         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispBndMpiBN(mshRankTrgt)
               memSize = vecNumBndMpiBN(mshRankTrgt)
               allocate(auxBndMpiBNMshRank_i8(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxBndMpiBNMshRank_i8,memSize,mpi_datatype_int8,mpiRankTrgt,target_displacement,memSize,mpi_datatype_int8,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numBndMpiNodeOrig
                  iNodeGSrl = bndMpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)
                  iPos = binarySearch_int_i8(auxBndMpiBNMshRank_i8,iNodeGSrl)
                  if(iPos.ne.0) then
                     bnd_numNodesToCommMshRank(iMshRank) = bnd_numNodesToCommMshRank(iMshRank) + 1
                     bnd_auxCommSchemeNumNodes(mshRankTrgt) = bnd_auxCommSchemeNumNodes(mshRankTrgt)+1
                  end if
               end do

               deallocate(auxBndMpiBNMshRank_i8)

            end if
         end do

         ! setting numRanksWithsComms
         bnd_numMshRanksWithComms(iMshRank) = 0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            if(bnd_auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               bnd_numMshRanksWithComms(iMshRank) = bnd_numMshRanksWithComms(iMshRank) + 1
            end if
         end do

         allocate(bnd_ranksToComm_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))
         allocate(bnd_commsMemSize_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))
         allocate(bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))
         allocate(bnd_commsMemPosInNgb_jv%vector(iMshRank)%elems(bnd_numMshRanksWithComms(iMshRank)))

         iAux=0
         memPos=1
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1
            bnd_auxCommSchemeMemPos(iMshRank,iMshRankTrgt)=memPos
            if(bnd_auxCommSchemeNumNodes(mshRankTrgt).ne.0) then
               iAux=iAux+1
               mshRankTrgt = iMshRankTrgt-1
               bnd_ranksToComm_jv%vector(iMshRank)%elems(iAux) = mshRankTrgt

               bnd_commsMemSize_jv%vector(iMshRank)%elems(iAux) = bnd_auxCommSchemeNumNodes(mshRankTrgt)
               bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems(iAux) = memPos
               memPos=memPos+bnd_auxCommSchemeNumNodes(mshRankTrgt)
            end if
         end do


         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_csNumNodes->',bnd_auxCommSchemeNumNodes(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_numNodesToComm',bnd_numNodesToCommMshRank(iMshRank)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_numRanksWC->',bnd_numMshRanksWithComms(iMshRank)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_ranksToComm->',bnd_ranksToComm_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_commsMemPosInLoc->',bnd_commsMemPosInLoc_jv%vector(iMshRank)%elems(:)
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'bnd_commsMemSize->',bnd_commsMemSize_jv%vector(iMshRank)%elems(:)

         allocate(bnd_nodesToComm_jv%vector(iMshRank)%elems(bnd_numNodesToCommMshRank(iMshRank)))
         iAux=0
         do iMshRankTrgt=1,numMshRanks2Part
            mshRankTrgt = iMshRankTrgt-1

            if(mshRankTrgt.ne.mshRankOrig) then
               mpiRankTrgt = mapMshRankToMpiRank(iMshRankTrgt)

               target_displacement = vecMemDispBndMpiBN(mshRankTrgt)
               memSize = vecNumBndMpiBN(mshRankTrgt)
               allocate(auxBndMpiBNMshRank_i8(memSize))

               call MPI_Win_lock(MPI_LOCK_SHARED,mpiRankTrgt,0,window_id,mpi_err)
               call MPI_Get(auxBndMpiBNMshRank_i8,memSize,mpi_datatype_int8,mpiRankTrgt,target_displacement,memSize,mpi_datatype_int8,window_id,mpi_err)
               call MPI_Win_unlock(mpiRankTrgt,window_id,mpi_err)

               do iNode=1,numBndMpiNodeOrig
                  iNodeGSrl = bndMpiBoundaryNodes_i8_jv%vector(iMshRank)%elems(iNode)
                  iPos = binarySearch_int_i8(auxBndMpiBNMshRank_i8,iNodeGSrl)
                  if(iPos.ne.0) then
                     iAux=iAux+1

                     iNodeLPos = binarySearch_int_i8(globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(:,1),iNodeGSrl)
                     iNodeL = globalIdSrlOrdered_i8_jm%matrix(iMshRank)%elems(iNodeLPos,2)

                     bnd_nodesToComm_jv%vector(iMshRank)%elems(iAux) = iNodeL
                  end if
               end do

               deallocate(auxBndMpiBNMshRank_i8)
            end if
         end do
      end do

      call MPI_Win_free(window_id,mpi_err)

      !---------------------------------------------------------------------------------------------------------

      allocate(bnd_auxCommSchemeMemPosAll(numMshRanks2Part**2))

      bnd_auxCommSchemeMemPosAll(:)=0
      do iMshRank=1,numMshRanksInMpiRank
         mshRank=mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,numMshRanks2Part
            iPos = mshRank*numMshRanks2Part+iMshRankTrgt
            bnd_auxCommSchemeMemPosAll(iPos) = bnd_auxCommSchemeMemPos(iMshRank,iMshRankTrgt)
         end do
      end do
      !write(*,*) '1.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      ! Create the window
      !--------------------------------------------------------------------------------------
      win_buffer_size = mpi_integer_size*(numMshRanks2Part**2)

      call MPI_Win_create(bnd_auxCommSchemeMemPosAll,win_buffer_size,mpi_integer_size,MPI_INFO_NULL,app_comm,window_id,mpi_err)
      call MPI_Win_fence(0,window_id,mpi_err)

      target_displacement=0
      do iMshRank=1,numMshRanks2Part
         mpiRank = mapMshRankToMpiRank(iMshRank)
         memSize = numMshRanks2Part
         memPos = (iMshRank-1)*numMshRanks2Part+1
         target_displacement = (iMshRank-1)*numMshRanks2Part

         call MPI_Get(bnd_auxCommSchemeMemPosAll(memPos),memSize,mpi_datatype_int,mpiRank,target_displacement,memSize,mpi_datatype_int,window_id,mpi_err)
      end do

      ! Wait for the MPI_Get issued to complete before going any further
      call MPI_Win_fence(0, window_id, mpi_err)
      call MPI_Win_free(window_id, mpi_err)
      !--------------------------------------------------------------------------------------
      !write(*,*) '2.[',mpi_rank,']auxCommSPAll',auxCommSchemeMemPosAll(:)

      do iMshRank=1,numMshRanksInMpiRank
         mshRankOrig = mshRanksInMpiRank(iMshRank)
         do iMshRankTrgt=1,bnd_numMshRanksWithComms(iMshRank)
            mshRankTrgt=bnd_ranksToComm_jv%vector(iMshRank)%elems(iMshRankTrgt)
            iPos = mshRankTrgt*numMshRanks2Part+mshRankOrig+1
            bnd_commsMemPosInNgb_jv%vector(iMshRank)%elems(iMshRankTrgt) = bnd_auxCommSchemeMemPosAll(iPos)
            !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'iMRT',iMshRankTrgt,'mRT',mshRankTrgt,'iPos',iPos
         end do
         !write(*,*) '[',mpi_rank,']mshRank',mshRankOrig,'commsMemPosInNgb->',commsMemPosInNgb_jv%vector(iMshRank)%elems(:)
      end do

      deallocate(bnd_auxCommSchemeMemPos)
      deallocate(bnd_auxCommSchemeMemPosAll)

   end subroutine generate_dof_and_boundary_mpi_comm_scheme_parallel

!------------------------------------------------------------------------------------------------------------------------------------

   subroutine set_nodesCoordinates(mnnode,mnpbou,mngaus,mshRank,isLinealOutput,numElemsInRank,numNodesInRank,globalIdSrlInRank_i8,listNodesInRank_i8,coordInRank,Ngp_l,connecParOrigMshRank,coordParMshRank,coordVTKMshRank)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou,mngaus
      logical,intent(in) :: isLinealOutput
      integer(4),intent(in) :: mshRank,numElemsInRank,numNodesInRank,connecParOrigMshRank(numElemsInRank,mnnode)
      integer(8),intent(in) :: globalIdSrlInRank_i8(numNodesInRank),listNodesInRank_i8(numNodesInRank)
      real(8),intent(in) :: coordInRank(numNodesInRank,3),Ngp_l(mngaus,mnnode)
      real(8),intent(out) :: coordParMshRank(numNodesInRank,3),coordVTKMshRank(numNodesInRank,3)
      integer(4) :: iPos,iNodeL
      integer(8) :: iNodeGsrl

      !if(mpi_rank.eq.0) write(*,*) ' # Creating nodes coordinates...'
      !------------------------------------------------------

      !- create the coordinate data for this process
      !!!$acc parallel loop
      do iNodeL=1,numNodesInRank
         iNodeGSrl=globalIdSrlInRank_i8(iNodeL)
         iPos = binarySearch_int_i8(listNodesInRank_i8,iNodeGSrl)
         !write (*,*) 'mshRank',mshRank,'iL', iNodeL,'iG', iNodeGSrl,'iPos',iPos,'[',coordInRank(iPos,1),']','[',coordInRank(iPos,2),']'
         coordParMshRank(iNodeL,1) = coordInRank(iPos,1)
         coordParMshRank(iNodeL,2) = coordInRank(iPos,2)
         coordParMshRank(iNodeL,3) = coordInRank(iPos,3)

         coordVTKMshRank(iNodeL,1) = coordInRank(iPos,1)
         coordVTKMshRank(iNodeL,2) = coordInRank(iPos,2)
         coordVTKMshRank(iNodeL,3) = coordInRank(iPos,3)
      end do
      !!!$acc end parallel loop

      call interpolateOriginalCoordinates(mnnode,mnpbou,mngaus,numElemsInRank,numNodesInRank,Ngp_l,connecParOrigMshRank,coordParMshRank)

      if(isLinealOutput) then
      !$acc kernels
      coordVTKMshRank(:,:) = coordParMshRank(:,:)
      !$acc end kernels   
      end if

   end subroutine set_nodesCoordinates

   subroutine interpolateOriginalCoordinates(mnnode,mnpbou,mngaus,numElemsInRank,numNodesInRank,Ngp_l,connecParOrigMshRank,coordParMshRank)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou,mngaus
      integer(4),intent(in) :: numElemsInRank,numNodesInRank,connecParOrigMshRank(numElemsInRank,mnnode)
      real(8),intent(in) :: Ngp_l(mngaus,mnnode)
      real(8),intent(inout) :: coordParMshRank(numNodesInRank,3)
      real(8), allocatable :: aux_1(:,:)
      integer(4) :: iElem,iNode,idime

      !if(mpi_rank.eq.0) write(*,*) "--| Interpolating nodes coordinates..."
      allocate(aux_1(numNodesInRank,ndime))
      aux_1(:,:) = coordParMshRank(:,:)
      do ielem = 1,numElemsInRank
         do inode = (2**ndime)+1,mnnode
            do idime = 1,ndime
               call var_interpolate_real8(mnnode,aux_1(connecParOrigMshRank(ielem,:),idime),Ngp_l(inode,:),coordParMshRank(connecParOrigMshRank(ielem,inode),idime))
            end do
         end do
      end do
      deallocate(aux_1)

   end subroutine interpolateOriginalCoordinates

   subroutine create_working_lists_parallel(mnnode,mnpbou,isPeriodic,mshRank,numElemsInRank,numNodesInRank,numBoundFacesInRank,connecParOrigInRank,&
            nPerInRank,masSlaInRank,connecParWorkInRank,connecOrigBoundFacesInRank,connecBoundFacesInRank,numWorkingNodesInRank,workingNodesInRank)
      implicit none
      integer(4),intent(in) :: mnnode,mnpbou
      logical,intent(in) :: isPeriodic
      integer,intent(in) :: mshRank,numElemsInRank,numNodesInRank,numBoundFacesInRank,connecParOrigInRank(numElemsInRank,mnnode),connecOrigBoundFacesInRank(numBoundFacesInRank,mnpbou)
      integer,intent(in) :: nPerInRank,masSlaInRank(nPerInRank,2)
      integer,intent(out) :: connecParWorkInRank(numElemsInRank,mnnode),numWorkingNodesInRank,connecBoundFacesInRank(numBoundFacesInRank,mnpbou)
      integer,intent(inout),allocatable :: workingNodesInRank(:)
      integer :: iNode,iNodeL,iElem,iBound,iPer,iAux,iPos
      integer :: iNodeL_Per,iNodeL_Per_Pair
      integer, allocatable :: aux_workingNodesInRank(:)

      !if(mpi_rank.eq.0) write(*,*) ' # Creating working lists...'

      !$acc kernels
      connecParWorkInRank(:,:)    = connecParOrigInRank(:,:)
      connecBoundFacesInRank(:,:) = connecOrigBoundFacesInRank(:,:)
      !$acc end kernels

      if(isPeriodic) then !do all stuff in case mesh is periodic
         !write(*,*) 'TO BE CHECK! subroutine create_working_list_parallel() for periodic cases'
         !----------------------------------------------------------------
         !-------------  CONNEC   ----------------------------------------

         do iElem = 1,numElemsInRank
            do iNode = 1,mnnode
               iNodeL = connecParWorkInRank(iElem,iNode)
               iPos = binarySearch_int_i4(masSlaInRank(:,2),iNodeL)
               if(iPos.ne.0) then !this node is slave, change!
                  iNodeL_Per_Pair = masSlaInRank(iPos,1)
                  connecParWorkInRank(iElem,iNode) = iNodeL_Per_Pair
               end if
            end do
         end do

         do iBound = 1,numBoundFacesInRank
            do iNode = 1,mnpbou
               iNodeL = connecBoundFacesInRank(iBound,iNode)
               iPos = binarySearch_int_i4(masSlaInRank(:,2),iNodeL)
               if(iPos.ne.0) then !this node is slave, change!
                  iNodeL_Per_Pair = masSlaInRank(iPos,1)
                  connecBoundFacesInRank(iBound,iNode) = iNodeL_Per_Pair
               end if
            end do
         end do

         !----------------------------------------------------------------
         !----------------WORKING NODES-----------------------------------

         numWorkingNodesInRank = numNodesInRank - nPerInRank

         allocate(aux_workingNodesInRank(numNodesInRank))
         allocate(workingNodesInRank(numWorkingNodesInRank))

         !$acc parallel loop
         do iNodeL = 1,numNodesInRank
            aux_workingNodesInRank(iNodeL) = iNodeL
         end do
         !$acc end parallel loop

         do iNodeL = 1,numNodesInRank
            iPos = binarySearch_int_i4(masSlaInRank(:,2),iNodeL)
            if(iPos.ne.0) then !this node is slave, change!
               aux_workingNodesInRank(iNodeL) = 0
            end if
         end do

         !TODO: GPU-PARALLELIZE THIS LOOP
         iAux = 0
         do iNodeL = 1,numNodesInRank
            if(aux_workingNodesInRank(iNodeL).ne.0) then
               iAux = iAux+1
               workingNodesInRank(iAux) = aux_workingNodesInRank(iNodeL)
               !write(*,*)'wNP(',iAux,')',workingNodesInRank(iAux)
            end if
         end do

         deallocate(aux_workingNodesInRank)

      else !non-periodic meshes
         numWorkingNodesInRank = numNodesInRank
         allocate(workingNodesInRank(numWorkingNodesInRank))

         !$acc parallel loop
         do iNodeL = 1,numWorkingNodesInRank
            workingNodesInRank(iNodeL) = iNodeL
         end do
         !$acc end parallel loop
      end if

   end subroutine create_working_lists_parallel

   subroutine evalShapeFunctions_Ngp_l(Ngp_l,dNgp,wgp,mporder,mnnode,mngaus,a2ijk)
      integer(4),intent(in) :: mporder,mnnode,mngaus
      integer(4),intent(in) :: a2ijk(mnnode)
      real(8),intent(out) :: Ngp_l(mngaus,mnnode), dNgp(ndime,mnnode,mngaus), wgp(mngaus)
      real(8) :: s, t, z
      integer(4) :: igaus
      real(8) :: xgp(mngaus,ndime)
      real(8) :: Ngp(mngaus,mnnode),dlxigp_ip(mngaus,ndime,mporder+1),dNgp_l(ndime,mnnode,mngaus)

      !*********************************************************

      call GaussLobattoLegendre_hex(mporder,mngaus,a2ijk,xgp,wgp)

      do igaus = 1,mngaus
         s = xgp(igaus,1)
         t = xgp(igaus,2)
         z = xgp(igaus,3)
         call hex_highorder(mporder,mnnode,s,t,z,a2ijk,Ngp(igaus,:),dNgp(:,:,igaus),Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
      end do

   end subroutine evalShapeFunctions_Ngp_l

   subroutine evalShapeFunctions_Ngp_l_2(Ngp_l,dNgp,wgp,mporder,mnnode,mngaus,a2ijk)
      integer(4),intent(in) :: mporder,mnnode,mngaus
      integer(4),intent(in) :: a2ijk(mnnode)
      real(8),intent(out) :: Ngp_l(mngaus,mnnode), dNgp(ndime,mnnode,mngaus), wgp(mngaus)
      real(8) :: s, t, z
      integer(4) :: igaus,orderGauss
      real(8) :: xgp(mngaus,ndime)
      real(8) :: Ngp(mngaus,mnnode),dlxigp_ip(mngaus,ndime,mporder+1),dNgp_l(ndime,mnnode,mngaus)

      !*********************************************************
      orderGauss = mngaus**(1.0_rp/3.0_rp)-1
      call GaussLobattoLegendre_hex_2(orderGauss,mngaus,xgp,wgp)

      do igaus = 1,mngaus
         s = xgp(igaus,1)
         t = xgp(igaus,2)
         z = xgp(igaus,3)
         call hex_highorder(mporder,mnnode,s,t,z,a2ijk,Ngp(igaus,:),dNgp(:,:,igaus),Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
      end do

   end subroutine evalShapeFunctions_Ngp_l_2

   subroutine get_listBoundaryElementsInRank_in_parallel(mporder,mnnode,mnpbou,gmsh2ijk,numElemsInRank,numNodesInRank,&
                              listNodesInRank_i8,connecInRank_i8,numBoundElemsRankPar,listBoundElemsRankPar)
      implicit none
      integer(4), intent(in) :: mporder,mnnode,mnpbou
      integer(4), intent(in) :: gmsh2ijk(mnnode)
      integer(4), intent(in) :: numElemsInRank,numNodesInRank!,numMasSlaNodesInRank
      integer(8), intent(in) :: listNodesInRank_i8(numNodesInRank),connecInRank_i8(numElemsInRank,mnnode)!,masSlaNodesInRank_i8(numMasSlaNodesInRank,2)

      integer(4), intent(out) :: numBoundElemsRankPar
      integer(4), allocatable, intent(inout) :: listBoundElemsRankPar(:)

      integer(4) :: nodeOwnedCnt(numNodesInRank)
      logical :: elemInBoundary(numElemsInRank)
      integer(4) :: ii,iElemL,indexNode
      integer(8) :: iNodeG
      integer(4),dimension(mnpbou) :: gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd
      !integer(4),parameter :: numQuadVert=4,numHexaVert=8
      !integer(4):: gmshQuadInnerVertInd(numQuadVert)
      !integer(4):: gmshHexVertInd(numHexaVert)
      !-------------------------------------------

      !$acc kernels
      nodeOwnedCnt(:)=0
      elemInBoundary(:)=.false.
      !$acc end kernels

      do iElemL=1,numElemsInRank
         do ii = 1,mnnode
            iNodeG = connecInRank_i8(iElemL,ii)
            indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
            nodeOwnedCnt(indexNode) = nodeOwnedCnt(indexNode) + 1
         end do
      end do

      !Ara ho faig boundaries i periodics junts, despres penso optimitzar
      !do ind=1,numMasSlaNodesInRank
      !   iNodeG = masSlaNodesInRank_i8(ind,1)
      !   indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
      !   nodeOwnedCnt(indexNode) = 2

      !   iNodeG = masSlaNodesInRank_i8(ind,2)
      !   indexNode = binarySearch_int_i8(listNodesInRank_i8,iNodeG)
      !   nodeOwnedCnt(indexNode) = 0
      !end do

      call get_gmshHexHOFacesIndex(mporder,mnpbou,gmshHexFaceFrontInd,gmshHexFaceBackInd,gmshHexFaceBottomInd,gmshHexFaceTopInd,gmshHexFaceLeftInd,gmshHexFaceRightInd)
      !call get_gmshHexHOVertIndex(mporder,gmshHexVertInd)
      !call get_gmshQuadHOInnerVertIndex(mporder,gmshQuadInnerVertInd)

      numBoundElemsRankPar = 0
      do iElemL=1,numElemsInRank
         call detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,gmshHexFaceFrontInd, listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
         call detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,gmshHexFaceLeftInd,  listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
         call detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,gmshHexFaceTopInd,   listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
         call detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,gmshHexFaceBackInd,  listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
         call detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,gmshHexFaceRightInd, listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
         call detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,gmshHexFaceBottomInd,listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
      end do

      do iElemL=1,numElemsInRank
         if(elemInBoundary(iElemL))  numBoundElemsRankPar = numBoundElemsRankPar + 1
      end do

      allocate(listBoundElemsRankPar(numBoundElemsRankPar))
      
      ii=0
      do iElemL=1,numElemsInRank
         if(elemInBoundary(iElemL)) then
           ii=ii+1
           listBoundElemsRankPar(ii) = iElemL
         end if
      end do

      !write(*,*) 'numBoundElemsRankPar',numBoundElemsRankPar,'listElems',listBoundElemsRankPar(:)

   end subroutine get_listBoundaryElementsInRank_in_parallel

   subroutine detect_is_elem_mpiOrMshBoundary(mporder,mnnode,mnpbou,numElemsInRank,numNodesInRank,iElemL,gmsh2ijk,face2ijk,listNodesInRank_i8,connecInRank_i8,nodeOwnedCnt,elemInBoundary)
      implicit none
      integer(4),intent(in) :: mporder,mnnode,mnpbou
      integer(4),intent(in) :: numElemsInRank,numNodesInRank
      integer(4),intent(in) :: iElemL,gmsh2ijk(mnnode),face2ijk(mnpbou),nodeOwnedCnt(numNodesInRank)
      integer(8),intent(in) :: listNodesInRank_i8(numNodesInRank),connecInRank_i8(numElemsInRank,mnnode)
      logical,intent(inout) :: elemInBoundary(numElemsInRank)
      integer(4) :: ind,indexNode_inFace
      integer(8) :: iNodeG_inFace
      logical :: isMshBound

      if(elemInBoundary(iElemL) .eqv. .true.) return !if already in boundary, go out!

      !# 1.Whatever face ------------------------------------------------------
      ind = face2ijk(mporder+3)
      iNodeG_inFace = connecInRank_i8(iElemL,ind)
      indexNode_inFace = binarySearch_int_i8(listNodesInRank_i8,iNodeG_inFace)
      if(nodeOwnedCnt(indexNode_inFace).eq.1) then !this face is boundary
         elemInBoundary(iElemL) = .true.
      end if

   end subroutine detect_is_elem_mpiOrMshBoundary



!---------------------------------------------------------------------------------------------------------------------------
!-------------------------DEPRECATED CODE SECTION---------------------------------------------------------------------------
#if 0

#if 0
BEGIN WHEN READING LINE BY LINE.. FUCKING SLOW
      !----------------------------------------------------------------------------------------------------------
      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      !call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write

      do iElem=1,numElemsInRank
         elemToRead=listElemsRank(iElem)
         ms_offset(2) = elemToRead-1
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         !call h5dread_f(dset_id,dtype,connecRank(iElem,:),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5dread_f(dset_id,dtype,connecRank(iElem,:),ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id)

         !connecRank(iElem,1:nnode) = auxConnec(:)
         !if(mpi_rank.eq.0)write(*,*) 'iElem',iElem,'e2r',elemToRead,'auxConnec',auxConnec(:)
         !if(mpi_rank.eq.0)write(*,*) 'connecRank',connecRank(iElem,:)
      end do

      !call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
      !----------------------------------------------------------------------------------------------------------
END WHEN READING LINE BY LINE.. FUCKING SLOW
#endif

#if 0
BEGIN WHEN READING LINE BY LINE.. FUCKING SLOW
      call h5dopen_f(gmsh_h5_fileId,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)!get filespace of the dataset
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)!get dimensions of the filespace
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err) ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      !call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err) ! Create property list for collective dataset write

      do iNode=1,numNodesInRank
         nodeToRead=listNodesRank(iNode)
         ms_offset(2) = nodeToRead-1
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
         !call h5dread_f(dset_id,dtype,auxCoords,fs_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)
         call h5dread_f(dset_id,dtype,auxCoords,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id)
         coordNodesRank(iNode,:) = auxCoords(:)
         !if(mpi_rank.eq.0)write(*,*) 'n2r',nodeToRead,'coords',coordNodesRank(iNode,:)
      end do

      !call h5pclose_f(plist_id,h5err)
      call h5sclose_f(mspace_id,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)
END WHEN READING LINE BY LINE.. FUCKING SLOW
#endif



   subroutine read_dims_file(file_path,file_name,numElems,numNodes,numBounds)
      implicit none
      character(500), intent(in) :: file_path, file_name
      integer(4), intent(out) :: numElems,numNodes,numBounds
      character(500) :: file_type, line

      write(file_type,*) ".dims.dat"

      if(mpi_rank.eq.0)write(*,*) "--| Reading dims.dat file..."
      open(unit=fileId_dims,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
      read(fileId_dims,*) line, numNodes
      if(mpi_rank.eq.0)write(*,*) "--| Nodes on mesh: ",numNodes
      read(fileId_dims,*) line, numElems
      if(mpi_rank.eq.0)write(*,*) "--| Elemensts on mesh: ",numElems
      read(fileId_dims,*) line, numBounds
      if(mpi_rank.eq.0)write(*,*) "--| Boundary faces on mesh : ",numBounds
      if(mpi_rank.eq.0)write(*,*) "--| End of dims.dat file!"
      close(unit=fileId_dims)

   end subroutine read_dims_file
!---------------------------------------------------------------------------------------------------------------------------
   subroutine open_geo_dat_file(file_path,file_name)
      implicit none
      character(500), intent(in) :: file_path, file_name
      character(128) :: file_type_geo

      write(file_type_geo,*) ".geo.dat"

      if(mpi_rank.eq.0) write(*,*) "--| Opening geo.dat file..."
      open(unit=fileId_geo,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type_geo)),status="old")
   end subroutine open_geo_dat_file

   subroutine close_geo_dat_file()
      implicit none

      if(mpi_rank.eq.0) write(*,*) "--| Closing geo.dat file"
      close(unit=fileId_geo)
   end subroutine close_geo_dat_file
!---------------------------------------------------------------------------------------------------------------------------
   subroutine open_fix_bou_file(file_path,file_name)
      implicit none
      character(500), intent(in) :: file_path, file_name
      character(128) :: file_type_bou

      write(file_type_bou,*) ".fix.bou"

      if(mpi_rank.eq.0) write(*,*) "--| Opening fix.bou file..."
      open(unit=fileId_bou,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type_bou)),status="old")
   end subroutine open_fix_bou_file

   subroutine close_fix_bou_file()
      implicit none

      if(mpi_rank.eq.0) write(*,*) "--| Closing fix.bou file"
      close(unit=fileId_bou)
   end subroutine close_fix_bou_file
!---------------------------------------------------------------------------------------------------------------------------

   subroutine read_element_section_geo_file(numElemsSrl)
      implicit none
      integer(4),intent(in) :: numElemsSrl
      integer(4) :: iline

      ! Nodes/element section, blank for now
      !------------------------------------------------------------------------------------------
      read(fileId_geo,*) ! Section header
      do iline = 1,numElemsSrl
         read(fileId_geo,*)
      end do
      read(fileId_geo,*) ! Section ender
      !------------------------------------------------------------------------------------------
   end subroutine read_element_section_geo_file

   subroutine read_geo_dat_file_in_parallel(file_path,file_name,isPeriodic,numElemsSrl,numNodesSrl,numBoundFacesSrl,&
                                            numElemsMpiRank,listElemsMpiRank,&
                                            numNodesMpiRank,listNodesMpiRank,numBoundsMpiRank,listElemsBoundsMpiRank,&
                                            boundFacesMpiRank,connecMpiRank,coordNodesMpiRank)
      implicit none
      character(500), intent(in) :: file_path, file_name
      logical, intent(in)        :: isPeriodic
      integer(4), intent(in)     :: numElemsSrl,numNodesSrl,numBoundFacesSrl
      integer(4), intent(inout)  :: numElemsMpiRank,numNodesMpiRank,numBoundsMpiRank
      integer, allocatable, intent(inout) :: listElemsMpiRank(:),listNodesMpiRank(:),listElemsBoundsMpiRank(:,:),connecMpiRank(:,:),boundFacesMpiRank(:,:)
      real(rp), allocatable, intent(inout) :: coordNodesMpiRank(:,:)
      integer(4)                 :: iElemMpiRankStart,iElemMpiRankEnd
      integer(4)                 :: auxCnt,iAux,jAux
      integer(4)                 :: nodeCnt
      integer(4)                 :: iline,iNode,iBound,iBoundBou,iElem,ind_gmsh,iNodeG,iElemG
      integer(4)                 :: iNodeG_inFace,iNodeG_inElem
      integer(4)                 :: iBoundNodes(npbou),aux(nnode+1),bou_aux(npbou+1),bou_code_aux
      character(2000)            :: line
      real(8)                    :: x,y,z
      logical                    :: vertexFound
      integer(4), allocatable    :: auxBoundFacesMpiRank(:,:)
      integer                    :: f_iNodeG(4),e_iNodeG(8) !for the vertex of squares of the face and the element

      !----------------------------------------------------------------------------------------------

      call do_element_partitioning_serial(numElemsSrl,iElemMpiRankStart,iElemMpiRankEnd,numElemsMpiRank)

      !----------------------------------------------------------------------------------------------

      ! Connectivity table section
      !------------------------------------------------------------------------------------------
      allocate(connecMpiRank(numElemsMpiRank,nnode))
      allocate(listElemsMpiRank(numElemsMpiRank))

      !$acc kernels
      connecMpiRank(:,:) = 0
      listElemsMpiRank(:) = 0
      !$acc end kernels

      auxCnt=0
      do iElem = iElemMpiRankStart,iElemMpiRankEnd
         auxCnt=auxCnt+1
         listElemsMpiRank(auxCnt) = iElem
      end do

      call open_geo_dat_file(file_path,file_name)
      call read_element_section_geo_file(numElemsSrl)
      call read_elem_connec_and_nodes_coords_from_geo_file(numElemsSrl,numNodesSrl,numElemsMpiRank,listElemsMpiRank,numNodesMpiRank,connecMpiRank,listNodesMpiRank,coordNodesMpiRank)

      ! Boundary nodes section
      !------------------------------------------------------------------------------------------

      if(isPeriodic) then
         if(mpi_rank.eq.0)write(*,*) "--| Shit for periodic...."
         write(*,*) 'ABORTING, SHIT FOR PERIODIC YET TO BE IMPLEMENTED!'
         call MPI_Abort(app_comm, -1, mpi_err)
         !call open_fix_bou_file(file_path,file_name)
         !call read_boundaries(numBoundFacesSrl,numElemsMpiRank,listElemsMpiRank,connecMpiRank,numBoundsMpiRank,listElemsBoundsMpiRank,boundFacesMpiRank)

         !call close_fix_bou_file()
      end if

      call close_geo_dat_file()

      !------------------------------------------------------------------------------------------

      if(mpi_rank.eq.0)write(*,*) "--| End of reading GEO.DAT file in parallel!"
      if(mpi_rank.eq.0)write(*,*) "--------------------------------------------"

   end subroutine read_geo_dat_file_in_parallel

   subroutine read_elem_connec_and_nodes_coords_from_geo_file(numElemsSrl,numNodesSrl,numElemsInRank,listElemsRank,numNodesInRank,connecRank,listNodesRank,coordNodesRank)
      implicit none
      integer,intent(in) :: numElemsSrl,numNodesSrl,numElemsInRank,listElemsRank(numElemsInRank)
      integer,intent(out) :: numNodesInRank,connecRank(numElemsInRank,nnode)
      integer,allocatable,intent(inout) :: listNodesRank(:)
      real(rp),allocatable,intent(inout) :: coordNodesRank(:,:)
      integer(4) :: auxConnec(nnode+1)
      integer(4) :: iLine,iElem,iNode,iNodeG,elemCnt,nodeCnt,iAux,jAux
      real(8) :: x,y,z
      character(2000) :: line
      integer(4) :: nextElemToRead,nextNodeToRead,prevNodeG

      integer,allocatable :: rawNodeListRank(:)

      if(mpi_rank.eq.0)write(*,*) "--| Reading element table..."

      !------------------------------------------------------------------------------------------
      read(fileId_geo,*) ! Section header
      elemCnt=0
      nextElemToRead=listElemsRank(elemCnt+1)
      do iElem=1,numElemsSrl
         !write(*,*) '(',mpi_rank,')iElem',iElem
         read(fileId_geo,'(a)') line
         if(iElem.eq.nextElemToRead) then
            read(line,*) (auxConnec(iNode), iNode=1,nnode+1)
            elemCnt=elemCnt+1

            !write(*,*) 'connecRank[',mpi_rank,'](',elemCnt,')'
            connecRank(elemCnt,1:nnode) = auxConnec(2:nnode+1)

            if(elemCnt.lt.numElemsInRank) nextElemToRead=listElemsRank(elemCnt+1)
            !if(mpi_rank.eq.0) write(*,*) iElem,connecRank(elemCnt,:)
         end if
      end do
      read(fileId_geo,*) ! Section ender

      !------------------------------------------------------------------------------------------

      allocate(rawNodeListRank(numElemsInRank*nnode))
      !$acc kernels
      rawNodeListRank(:) = 0
      !$acc end kernels

      ! add all the iNodeGs of this mpirank
      nodeCnt=0
      do iAux = 1,numElemsInRank
         do jAux = 1,nnode
            iNodeG = connecRank(iAux,jAux)
            nodeCnt=nodeCnt+1
            rawNodeListRank(nodeCnt) = iNodeG
         end do
      end do

      !sorting the nodes id
      call quicksort_array_int(rawNodeListRank)

      prevNodeG=0
      nodeCnt=0
      do iAux = 1,(numElemsInRank*nnode)
         if((rawNodeListRank(iAux)).ne.(prevNodeG)) then
            nodeCnt=nodeCnt+1
            prevNodeG = rawNodeListRank(iAux)
         end if
      end do
      numNodesInRank = nodeCnt
      !write(*,*) 'newMethod.[',mpi_rank,'] numNodesInRank ',numNodesInRank

      allocate(listNodesRank(numNodesInRank))
      allocate(coordNodesRank(numNodesInRank,ndime))

      !$acc kernels
      coordNodesRank(:,:) = 0
      listNodesRank(:) = 0
      !$acc end kernels

      prevNodeG=0
      nodeCnt=0
      do iAux = 1,(numElemsInRank*nnode)
         if((rawNodeListRank(iAux)).ne.(prevNodeG)) then
            nodeCnt=nodeCnt+1
            listNodesRank(nodeCnt) = rawNodeListRank(iAux)
            prevNodeG = rawNodeListRank(iAux)
         end if
      end do

      ! Nodal coordinates section
      !------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0)write(*,*) "--| Reading coordinates..."
      if (ndime .ne. 3) then
         write(*,*) 'Crashing in read_elem_connec_and_nodes_coords_from_geo_file! ndime is not 3! EXIT!'
         call MPI_Abort(app_comm,-1,mpi_err)
      end if

      nodeCnt=0
      nextNodeToRead=listNodesRank(nodeCnt+1)
      read(fileId_geo,*) ! Section header
      do iLine =1,numNodesSrl
         read(fileId_geo,*) iNode,x,y,z !iLine is equal to iNode (or should be!)
         !write(*,*) 'iNode',iNode,x,y,z
         if(iNode.eq.nextNodeToRead) then
            nodeCnt=nodeCnt+1
            coordNodesRank(nodeCnt,1)=real(x,rp)
            coordNodesRank(nodeCnt,2)=real(y,rp)
            coordNodesRank(nodeCnt,3)=real(z,rp)

            if(nodeCnt.lt.numNodesInRank) nextNodeToRead=listNodesRank(nodeCnt+1)
         end if
      end do
      read(fileId_geo,*) ! Section ender

      deallocate(rawNodeListRank)

   end subroutine read_elem_connec_and_nodes_coords_from_geo_file

   subroutine read_boundaries(numBoundFacesSrl,numElemsInRank,listElemsInRank,connecInRank,numBoundsInRank,maxBoundCodeInRank,listElemsBoundsInRank,listBoundFacesInRank,boundFacesCodesInRank,boundFacesInRank)
      implicit none
      integer,intent(in) :: numBoundFacesSrl,numElemsInRank,listElemsInRank(numElemsInRank),connecInRank(numElemsInRank,nnode)
      integer(4), intent(out) :: numBoundsInRank,maxBoundCodeInRank
      integer, allocatable, intent(inout) :: listElemsBoundsInRank(:,:),listBoundFacesInRank(:),boundFacesCodesInRank(:),boundFacesInRank(:,:)
      integer(4) :: iBoundNodes(npbou),bou_aux(npbou+1),bou_code_aux
      character(2000) :: line
      integer(4),allocatable :: auxBoundFacesInRank(:,:)
      integer(4) :: iLine,iNode,iBound,iBoundBou,iElem,ind_gmsh,iElemG,iAux,jAux
      integer(4) :: iNodeG_inFace,iNodeG_inElem,nodeCnt
      integer :: f_iNodeG(4),e_iNodeG(8) !for the vertex of squares of the face and the element
      logical :: vertexFound

      read(fileId_geo,*) line! Section header
      read(fileId_bou,*) line! Section header

      !for the moment I will use this auxiliar array, but if we see that is allocating too much memory
      !maybe we can try another strategy, without allocating an array of numBoundFaceSrl
      !but doing a double loop
      !1st loop: determine the bounds that are in the msh rank
      !2nd loop: store only the nodes in the rank
      allocate(auxBoundFacesInRank(numBoundFacesSrl,2+npbou))
      auxBoundFacesInRank(:,:) = 0
      allocate(listElemsBoundsInRank(numElemsInRank,maxBoundsPerElem))
      listElemsBoundsInRank(:,:) = 0

      numBoundsInRank=0
      do iLine = 1,numBoundFacesSrl
         !reading boundary nodes-------------------------------
         read(fileId_geo,'(a)') line
         read(line,*) (bou_aux(iNode), iNode=1,npbou+1)
         iBound = bou_aux(1) !iLine & iBound must be the same
         iBoundNodes = bou_aux(2:npbou+1)
         !-----------------------------------------------------

         !reading boundary code--------------------------------
         read(fileId_bou,'(a)') line
         read(line,*) iBoundBou,bou_code_aux
         !-----------------------------------------------------

         if(iBound.ne.iBoundBou) then
            write(*,*) 'Error in read_boundaries: iBound',iBound,' not equal to iBoundBou',iBoundBou
            call MPI_Abort(app_comm, -1, mpi_err)
         end if
         !write(*,*) 'iBound',iBound,'iBoundBou',iBoundBou,'bouCode',bou_code_aux,'iBoundNodes',iBoundNodes

         !fill the corners of the face to check
         do iAux=1,4
            ind_gmsh = gmsh2ij(posFaceVertices(iAux))
            f_iNodeG(iAux) = iBoundNodes(ind_gmsh)
         end do
         !write(*,*) '[',mpi_rank,']iBound',iBound,' -> f_iNodeG ',f_iNodeG(:)

         !Search to which element this face belongs
         do iElem=1,numElemsInRank
            nodeCnt=0
            iElemG=listElemsInRank(iElem)
            !fill the corners of the element
            do iAux=1,8
               ind_gmsh = gmsh2ijk(posElemVertices(iAux))
               e_iNodeG(iAux) = connecInRank(iElem,ind_gmsh)
            end do

            !isBoundInElem=.false.
            fLoop: do iAux=1,4
               vertexFound=.false.
               iNodeG_inFace = f_iNodeG(iAux)
               eLoop: do jAux=1,8
                  iNodeG_inElem = e_iNodeG(jAux)
                  if(iNodeG_inFace .eq. iNodeG_inElem) then
                     nodeCnt=nodeCnt+1
                     vertexFound=.true.
                     exit eLoop
                  end if
               end do eLoop
               if(.not.(vertexFound)) exit fLoop
            end do fLoop
            if(nodeCnt.ge.4) then
               !isBoundInElem=.true.
               numBoundsInRank=numBoundsInRank+1
               !write(*,*) '[',mpi_rank,']iBound',iBound,' -> elem ',iElemG
               auxBoundFacesInRank(numBoundsInRank,1)=iBound
               auxBoundFacesInRank(numBoundsInRank,2)=bou_code_aux
               auxBoundFacesInRank(numBoundsInRank,3:npbou+2)=iBoundNodes(1:npbou)

               addLoop : do iAux = 1,maxBoundsPerElem
                  if(listElemsBoundsInRank(iElem,iAux).eq.0) then
                     !write(*,*) '[',mpi_rank,'] adding iBound',iBound,'to elem',iElemG
                     listElemsBoundsInRank(iElem,iAux) = iBound

                     exit addLoop
                  end if
               end do addLoop
            end if

         end do
      end do

      read(fileId_geo,*) ! Section ender
      read(fileId_bou,*) ! Section ender

      !write(*,*) '[',mpi_rank,'] numBoundsInRank',numBoundsInRank
      !do iElem=1,numElemsInRank
      !   if(listElemsBoundsInRank(iElem,1).ne.0) then
      !      write(*,*) 'i',iElem,'iG',listElemsInRank(iElem),'bounds',listElemsBoundsInRank(iElem,:)
      !   end if
      !end do

      allocate(listBoundFacesInRank(numBoundsInRank))
      allocate(boundFacesCodesInRank(numBoundsInRank))
      allocate(boundFacesInRank(numBoundsInRank,npbou))

      maxBoundCodeInRank = 0
      do iBound=1,numBoundsInRank
         listBoundFacesInRank(iBound)  = auxBoundFacesInRank(iBound,1)
         boundFacesCodesInRank(iBound) = auxBoundFacesInRank(iBound,2)
         boundFacesInRank(iBound,:)    = auxBoundFacesInRank(iBound,3:npbou+2)

         maxBoundCodeInRank = max(maxBoundCodeInRank,boundFacesCodesInRank(iBound))
         !write(*,*) '[',mpi_rank,']bfmpirank(',iBound,')',boundFacesInRank(iBound,:)
      end do

      deallocate(auxBoundFacesInRank)

   end subroutine read_boundaries

#endif
!---------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------


end module mod_meshConversorTool
