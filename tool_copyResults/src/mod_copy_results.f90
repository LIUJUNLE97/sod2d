module mod_copy_results
   use mod_mpi
   use mod_ijk_indices
   use mod_mpi_mesh
   use mod_comms
   use mod_hdf5
   use mod_partition_utils

   implicit none

contains

   subroutine copy_results_same_mesh_Npartitions(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName,target_Nprocs,&
                                                results_first,results_last,results_step)
      implicit none
      character(len=*), intent(in)  :: mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName
      integer(4),intent(in) :: target_Nprocs,results_first,results_last,results_step
      character(512) :: source_meshFile_h5_full_name,target_meshFile_h5_full_name
      character(512) :: source_base_resultsFile_h5,target_base_resultsFile_h5
      character(512) :: source_full_resultsFile_h5,target_full_resultsFile_h5
      integer(4) :: mporder,mnnode,mngaus,mnpbou
      integer(4) :: res_inst,targetRank,target_numNodesRankPar
      integer(4) :: iNodeTrgt,iNodeMpi,dsetScaCnt,dsetVecCnt
      logical :: evalMeshQuality=.false.
      integer(8),dimension(0:target_Nprocs-1) :: iNodeStartPar_i8

      integer(4) :: iTrgtRank,trgtRank
      integer(4) :: trgtRankInMpiRankStart,trgtRankInMpiRankEnd,numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),allocatable :: trgtRanksInMpiRank(:),mapTrgtRankToMpiRank(:)

      integer(8) :: numNodesTrgtTotal_i8
      integer(8),allocatable :: trgtRankNodeStart_i8(:),trgtRankNodeEnd_i8(:)
      integer(4),allocatable :: mpiRankElemStart(:),mpiRankElemEnd(:),trgtRankElemStart(:),trgtRankElemEnd(:)
      integer(4),allocatable :: numElemsTrgtRank(:),numNodesTrgtRank(:)
      integer(4),allocatable :: numElemsVTKTrgtRank(:),sizeConnecVTKTrgtRank(:)
      type(jagged_vector_int4) :: listNodesTrgtRankMpiId_jv,mapNodeIdTrgtToMpi_jv,connecVTKTrgtRank_jv
      type(jagged_matrix_int4) :: connecTrgtRankMpiId_jm,connecParOrigTrgtRank_jm
      type(jagged_matrix_rp) :: coordTrgtRank_jm,coordVTKTrgtRank_jm
      type(jagged_vector_rp) :: quality_jv

      character(len=256),dimension(max_num_saved_fields) :: datasetsScalarFields,datasetsVectorFields
      logical,dimension(max_num_saved_fields) :: datasetsScaExist,datasetsVecExist


      integer(hid_t) :: targetRes_hdf5_file_id,sourceRes_hdf5_file_id,targetMesh_hdf5_file_id
      integer(hsize_t),dimension(1) :: ms_dims
      integer(hssize_t),dimension(1) :: ms_offset
      integer(4),allocatable :: aux_array_i4(:)
      integer(4) :: target_vecNumNodesRankPar(target_Nprocs)

      !----------------------------------------------------------------------------------------------
      call init_hdf5_interface()
      !----------------------------------------------------------------------------------------------

      call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName,mpi_size,source_meshFile_h5_full_name)
      call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName,target_Nprocs,target_meshFile_h5_full_name)
      call set_hdf5_resultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,mpi_size,source_base_resultsFile_h5)
      call set_hdf5_resultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,target_Nprocs,target_base_resultsFile_h5)

      if(mpi_rank.eq.0) then
         write(*,*) 'Source mesh: ',trim(source_meshFile_h5_full_name)
         write(*,*) 'Target mesh: ',trim(target_meshFile_h5_full_name)
      end if

      call get_mesh_porder_from_hdf5(source_meshFile_h5_full_name,mporder)
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)

      !Loading original mesh file
      call load_hdf5_meshFile(source_meshFile_h5_full_name)

      if(mpi_rank.eq.0) then
         write(*,*) '# Generating new mesh file:',trim(target_meshFile_h5_full_name)
      end if

      !Distributing the target_Nprocs in the current Nprocs
      call distribute_ranks2Part_in_mpiRank(target_Nprocs,trgtRankInMpiRankStart,trgtRankInMpiRankEnd,numTrgtRanksInMpiRank,maxNumTrgtRanks,trgtRanksInMpiRank,mapTrgtRankToMpiRank)


      allocate(   numElemsTrgtRank(numTrgtRanksInMpiRank))
      allocate(  trgtRankElemStart(numTrgtRanksInMpiRank))
      allocate(    trgtRankElemEnd(numTrgtRanksInMpiRank))
      allocate(   mpiRankElemStart(numTrgtRanksInMpiRank))
      allocate(     mpiRankElemEnd(numTrgtRanksInMpiRank))
      allocate(   numNodesTrgtRank(numTrgtRanksInMpiRank))

      allocate(listNodesTrgtRankMpiId_jv%vector(numTrgtRanksInMpiRank))
      allocate(   connecTrgtRankMpiId_jm%matrix(numTrgtRanksInMpiRank))
      allocate(         coordTrgtRank_jm%matrix(numTrgtRanksInMpiRank))
      allocate(               quality_jv%vector(numTrgtRanksInMpiRank))


      if(mpi_rank.eq.0) write(*,*) ' 1. Generating new partitioning and connec...'

      do iTrgtRank=1,numTrgtRanksInMpiRank
         trgtRank = trgtRanksInMpiRank(iTrgtRank)

         call do_element_partitioning_serial(numElemsRankPar,mpiRankElemStart(iTrgtRank),mpiRankElemEnd(iTrgtRank),&
                                             numElemsTrgtRank(iTrgtRank),(iTrgtRank-1),numTrgtRanksInMpiRank)


         call generate_connec_and_coords_for_targetRank(mnnode,trgtRank,&
                  mpiRankElemStart(iTrgtRank),mpiRankElemEnd(iTrgtRank),numElemsTrgtRank(iTrgtRank),&
                  numNodesTrgtRank(iTrgtRank),connecTrgtRankMpiId_jm%matrix(iTrgtRank)%elems,&
                  listNodesTrgtRankMpiId_jv%vector(iTrgtRank)%elems,coordTrgtRank_jm%matrix(iTrgtRank)%elems)

         !modify with the rankElemStart as offset
         trgtRankElemStart(iTrgtRank) = mpiRankElemStart(iTrgtRank) + (rankElemStart-1)
         trgtRankElemEnd(iTrgtRank)   = mpiRankElemEnd(iTrgtRank)   + (rankElemStart-1)

      end do
      
      allocate(trgtRankNodeStart_i8(numTrgtRanksInMpiRank))
      allocate(trgtRankNodeEnd_i8(numTrgtRanksInMpiRank))

      trgtRankNodeStart_i8(:) = 0
      trgtRankNodeEnd_i8(:) = 0

      call define_parallelNodePartitioning(target_Nprocs,numTrgtRanksInMpiRank,numNodesTrgtRank,trgtRanksInMpiRank,mapTrgtRankToMpiRank,&
                                          trgtRankNodeStart_i8,trgtRankNodeEnd_i8,iNodeStartPar_i8,numNodesTrgtTotal_i8)

      allocate(mapNodeIdTrgtToMpi_jv%vector(numTrgtRanksInMpiRank))
      allocate(coordVTKTrgtRank_jm%matrix(numTrgtRanksInMpiRank))
      allocate(connecVTKTrgtRank_jv%vector(numTrgtRanksInMpiRank))
      allocate(connecParOrigTrgtRank_jm%matrix(numTrgtRanksInMpiRank))

      allocate(numElemsVTKTrgtRank(numTrgtRanksInMpiRank))
      allocate(sizeConnecVTKTrgtRank(numTrgtRanksInMpiRank))

      if(mpi_rank.eq.0) write(*,*) ' 2. Generating new node ordering and mapping between meshes...'

      do iTrgtRank=1,numTrgtRanksInMpiRank
         trgtRank = trgtRanksInMpiRank(iTrgtRank)

         numElemsVTKTrgtRank(iTrgtRank)   = numElemsTrgtRank(iTrgtRank)*mesh_numVTKElemsPerMshElem
         sizeConnecVTKTrgtRank(iTrgtRank) = numElemsVTKTrgtRank(iTrgtRank)*mesh_VTKnnode

         !write(*,*) 'numElemsVTK(',trgtRank,')',numElemsVTKTrgtRank(iTrgtRank)
         !write(*,*) 'sizeConnec(',trgtRank,')',sizeConnecVTKTrgtRank(iTrgtRank)

         allocate(mapNodeIdTrgtToMpi_jv%vector(iTrgtRank)%elems(numNodesTrgtRank(iTrgtRank)))
         allocate(coordVTKTrgtRank_jm%matrix(iTrgtRank)%elems(numNodesTrgtRank(iTrgtRank),3))
         allocate(connecVTKTrgtRank_jv%vector(iTrgtRank)%elems(sizeConnecVTKTrgtRank(iTrgtRank)))
         allocate(connecParOrigTrgtRank_jm%matrix(iTrgtRank)%elems(numElemsTrgtRank(iTrgtRank),mnnode))


         call reorder_nodes_in_trgtRank(mporder,mnnode,trgtRank,target_Nprocs,&
               numElemsTrgtRank(iTrgtRank),numNodesTrgtRank(iTrgtRank),sizeConnecVTKTrgtRank(iTrgtRank),isMeshLinealOutput,&
               listNodesTrgtRankMpiId_jv%vector(iTrgtRank)%elems,connecTrgtRankMpiId_jm%matrix(iTrgtRank)%elems,&
               iNodeStartPar_i8(iTrgtRank),mapNodeIdTrgtToMpi_jv%vector(iTrgtRank)%elems,&
               connecVTKTrgtRank_jv%vector(iTrgtRank)%elems,connecParOrigTrgtRank_jm%matrix(iTrgtRank)%elems,&
               coordVTKTrgtRank_jm%matrix(iTrgtRank)%elems,coordTrgtRank_jm%matrix(iTrgtRank)%elems)

      end do


      !Writing new results file in target N procs
      !---------------------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) ' 3. Generating new mesh partitioned in',target_Nprocs,'N procs'

      !call create_hdf5_file(target_full_resultsFile_h5,targetRes_hdf5_file_id)
      call create_hdf5_file(target_meshFile_h5_full_name,targetMesh_hdf5_file_id)
      
      call create_groups_datasets_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,targetMesh_hdf5_file_id,isMeshLinealOutput,.false.,&
                                                target_Nprocs,totalNumElements,numNodesTrgtTotal_i8,mesh_VTKnnode,mesh_numVTKElemsPerMshElem)

      do iTrgtRank=1,numTrgtRanksInMpiRank
         trgtRank = trgtRanksInMpiRank(iTrgtRank)

        call write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,targetMesh_hdf5_file_id,evalMeshQuality,trgtRank,target_Nprocs,&
            numElemsTrgtRank(iTrgtRank),numElemsVTKTrgtRank(iTrgtRank),sizeConnecVTKTrgtRank(iTrgtRank),mesh_VTKnnode,mesh_numVTKElemsPerMshElem,&
            trgtRankElemStart(iTrgtRank),trgtRankElemEnd(iTrgtRank),trgtRankNodeStart_i8(iTrgtRank),trgtRankNodeEnd_i8(iTrgtRank),numNodesTrgtRank(iTrgtRank),&
            coordVTKTrgtRank_jm%matrix(iTrgtRank)%elems,connecVTKTrgtRank_jv%vector(iTrgtRank)%elems,quality_jv%vector(iTrgtRank)%elems)
      end do

      do iTrgtRank=(numTrgtRanksInMpiRank+1),maxNumTrgtRanks
        call dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(targetMesh_hdf5_file_id,evalMeshQuality,target_Nprocs)
      end do

      call close_hdf5_file(targetMesh_hdf5_file_id)
      !----------------------------------------------------------------------------------------------------------------------------------     

      do res_inst=results_first,results_last,results_step

         call set_hdf5_resultsFile_name(source_base_resultsFile_h5,res_inst,source_full_resultsFile_h5)
         call set_hdf5_resultsFile_name(target_base_resultsFile_h5,res_inst,target_full_resultsFile_h5)

         if(mpi_rank.eq.0) then
            write(*,*) '# Doing results inst:',res_inst
            write(*,*) ' - Source resFile: ',trim(source_full_resultsFile_h5)
            write(*,*) ' - Target resFile: ',trim(target_full_resultsFile_h5)
         end if

         call open_hdf5_file(source_full_resultsFile_h5,sourceRes_hdf5_file_id)

         call check_datasets_in_source_results_file(sourceRes_hdf5_file_id,dsetScaCnt,dsetVecCnt,&
               datasetsScalarFields,datasetsVectorFields,datasetsScaExist,datasetsVecExist)

         call create_hdf5_file(target_full_resultsFile_h5,targetRes_hdf5_file_id)

         call create_vtkhdf_unstructuredGrid_struct_for_resultsFile(target_meshFile_h5_full_name,targetRes_hdf5_file_id)

         call create_fields_datasets_in_new_results_file(targetRes_hdf5_file_id,numNodesTrgtTotal_i8,&
               dsetScaCnt,dsetVecCnt,datasetsScalarFields,datasetsVectorFields,datasetsScaExist,datasetsVecExist)

         call copy_datasets_results_from_source_to_target(sourceRes_hdf5_file_id,targetRes_hdf5_file_id,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                  numNodesTrgtRank,trgtRankNodeStart_i8,numNodesTrgtTotal_i8,mapNodeIdTrgtToMpi_jv,&
                  dsetScaCnt,dsetVecCnt,datasetsScalarFields,datasetsVectorFields,datasetsScaExist,datasetsVecExist)

         call close_hdf5_file(targetRes_hdf5_file_id)
         if(mpi_rank.eq.0) write(*,*) ' # New results file ',trim(target_full_resultsFile_h5),' succesfully generated!'

         call close_hdf5_file(sourceRes_hdf5_file_id)

      end do

      !---------------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------
      call end_hdf5_interface()
      !----------------------------------------------------------------------------------------------

   end subroutine copy_results_same_mesh_Npartitions

   subroutine get_mesh_porder_from_hdf5(meshFile_h5_full_name,mporder)
      implicit none
      character(512),intent(in) :: meshFile_h5_full_name
      integer(4),intent(out) :: mporder
      integer(hid_t) :: hdf5_file_id
      character(128) :: dsetname
      integer(hsize_t),dimension(1) :: ms_dims
      integer(hssize_t),dimension(1) :: ms_offset
      integer(4),dimension(1) :: aux_array

      call open_hdf5_file(meshFile_h5_full_name,hdf5_file_id)

      !------------------------------------------------------------------------------------------------
      ms_dims(1) = 1
      ms_offset(1) = 0

      dsetname = '/order/porder'
      call read_dataspace_1d_int4_hyperslab_parallel(hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array)
      mporder = aux_array(1)

      call close_hdf5_file(hdf5_file_id)

   end subroutine get_mesh_porder_from_hdf5

   subroutine generate_connec_and_coords_for_targetRank(mnnode,trgtRank,trgtRankElemStart,trgtRankElemEnd,numElemsTrgtRank,&
            numNodesTrgtRank,connecTrgtRank,listNodesTrgtRank,coordNodesTrgtRank)
      implicit none
      integer(4),intent(in) :: mnnode,trgtRank,trgtRankElemStart,trgtRankElemEnd,numElemsTrgtRank
      integer(4),intent(out) :: numNodesTrgtRank
      integer(4),allocatable,intent(inout) :: connecTrgtRank(:,:),listNodesTrgtRank(:)
      real(rp),allocatable,intent(inout) :: coordNodesTrgtRank(:,:)

      integer(4) :: iElemMpiRank,iElemTrgtRank,elemOffset,iNodeMpi,iAux,jAux,nodeCnt
      integer(4),allocatable :: rawNodeListTrgtRank(:)

      allocate(connecTrgtRank(numElemsTrgtRank,mnnode))

      elemOffset = (trgtRankElemStart-1)

      do iElemTrgtRank=1,numElemsTrgtRank
         iElemMpiRank = iElemTrgtRank + elemOffset 
         connecTrgtRank(iElemTrgtRank,:) = connecParOrig(iElemMpiRank,:)
         !write(*,*) 'connecTrgt',connecTrgtRank(iElemTrgtRank,1:16)
      end do

      allocate(rawNodeListTrgtRank(numElemsTrgtRank*mnnode))
      !$acc kernels
      rawNodeListTrgtRank(:) = 0
      !$acc end kernels

      ! add all the iNodes of this mpirank
      nodeCnt=0
      do iAux = 1,numElemsTrgtRank
         do jAux = 1,mnnode
            iNodeMpi = connecTrgtRank(iAux,jAux)
            nodeCnt=nodeCnt+1
            rawNodeListTrgtRank(nodeCnt) = iNodeMpi
         end do
      end do

      !sorting the nodes id
      call quicksort_array_int4(rawNodeListTrgtRank)

      jAux=0
      nodeCnt=0
      do iAux = 1,(numElemsTrgtRank*mnnode)
         if((rawNodeListTrgtRank(iAux)).ne.(jAux)) then
            nodeCnt=nodeCnt+1
            jAux = rawNodeListTrgtRank(iAux)
         end if
      end do
      numNodesTrgtRank = nodeCnt
      !write(*,*) 'newMethod.[',mpi_rank,'] numNodesTrgtRank ',numNodesTrgtRank

      allocate(listNodesTrgtRank(numNodesTrgtRank))
      allocate(coordNodesTrgtRank(numNodesTrgtRank,ndime))

      !$acc kernels
      coordNodesTrgtRank(:,:) = 0
      listNodesTrgtRank(:) = 0
      !$acc end kernels

      jAux=0
      nodeCnt=0
      do iAux = 1,(numElemsTrgtRank*mnnode)
         if((rawNodeListTrgtRank(iAux)).ne.(jAux)) then
            nodeCnt=nodeCnt+1
            iNodeMpi = rawNodeListTrgtRank(iAux)
            listNodesTrgtRank(nodeCnt)    = iNodeMpi
            coordNodesTrgtRank(nodeCnt,:) = coordPar(iNodeMpi,:)
            jAux = rawNodeListTrgtRank(iAux)
         end if
      end do

      deallocate(rawNodeListTrgtRank)
      !write(*,*) 'listNodesTrgtRank',listNodesTrgtRank(:)

   end subroutine generate_connec_and_coords_for_targetRank


   subroutine reorder_nodes_in_trgtRank(mporder,mnnode,trgtRank,target_Nprocs,&
               numElemsTrgtRank,numNodesTrgtRank,sizeConnecVTKTrgtRank,linOutput,listNodesTrgtRankMpi,connecTrgtRankMpi,iNodeStartPar_i8,&
               mapNodeIdNewToOld,connecVTKTrgtRank,connecParOrigTrgtRank,coordVTK,coord)
      implicit none
      integer(4),intent(in) :: mporder,mnnode
      integer(4),intent(in) :: trgtRank,target_Nprocs,numElemsTrgtRank,numNodesTrgtRank,sizeConnecVTKTrgtRank
      logical,intent(in)    :: linOutput
      integer(4),intent(in) :: listNodesTrgtRankMpi(numNodesTrgtRank),connecTrgtRankMpi(numElemsTrgtRank,mnnode)
      integer(8),dimension(0:target_Nprocs-1),intent(in) :: iNodeStartPar_i8

      integer(4),intent(out) :: mapNodeIdNewToOld(numNodesTrgtRank)
      integer(4),intent(out) :: connecVTKTrgtRank(sizeConnecVTKTrgtRank),connecParOrigTrgtRank(numElemsTrgtRank,mnnode)      
      real(rp),intent(out) :: coordVTK(numNodesTrgtRank,ndime),coord(numNodesTrgtRank,ndime)
      !integer(8),intent(out) :: globalIdSrlTrgtRank_i8(numNodesTrgtRank),globalIdParTrgtRank_i8(numNodesTrgtRank)
      !integer(8),intent(out),dimension(numNodesTrgtRank,2) :: globalIdSrlOrderedTrgtRank_i8

      integer(4),dimension(numNodesTrgtRank) :: isNodeAdded
      integer(4) :: m,indConn,indexIJK,indexGMSH,indexNew,nodeIndexCnt,indPosListNodes
      integer(4) :: iNodeL,iNodeMpi,iElem
      integer(8) :: iNodeGsrl,iNodeGpar
      integer(4) :: ii,jj,kk,igp,jgp,kgp
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)

      !$acc kernels
      isNodeAdded(:)=-1
      !globalIdSrlOrderedTrgtRank_i8(:,:)=-1
      !$acc end kernels
      nodeIndexCnt = 0
      indConn = -1

      !----------------------------------------------------------------------------------------------------------

      do iElem=1,numElemsTrgtRank
         !write(*,*) '# iElem',iElem
         do m=1,mnnode

            iNodeMpi = connecTrgtRankMpi(iElem,m)
            indPosListNodes = binarySearch_int_i4(listNodesTrgtRankMpi,iNodeMpi)

            if(isNodeAdded(indPosListNodes) < 0) then !node not added put it in the list
               nodeIndexCnt=nodeIndexCnt+1
               iNodeL = nodeIndexCnt

               isNodeAdded(indPosListNodes)  = iNodeL

               iNodeGPar = int(iNodeL,8) + iNodeStartPar_i8(trgtRank) - 1

               !globalIdSrlTrgtRank_i8(iNodeL) = iNodeGsrl
               !globalIdParTrgtRank_i8(iNodeL) = iNodeGPar

               !globalIdSrlOrderedTrgtRank_i8(iNodeL,1) = iNodeGsrl
               !globalIdSrlOrderedTrgtRank_i8(iNodeL,2) = iNodeL

               coordVTK(iNodeL,:) = coord(indPosListNodes,:)

            else
               iNodeL = isNodeAdded(indPosListNodes)
            endif

            mapNodeIdNewToOld(iNodeL)      = iNodeMpi
            connecParOrigTrgtRank(iElem,m) = iNodeL

            !write(*,*) ' -inode',m,'iNode',iNodeL,'iNodeMpi',iNodeMpi

            if(.not.(linOutput)) then
               indConn = (iElem-1)*mnnode + mesh_vtk2ijk(m)
               connecVTKTrgtRank(indConn) = iNodeL
            end if

         end do

         if(linOutput) then
            indConn = (iElem-1)*((mporder**3)*8)
            !write(*,*) 'connecVTKRANG',indConn,'-',indConn+((mporder**3)*8)
            do ii=0,(mporder-1)
               do jj=0,(mporder-1)
                  do kk=0,(mporder-1)
                     !1.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElem,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL,'indConn',indConn
                     !2.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !3.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !4.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !5.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !6.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !7.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii+1)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                     !8.-------------------------------
                     igp=ijk_gmsh_to_sod2d(ii)-1
                     jgp=ijk_gmsh_to_sod2d(jj+1)-1
                     kgp=ijk_gmsh_to_sod2d(kk+1)-1
                     indexIJK = get_indexIJK_sod2d(mporder,igp,jgp,kgp)
                     indexNew = mesh_gmsh2ijk(indexIJK)
                     iNodeL = connecParOrigTrgtRank(iElem,indexNew)
                     indConn = indConn+1
                     connecVTKTrgtRank(indConn) = iNodeL
                     !write(*,*) 'iE',iElemL,'ijk',indexIJK,'new',indexNew,'iNodeL',iNodeL!,'indConn',indConn
                  end do
               end do
            end do
         end if

      end do

      deallocate(ijk_gmsh_to_sod2d,ijk_sod2d_to_gmsh)

   end subroutine reorder_nodes_in_trgtRank

   subroutine copy_datasets_results_from_source_to_target(source_file_id,target_file_id,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                                                         numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv,&
                                                         dsetScaCnt,dsetVecCnt,datasetsScalarFields,datasetsVectorFields,datasetsScaExist,datasetsVecExist)
      implicit none
      integer(hid_t),intent(in) :: source_file_id,target_file_id
      integer(4),intent(in) :: numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),intent(in) :: numNodesTrgtRank(numTrgtRanksInMpiRank)
      integer(8),intent(in) :: trgtRankNodeStart(numTrgtRanksInMpiRank),numNodesTrgtTotal
      type(jagged_vector_int4),intent(in) :: mapNodeIdTrgtToMpi_jv
      integer(4),intent(in) :: dsetScaCnt,dsetVecCnt
      character(256),intent(in) :: datasetsScalarFields(max_num_saved_fields),datasetsVectorFields(max_num_saved_fields)
      logical,intent(in) :: datasetsScaExist(max_num_saved_fields),datasetsVecExist(max_num_saved_fields)

      real(rp_vtk) :: sourceScalarField(numNodesRankPar),sourceVectorField(numNodesRankPar,ndime)
      integer(4) :: iTrgtRank,iSca,iVec
      character(512) :: dsetname

      do iSca=1,dsetScaCnt
         dsetname = datasetsScalarFields(iSca)
         if(datasetsScaExist(iSca)) then 
            if(mpi_rank.eq.0) write(*,*) '  - Copying scalar field',trim(adjustl(dsetname)),'(id',iSca,')...'
            call read_and_load_source_scalarfield(dsetname,source_file_id,sourceScalarField)
            call copy_scalarfield_dataset_results_in_target(dsetname,sourceScalarField,target_file_id,numTrgtRanksInMpiRank,&
                     maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
         end if
      end do

      do iVec=1,dsetVecCnt
         dsetname = datasetsVectorFields(iVec)
         if(datasetsVecExist(iVec)) then 
            if(mpi_rank.eq.0) write(*,*) '  - Copying vector field',trim(adjustl(dsetname)),'(id',iVec,')...'
            call read_and_load_source_vectorfield(dsetname,source_file_id,sourceVectorField)
            call copy_vectorfield_dataset_results_in_target(dsetname,sourceVectorField,target_file_id,numTrgtRanksInMpiRank,&
                     maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
         end if
      end do

   end subroutine copy_datasets_results_from_source_to_target

   subroutine copy_scalarfield_dataset_results_in_target(dsetname,sourceScalarField,target_file_id,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                                             numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id
      integer(4),intent(in) :: numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),intent(in) :: numNodesTrgtRank(numTrgtRanksInMpiRank)
      integer(8),intent(in) :: trgtRankNodeStart(numTrgtRanksInMpiRank),numNodesTrgtTotal
      type(jagged_vector_int4),intent(in) :: mapNodeIdTrgtToMpi_jv
      real(rp_vtk),intent(in) :: sourceScalarField(numNodesRankPar)

      integer(4) :: iTrgtRank

	   !--------------------------------------------------------------------------------
      do iTrgtRank=1,numTrgtRanksInMpiRank
         call copy_scalarfield_result_in_trgtRank(dsetname,sourceScalarField,target_file_id,numNodesTrgtRank(iTrgtRank),&
               trgtRankNodeStart(iTrgtRank),numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv%vector(iTrgtRank)%elems)
      end do

      do iTrgtRank=(numTrgtRanksInMpiRank+1),maxNumTrgtRanks
         call dummy_copy_scalarfield_result_in_trgtRank(dsetname,target_file_id)
      end do
      !----------------------------------------------------------------------

   end subroutine copy_scalarfield_dataset_results_in_target

   subroutine copy_vectorfield_dataset_results_in_target(dsetname,sourceVectorField,target_file_id,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                                             numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id
      integer(4),intent(in) :: numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),intent(in) :: numNodesTrgtRank(numTrgtRanksInMpiRank)
      integer(8),intent(in) :: trgtRankNodeStart(numTrgtRanksInMpiRank),numNodesTrgtTotal
      type(jagged_vector_int4),intent(in) :: mapNodeIdTrgtToMpi_jv
      real(rp_vtk),intent(in) :: sourceVectorField(numNodesRankPar,ndime)

      integer(4) :: iTrgtRank

	   !--------------------------------------------------------------------------------
      do iTrgtRank=1,numTrgtRanksInMpiRank
         call copy_vectorfield_result_in_trgtRank(dsetname,sourceVectorField,target_file_id,numNodesTrgtRank(iTrgtRank),&
               trgtRankNodeStart(iTrgtRank),numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv%vector(iTrgtRank)%elems)
      end do

      do iTrgtRank=(numTrgtRanksInMpiRank+1),maxNumTrgtRanks
         call dummy_copy_vectorfield_result_in_trgtRank(dsetname,target_file_id)
      end do
      !----------------------------------------------------------------------

   end subroutine copy_vectorfield_dataset_results_in_target


   subroutine read_and_load_source_scalarfield(dsetname,source_file_id,sourceScalarField)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: source_file_id
      real(rp_vtk),intent(out) :: sourceScalarField(numNodesRankPar)

      integer(hsize_t) :: src_ds_dims(1),src_ms_dims(1)
      integer(hssize_t) :: src_ms_offset(1)

      !---------------------------------------------------------------------------------------------------------------------------

      src_ds_dims(1) = int(totalNumNodesPar,hsize_t)
      src_ms_dims(1) = int(numNodesRankPar,hsize_t)
      src_ms_offset(1) = int(rankNodeStart,hssize_t)-1

      call read_dataspace_1d_real_rp_vtk_hyperslab_parallel(source_file_id,dsetname,src_ms_dims,src_ms_offset,sourceScalarField)

   end subroutine read_and_load_source_scalarfield

   subroutine read_and_load_source_vectorfield(dsetname,source_file_id,sourceVectorField)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: source_file_id
      real(rp_vtk),intent(out) :: sourceVectorField(numNodesRankPar,ndime)

      integer(hsize_t) :: src_ds_dims2d(2),src_ms_dims2d(2)
      integer(hssize_t) :: src_ms_offset2d(2)

      !---------------------------------------------------------------------------------------------------------------------------

      src_ds_dims2d(1) = int(ndime,hsize_t)
      src_ds_dims2d(2) = int(totalNumNodesPar,hsize_t)
      src_ms_dims2d(1) = int(ndime,hsize_t)
      src_ms_dims2d(2) = int(numNodesRankPar,hsize_t)
      src_ms_offset2d(1) = 0
      src_ms_offset2d(2) = int(rankNodeStart,hssize_t)-1

      call read_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(source_file_id,dsetname,src_ms_dims2d,src_ms_offset2d,sourceVectorField)

   end subroutine read_and_load_source_vectorfield

   subroutine copy_scalarfield_result_in_trgtRank(dsetname,sourceScalarField,target_file_id,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdNewToOld)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id
      integer(4),intent(in) :: numNodesTrgtRank
      integer(8),intent(in) :: trgtRankNodeStart,numNodesTrgtTotal
      integer(4),intent(in) :: mapNodeIdNewToOld(numNodesTrgtRank)
      real(rp_vtk),intent(in) :: sourceScalarField(numNodesRankPar)

      integer(4) :: iNodeTrgt,iNodeSrc
      integer(hsize_t) :: trgt_ds_dims(1),trgt_ms_dims(1)
      integer(hssize_t) :: trgt_ms_offset(1)
      real(rp_vtk) :: targetScalarField(numNodesTrgtRank)

      !---------------------------------------------------------------------------------------------------------------------------

      trgt_ds_dims(1) = int(numNodesTrgtTotal,hsize_t)
      trgt_ms_dims(1) = int(numNodesTrgtRank,hsize_t)
      trgt_ms_offset(1) = int(trgtRankNodeStart,hssize_t)-1

      do iNodeTrgt=1,numNodesTrgtRank
         iNodeSrc = mapNodeIdNewToOld(iNodeTrgt)
         targetScalarField(iNodeTrgt) = sourceScalarField(iNodeSrc)
      end do

      call write_dataspace_1d_real_rp_vtk_hyperslab_parallel(target_file_id,dsetname,trgt_ms_dims,trgt_ms_offset,targetScalarField)

   end subroutine copy_scalarfield_result_in_trgtRank

   subroutine copy_vectorfield_result_in_trgtRank(dsetname,sourceVectorField,target_file_id,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdNewToOld)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id
      integer(4),intent(in) :: numNodesTrgtRank
      integer(8),intent(in) :: trgtRankNodeStart,numNodesTrgtTotal
      integer(4),intent(in) :: mapNodeIdNewToOld(numNodesTrgtRank)
      real(rp_vtk),intent(in) :: sourceVectorField(numNodesRankPar,ndime)

      integer(4) :: iNodeTrgt,iNodeSrc
      integer(hsize_t) :: trgt_ds_dims2d(2),trgt_ms_dims2d(2)
      integer(hssize_t) :: trgt_ms_offset2d(2)
      real(rp_vtk) :: targetVectorField(numNodesTrgtRank,ndime)

      !---------------------------------------------------------------------------------------------------------------------------

      trgt_ds_dims2d(1) = int(ndime,hsize_t)
      trgt_ds_dims2d(2) = int(numNodesTrgtTotal,hsize_t)
      trgt_ms_dims2d(1) = int(ndime,hsize_t)
      trgt_ms_dims2d(2) = int(numNodesTrgtRank,hsize_t)
      trgt_ms_offset2d(1) = 0
      trgt_ms_offset2d(2) = int(trgtRankNodeStart,hssize_t)-1

      do iNodeTrgt=1,numNodesTrgtRank
         iNodeSrc = mapNodeIdNewToOld(iNodeTrgt)
         targetVectorField(iNodeTrgt,:) = sourceVectorField(iNodeSrc,:)
      end do

      call write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(target_file_id,dsetname,trgt_ms_dims2d,trgt_ms_offset2d,targetVectorField)

   end subroutine copy_vectorfield_result_in_trgtRank

   subroutine dummy_copy_scalarfield_result_in_trgtRank(dsetname,target_file_id)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id

      integer(hsize_t) :: trgt_ds_dims(1),trgt_ms_dims(1)
      integer(hssize_t) :: trgt_ms_offset(1)
      real(rp_vtk),allocatable :: targetScalarField(:)
      !---------------------------------------------------------------------------------------------------------------------------
      
      allocate(targetScalarField(0))

      trgt_ds_dims(1) = 0
      trgt_ms_dims(1) = 0
      trgt_ms_offset(1) = 0

      call write_dataspace_1d_real_rp_vtk_hyperslab_parallel(target_file_id,dsetname,trgt_ms_dims,trgt_ms_offset,targetScalarField)

      deallocate(targetScalarField)

   end subroutine dummy_copy_scalarfield_result_in_trgtRank

   subroutine dummy_copy_vectorfield_result_in_trgtRank(dsetname,target_file_id)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id

      integer(hsize_t) :: trgt_ds_dims2d(2),trgt_ms_dims2d(2)
      integer(hssize_t) :: trgt_ms_offset2d(2)
      real(rp_vtk),allocatable :: targetVectorField(:,:)
      !---------------------------------------------------------------------------------------------------------------------------

      allocate(targetVectorField(0,0))

      trgt_ds_dims2d(1) = 0
      trgt_ds_dims2d(2) = 0
      trgt_ms_dims2d(1) = 0
      trgt_ms_dims2d(2) = 0
      trgt_ms_offset2d(1) = 0
      trgt_ms_offset2d(2) = 0

      call write_dataspace_2d_tr_real_rp_vtk_hyperslab_parallel(target_file_id,dsetname,trgt_ms_dims2d,trgt_ms_offset2d,targetVectorField)

      deallocate(targetVectorField)

   end subroutine dummy_copy_vectorfield_result_in_trgtRank

   subroutine check_datasets_in_source_results_file(src_file_id,dsetScaCnt,dsetVecCnt,datasetsScalarFields,datasetsVectorFields,&
                                                   datasetsScaExist,datasetsVecExist)
      implicit none
      integer(hid_t),intent(in) :: src_file_id
      integer(4),intent(inout) :: dsetScaCnt,dsetVecCnt
      character(256),intent(out) :: datasetsScalarFields(max_num_saved_fields),datasetsVectorFields(max_num_saved_fields)
      logical,intent(out) :: datasetsScaExist(max_num_saved_fields),datasetsVecExist(max_num_saved_fields)
      integer(4) :: h5err
      character(512) :: groupname,dsetname
      logical :: dset_exists

      integer(4) :: iSca,iVec
	   !--------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------------------------------
      groupname = '/VTKHDF/PointData/'
      dsetScaCnt = 0
      dsetVecCnt = 0

      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'rho'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'mu_fluid'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'pr'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'E'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'eta'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'csound'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'machno'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'divU'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'qcrit'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'mut'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'mue'
      dsetScaCnt=dsetScaCnt+1
      datasetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'u'
      dsetVecCnt=dsetVecCnt+1
      datasetsVectorFields(dsetVecCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'gradRho'
      dsetVecCnt=dsetVecCnt+1
      datasetsVectorFields(dsetVecCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'curlU'
      dsetVecCnt=dsetVecCnt+1
      datasetsVectorFields(dsetVecCnt) = dsetname
      !--------------------------------------------------

      !--------------------------------------------------------------------------------
      datasetsScaExist(:) = .false.
      datasetsVecExist(:) = .false.

      do iSca=1,dsetScaCnt
         dsetname = datasetsScalarFields(iSca)
         call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
         datasetsScaExist(iSca) = dset_exists
         !write(*,*) 'checking if scaf',dsetname,'exists:',dset_exists
      end do

      do iVec=1,dsetVecCnt
         dsetname = datasetsVectorFields(iVec)
         call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
         datasetsVecExist(iVec) = dset_exists
         !write(*,*) 'checking if vec',dsetname,'exists:',dset_exists
      end do

   end subroutine check_datasets_in_source_results_file

   subroutine create_fields_datasets_in_new_results_file(trgt_file_id,numNodesParTotal_i8,dsetScaCnt,dsetVecCnt,&
                                    datasetsScalarFields,datasetsVectorFields,datasetsScaExist,datasetsVecExist)
      implicit none
      integer(hid_t),intent(in) :: trgt_file_id
      integer(8),intent(in) :: numNodesParTotal_i8
      integer(4),intent(in) :: dsetScaCnt,dsetVecCnt
      character(256),intent(in) :: datasetsScalarFields(max_num_saved_fields),datasetsVectorFields(max_num_saved_fields)
      logical,intent(in) :: datasetsScaExist(max_num_saved_fields),datasetsVecExist(max_num_saved_fields)
      integer(hid_t) :: dtype
      integer(hsize_t) :: ds_dims(1),ds_dims2d(2),aux_ds_dims
      integer(4) :: ds_rank,h5err,iSca,iVec
      character(512) :: dsetname
      !--------------------------------------------------------------------------------

      call select_dtype_rp_vtk(dtype)

      !----------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numNodesParTotal_i8

      do iSca=1,dsetScaCnt
         dsetname = datasetsScalarFields(iSca)
         if(datasetsScaExist(iSca)) then 
            !write(*,*) 'creating datascavec',dsetname
            call create_dataspace_hdf5(trgt_file_id,dsetname,ds_rank,ds_dims,dtype)
         end if
      end do
      !----------------------------------------------------------------------
      ds_rank = 2
      ds_dims2d(1) = ndime
      ds_dims2d(2) = numNodesParTotal_i8

      do iVec=1,dsetVecCnt
         dsetname = datasetsVectorFields(iVec)
         if(datasetsVecExist(iVec)) then 
            !write(*,*) 'creating datasetvec',dsetname
            call create_dataspace_hdf5(trgt_file_id,dsetname,ds_rank,ds_dims2d,dtype)
         end if
      end do
      !----------------------------------------------------------------------

   end subroutine create_fields_datasets_in_new_results_file

end module mod_copy_results
