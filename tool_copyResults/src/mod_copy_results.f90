module mod_copy_results
   use mod_mpi
   use mod_ijk_indices
   use mod_mpi_mesh
   use mod_comms
   use mod_hdf5
   use mod_saveFields
   use mod_partition_utils

   implicit none

contains

   subroutine copy_results_same_mesh_Npartitions(mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName,target_Nprocs,&
                                                type_resultsFile,generateMesh,results_first,results_last,results_step)
      implicit none
      character(len=*), intent(in)  :: mesh_h5_filePath,mesh_h5_fileName,results_h5_filePath,results_h5_fileName
      integer(4),intent(in) :: target_Nprocs,type_resultsFile,results_first,results_last,results_step
      logical,intent(in) :: generateMesh
      character(512) :: source_meshFile_h5_full_name,target_meshFile_h5_full_name
      character(512) :: source_base_resultsFile_h5,target_base_resultsFile_h5
      character(512) :: source_full_resultsFile_h5,target_full_resultsFile_h5
      character(512) :: mappingFile_h5_full_name
      integer(4) :: mporder,mnnode,mngaus,mnpbou,res_ii,res_step,restart_file,load_step
      real(rp) :: time
      character(128) :: dsetname

      integer(4) :: numTrgtRanksInMpiRank,maxNumTrgtRanks

      integer(8) :: numNodesTrgtTotal_i8
      integer(8),allocatable :: trgtRankNodeStart_i8(:),trgtRankNodeEnd_i8(:)
      integer(4),allocatable :: mapNodeRankTrgt(:,:),numNodesTrgtRank(:)
      type(jagged_vector_int4) :: mapNodeIdTrgtToMpi_jv

      integer(hid_t) :: targetRes_hdf5_file_id,sourceRes_hdf5_file_id

      integer(4) :: numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S
      character(len=256),allocatable :: dsetsScalarFieldsOrig(:),dsetsVectorFieldsOrig(:),dsetsV2SFieldsOrig(:)
      character(len=256),allocatable :: dsetsScalarFieldsTrgt(:),dsetsVectorFieldsTrgt(:),dsetsV2SFieldsTrgt(:,:)

      integer(4) :: iRank,app_comm_group,smallWorld_group,smallWorld_comm
      integer(4),allocatable :: smallWorldMpiRanks(:)


      !----------------------------------------------------------------------------------------------
      call init_hdf5_interface()
      !----------------------------------------------------------------------------------------------

      !--- inital check -------------
      if(generateMesh) then
         if(mpi_size.gt.target_Nprocs) then
            write(*,*) "If generateMesh=.true. then must be (mpi_size<=target_Nprocs).",&
                       "Fix current values mpi_size",mpi_size,"target_Nprocs",target_Nprocs,&
                       "Aborting!"
       	   call MPI_Abort(app_comm,-1,mpi_err)
         endif
      else
         !if(target_Nprocs.gt.mpi_size) then
         !   write(*,*) "If generateMesh=.false. then must be (target_Nprocs<=mpi_size).",&
         !              "Fix current values mpi_size",mpi_size,"target_Nprocs",target_Nprocs,&
         !              "Aborting!"
         !  call MPI_Abort(app_comm,-1,mpi_err)
         !endif
      endif

      call init_saveFields()

      call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName,mpi_size,source_meshFile_h5_full_name)
      call set_hdf5_meshFile_name(mesh_h5_filePath,mesh_h5_fileName,target_Nprocs,target_meshFile_h5_full_name)

      if(mpi_rank.eq.0) then
         write(*,*) 'Source mesh: ',trim(source_meshFile_h5_full_name)
         write(*,*) 'Target mesh: ',trim(target_meshFile_h5_full_name)
      end if

      call get_mesh_porder_from_hdf5(source_meshFile_h5_full_name,mporder)
      call get_porder_values(mporder,mnnode,mngaus,mnpbou)

      !Loading source mesh file
      call load_hdf5_meshFile(source_meshFile_h5_full_name)

      if(generateMesh) then

         if(mpi_size.eq.target_Nprocs) then
            if(mpi_rank.eq.0) write(*,*) "For generateMesh .true. is not allowed mpi_size == target_Nranks! Aborting!"
       	   call MPI_Abort(app_comm,-1,mpi_err)
         end if

         call generate_new_mesh_for_viz(target_meshFile_h5_full_name,mporder,mnnode,target_Nprocs,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                                       numNodesTrgtRank,numNodesTrgtTotal_i8,trgtRankNodeStart_i8,trgtRankNodeEnd_i8,mapNodeIdTrgtToMpi_jv)

      else
         if(mpi_size .eq. target_Nprocs) then
            call do_naive_one_to_one_mapping(numTrgtRanksInMpiRank,maxNumTrgtRanks,numNodesTrgtRank,&
                                       numNodesTrgtTotal_i8,trgtRankNodeStart_i8,trgtRankNodeEnd_i8,mapNodeIdTrgtToMpi_jv)
         else

            allocate(mapNodeRankTrgt(numNodesRankPar,2))
            allocate(numNodesTrgtRank(0:target_Nprocs-1))
            allocate(trgtRankNodeStart_i8(target_Nprocs))
            allocate(trgtRankNodeEnd_i8(target_Nprocs))

            call open_target_mesh_and_do_mapping(target_meshFile_h5_full_name,target_Nprocs,mapNodeRankTrgt,&
                                                numNodesTrgtRank,numNodesTrgtTotal_i8,trgtRankNodeStart_i8,trgtRankNodeEnd_i8)

            if(type_resultsFile.eq.5) then

               call set_hdf5_mappingFile_name(mesh_h5_filePath,mesh_h5_fileName,mpi_size,target_Nprocs,mappingFile_h5_full_name)

               call save_mapping_file(mappingFile_h5_full_name,mapNodeRankTrgt)

               call end_hdf5_interface()
               deallocate(mapNodeRankTrgt,numNodesTrgtRank)
               return
            end if
         end if
      end if
      
      !-----------------------------------------------------------------------------------------------
      !   Generacio de fitxers resultats!

      !if((.not.(generateMesh)).and.(type_resultsFile.lt.4)) then
      !   if(mpi_rank.eq.0) write(*,*) "Generate mesh and type_resultsFile<4 not yet supported! Aborting!"
    	!   call MPI_Abort(app_comm,-1,mpi_err)
      !endif

      if(type_resultsFile.eq.4) then
         if(generateMesh) then
            if(mpi_rank.eq.0) write(*,*) "For type_resultsFile=4 only generateMesh=.false. supported! Aborting!"
       	   call MPI_Abort(app_comm,-1,mpi_err)
         end if
         if(mpi_size.ne.target_Nprocs) then
            if(mpi_rank.eq.0) write(*,*) "For type_resultsFile=4 only mpi_size == target_Nranks supported! Aborting!"
       	   call MPI_Abort(app_comm,-1,mpi_err)
         end if
      endif

      res_step = results_step

      if(type_resultsFile .eq. 1) then
         call set_hdf5_resultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,mpi_size,source_base_resultsFile_h5)
         call set_hdf5_resultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,target_Nprocs,target_base_resultsFile_h5)
      else if(type_resultsFile .eq. 2) then
         call set_hdf5_avgResultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,mpi_size,source_base_resultsFile_h5)
         call set_hdf5_avgResultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,target_Nprocs,target_base_resultsFile_h5)
      else if(type_resultsFile .eq. 3) then
         call set_hdf5_restartFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,mpi_size,source_base_resultsFile_h5)
         call set_hdf5_restartFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,target_Nprocs,target_base_resultsFile_h5)
      else if(type_resultsFile .eq. 4) then
         call set_hdf5_resultsFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,mpi_size,source_base_resultsFile_h5)
         call set_hdf5_restartFile_baseName(results_h5_filePath,results_h5_fileName,mesh_h5_fileName,target_Nprocs,target_base_resultsFile_h5)
         restart_file = results_step
         res_step = 1
      else
         write(*,*) "Wrong type_resultsFile! Must be 1,2,3, 4 or 5 (1:inst, 2:avg, 3:restart, 4:inst_2_restart, 5:mapping)! Aborting!"
       	call MPI_Abort(app_comm,-1,mpi_err)
      end if

      do res_ii=results_first,results_last,res_step

         if(type_resultsFile .eq. 1) then
            ! INST RESULTS
            call set_hdf5_resultsFile_name(source_base_resultsFile_h5,res_ii,source_full_resultsFile_h5)
            call set_hdf5_resultsFile_name(target_base_resultsFile_h5,res_ii,target_full_resultsFile_h5)
         else if(type_resultsFile .eq. 2) then
            ! AVG RESULTS
            call set_hdf5_avgResultsFile_name(source_base_resultsFile_h5,res_ii,source_full_resultsFile_h5)
            call set_hdf5_avgResultsFile_name(target_base_resultsFile_h5,res_ii,target_full_resultsFile_h5)
         else if(type_resultsFile .eq. 3) then
            ! RESTART FILES
            call set_hdf5_restartFile_name(source_base_resultsFile_h5,res_ii,source_full_resultsFile_h5)
            call set_hdf5_restartFile_name(target_base_resultsFile_h5,res_ii,target_full_resultsFile_h5)
         else if(type_resultsFile .eq. 4) then
            ! INST 2 RESTART
            call set_hdf5_resultsFile_name(source_base_resultsFile_h5,res_ii,source_full_resultsFile_h5)
            call set_hdf5_restartFile_name(target_base_resultsFile_h5,restart_file,target_full_resultsFile_h5)
         end if

         if(mpi_rank.eq.0) then
            write(*,*) '# Doing results ii:',res_ii
            write(*,*) ' - Source resFile: ',trim(source_full_resultsFile_h5)
            write(*,*) ' - Target resFile: ',trim(target_full_resultsFile_h5)
         end if

         call open_hdf5_file(source_full_resultsFile_h5,sourceRes_hdf5_file_id)

         call generate_dsets_to_copy(type_resultsFile,sourceRes_hdf5_file_id,numDsetSca,numDsetVec,numDsetV2S,&
                                     maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
                                     dsetsScalarFieldsOrig,dsetsVectorFieldsOrig,dsetsV2SFieldsOrig,&
                                     dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)

         call create_hdf5_file(target_full_resultsFile_h5,targetRes_hdf5_file_id)

         if(type_resultsFile .le. 2) then
            call create_vtkhdf_unstructuredGrid_struct_for_resultsFile(target_meshFile_h5_full_name,targetRes_hdf5_file_id)
         end if

         call create_fields_dsets_in_new_results_file(targetRes_hdf5_file_id,numNodesTrgtTotal_i8,&
               numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
               dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)

         if(mpi_size.le.target_Nprocs) then
            write(*,*) 'classical copy results...'
            call copy_dsets_results_for_generated_mesh(sourceRes_hdf5_file_id,targetRes_hdf5_file_id,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                     numNodesTrgtRank,trgtRankNodeStart_i8,numNodesTrgtTotal_i8,mapNodeIdTrgtToMpi_jv,&
                     numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
                     dsetsScalarFieldsOrig,dsetsVectorFieldsOrig,dsetsV2SFieldsOrig,&
                     dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)
         else
            write(*,*) 'new feature copy results... implementing!'
            call copy_dsets_results_for_loaded_mesh(sourceRes_hdf5_file_id,targetRes_hdf5_file_id,target_Nprocs,&
                     numNodesTrgtRank,trgtRankNodeStart_i8,numNodesTrgtTotal_i8,mapNodeRankTrgt,&
                     numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
                     dsetsScalarFieldsOrig,dsetsVectorFieldsOrig,dsetsV2SFieldsOrig,&
                     dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)

#if 0
            !idea boja....
            !if(mpi_rank.lt.target_Nprocs) then
            allocate(smallWorldMpiRanks(0:target_Nprocs-1))
            call MPI_Comm_group(app_comm, app_comm_group, mpi_err)
            do iRank=0,(target_Nprocs-1)
               smallWorldMpiRanks(iRank) = iRank
            end do
            call MPI_Group_incl(app_comm_group, target_Nprocs, smallWorldMpiRanks, smallWorld_group, mpi_err)
            call MPI_Comm_create(app_comm, smallWorld_group, smallWorld_comm, mpi_err)
            !call MPI_Comm_rank(smallWorld_comm, mpi_rank, mpi_err)
            !call MPI_Comm_size(smallWorld_comm, mpi_size, mpi_err)

               !call load_hdf5_meshFile(target_meshFile_h5_full_name)
            !endif 

            !restore old mpi_rank and mpi_size
            !call MPI_Comm_rank(app_comm, mpi_rank, mpi_err)
            !call MPI_Comm_size(app_comm, mpi_size, mpi_err)

            !call load_hdf5_meshFile(source_meshFile_h5_full_name)
#endif
         end if

         deallocate(dsetsScalarFieldsOrig,dsetsVectorFieldsOrig)
         deallocate(dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt)
         deallocate(dsetsV2SFieldsOrig,dsetsV2SFieldsTrgt)

         if(type_resultsFile .eq. 1) then
            ! INST RESULTS
            dsetname = 'time'
            call read_real_rp_vtk_in_rp_dataset_hdf5_file(sourceRes_hdf5_file_id,dsetname,time)
            call save_real_rp_in_rp_vtk_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,time)
         else if(type_resultsFile .eq. 2) then
            dsetname = 'elapsed_avgTime'
            call read_real_rp_vtk_in_rp_dataset_hdf5_file(sourceRes_hdf5_file_id,dsetname,time)
            call save_real_rp_in_rp_vtk_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,time)

            dsetname = 'initial_avgTime'
            call read_real_rp_vtk_in_rp_dataset_hdf5_file(sourceRes_hdf5_file_id,dsetname,time)
            call save_real_rp_in_rp_vtk_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,time)
         else if(type_resultsFile .eq. 3) then
            dsetname = 'time'
            call read_real_rp_vtk_in_rp_dataset_hdf5_file(sourceRes_hdf5_file_id,dsetname,time)
            call save_real_rp_in_rp_vtk_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,time)

            dsetname = 'istep'
            call read_int4_in_dataset_hdf5_file(sourceRes_hdf5_file_id,dsetname,load_step)
            call save_int4_in_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,load_step)
         else if(type_resultsFile .eq. 4) then
            dsetname = 'time'
            call read_real_rp_vtk_in_rp_dataset_hdf5_file(sourceRes_hdf5_file_id,dsetname,time)
            call save_real_rp_in_rp_vtk_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,time)

            dsetname = 'istep'
            call save_int4_in_dataset_hdf5_file(targetRes_hdf5_file_id,dsetname,res_ii)
         end if

         call close_hdf5_file(targetRes_hdf5_file_id)
         if(mpi_rank.eq.0) write(*,*) ' # New results file ',trim(target_full_resultsFile_h5),' succesfully generated!'

         call close_hdf5_file(sourceRes_hdf5_file_id)

      end do

      !---------------------------------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------
      call end_hdf5_interface()
      !----------------------------------------------------------------------------------------------

   end subroutine copy_results_same_mesh_Npartitions

   subroutine do_naive_one_to_one_mapping(numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                                       numNodesTrgtRank,numNodesTrgtTotal_i8,trgtRankNodeStart_i8,trgtRankNodeEnd_i8,mapNodeIdTrgtToMpi_jv)
      implicit none
      integer(4),intent(inout) :: numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),intent(inout),allocatable :: numNodesTrgtRank(:)
      integer(8),intent(inout) :: numNodesTrgtTotal_i8
      integer(8),intent(inout),allocatable :: trgtRankNodeStart_i8(:),trgtRankNodeEnd_i8(:)
      type(jagged_vector_int4),intent(inout) :: mapNodeIdTrgtToMpi_jv
      integer(4) :: iNode

      numTrgtRanksInMpiRank = 1
      maxNumTrgtRanks = 1

      allocate(numNodesTrgtRank(1))
      allocate(trgtRankNodeStart_i8(1))
      allocate(trgtRankNodeEnd_i8(1))
      numNodesTrgtRank(1)     = numNodesRankPar
      trgtRankNodeStart_i8(1) = rankNodeStart
      trgtRankNodeEnd_i8(1)   = rankNodeEnd

      numNodesTrgtTotal_i8 = totalNumNodesPar

      allocate(mapNodeIdTrgtToMpi_jv%vector(1))
      allocate(mapNodeIdTrgtToMpi_jv%vector(1)%elems(numNodesRankPar))

      do iNode=1,numNodesRankPar
         mapNodeIdTrgtToMpi_jv%vector(1)%elems(iNode) = iNode
      end do

   end subroutine

   subroutine generate_new_mesh_for_viz(target_meshFile_h5_full_name,mporder,mnnode,target_Nprocs,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                                       numNodesTrgtRank,numNodesTrgtTotal_i8,trgtRankNodeStart_i8,trgtRankNodeEnd_i8,mapNodeIdTrgtToMpi_jv)
      implicit none
      character(*),intent(in) :: target_meshFile_h5_full_name
      integer(4),intent(in) :: mporder,mnnode,target_Nprocs
      integer(4),intent(inout) :: numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),intent(inout),allocatable :: numNodesTrgtRank(:)
      integer(8),intent(inout) :: numNodesTrgtTotal_i8
      integer(8),intent(inout),allocatable :: trgtRankNodeStart_i8(:),trgtRankNodeEnd_i8(:)
      type(jagged_vector_int4),intent(inout) :: mapNodeIdTrgtToMpi_jv

      integer(4) :: trgtRankInMpiRankStart,trgtRankInMpiRankEnd
      integer(4),allocatable :: trgtRanksInMpiRank(:),mapTrgtRankToMpiRank(:)
      integer(4),allocatable :: mpiRankElemStart(:),mpiRankElemEnd(:),trgtRankElemStart(:),trgtRankElemEnd(:)
      integer(4),allocatable :: numElemsTrgtRank(:)
      integer(4),allocatable :: numElemsVTKTrgtRank(:),sizeConnecVTKTrgtRank(:)
      type(jagged_vector_int4) :: listNodesTrgtRankMpiId_jv,connecVTKTrgtRank_jv
      type(jagged_matrix_int4) :: connecTrgtRankMpiId_jm,connecParOrigTrgtRank_jm
      type(jagged_matrix_real8) :: coordTrgtRank_jm,coordVTKTrgtRank_jm
      type(jagged_matrix_real8) :: quality_jm

      integer(8),dimension(0:target_Nprocs-1) :: iNodeStartPar_i8
      integer(4) :: connecChunkSize = 10000000
      logical :: evalMeshQuality=.false.
      integer(4) :: iTrgtRank,trgtRank
      integer(hid_t) :: targetMesh_hdf5_file_id

      !---------------------------------------------------------------------------------------------
      ! Generacio malla dummy, en un futur fer possible que lusuari introdeuxi una ja feta!
      if(mpi_rank.eq.0) write(*,*) '# Generating new target mesh file:',trim(target_meshFile_h5_full_name)
      if(mpi_rank.eq.0) write(*,*) ' 1. Generating new partitioning and connec...'

      !Distributing the target_Nprocs in the current Nprocs
      call distribute_ranks2Part_in_mpiRank(target_Nprocs,trgtRankInMpiRankStart,trgtRankInMpiRankEnd,&
                        numTrgtRanksInMpiRank,maxNumTrgtRanks,trgtRanksInMpiRank,mapTrgtRankToMpiRank)

      !---------------------------------------------------------
      allocate(    numNodesTrgtRank(numTrgtRanksInMpiRank))
      allocate(trgtRankNodeStart_i8(numTrgtRanksInMpiRank))
      allocate(  trgtRankNodeEnd_i8(numTrgtRanksInMpiRank))
      numNodesTrgtRank(:)     = 0
      trgtRankNodeStart_i8(:) = 0
      trgtRankNodeEnd_i8(:)   = 0

      allocate(mapNodeIdTrgtToMpi_jv%vector(numTrgtRanksInMpiRank))
      !----------------------------------------------------------------------------

      allocate(    numElemsTrgtRank(numTrgtRanksInMpiRank))
      allocate(   trgtRankElemStart(numTrgtRanksInMpiRank))
      allocate(     trgtRankElemEnd(numTrgtRanksInMpiRank))
      allocate(    mpiRankElemStart(numTrgtRanksInMpiRank))
      allocate(      mpiRankElemEnd(numTrgtRanksInMpiRank))
      numElemsTrgtRank(:)     = 0
      trgtRankElemStart(:)    = 0
      trgtRankElemEnd(:)      = 0
      mpiRankElemStart(:)     = 0
      mpiRankElemEnd(:)       = 0

      allocate(listNodesTrgtRankMpiId_jv%vector(numTrgtRanksInMpiRank))
      allocate(   connecTrgtRankMpiId_jm%matrix(numTrgtRanksInMpiRank))
      allocate(         coordTrgtRank_jm%matrix(numTrgtRanksInMpiRank))
      allocate(               quality_jm%matrix(numTrgtRanksInMpiRank))

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
      
      call define_parallelNodePartitioning(target_Nprocs,numTrgtRanksInMpiRank,numNodesTrgtRank,trgtRanksInMpiRank,mapTrgtRankToMpiRank,&
                                          trgtRankNodeStart_i8,trgtRankNodeEnd_i8,iNodeStartPar_i8,numNodesTrgtTotal_i8)

      allocate(coordVTKTrgtRank_jm%matrix(numTrgtRanksInMpiRank))
      allocate(connecVTKTrgtRank_jv%vector(numTrgtRanksInMpiRank))
      allocate(connecParOrigTrgtRank_jm%matrix(numTrgtRanksInMpiRank))

      allocate(numElemsVTKTrgtRank(numTrgtRanksInMpiRank))
      allocate(sizeConnecVTKTrgtRank(numTrgtRanksInMpiRank))
      numElemsVTKTrgtRank(:)   = 0
      sizeConnecVTKTrgtRank(:) = 0

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
               coordTrgtRank_jm%matrix(iTrgtRank)%elems,coordVTKTrgtRank_jm%matrix(iTrgtRank)%elems)

      end do

      !---------------------------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) write(*,*) ' 3. Generating new mesh partitioned in',target_Nprocs,'N procs'

      call create_hdf5_file(target_meshFile_h5_full_name,targetMesh_hdf5_file_id)
      
      call create_groups_datasets_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,targetMesh_hdf5_file_id,isMeshLinealOutput,.false.,&
                                                target_Nprocs,totalNumElements,numNodesTrgtTotal_i8,mesh_VTKnnode,mesh_numVTKElemsPerMshElem)

      do iTrgtRank=1,numTrgtRanksInMpiRank
         trgtRank = trgtRanksInMpiRank(iTrgtRank)

        call write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(mporder,mnnode,targetMesh_hdf5_file_id,evalMeshQuality,trgtRank,target_Nprocs,&
            numElemsTrgtRank(iTrgtRank),numElemsVTKTrgtRank(iTrgtRank),sizeConnecVTKTrgtRank(iTrgtRank),mesh_VTKnnode,mesh_numVTKElemsPerMshElem,&
            trgtRankElemStart(iTrgtRank),trgtRankElemEnd(iTrgtRank),trgtRankNodeStart_i8(iTrgtRank),trgtRankNodeEnd_i8(iTrgtRank),numNodesTrgtRank(iTrgtRank),&
            coordVTKTrgtRank_jm%matrix(iTrgtRank)%elems,connecVTKTrgtRank_jv%vector(iTrgtRank)%elems,quality_jm%matrix(iTrgtRank)%elems,connecChunkSize)
      end do

      do iTrgtRank=(numTrgtRanksInMpiRank+1),maxNumTrgtRanks
        call dummy_write_mshRank_data_vtkhdf_unstructuredGrid_meshFile(targetMesh_hdf5_file_id,evalMeshQuality,target_Nprocs,connecChunkSize)
      end do

      call close_hdf5_file(targetMesh_hdf5_file_id)
      !----------------------------------------------------------------------------------------------------------------------------------     

   end subroutine generate_new_mesh_for_viz

   subroutine open_target_mesh_and_do_mapping(target_meshFile_h5_full_name,target_Nprocs,mapNodeRankTrgt,numNodesMapTrgtRank,&
                                             numNodesTrgtTotal_i8,trgtRankNodeStart_i8,trgtRankNodeEnd_i8)
      implicit none
      character(*),intent(in) :: target_meshFile_h5_full_name
      integer(4),intent(in) :: target_Nprocs
      integer(4),intent(inout) :: mapNodeRankTrgt(numNodesRankPar,2),numNodesMapTrgtRank(0:target_Nprocs-1)
      integer(8),intent(inout) :: numNodesTrgtTotal_i8,trgtRankNodeStart_i8(target_Nprocs),trgtRankNodeEnd_i8(target_Nprocs)
      !-------------------------------------------------------
      integer(4) :: numNodesTrgtRank(target_Nprocs),iNodeL
      integer(8) :: iNodeG

      integer(4) :: iTrgtRank,trgtRank,iPos,h5err
      integer(4),allocatable :: aux_array_i4(:)
      integer(8),allocatable :: aux_array_i8(:),globalIdSrlOrdered_i8(:,:)

      character(128) :: dsetname
      integer(hid_t) :: targetMesh_hdf5_file_id,dset_id,fspace_id
      integer(hsize_t), dimension(1) :: ms_dims,fs_dims,fs_maxdims
      integer(hssize_t), dimension(1) :: ms_offset

      if(mpi_rank.eq.0) write(*,*) '# Opening target mesh file:',trim(target_meshFile_h5_full_name)

      mapNodeRankTrgt(:,:) = 0
      numNodesMapTrgtRank(:) = 0

      call open_hdf5_file(target_meshFile_h5_full_name,targetMesh_hdf5_file_id)

      !------------------------------------------------------------------------
      dsetname = '/globalIds/globalIdSrl'
      call h5dopen_f(targetMesh_hdf5_file_id,dsetname,dset_id,h5err)
      call h5dget_space_f(dset_id,fspace_id,h5err)
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)
      call h5sclose_f(fspace_id,h5err)
      call h5dclose_f(dset_id,h5err)

      numNodesTrgtTotal_i8 = fs_dims(1)
      !------------------------------------------------------------------------

      do iTrgtRank=1,target_Nprocs
         trgtRank = (iTrgtRank-1)

         ms_dims(1) = 1
         ms_offset(1) = int(trgtRank,hssize_t)

         allocate(aux_array_i8(1))

         dsetname = '/Parallel_data/rankNodeStart'
         call read_dataspace_1d_int8_hyperslab_parallel(targetMesh_hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
         trgtRankNodeStart_i8(iTrgtRank)=aux_array_i8(1)

         dsetname = '/Parallel_data/rankNodeEnd'
         call read_dataspace_1d_int8_hyperslab_parallel(targetMesh_hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)
         trgtRankNodeEnd_i8(iTrgtRank)=aux_array_i8(1)

         numNodesTrgtRank(iTrgtRank) = trgtRankNodeEnd_i8(iTrgtRank) - trgtRankNodeStart_i8(iTrgtRank) + 1

         deallocate(aux_array_i8)

         !---------------------------------------------------------------------------------------------------------------

         ms_dims(1) = int(numNodesTrgtRank(iTrgtRank),hsize_t)
         ms_offset(1) = int(trgtRankNodeStart_i8(iTrgtRank),hssize_t)-1

         allocate(aux_array_i8(numNodesTrgtRank(iTrgtRank)))
         allocate(globalIdSrlOrdered_i8(numNodesTrgtRank(iTrgtRank),2))

         dsetname = '/globalIds/globalIdSrl'
         call read_dataspace_1d_int8_hyperslab_parallel(targetMesh_hdf5_file_id,dsetname,ms_dims,ms_offset,aux_array_i8)

         do iNodeL=1,numNodesTrgtRank(iTrgtRank)
            globalIdSrlOrdered_i8(iNodeL,1) = aux_array_i8(iNodeL)
            globalIdSrlOrdered_i8(iNodeL,2) = iNodeL
         end do

         call quicksort_matrix_int8(globalIdSrlOrdered_i8,1)

         deallocate(aux_array_i8)

         do iNodeL=1,numNodesRankPar
            iNodeG = globalIdSrl(iNodeL)

            iPos = binarySearch_int_i8(globalIdSrlOrdered_i8(:,1),iNodeG)
            if(iPos.ne.0) then
               mapNodeRankTrgt(iNodeL,1) = globalIdSrlOrdered_i8(iPos,2)
               mapNodeRankTrgt(iNodeL,2) = trgtRank
               !if(mpi_rank.eq.0) write(*,*) '[',mpi_rank,']iNodeL',iNodeL,'iNodeG',iNodeG,'iPos',iPos,'trgtRank',trgtRank,'numNiTR',numNodesMapTrgtRank(iTrgtRank)
            end if
         end do

         deallocate(globalIdSrlOrdered_i8)

      end do

      do iNodeL=1,numNodesRankPar
         trgtRank = mapNodeRankTrgt(iNodeL,2)
         numNodesMapTrgtRank(trgtRank) = numNodesMapTrgtRank(trgtRank) + 1
      end do

      !write(*,*) '[',mpi_rank,']numNodesMapTrgtRank:',numNodesMapTrgtRank(:),'numNodesRankPar',numNodesRankPar
      !do iNodeL=1,numNodesRankPar
      !   if(mpi_rank.eq.0) write(*,*) '[',mpi_rank,']iNodeL',iNodeL,'->iNL',mapNodeRankTrgt(iNodeL,1),'iTR',mapNodeRankTrgt(iNodeL,2)
      !end do

      call close_hdf5_file(targetMesh_hdf5_file_id)

   end subroutine open_target_mesh_and_do_mapping

   subroutine save_mapping_file(mappingFile_h5_full_name,mapNodeRankTrgt)
      implicit none
      character(512),intent(in) :: mappingFile_h5_full_name
      integer(4),intent(in) :: mapNodeRankTrgt(numNodesRankPar,2)
      integer(4),parameter :: ds_rank=1
      integer(hid_t) :: mapping_hdf5_file_id
      character(128) :: dsetname
      integer(hid_t) :: dtype
      integer(hsize_t),dimension(ds_rank) :: ms_dims,ds_dims
      integer(hssize_t),dimension(ds_rank) :: ms_offset
      !------------------------------------------------------------------------------------------------

      if(mpi_rank.eq.0) write(*,*) 'Saving mapping in: ',trim(mappingFile_h5_full_name)

      call create_hdf5_file(mappingFile_h5_full_name,mapping_hdf5_file_id)

      !-----------------------------------------------------------------------------------------------
      dtype = h5_datatype_int4
      ds_dims(1) = int(totalNumNodesPar,hsize_t)
      ms_dims(1) = int(numNodesRankPar,hsize_t)
      ms_offset(1) = int(rankNodeStart,hssize_t)-1
      !-----------------------------------------------------------------------------------------------

      dsetname = '/mappingNode'
      call create_dataspace_hdf5(mapping_hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(mapping_hdf5_file_id,dsetname,ms_dims,ms_offset,mapNodeRankTrgt(:,1))

      dsetname = '/mappingRank'
      call create_dataspace_hdf5(mapping_hdf5_file_id,dsetname,ds_rank,ds_dims,dtype)
      call write_dataspace_1d_int4_hyperslab_parallel(mapping_hdf5_file_id,dsetname,ms_dims,ms_offset,mapNodeRankTrgt(:,2))

      call close_hdf5_file(mapping_hdf5_file_id)

   end subroutine

   subroutine set_hdf5_mappingFile_name(file_path,file_name,numSrcRanks,numTrgtRanks,mappingFile_h5_full_name)
      implicit none
      character(len=*),intent(in) :: file_path,file_name
      integer,intent(in) :: numSrcRanks,numTrgtRanks
      character(len=*),intent(out) :: mappingFile_h5_full_name
      character(len=12) :: aux_srcRanks,aux_trgtRanks

      write(aux_srcRanks,'(I0)') numSrcRanks
      write(aux_trgtRanks,'(I0)') numTrgtRanks
      mappingFile_h5_full_name = trim(adjustl(file_path))//trim(adjustl(file_name))//'-mapping-'//trim(aux_srcRanks)//'-to-'//trim(aux_trgtRanks)//'.hdf'
   end subroutine set_hdf5_mappingFile_name

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
      real(8),allocatable,intent(inout) :: coordNodesTrgtRank(:,:)

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
            coordNodesTrgtRank(nodeCnt,:) = real(coordPar(iNodeMpi,:),8)
            jAux = rawNodeListTrgtRank(iAux)
         end if
      end do

      deallocate(rawNodeListTrgtRank)
      !write(*,*) 'listNodesTrgtRank',listNodesTrgtRank(:)

   end subroutine generate_connec_and_coords_for_targetRank


   subroutine reorder_nodes_in_trgtRank(mporder,mnnode,trgtRank,target_Nprocs,&
               numElemsTrgtRank,numNodesTrgtRank,sizeConnecVTKTrgtRank,linOutput,listNodesTrgtRankMpi,connecTrgtRankMpi,iNodeStartPar_i8,&
               mapNodeIdNewToOld,connecVTKTrgtRank,connecParOrigTrgtRank,coord,coordVTK)
      implicit none
      integer(4),intent(in) :: mporder,mnnode,trgtRank,target_Nprocs,numElemsTrgtRank,numNodesTrgtRank,sizeConnecVTKTrgtRank
      logical,intent(in) :: linOutput
      integer(4),intent(in) :: listNodesTrgtRankMpi(numNodesTrgtRank),connecTrgtRankMpi(numElemsTrgtRank,mnnode)
      integer(8),dimension(0:target_Nprocs-1),intent(in) :: iNodeStartPar_i8

      integer(4),intent(out) :: mapNodeIdNewToOld(numNodesTrgtRank)
      integer(4),intent(out) :: connecVTKTrgtRank(sizeConnecVTKTrgtRank),connecParOrigTrgtRank(numElemsTrgtRank,mnnode)      
      real(8),intent(in) :: coord(numNodesTrgtRank,ndime)
      real(8),intent(out) :: coordVTK(numNodesTrgtRank,ndime)


      integer(4) :: isNodeAdded(numNodesTrgtRank),auxVTKorder(mnnode)
      integer(4) :: m,indConn,indexIJK,indexGMSH,indexNew,nodeIndexCnt,indPosListNodes
      integer(4) :: iNodeL,iNodeMpi,iElem
      integer(8) :: iNodeGsrl,iNodeGpar
      integer(4) :: ii,jj,kk,igp,jgp,kgp,indexVTK
      integer(4),allocatable :: ijk_sod2d_to_gmsh(:),ijk_gmsh_to_sod2d(:)

      call set_allocate_array_ijk_sod2d_criteria(mporder,ijk_sod2d_to_gmsh,ijk_gmsh_to_sod2d)

      !$acc kernels
      isNodeAdded(:)=-1
      !globalIdSrlOrderedTrgtRank_i8(:,:)=-1
      !$acc end kernels
      nodeIndexCnt = 0
      indConn = -1

      !----------------------------------------------------------------------------------------------------------
      do kk = 0,mporder
         do ii = 0,mporder
            do jj = 0,mporder

               indexIJK = ((mporder+1)**2)*kk+(mporder+1)*ii+jj+1
               auxVTKorder(mesh_a2ijk(indexIJK)) = mesh_vtk2ijk(indexIJK)

            end do
         end do
      end do

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
               indConn = (iElem-1)*mnnode + auxVTKorder(m)
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

   subroutine copy_dsets_results_for_generated_mesh(source_file_id,target_file_id,numTrgtRanksInMpiRank,maxNumTrgtRanks,&
                           numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv,&
                           numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
                           dsetsScalarFieldsOrig,dsetsVectorFieldsOrig,dsetsV2SFieldsOrig,dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)
      implicit none
      integer(hid_t),intent(in) :: source_file_id,target_file_id
      integer(4),intent(in) :: numTrgtRanksInMpiRank,maxNumTrgtRanks
      integer(4),intent(in) :: numNodesTrgtRank(numTrgtRanksInMpiRank)
      integer(8),intent(in) :: trgtRankNodeStart(numTrgtRanksInMpiRank),numNodesTrgtTotal
      type(jagged_vector_int4),intent(in) :: mapNodeIdTrgtToMpi_jv
      integer(4),intent(in) :: numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S
      character(256),intent(in) :: dsetsScalarFieldsOrig(maxNumDsetSca),dsetsVectorFieldsOrig(maxNumDsetVec),dsetsV2SFieldsOrig(maxNumDsetV2S)
      character(256),intent(in) :: dsetsScalarFieldsTrgt(maxNumDsetSca),dsetsVectorFieldsTrgt(maxNumDsetVec),dsetsV2SFieldsTrgt(maxNumDsetV2S,ndime)

      real(rp_vtk) :: sourceScalarField(numNodesRankPar),sourceVectorField(numNodesRankPar,ndime)
      integer(4) :: iTrgtRank,iSca,iVec,iDim
      character(512) :: dsetnameOrig,dsetnameTrgt

      do iSca=1,numDsetSca
         dsetnameOrig = dsetsScalarFieldsOrig(iSca)
         dsetnameTrgt = dsetsScalarFieldsTrgt(iSca)
         if(mpi_rank.eq.0) write(*,*) ' # Copying scalar field ',trim(adjustl(dsetnameOrig)),' (id',iSca,')...'
         call read_and_load_source_scalarfield(dsetnameOrig,source_file_id,sourceScalarField)
         call copy_scalarfield_dataset_results_in_target(dsetnameTrgt,sourceScalarField,target_file_id,numTrgtRanksInMpiRank,&
                  maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
      end do

      do iVec=1,numDsetVec
         dsetnameOrig = dsetsVectorFieldsOrig(iVec)
         dsetnameTrgt = dsetsVectorFieldsTrgt(iVec)
         if(mpi_rank.eq.0) write(*,*) ' # Copying vector field ',trim(adjustl(dsetnameOrig)),' (id',iVec,')...'
         call read_and_load_source_vectorfield(dsetnameOrig,source_file_id,sourceVectorField)
         call copy_vectorfield_dataset_results_in_target(dsetnameTrgt,sourceVectorField,target_file_id,numTrgtRanksInMpiRank,&
                  maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
      end do

      do iVec=1,numDsetV2S
         dsetnameOrig = dsetsV2SFieldsOrig(iVec)
         if(mpi_rank.eq.0) write(*,*) ' # Copying vector field ',trim(adjustl(dsetnameOrig)),' (id',iVec,')...'
         call read_and_load_source_vectorfield(dsetnameOrig,source_file_id,sourceVectorField)
         do iDim=1,ndime
            dsetnameTrgt = dsetsV2SFieldsTrgt(iVec,iDim)
            if(mpi_rank.eq.0) write(*,*) '  - in scalar field ',trim(adjustl(dsetnameOrig)),' (dim',iDim,')...'
            call copy_scalarfield_dataset_results_in_target(dsetnameTrgt,sourceVectorField(:,iDim),target_file_id,numTrgtRanksInMpiRank,&
                  maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
         end do
      end do

   end subroutine copy_dsets_results_for_generated_mesh

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

   subroutine copy_dsets_results_for_loaded_mesh(source_file_id,target_file_id,target_Nprocs,&
                           numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeRankTrgt,&
                           numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
                           dsetsScalarFieldsOrig,dsetsVectorFieldsOrig,dsetsV2SFieldsOrig,dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)
      implicit none
      integer(hid_t),intent(in) :: source_file_id,target_file_id
      integer(4),intent(in) :: target_Nprocs,numNodesTrgtRank(0:target_Nprocs-1),mapNodeRankTrgt(numNodesRankPar,2)
      integer(8),intent(in) :: trgtRankNodeStart(target_Nprocs),numNodesTrgtTotal
      integer(4),intent(in) :: numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S
      character(256),intent(in) :: dsetsScalarFieldsOrig(maxNumDsetSca),dsetsVectorFieldsOrig(maxNumDsetVec),dsetsV2SFieldsOrig(maxNumDsetV2S)
      character(256),intent(in) :: dsetsScalarFieldsTrgt(maxNumDsetSca),dsetsVectorFieldsTrgt(maxNumDsetVec),dsetsV2SFieldsTrgt(maxNumDsetV2S,ndime)

      real(rp_vtk) :: sourceScalarField(numNodesRankPar),sourceVectorField(numNodesRankPar,ndime)
      integer(4) :: iTrgtRank,iSca,iVec,iDim
      character(512) :: dsetnameOrig,dsetnameTrgt

      do iSca=1,numDsetSca
         dsetnameOrig = dsetsScalarFieldsOrig(iSca)
         dsetnameTrgt = dsetsScalarFieldsTrgt(iSca)
         if(mpi_rank.eq.0) write(*,*) ' # Copying scalar field ',trim(adjustl(dsetnameOrig)),' (id',iSca,')...'
         call read_and_load_source_scalarfield(dsetnameOrig,source_file_id,sourceScalarField)
         call copy_scalarfield_dataset_results_for_loaded_mesh(dsetnameTrgt,sourceScalarField,target_file_id,&
                  target_Nprocs,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeRankTrgt)
      end do

      do iVec=1,numDsetVec
         dsetnameOrig = dsetsVectorFieldsOrig(iVec)
         dsetnameTrgt = dsetsVectorFieldsTrgt(iVec)
         if(mpi_rank.eq.0) write(*,*) ' # Copying vector field ',trim(adjustl(dsetnameOrig)),' (id',iVec,')...'
         call read_and_load_source_vectorfield(dsetnameOrig,source_file_id,sourceVectorField)
         !call copy_vectorfield_dataset_results_in_target(dsetnameTrgt,sourceVectorField,target_file_id,numTrgtRanksInMpiRank,&
         !         maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
      end do

      do iVec=1,numDsetV2S
         dsetnameOrig = dsetsV2SFieldsOrig(iVec)
         if(mpi_rank.eq.0) write(*,*) ' # Copying vector field ',trim(adjustl(dsetnameOrig)),' (id',iVec,')...'
         call read_and_load_source_vectorfield(dsetnameOrig,source_file_id,sourceVectorField)
         !do iDim=1,ndime
         !   dsetnameTrgt = dsetsV2SFieldsTrgt(iVec,iDim)
         !   if(mpi_rank.eq.0) write(*,*) '  - in scalar field ',trim(adjustl(dsetnameOrig)),' (dim',iDim,')...'
         !   call copy_scalarfield_dataset_results_in_target(dsetnameTrgt,sourceVectorField(:,iDim),target_file_id,numTrgtRanksInMpiRank,&
         !         maxNumTrgtRanks,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv)
         !end do
      end do

   end subroutine copy_dsets_results_for_loaded_mesh

   subroutine copy_scalarfield_dataset_results_for_loaded_mesh(dsetname,sourceScalarField,target_file_id,&
                                             target_Nprocs,numNodesTrgtRank,trgtRankNodeStart,numNodesTrgtTotal,mapNodeRankTrgt)
      implicit none
      character(len=*),intent(in) :: dsetname
      integer(hid_t),intent(in) :: target_file_id
      integer(4),intent(in) :: target_Nprocs,numNodesTrgtRank(0:target_Nprocs-1),mapNodeRankTrgt(numNodesRankPar,2)
      integer(8),intent(in) :: trgtRankNodeStart(target_Nprocs),numNodesTrgtTotal
      real(rp_vtk),intent(in) :: sourceScalarField(numNodesRankPar)

      integer(4) :: iNodeSrc,iNodeTrgt,iRankTrgt,iNodeCnt,trgtRank,numNodesInTrgt,minNumNodesInTrgt
      integer(8) :: trgtRankOffset

      real(rp_vtk),allocatable :: targetScalarField(:)

      integer(hsize_t),allocatable :: ms_coords(:,:)
      integer(4),parameter :: ms_rank=1
      integer(4) :: h5err
      integer(hid_t) :: dset_id,fspace_id,mspace_id,plist_id,dtype
      integer(hsize_t) :: ms_numElems,fs_dims(ms_rank),fs_maxdims(ms_rank),ms_dims(ms_rank)
      integer(hssize_t) :: ms_offset(ms_rank)
	   !--------------------------------------------------------------------------------

      do trgtRank=0,target_Nprocs-1

         numNodesInTrgt = numNodesTrgtRank(trgtRank)
         !numNodesInTrgt = max(1,numNodesTrgtRank(trgtRank))
         minNumNodesInTrgt = max(1,numNodesInTrgt)

         write(*,*) '[',mpi_rank,']trgtRank',trgtRank,'numNodesInTrgt',numNodesInTrgt,'minNumNodesInTrgt',minNumNodesInTrgt
         allocate(targetScalarField(numNodesInTrgt))
         allocate(ms_coords(ms_rank,numNodesInTrgt))

         ms_numElems = int(numNodesInTrgt,hsize_t)
         ms_dims(1) = int(numNodesInTrgt,hsize_t)

         targetScalarField(:) = 0
         ms_coords(:,:) = 1

         iNodeCnt = 0
         do iNodeSrc=1,numNodesRankPar
            if(mapNodeRankTrgt(iNodeSrc,2).eq.trgtRank) then
               iNodeTrgt = mapNodeRankTrgt(iNodeSrc,1)
               iRankTrgt = mapNodeRankTrgt(iNodeSrc,2)
               trgtRankOffset = trgtRankNodeStart(iRankTrgt+1) - 1
               
               iNodeCnt = iNodeCnt + 1

               targetScalarField(iNodeCnt) = sourceScalarField(iNodeSrc)

               ms_coords(1,iNodeCnt) = iNodeTrgt + trgtRankOffset !iNodeTrgt
               !if(iNodeSrc.ge.1000 .and. iNodeSrc.le.(numNodesRankPar-1000)) then
               !   ms_coords(1,iNodeCnt) = iNodeTrgt + trgtRankOffset !iNodeTrgt
               !endif
               !if(iNodeTrgt.le.0) write(*,*) '[',mpi_rank,']iNodeTrgt',iNodeTrgt,'iNodeSrc',iNodeSrc,'trgtRank',trgtRank
               !if(iNodeTrgt.gt.7000)write(*,*) '[',mpi_rank,']iNodeTrgt',iNodeTrgt,'iNodeSrc',iNodeSrc,'trgtRank',trgtRank
               !if(mpi_rank.eq.2)write(*,*) '[',mpi_rank,']iNodeTrgt',iNodeTrgt,'iNodeSrc',iNodeSrc,'trgtRank',trgtRank
            end if
         end do
         
      !--------------------------------------------------------------------

      call select_dtype_rp_vtk(dtype)

      !------------------------------------------------------------------------------------------------------
      !call open_create_dataspace_hyperslab_parallel(file_id,dsetname,ms_rank,ms_dims,ms_offset,&
      !                                   dset_id,fspace_id,mspace_id,plist_id,fs_dims,fs_maxdims)
      call h5dopen_f(target_file_id, dsetname, dset_id, h5err)

      !get filespace of the dataset
      call h5dget_space_f(dset_id, fspace_id, h5err)

      !get dimensions of the filespace
      call h5sget_simple_extent_dims_f(fspace_id,fs_dims,fs_maxdims,h5err)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5screate_simple_f(ms_rank,ms_dims,mspace_id,h5err)

      if(numNodesInTrgt.gt.0) then
         call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err)
      else
         ms_offset(1) = 0
         call h5sselect_hyperslab_f(fspace_id,H5S_SELECT_SET_F,ms_offset,ms_dims,h5err)
      endif

      ! Create property list for collective dataset write
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5err)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F,h5err)
      !------------------------------------------------------------------------------------------------------


      call h5dwrite_f(dset_id,dtype,targetScalarField,ms_dims,h5err,file_space_id=fspace_id,mem_space_id=mspace_id,xfer_prp=plist_id)

      call close_dataspace_hyperslab_parallel(dset_id,fspace_id,mspace_id,plist_id)


!--------------------




         !call h5sselect_elements_f(fspace_id,H5S_SELECT_SET_F,ms_rank,ms_numElems,ms_coords,h5err)
         !call h5dwrite_f(dset_id, dtype, targetScalarField, h5err, mem_space_id=mspace_id, file_space_id=fspace_id, xfer_prp=plist_id)


         !if(mpi_rank.eq.0) write(*,*) '[',mpi_rank,']numNodesInTrgt(',trgtRank,')',numNodesInTrgt,'iNodeCnt',iNodeCnt,'numNodesRankPar',numNodesRankPar

      !   call copy_scalarfield_result_in_trgtRank(dsetname,sourceScalarField,target_file_id,numNodesTrgtRank(iTrgtRank),&
      !         trgtRankNodeStart(iTrgtRank),numNodesTrgtTotal,mapNodeIdTrgtToMpi_jv%vector(iTrgtRank)%elems)

         deallocate(targetScalarField,ms_coords)
      end do


   end subroutine copy_scalarfield_dataset_results_for_loaded_mesh

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

   subroutine generate_dsets_to_copy(type_resultsFile,src_file_id,numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,&
            dsetsScalarFieldsOrig,dsetsVectorFieldsOrig,dsetsV2SFieldsOrig,dsetsScalarFieldsTrgt,dsetsVectorFieldsTrgt,dsetsV2SFieldsTrgt)
      implicit none
      integer(4),intent(in) :: type_resultsFile
      integer(hid_t),intent(in) :: src_file_id
      integer(4),intent(inout) :: numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S
      character(256),intent(inout),allocatable :: dsetsScalarFieldsOrig(:),dsetsVectorFieldsOrig(:),dsetsV2SFieldsOrig(:)
      character(256),intent(inout),allocatable :: dsetsScalarFieldsTrgt(:),dsetsVectorFieldsTrgt(:),dsetsV2SFieldsTrgt(:,:)
      character(512) :: groupname,dsetname
      integer(4) :: ii,h5err
      logical :: dset_exists
   
      groupname = '/VTKHDF/PointData/'
      numDsetSca = 0
      numDsetVec = 0
      numDsetV2S = 0
      maxNumDsetSca = 0
      maxNumDsetVec = 0
      maxNumDsetV2S = 0

      if(type_resultsFile .eq. 1) then
         ! INST RESULTS ------------------------------------------------------------------------------
         maxNumDsetSca = numNodeScalarFields+numElGPScalarFields
         maxNumDsetVec = numNodeVectorFields

         allocate(dsetsScalarFieldsOrig(maxNumDsetSca))
         allocate(dsetsVectorFieldsOrig(maxNumDsetVec))
         allocate(dsetsScalarFieldsTrgt(maxNumDsetSca))
         allocate(dsetsVectorFieldsTrgt(maxNumDsetVec))
         allocate(dsetsV2SFieldsOrig(maxNumDsetV2S))
         allocate(dsetsV2SFieldsTrgt(maxNumDsetV2S,ndime))

         do ii=1,numNodeScalarFields
            dsetname = trim(adjustl(groupname))//trim(adjustl(nodeScalarNameFields(ii)))
            call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
            if(dset_exists) then
               numDsetSca = numDsetSca + 1
               dsetsScalarFieldsOrig(numDsetSca) = dsetname
               dsetsScalarFieldsTrgt(numDsetSca) = dsetname
            end if
         end do

         do ii=1,numElGPScalarFields
            dsetname = trim(adjustl(groupname))//trim(adjustl(elGPScalarNameFields(ii)))
            call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
            if(dset_exists) then
               numDsetSca = numDsetSca + 1
               dsetsScalarFieldsOrig(numDsetSca) = dsetname
               dsetsScalarFieldsTrgt(numDsetSca) = dsetname
            end if
         end do

         do ii=1,numNodeVectorFields
            dsetname = trim(adjustl(groupname))//trim(adjustl(nodeVectorNameFields(ii)))
            call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
            if(dset_exists) then
               numDsetVec = numDsetVec + 1
               dsetsVectorFieldsOrig(numDsetVec) = dsetname
               dsetsVectorFieldsTrgt(numDsetVec) = dsetname
            end if
         end do
         !--------------------------------------------------------------------------------------------
      else if(type_resultsFile .eq. 2) then
         ! AVG RESULTS -------------------------------------------------------------------------------
         maxNumDsetSca = numAvgNodeScalarFields
         maxNumDsetVec = numAvgNodeVectorFields

         allocate(dsetsScalarFieldsOrig(maxNumDsetSca))
         allocate(dsetsScalarFieldsTrgt(maxNumDsetSca))
         allocate(dsetsVectorFieldsOrig(maxNumDsetVec))
         allocate(dsetsVectorFieldsTrgt(maxNumDsetVec))
         allocate(dsetsV2SFieldsOrig(maxNumDsetV2S))
         allocate(dsetsV2SFieldsTrgt(maxNumDsetV2S,ndime))

         do ii=1,numAvgNodeScalarFields
            dsetname = trim(adjustl(groupname))//trim(adjustl(avgNodeScalarNameFields(ii)))
            call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
            if(dset_exists) then
               numDsetSca = numDsetSca + 1
               dsetsScalarFieldsOrig(numDsetSca) = dsetname
               dsetsScalarFieldsTrgt(numDsetSca) = dsetname
            end if
         end do

         do ii=1,numAvgNodeVectorFields
            dsetname = trim(adjustl(groupname))//trim(adjustl(avgNodeVectorNameFields(ii)))
            call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
            if(dset_exists) then
               numDsetVec = numDsetVec + 1
               dsetsVectorFieldsOrig(numDsetVec) = dsetname
               dsetsVectorFieldsTrgt(numDsetVec) = dsetname
            end if
         end do
         !--------------------------------------------------------------------------------------------

      else if(type_resultsFile .eq. 3) then
         ! RESTART FILES
         maxNumDsetSca = numRestartFields

         allocate(dsetsScalarFieldsOrig(maxNumDsetSca))
         allocate(dsetsScalarFieldsTrgt(maxNumDsetSca))
         allocate(dsetsVectorFieldsOrig(maxNumDsetVec))
         allocate(dsetsVectorFieldsTrgt(maxNumDsetVec))
         allocate(dsetsV2SFieldsOrig(maxNumDsetV2S))
         allocate(dsetsV2SFieldsTrgt(maxNumDsetV2S,ndime))

         do ii=1,numRestartFields
            dsetname = trim(adjustl(restartNameFields(ii)))
            call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
            if(dset_exists) then
               numDsetSca = numDsetSca + 1
               dsetsScalarFieldsOrig(numDsetSca) = dsetname
               dsetsScalarFieldsTrgt(numDsetSca) = dsetname
            end if
         end do
      else if(type_resultsFile .eq. 4) then
         ! INST TO RESTART FILES
         maxNumDsetSca = 5 !rho,pr,ener,mut,mue
         maxNumDsetV2S = 1 !u_x,u_y,u_z

         allocate(dsetsScalarFieldsOrig(maxNumDsetSca))
         allocate(dsetsScalarFieldsTrgt(maxNumDsetSca))
         allocate(dsetsVectorFieldsOrig(maxNumDsetVec))
         allocate(dsetsVectorFieldsTrgt(maxNumDsetVec))
         allocate(dsetsV2SFieldsOrig(maxNumDsetV2S))
         allocate(dsetsV2SFieldsTrgt(maxNumDsetV2S,ndime))

         !RHO
         numDsetSca = numDsetSca + 1
         dsetname = trim(adjustl(groupname))//trim(adjustl(nodeScalarNameFields(indNS_rho)))
         dsetsScalarFieldsOrig(numDsetSca) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_rho)))
         dsetsScalarFieldsTrgt(numDsetSca) = dsetname

         !PR
         numDsetSca = numDsetSca + 1
         dsetname = trim(adjustl(groupname))//trim(adjustl(nodeScalarNameFields(indNS_pr)))
         dsetsScalarFieldsOrig(numDsetSca) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_pr)))
         dsetsScalarFieldsTrgt(numDsetSca) = dsetname

         !ENER
         numDsetSca = numDsetSca + 1
         dsetname = trim(adjustl(groupname))//trim(adjustl(nodeScalarNameFields(indNS_ener)))
         dsetsScalarFieldsOrig(numDsetSca) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_ener)))
         dsetsScalarFieldsTrgt(numDsetSca) = dsetname

         !MUT
         numDsetSca = numDsetSca + 1
         dsetname = trim(adjustl(groupname))//trim(adjustl(elGPScalarNameFields(indES_mut)))
         dsetsScalarFieldsOrig(numDsetSca) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_mut)))
         dsetsScalarFieldsTrgt(numDsetSca) = dsetname

         !MUE
         numDsetSca = numDsetSca + 1
         dsetname = trim(adjustl(groupname))//trim(adjustl(elGPScalarNameFields(indES_mue)))
         dsetsScalarFieldsOrig(numDsetSca) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_mue)))
         dsetsScalarFieldsTrgt(numDsetSca) = dsetname

         !VEL (ux,uy,uz)
         numDsetV2S = numDsetV2S + 1
         dsetname = trim(adjustl(groupname))//trim(adjustl(nodeVectorNameFields(indNV_vel)))
         dsetsV2SFieldsOrig(numDsetV2S) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_ux)))
         dsetsV2SFieldsTrgt(numDsetV2S,1) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_uy)))
         dsetsV2SFieldsTrgt(numDsetV2S,2) = dsetname
         dsetname = trim(adjustl(restartNameFields(indRF_uz)))
         dsetsV2SFieldsTrgt(numDsetV2S,3) = dsetname

         !!!do ii=1,numRestartFields
         !!!   !dsetname = trim(adjustl(groupname))//trim(adjustl(nodeScalarNameFields(ii)))
         !!!   !dsetname = trim(adjustl(restartNameFields(ii)))
         !!!   call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
         !!!   if(dset_exists) then
         !!!      numDsetSca = numDsetSca + 1
         !!!      dsetsScalarFieldsOrig(numDsetSca) = dsetname
         !!!      dsetsScalarFieldsTrgt(numDsetSca) = dsetname
         !!!   end if
         !!!end do

      end if

   end subroutine generate_dsets_to_copy

#if 0
   subroutine check_dsets_in_source_results_file(src_file_id,dsetScaCnt,dsetVecCnt,dsetsScalarFields,dsetsVectorFields,&
                                                   dsetsScaExist,dsetsVecExist)
      implicit none
      integer(hid_t),intent(in) :: src_file_id
      integer(4),intent(in) :: dsetScaCnt,dsetVecCnt
      character(256),intent(out) :: dsetsScalarFields(max_num_saved_fields),dsetsVectorFields(max_num_saved_fields)
      logical,intent(out) :: dsetsScaExist(max_num_saved_fields),dsetsVecExist(max_num_saved_fields)
      integer(4) :: h5err
      character(512) :: groupname,dsetname
      logical :: dset_exists

      integer(4) :: iSca,iVec
	   !--------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------------------------------------------
#if 0


      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'rho'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'mu_fluid'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'pr'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'E'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'eta'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'csound'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'machno'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'divU'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'qcrit'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'mut'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'mue'
      dsetScaCnt=dsetScaCnt+1
      dsetsScalarFields(dsetScaCnt) = dsetname
      !--------------------------------------------------
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'u'
      dsetVecCnt=dsetVecCnt+1
      dsetsVectorFields(dsetVecCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'gradRho'
      dsetVecCnt=dsetVecCnt+1
      dsetsVectorFields(dsetVecCnt) = dsetname
      !--------------------------------------------------
      dsetname = trim(adjustl(groupname))//'curlU'
      dsetVecCnt=dsetVecCnt+1
      dsetsVectorFields(dsetVecCnt) = dsetname
      !--------------------------------------------------
#endif
      !--------------------------------------------------------------------------------
      dsetsScaExist(:) = .false.
      dsetsVecExist(:) = .false.

      do iSca=1,dsetScaCnt
         dsetname = dsetsScalarFields(iSca)
         call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
         dsetsScaExist(iSca) = dset_exists
         !write(*,*) 'checking if scaf',dsetname,'exists:',dset_exists
      end do

      do iVec=1,dsetVecCnt
         dsetname = dsetsVectorFields(iVec)
         call h5lexists_f(src_file_id,dsetname,dset_exists,h5err)
         dsetsVecExist(iVec) = dset_exists
         !write(*,*) 'checking if vec',dsetname,'exists:',dset_exists
      end do

   end subroutine check_dsets_in_source_results_file
#endif

   subroutine create_fields_dsets_in_new_results_file(trgt_file_id,numNodesParTotal_i8,numDsetSca,numDsetVec,numDsetV2S,&
               maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S,dsetsScalarFields,dsetsVectorFields,dsetsV2SFields)
      implicit none
      integer(hid_t),intent(in) :: trgt_file_id
      integer(8),intent(in) :: numNodesParTotal_i8
      integer(4),intent(in) :: numDsetSca,numDsetVec,numDsetV2S,maxNumDsetSca,maxNumDsetVec,maxNumDsetV2S
      character(256),intent(in) :: dsetsScalarFields(maxNumDsetSca),dsetsVectorFields(maxNumDsetVec),dsetsV2SFields(maxNumDsetV2S,ndime)
      integer(hid_t) :: dtype
      integer(hsize_t) :: ds_dims(1),ds_dims2d(2),aux_ds_dims
      integer(4) :: ds_rank,h5err,iSca,iVec,iDim
      character(512) :: dsetname
      !--------------------------------------------------------------------------------

      call select_dtype_rp_vtk(dtype)

      !----------------------------------------------------------------------
      ds_rank = 1
      ds_dims(1) = numNodesParTotal_i8

      do iSca=1,numDsetSca
         dsetname = dsetsScalarFields(iSca)
         call create_dataspace_hdf5(trgt_file_id,dsetname,ds_rank,ds_dims,dtype)
      end do

      do iSca=1,numDsetV2S
         do iDim=1,ndime
            dsetname = dsetsV2SFields(iSca,iDim)
            call create_dataspace_hdf5(trgt_file_id,dsetname,ds_rank,ds_dims,dtype)
         end do
      end do

      !----------------------------------------------------------------------
      ds_rank = 2
      ds_dims2d(1) = ndime
      ds_dims2d(2) = numNodesParTotal_i8

      do iVec=1,numDsetVec
         dsetname = dsetsVectorFields(iVec)
         call create_dataspace_hdf5(trgt_file_id,dsetname,ds_rank,ds_dims2d,dtype)
      end do
      !----------------------------------------------------------------------

   end subroutine create_fields_dsets_in_new_results_file

end module mod_copy_results
