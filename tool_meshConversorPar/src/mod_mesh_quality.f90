module mod_mesh_quality
   use mod_constants
   use mod_mpi
   use elem_hex
   use quadrature_rules
   use jacobian_oper
   use mod_custom_types
   implicit none

contains
    subroutine ideal_hexa(mnnode, nelem, npoin, ielem, coord, connec, idealJ)
        implicit none
        integer(4),parameter  :: nedge = 12
        integer(4),intent(in) :: mnnode, ielem, npoin, nelem
        integer(4),intent(in) :: connec(nelem,mnnode)
        real(8), intent(in)   :: coord(npoin,ndime)
        real(8), intent(out)  :: idealJ(ndime,ndime)
        real(8)               :: dist(nedge,ndime),dist2(nedge),h1,h2,h3

        call get_hexa_edges_dist_r8(mnnode,ielem,nelem,npoin,connec,coord,dist)
        dist2(:) = sqrt(dist(:,1)*dist(:,1)+dist(:,2)*dist(:,2)+dist(:,3)*dist(:,3))
        idealJ(:,:) = 0.0d0
        idealJ(1,1) = 1.0d0/((dist2(1)+dist2(3)+dist2(5)+dist2(7))/4.0d0)
        idealJ(2,2) = 1.0d0/((dist2(2)+dist2(4)+dist2(6)+dist2(8))/4.0d0)
        idealJ(3,3) = 1.0d0/((dist2(9)+dist2(10)+dist2(11)+dist2(12))/4.0d0)
    end subroutine

    subroutine shape_measure(elemJ, idealJ, eta, isInvalid)
        implicit none
        real(8),intent(in)  :: elemJ(ndime,ndime), idealJ(ndime,ndime)
        real(8),intent(out) :: eta
        logical, intent(out) :: isInvalid
        real(8) :: S(ndime,ndime), S2(ndime,ndime), sigma, Sf, detS
        real(8),parameter :: d=3.0d0
        real(8)                :: Saux(ndime,ndime),detaux
        S     = matmul(elemJ, idealJ)

        detS  = det_3x3(S)

        sigma = (detS + abs(detS))/2
        S2    = matmul(transpose(S), S)
        Sf    = S2(1,1) + S2(2,2) + S2(3,3)
        eta   = Sf/(d*sigma**(2.0d0/d)) 
        
        isInvalid = detS < tiny(1.0d0)
    end subroutine

   subroutine eval_meshQuality_and_writeInFile(mesh_fileName,numMshRanks2Part,typeCoords,mnnode,mngaus,numMshRanksInMpiRank,numVTKElemsPerMshElem,&
                                             numElemsVTKMshRank,numElemsMshRank,numNodesMshRank,mshRanksInMpiRank,&
                                             dNgp,wgp,connecParOrig_jm,elemGidMshRank_jv,coordPar_jm,quality_jm)
      implicit none
      character(len=*),intent(in) :: mesh_fileName
      integer(4),intent(in) :: numMshRanks2Part,typeCoords,mnnode,mngaus,numMshRanksInMpiRank,numVTKElemsPerMshElem
      integer(4),intent(in),dimension(numMshRanksInMpiRank) :: numElemsVTKMshRank,numElemsMshRank,numNodesMshRank,mshRanksInMpiRank
      real(8),intent(in) :: dNgp(mngaus,mnnode),wgp(mngaus) 
      type(jagged_matrix_int4),intent(in) :: connecParOrig_jm
      type(jagged_vector_int4),intent(in) :: elemGidMshRank_jv
      type(jagged_matrix_real8),intent(in) :: coordPar_jm
      type(jagged_matrix_real8),intent(inout) :: quality_jm
      !--------------------------------------------------------------------------------
      real(8) :: maxQuality(2),minQuality(2),avgQuality(2)
      real(8) :: rankMaxQuality(2),rankMinQuality(2),rankAvgQuality(2)
      integer(4) :: iMshRank,numTangledElemsMpiRank,numTangledElemsTotal
      !--------------------------------------------------------------------------------
      character(len=12) :: aux_numRanks,aux_mpiRank,aux_typeCoord
      character(len=512) :: meshQualityBase,meshQuality_filename,meshQualitySum_filename
      !--------------------------------------------------------------------------------

      meshQualityBase = 'meshQuality'
      if(typeCoords==1) then
         aux_typeCoord = '_gll_'
      elseif(typeCoords==2) then
         aux_typeCoord = '_equi_'
      else
         aux_typeCoord = '_error_'
      endif

      if(numMshRanksInMpiRank .gt. 0) then
         ! Open a file for outputting realted information
         write(aux_mpiRank,'(I0)') mpi_rank
         write(aux_numRanks,'(I0)') numMshRanks2Part
         !-----------------------------------------------------------------------------------------------
         meshQuality_filename = trim(adjustl(meshQualityBase))//trim(adjustl(aux_typeCoord))//trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'_rank'//trim(aux_mpiRank)//'.dat'
         open(unit=555, file=meshQuality_filename, status="unknown", action="write", form="formatted")
         write(555,*) "--| Evaluating mesh quality"
         write(555,*) "----| List of tangled elements (GMSH global numeration)"
         call flush(555)
      end if
      !-----------------------------------------------------------------------------------------------
      if(mpi_rank.eq.0) then
         meshQualitySum_filename = trim(adjustl(meshQualityBase))//trim(adjustl(aux_typeCoord))//"summary_"//trim(adjustl(mesh_fileName))//'-'//trim(aux_numRanks)//'.dat'
         open(unit=554, file=meshQualitySum_filename, status="unknown", action="write", form="formatted")
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

      call eval_MeshQuality(numMshRanksInMpiRank,numElemsVTKMshRank,numElemsMshRank,numNodesMshRank,mnnode,mngaus,coordPar_jm,&
                           connecParOrig_jm,dNgp,wgp,mshRanksInMpiRank,elemGidMshRank_jv,numVTKElemsPerMshElem,quality_jm,&
                           rankMaxQuality, rankMinQuality, rankAvgQuality)

      ! Compute max, min and avg across all ranks
      call MPI_Reduce(numTangledElemsMpiRank,numTangledElemsTotal,1,mpi_datatype_int4,MPI_SUM,0,app_comm,mpi_err)
      call MPI_Reduce(rankMaxQuality,maxQuality,2,mpi_datatype_real8,MPI_MAX,0,app_comm,mpi_err)
      call MPI_Reduce(rankMinQuality,minQuality,2,mpi_datatype_real8,MPI_MIN,0,app_comm,mpi_err)
      call MPI_Reduce(rankAvgQuality,avgQuality,2,mpi_datatype_real8,MPI_SUM,0,app_comm,mpi_err)
      avgQuality(:) = avgQuality(:) / min(numMshRanks2Part,mpi_size)
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

   end subroutine eval_meshQuality_and_writeInFile

   subroutine eval_MeshQuality(numMshRanksInMpiRank,numElemsVTKMshRank,numElemsMshRank,numNodesMshRank,mnnode,mngaus,coordPar_jm,&
                               connecParOrig_jm,dNgp,wgp,mshRanksInMpiRank,elemGidMshRank_jv,numVTKElemsPerMshElem,quality_jm,&
                               rankMaxQuality,rankMinQuality,rankAvgQuality)
      implicit none
      integer(4), intent(in) :: numMshRanksInMpiRank            
      integer(4), intent(in) :: numElemsVTKMshRank(:)           
      integer(4), intent(in) :: numElemsMshRank(:)              
      integer(4), intent(in) :: numNodesMshRank(:)              
      integer(4), intent(in) :: mnnode, mngaus                  
      type(jagged_matrix_real8), intent(in) :: coordPar_jm
      type(jagged_matrix_int4), intent(in) :: connecParOrig_jm 
      real(8), intent(in) :: dNgp(mngaus, mnnode)               
      real(8), intent(in) :: wgp(mngaus)                        
      integer(4), intent(in) :: mshRanksInMpiRank(:)            
      type(jagged_vector_int4), intent(in) :: elemGidMshRank_jv 
      integer(4), intent(in) :: numVTKElemsPerMshElem           
      type(jagged_matrix_real8), intent(inout) :: quality_jm
      real(8) :: quality_elem, auxAvg(2), elemVertOrig(3)
      real(8) :: quality_vec(2)
      integer(4) :: iElem, iNode, iElemVTK, iElemGid, iMshRank, mshRank
      integer(4) :: numTangledElemsMpiRank, numTangledElemsMshRank
      real(8), intent(out) :: rankMaxQuality(2), rankMinQuality(2), rankAvgQuality(2)
   
      ! Initialize vars to get high-level info about mesh quality
      rankMaxQuality(:) = 0.0d0          ! Max absolute
      rankMinQuality(:) = 100000.0d0     ! Min absolute
      rankAvgQuality(:) = 0.0d0          ! Avg across all ranks
      numTangledElemsMpiRank=0

      do iMshRank=1,numMshRanksInMpiRank
          allocate(quality_jm%matrix(iMshRank)%elems(numElemsVTKMshRank(iMshRank),2))
          auxAvg(:) = 0.0d0
          numTangledElemsMshRank=0
          do iElem = 1, numElemsMshRank(iMshRank)
              call eval_ElemQuality(mnnode,mngaus,numNodesMshRank(iMshRank),numElemsMshRank(iMshRank),iElem,coordPar_jm%matrix(iMshRank)%elems,connecParOrig_jm%matrix(iMshRank)%elems, dNgp, wgp, quality_vec)
              ! Compute rank max, min and avg
              quality_elem = quality_vec(1)
              rankMaxQuality(1) = max(rankMaxQuality(1), quality_vec(1))
              rankMinQuality(1) = min(rankMinQuality(1), quality_vec(1))
              auxAvg(1) = auxAvg(1) + quality_vec(1)
              rankMaxQuality(2) = max(rankMaxQuality(2), quality_vec(2))
              rankMinQuality(2) = min(rankMinQuality(2), quality_vec(2))
              auxAvg(2) = auxAvg(2) + quality_vec(2)
              if (any(quality_vec < 0)) then
                  numTangledElemsMshRank=numTangledElemsMshRank+1
                  ! Write data to file
                  mshRank = mshRanksInMpiRank(iMshRank)
                  iElemGid = elemGidMshRank_jv%vector(iMshRank)%elems(iElem)

                  iNode = connecParOrig_jm%matrix(iMshRank)%elems(iElem,1)
                  elemVertOrig(:)= coordPar_jm%matrix(iMshRank)%elems(iNode,:)

                  write(555,*) "# Tangled element",numTangledElemsMshRank," found in mshRank",mshRank
                  write(555,*) "  elem(lid)",iElem
                  write(555,*) "  elem(gid)",iElemGid
                  write(555,*) "  anisotropic quality",quality_vec(1)
                  write(555,*) "  isotropic quality",quality_vec(2)
                  write(555,*) "  vertex0",elemVertOrig(:)
              end if
              do iElemVTK = 1, numVTKElemsPerMshElem
                  quality_jm%matrix(iMshRank)%elems((iElem-1)*numVTKElemsPerMshElem + iElemVTK,:) = quality_vec
              end do
          end do
          auxAvg(:) = auxAvg(:) / numElemsMshRank(iMshRank)
          rankAvgQuality(:) = rankAvgQuality(:) + auxAvg(:)
          numTangledElemsMpiRank=numTangledElemsMpiRank+numTangledElemsMshRank
      end do
      if(numMshRanksInMpiRank.gt.0) then
         rankAvgQuality(:) = rankAvgQuality(:) / numMshRanksInMpiRank
      end if

      !Final writing in the output files of each rank
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

   end subroutine eval_MeshQuality

    subroutine eval_ElemQuality(mnnode,mngaus,npoin, nelem, ielem, coordPar, connecParOrig, dNgp, wgp, quality_vec,qual_gauss)
        implicit none
        integer(4), intent(in) :: mnnode,mngaus,npoin,nelem,ielem,connecParOrig(nelem,mnnode)
        real(8), intent(in) :: coordPar(npoin,ndime),dNgp(ndime,mnnode,mngaus),wgp(mngaus)
        real(8), intent(out) :: quality_vec(2)
        real(8), optional, intent(inout) :: qual_gauss(mngaus)
        integer(4) :: igaus, idime
        real(8) :: elemJ(ndime,ndime),idealJ(ndime,ndime),gpvol,gpvolIdeal
        real(8) :: eta,volume,volume_cube,modulus
        real(8) :: eta_elem, eta_cube, quality
        real(8) :: idealCubeJ(ndime,ndime)
        real(8) :: detIdeal,detIdealCube
        logical :: isInvalid

        idealCubeJ = 0.0d0
        do idime = 1,ndime
            idealCubeJ(idime,idime) = 1.0d0
        end do
        detIdealCube = 1.0d0
        eta_elem = 0.0d0
        eta_cube = 0.0d0
        volume = 0.0d0
        volume_cube = 0.0d0
        call ideal_hexa(mnnode,nelem,npoin,ielem,coordPar,connecParOrig,idealJ) !Assumim que el jacobià de l'element ideal és constant
        detIdeal = det_3x3(idealJ)
        detIdeal = 1.0_rp/detIdeal
        do igaus = 1,mngaus
            call compute_jacobian(mnnode,mngaus,nelem,npoin,ielem,igaus,dNgp,wgp(igaus),coordPar,connecParOrig,elemJ,gpvol)
            elemJ = transpose(elemJ)
            gpvolIdeal = detIdeal*wgp(igaus)
            call shape_measure(elemJ, idealJ, eta, isInvalid)
            eta_elem = eta_elem + eta*eta*gpvolIdeal
            volume = volume + 1*1*gpvolIdeal
            
            gpvolIdeal = detIdealCube*wgp(igaus)
            call shape_measure(elemJ, idealCubeJ, eta, isInvalid)
            eta_cube = eta_cube + eta*eta*gpvolIdeal
            volume_cube = volume_cube + 1*1*gpvolIdeal

!             if(present(qual_gauss)) then
!               if(isInvalid) then
!                 qual_gauss(igaus) = huge(1.0d0)
!               else
!                 qual_gauss(igaus) = eta!*gpvolIdeal
!               end if
!             end if

            if(eta>100000) then
                print*,igaus, ' ',eta
            end if
        end do

        eta_elem = sqrt(eta_elem)/sqrt(volume)
        quality = 1.0d0/eta_elem
        modulus = modulo(quality, 1.0d0)
        if (int(modulus) .ne. 0) then
            quality = -1.0d0
        end if
        quality_vec(1) = quality
        
        if (present(qual_gauss)) then
          qual_gauss(:) = quality  ! all gauss points the elemental value -> here with respect to cube (isotropic)
          if (quality<0.01) then
              qual_gauss(:) = -10000.0d0
          end if
        end if
        
        eta_elem = eta_cube
        eta_elem = sqrt(eta_elem)/sqrt(volume_cube)
        quality = 1.0_rp/eta_elem
        modulus = modulo(quality, 1.0d0)
        if (int(modulus) .ne. 0) then
            quality = -1.0d0
        end if
        quality_vec(2) = quality
        
    end subroutine eval_ElemQuality

    function det_3x3(A) result(det)
        implicit none
        real(8), intent(in) :: A(:,:)         ! Input 3x3 matrix
        real(8) :: det           ! Determinant of the matrix
    
        det = A(1, 1)*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2)) &
            - A(1, 2)*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1)) &
            + A(1, 3)*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
    end function det_3x3

     subroutine inverse_matrix_3x3(A, A_inv, det, success)
        implicit none
        ! Input
        real(8), intent(in) :: A(:,:)         ! Input 3x3 matrix
        ! Output
        real(rp), intent(out) :: A_inv(3,3)    ! Inverse of the matrix
        real(rp), intent(out) :: det           ! Determinant of the matrix
        logical, intent(out) :: success        ! True if inversion is successful
    
        ! Local variables
        real(8) :: cof(3, 3)
        integer :: i, j
    
        ! Check matrix dimensions
        if (size(A, 1) /= 3 .or. size(A, 2) /= 3) then
            print *, "Error: Input matrix must be 3x3."
            success = .false.
            return
        end if
    
        ! Calculate the determinant
        det = A(1, 1)*(A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2)) &
            - A(1, 2)*(A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1)) &
            + A(1, 3)*(A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
    
        ! Check if the matrix is invertible
        if (abs(det) < 1.0E-10_rp) then
            success = .false.
            return
        end if
    
        ! Calculate the adjugate (transpose of the cofactor matrix)
        cof(1, 1) =   (A(2, 2)*A(3, 3) - A(2, 3)*A(3, 2))
        cof(1, 2) = - (A(1, 2)*A(3, 3) - A(1, 3)*A(3, 2))
        cof(1, 3) =   (A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))
    
        cof(2, 1) = - (A(2, 1)*A(3, 3) - A(2, 3)*A(3, 1))
        cof(2, 2) =   (A(1, 1)*A(3, 3) - A(1, 3)*A(3, 1))
        cof(2, 3) = - (A(1, 1)*A(2, 3) - A(1, 3)*A(2, 1))
    
        cof(3, 1) =   (A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))
        cof(3, 2) = - (A(1, 1)*A(3, 2) - A(1, 2)*A(3, 1))
        cof(3, 3) =   (A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
    
        ! Compute the inverse
        do i = 1, 3
            do j = 1, 3
                A_inv(i, j) = cof(i, j) / det
            end do
        end do
    
        success = .true.
    end subroutine inverse_matrix_3x3
    

end module
