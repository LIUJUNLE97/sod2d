module mod_quality
    use mod_constants
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


    subroutine eval_ElemQuality_simple(mnnode,mngaus,npoin, nelem, ielem, coordPar, connecParOrig, dNgp, wgp, quality_vec,qual_gauss)
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
        
    end subroutine eval_ElemQuality_simple

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
    

end module mod_quality
