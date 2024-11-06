module mod_mesh_quality
    use mod_constants
    use elem_hex
    use quadrature_rules
    use jacobian_oper
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

    subroutine shape_measure(elemJ, idealJ, eta)
        implicit none
		real(8),intent(in)  :: elemJ(ndime,ndime), idealJ(ndime,ndime)
		real(8),intent(out) :: eta
		real(8) :: S(ndime,ndime), S2(ndime,ndime), sigma, Sf, detS
        real(8),parameter :: d=3.0d0
        S     = matmul(elemJ, idealJ)
        detS  = (S(1,1)*S(2,2)*S(3,3)+S(1,2)*S(2,3)*S(3,1)+S(1,3)*S(2,1)*S(3,2)-S(3,1)*S(2,2)*S(1,3)-S(3,2)*S(2,3)*S(1,1)-S(3,3)*S(2,1)*S(1,2))
        sigma = (detS + abs(detS))/2
        S2    = matmul(transpose(S), S)
        Sf    = S2(1,1) + S2(2,2) + S2(3,3)
        eta   = Sf/(d*sigma**(2.0d0/d)) 
    end subroutine

    subroutine eval_MeshQuality(mnnode,mngaus,npoin, nelem, ielem, coordPar, connecParOrig, dNgp, wgp, quality)
        integer(4), intent(in) :: mnnode,mngaus,npoin,nelem,ielem,connecParOrig(nelem,mnnode)
        real(8), intent(in) :: coordPar(npoin,ndime),dNgp(ndime,mnnode,mngaus),wgp(mngaus)
        real(8), intent(out) :: quality
        integer(4) :: igaus
        real(8) :: elemJ(ndime,ndime),idealJ(ndime,ndime),gpvol
        real(8) :: eta,volume,modulus
        real(8) :: eta_elem

        eta_elem = 0.0d0
        volume = 0.0d0
        call ideal_hexa(mnnode,nelem,npoin,ielem,coordPar,connecParOrig,idealJ) !Assumim que el jacobià de l'element ideal és constant
        do igaus = 1,mngaus
            call compute_jacobian(mnnode,mngaus,nelem,npoin,ielem,igaus,dNgp,wgp(igaus),coordPar,connecParOrig,elemJ,gpvol)
            elemJ = transpose(elemJ)
            call shape_measure(elemJ, idealJ, eta)
            eta_elem = eta_elem + eta*eta*gpvol
            volume = volume + 1*1*gpvol
        end do
        eta_elem = sqrt(eta_elem)/sqrt(volume)
        quality = 1.0d0/eta_elem
        modulus = modulo(quality, 1.0d0)
        if (int(modulus) .ne. 0) then
            quality = -1.0d0
        end if
  
     end subroutine eval_MeshQuality
end module
