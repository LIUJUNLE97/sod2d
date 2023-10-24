module mod_mesh_quality
    use mod_constants
    use elem_hex
    implicit none
contains
    subroutine ideal_hexa(mnnode, nelem, npoin, ielem, coord, connec, idealJ)
        implicit none
        integer(4), intent(in)  :: mnnode, ielem, npoin, nelem
        integer(4), intent(in)  :: connec(nelem,nnode)
        real(rp),   intent(in)  :: coord(npoin,ndime)
        real(rp),   intent(out) :: idealJ(ndime,ndime)
        integer(4)              :: ncorner, nedge
        real(rp)                :: dist(12,ndime), dist2(12), h1,h2,h3

        call hexa_edges(mnnode,ielem,nelem,npoin,connec,coord,ncorner,nedge,dist)
        dist2(:) = sqrt(dist(:,1)*dist(:,1)+dist(:,2)*dist(:,2)+dist(:,3)*dist(:,3))
        idealJ(:,:) = 0.0_rp
        idealJ(1,1) = 1.0_rp/((dist2(1)+dist2(3)+dist2(5)+dist2(7))/4.0_rp)
        idealJ(2,2) = 1.0_rp/((dist2(2)+dist2(4)+dist2(6)+dist2(8))/4.0_rp)
        idealJ(3,3) = 1.0_rp/((dist2(9)+dist2(10)+dist2(11)+dist2(12))/4.0_rp)
    end subroutine

    subroutine shape_measure(elemJ, idealJ, eta)
        implicit none
		real(rp),   intent(in)  :: elemJ(ndime,ndime), idealJ(ndime,ndime)
		real(rp),   intent(out) :: eta
		real(rp)                :: S(ndime,ndime), S2(ndime,ndime), sigma, Sf, detS
        real(rp)             :: d=3.0_rp
        S     = matmul(elemJ, idealJ)
        detS  = (S(1,1)*S(2,2)*S(3,3)+S(1,2)*S(2,3)*S(3,1)+S(1,3)*S(2,1)*S(3,2)-S(3,1)*S(2,2)*S(1,3)-S(3,2)*S(2,3)*S(1,1)-S(3,3)*S(2,1)*S(1,2))
        sigma = (detS + abs(detS))/2
        S2    = matmul(transpose(S), S)
        Sf    = S2(1,1) + S2(2,2) + S2(3,3)
        eta   = Sf/(d*sigma**(2.0_rp/d)) 
    end subroutine
end module
