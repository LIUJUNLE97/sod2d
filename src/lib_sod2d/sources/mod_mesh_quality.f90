module mod_mesh_quality
    use mod_constants
    implicit none
contains
    subroutine ideal_hexa(nelem, npoin, ielem, coord, connec, idealJ)
        implicit none
        integer(4), intent(in)  :: ielem, npoin, nelem
		integer(4), intent(in)  :: connec(nelem,nnode)
		real(rp),   intent(in)  :: coord(npoin,ndime)
		real(rp),   intent(out) :: idealJ(ndime,ndime)
        idealJ(:,:) = 0.0_rp
        idealJ(1,1) = 1.0_rp
        idealJ(2,2) = 1.0_rp
        idealJ(3,3) = 1.0_rp
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
