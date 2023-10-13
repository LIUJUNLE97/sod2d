module mod_mesh_quality
    use mod_constants
    implict none
    contains
        subroutine ideal_hexa()
            implicit none
        end subroutine

        subroutine shape_measure(elemJ, idealJ, eta)
            implicit none
			real(rp),   intent(in)  :: elemJ(ndime,ndime), idealJ(ndime)
			real(rp),   intent(out) :: eta
			real(rp)                :: S(ndime,ndime), S2(ndime,ndime), sigma, Sf, detS

            S     = matmul(elemJ, idealJ)
            detS  = (S(1,1)*S(2,2)*S(3,3)+S(1,2)*S(2,3)*S(3,1)+S(1,3)*S(2,1)*S(3,2)-S(3,1)*S(2,2)*S(1,3)-S(3,2)*S(2,3)*S(1,1)-S(3,3)*S(2,1)*S(1,2))
            sigma = (detS + abs(detS))/2
            S2    = matmul(transpose(S), S)
            Sf    = S2(1,1) + S2(2,2) + S2(3,3)
            eta   = Sf*Sf/(2*sigma) 

        end subroutine
end module