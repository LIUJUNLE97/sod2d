module mod_bc_routines

   use mod_constants

   implicit none

      contains

         subroutine temporary_bc_routine(npoin,nboun,bou_codes,bound,nbnodes,lbnodes,aux_rho,aux_q)

            implicit none

            integer(4), intent(in)    :: npoin, nboun, bou_codes(nboun,2), bound(nboun,npbou)
            integer(4), intent(in)    :: nbnodes, lbnodes(nbnodes)
            real(8),    intent(inout) :: aux_rho(npoin), aux_q(npoin,ndime)
            integer(4)                :: iboun, bcode, ipbou

            if (ndime == 3) then
               !
               ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
               ! Nodes belonging to both codes will be zeroed on both directions.
               ! Like this, there's no need to fnd intersections.
               !
               !$acc parallel loop gang
               do iboun = 1,nboun
                  bcode = bou_codes(iboun,2) ! Boundary element code
                 ! if (bcode == 3) then ! inlet
                 !    !$acc loop vector
                 !    do ipbou = 1,npbou
                 !       aux_rho(bound(iboun,ipbou)) = 1.0d0
                 !    end do
                 ! end if
               end do
               !$acc end parallel loop
            end if
            !
            ! Janky boundary conditions. TODO: Fix this shite...
            !
            if (ndime == 2) then
               !$acc kernels
               aux_q(lbnodes,2) = 0.0d0
               !$acc end kernels
            else if (ndime == 3) then
               !
               ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
               ! Nodes belonging to both codes will be zeroed on both directions.
               ! Like this, there's no need to fnd intersections.
               !
               !$acc parallel loop gang
               do iboun = 1,nboun
                  bcode = bou_codes(iboun,2) ! Boundary element code
                  if (bcode == 1) then
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),2) = 0.0d0
                     end do
                  else if (bcode == 2) then
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),3) = 0.0d0
                     end do
                  else if (bcode == 3) then ! non_slip wall
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),1) = 0.0d0
                        aux_q(bound(iboun,ipbou),2) = 0.0d0
                        aux_q(bound(iboun,ipbou),3) = 0.0d0
                     end do
                  end if
               end do
               !$acc end parallel loop
            end if
         end subroutine temporary_bc_routine
end module mod_bc_routines
