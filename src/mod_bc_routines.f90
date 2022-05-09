module mod_bc_routines
   implicit none
      contains
         subroutine temporary_bc_routine()
            ! TODO: This wont compile yet
            implicit none
               if (nboun .ne. 0) then
                  if (ndime == 3) then
                     !
                     ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                     ! Nodes belonging to both codes will be zeroed on both directions.
                     ! Like this, there's no need to fnd intersections.
                     !
                     !$acc parallel loop gang
                     do iboun = 1,nboun
                        bcode = bou_codes(iboun,2) ! Boundary element code
                        if (bcode == 3) then ! inlet
                           !$acc loop vector
                           do ipbou = 1,npbou
                              rho_1(bound(iboun,ipbou)) = 1.0d0
                           end do
                        end if
                     end do
                     !$acc end parallel loop
                  end if
                end if
                !
                ! Janky boundary conditions. TODO: Fix this shite...
                !
                if (nboun .ne. 0) then
                   if (ndime == 2) then
                      !$acc kernels
                      q_1(lbnodes,2) = 0.0d0
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
                               q_1(bound(iboun,ipbou),2) = 0.0d0
                            end do
                         else if (bcode == 2) then
                            !$acc loop vector
                            do ipbou = 1,npbou
                               q_1(bound(iboun,ipbou),3) = 0.0d0
                            end do
                         else if (bcode == 3) then ! non_slip wall
                            !$acc loop vector
                            do ipbou = 1,npbou
                               q_1(bound(iboun,ipbou),1) = 0.0d0
                               q_1(bound(iboun,ipbou),2) = 0.0d0
                               q_1(bound(iboun,ipbou),3) = 0.0d0
                            end do
                         end if
                      end do
                      !$acc end parallel loop
                   end if
                end if
         end subroutine temporary_bc_routine
end module mod_bc_routines