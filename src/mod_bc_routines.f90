module mod_bc_routines

   use mod_constants

   implicit none

      contains

         subroutine temporary_bc_routine(npoin,nboun,bou_codes,bound,nbnodes,lbnodes,aux_rho,aux_q,aux_p,aux_c,aux_u,dt,bcc)

            implicit none

            integer(4), intent(in)    :: npoin, nboun, bou_codes(nboun,2), bound(nboun,npbou)
            integer(4), intent(in)    :: nbnodes, lbnodes(nbnodes),bcc
            real(8),    intent(inout) :: aux_rho(npoin),aux_q(npoin,ndime),aux_p(npoin),aux_c(npoin),aux_u(npoin,ndime),dt
            real(8)                   :: ubc,L5,L1
            integer(4)                :: iboun, bcode, ipbou

            if(bcc>0) then
               if (ndime == 3) then
                  !
                  ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                  ! Nodes belonging to both codes will be zeroed on both directions.
                  ! Like this, there's no need to fnd intersections.
                  !
                  !$acc parallel loop gang
                  do iboun = 1,nboun
                     bcode = bou_codes(iboun,2) ! Boundary element code
                     if (bcode == 3) then ! wall
                        !$acc loop vector
                        do ipbou = 1,npbou
                           aux_rho(bound(iboun,ipbou)) = 1.0d0
                        end do
                     end if
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
            else
               if (ndime == 3) then
                  !
                  ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
                  ! Nodes belonging to both codes will be zeroed on both directions.
                  ! Like this, there's no need to fnd intersections.
                  !
                  !$acc parallel loop gang
                  do iboun = 1,nboun
                     bcode = bou_codes(iboun,2) ! Boundary element code
                     if (bcode == 4) then ! inlet just for aligened inlets with x
                        !$acc loop vector
                        do ipbou = 1,npbou
                           L1=nscbc_sigma*aux_c(bound(iboun,ipbou))*(aux_p(bound(iboun,ipbou))-nscbc_p_inf)/nscbc_L
                           L5=nscbc_sigma*aux_c(bound(iboun,ipbou))*aux_c(bound(iboun,ipbou))*aux_rho(bound(iboun,ipbou))*(aux_u(bound(iboun,ipbou),1)-nscbc_ut)/nscbc_L 

                           ubc = aux_u(bound(iboun,ipbou),1) - dt*(L5-L1)/(2.0d0*aux_c(bound(iboun,ipbou))*aux_rho(bound(iboun,ipbou)))
                           aux_q(bound(iboun,ipbou),1) = ubc*aux_rho(bound(iboun,ipbou))
                           aux_q(bound(iboun,ipbou),2) = 0.0d0
                           aux_q(bound(iboun,ipbou),3) = 0.0d0
                           aux_u(bound(iboun,ipbou),1) = ubc
                           aux_u(bound(iboun,ipbou),2) = 0.0d0
                           aux_u(bound(iboun,ipbou),3) = 0.0d0

                           !ubc = aux_p(bound(iboun,ipbou)) - dt*(L5+L1)/2.0d0
                           !aux_p(bound(iboun,ipbou)) = ubc
                        end do
                     else if (bcode == 5) then ! outlet
                        ! we assume neuman
                        !$acc loop vector
                        do ipbou = 1,npbou
                           L1=nscbc_sigma*aux_c(bound(iboun,ipbou))*(aux_p(bound(iboun,ipbou))-nscbc_p_inf)/nscbc_L
                           L5=nscbc_sigma*aux_c(bound(iboun,ipbou))*aux_c(bound(iboun,ipbou))*aux_rho(bound(iboun,ipbou))*(aux_u(bound(iboun,ipbou),1)-nscbc_ut)/nscbc_L 
                           
                           ubc = aux_p(bound(iboun,ipbou)) - dt*(L1)/2.0d0
                           !ubc = aux_p(bound(iboun,ipbou)) - dt*(L5+L1)/2.0d0
                           aux_p(bound(iboun,ipbou)) = ubc

                           !ubc = aux_u(bound(iboun,ipbou),1) - dt*(L5-L1)/(2.0d0*aux_c(bound(iboun,ipbou))*aux_rho(bound(iboun,ipbou)))
                           !aux_q(bound(iboun,ipbou),1) = ubc*aux_rho(bound(iboun,ipbou))
                           !aux_q(bound(iboun,ipbou),2) = 0.0d0
                           !aux_q(bound(iboun,ipbou),3) = 0.0d0
                           !aux_u(bound(iboun,ipbou),1) = ubc
                           !aux_u(bound(iboun,ipbou),2) = 0.0d0
                           !aux_u(bound(iboun,ipbou),3) = 0.0d0
                        end do
                     end if
                  end do
                  !$acc end parallel loop
               end if
            end if
         end subroutine temporary_bc_routine
end module mod_bc_routines
