module mod_bc_routines

   use mod_constants

   implicit none

      contains

         subroutine temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bound,nbnodes,lbnodes,lnbn,aux_rho,aux_q,aux_u,aux_p,aux_E)

            implicit none

            integer(4), intent(in)    :: npoin, nboun, bou_codes(nboun,2), bound(nboun,npbou)
            integer(4), intent(in)    :: nbnodes, lbnodes(nbnodes),lnbn(nboun,npbou)
            real(8),    intent(inout) :: aux_rho(npoin),aux_q(npoin,ndime),aux_u(npoin,ndime),aux_p(npoin),aux_E(npoin)
            integer(4)                :: iboun, bcode, ipbou
            real(8)                   :: cin,R_plus,R_minus,v_b,c_b,s_b,rho_b,p_b,rl,rr

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
                  if (bcode == 4) then ! inlet just for aligened inlets with x
                     !$acc loop vector
                     do ipbou = 1,npbou
                        cin = sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou)))
                        R_plus  = aux_q(lnbn(iboun,ipbou),1)/aux_rho(lnbn(iboun,ipbou)) + 2.0d0*cin/(nscbc_gamma_inf-1.0d0)
                        R_minus = nscbc_u_inf - 2.0d0*nscbc_c_inf/(nscbc_gamma_inf-1.0d0)
                        v_b = (R_plus+R_minus)*0.5d0
                        c_b  = (nscbc_gamma_inf-1.0d0)*(R_plus-R_minus)/4.0d0
                        s_b = nscbc_c_inf**2/(nscbc_gamma_inf*nscbc_rho_inf**(nscbc_gamma_inf-1.0d0))
                        rho_b = (c_b**2/(nscbc_gamma_inf*s_b))**(1.0/(nscbc_gamma_inf-1.0d0))
                        p_b = rho_b*c_b**2/nscbc_gamma_inf

                        !write(*,*) "v_b ",v_b," rb ",rho_b

                       !aux_q(bound(iboun,ipbou),1) = v_b*rho_b
                       !aux_q(bound(iboun,ipbou),2) = 0.0d0
                       !aux_q(bound(iboun,ipbou),3) = 0.0d0

                       !aux_rho(bound(iboun,ipbou)) = rho_b

                       !aux_E(bound(iboun,ipbou)) = rho_b*0.5d0*v_b**2 + p_b/(nscbc_gamma_inf-1.0d0)

                       !aux_p(bound(iboun,ipbou)) = p_b

                       !aux_u(bound(iboun,ipbou),1) = v_b
                       !aux_u(bound(iboun,ipbou),2) = 0.0d0
                       !aux_u(bound(iboun,ipbou),3) = 0.0d0

                        rl = sqrt(nscbc_rho_inf)
                        rr = sqrt(aux_rho(lnbn(iboun,ipbou)))

                        aux_q(bound(iboun,ipbou),1) = (rl*nscbc_u_inf*nscbc_rho_inf + rr*aux_q(lnbn(iboun,ipbou),1))/(rr+rl) 
                        aux_q(bound(iboun,ipbou),2) = 0.0d0
                        aux_q(bound(iboun,ipbou),3) = 0.0d0

                        aux_rho(bound(iboun,ipbou)) = sqrt(aux_rho(lnbn(iboun,ipbou))*nscbc_rho_inf)

                        aux_E(bound(iboun,ipbou)) = (aux_E(lnbn(iboun,ipbou))*rr + &
                                                    rl*(nscbc_rho_inf*0.5d0*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0d0)))/(rr+rl)

                        aux_u(bound(iboun,ipbou),1) = (rr*aux_u(lnbn(iboun,ipbou),1)+rl*nscbc_u_inf)/(rr+rl)
                        aux_u(bound(iboun,ipbou),2) = 0.0d0
                        aux_u(bound(iboun,ipbou),3) = 0.0d0

                        aux_p(bound(iboun,ipbou)) = (rr*aux_p(lnbn(iboun,ipbou))+rl*nscbc_p_inf)/(rr+rl)
                     end do
                  else if (bcode == 5) then ! inlet just for aligened inlets with x
                     !$acc loop vector
                     do ipbou = 1,npbou
                        cin = sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou)))
                        R_plus  = aux_q(lnbn(iboun,ipbou),1)/aux_rho(lnbn(iboun,ipbou)) + 2.0d0*cin/(nscbc_gamma_inf-1.0d0)
                        R_minus = nscbc_u_inf - 2.0d0*nscbc_c_inf/(nscbc_gamma_inf-1.0d0)
                        v_b = (R_plus+R_minus)*0.5d0
                        c_b  = (nscbc_gamma_inf-1.0d0)*(R_plus-R_minus)/4.0d0
                        s_b = cin**2/(nscbc_gamma_inf*aux_rho(lnbn(iboun,ipbou))**(nscbc_gamma_inf-1.0d0))
                        rho_b = (c_b**2/(nscbc_gamma_inf*s_b))**(1.0/(nscbc_gamma_inf-1.0d0))
                        p_b = rho_b*c_b**2/nscbc_gamma_inf

                       !aux_q(bound(iboun,ipbou),1) = v_b*rho_b
                       !aux_q(bound(iboun,ipbou),2) = 0.0d0
                       !aux_q(bound(iboun,ipbou),3) = 0.0d0

                       !aux_rho(bound(iboun,ipbou)) = rho_b

                       !aux_E(bound(iboun,ipbou)) = rho_b*0.5d0*v_b**2 + p_b/(nscbc_gamma_inf-1.0d0)

                       !aux_p(bound(iboun,ipbou)) = p_b

                       !aux_u(bound(iboun,ipbou),1) = v_b
                       !aux_u(bound(iboun,ipbou),2) = 0.0d0
                       !aux_u(bound(iboun,ipbou),3) = 0.0d0

                        rl = sqrt(nscbc_rho_inf)
                        rr = sqrt(aux_rho(lnbn(iboun,ipbou)))

                        aux_q(bound(iboun,ipbou),1) = (rl*nscbc_u_inf*nscbc_rho_inf + rr*aux_q(lnbn(iboun,ipbou),1))/(rr+rl) 
                        aux_q(bound(iboun,ipbou),2) = 0.0d0
                        aux_q(bound(iboun,ipbou),3) = 0.0d0

                        aux_rho(bound(iboun,ipbou)) = sqrt(aux_rho(lnbn(iboun,ipbou))*nscbc_rho_inf)

                        aux_E(bound(iboun,ipbou)) = (aux_E(lnbn(iboun,ipbou))*rr + &
                                                    rl*(nscbc_rho_inf*0.5d0*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0d0)))/(rr+rl)

                        aux_u(bound(iboun,ipbou),1) = (rr*aux_u(lnbn(iboun,ipbou),1)+rl*nscbc_u_inf)/(rr+rl)
                        aux_u(bound(iboun,ipbou),2) = 0.0d0
                        aux_u(bound(iboun,ipbou),3) = 0.0d0

                        aux_p(bound(iboun,ipbou)) = (rr*aux_p(lnbn(iboun,ipbou))+rl*nscbc_p_inf)/(rr+rl)
                     end do
                  else if (bcode == 1) then
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),2) = 0.0d0
                        aux_u(bound(iboun,ipbou),2) = 0.0d0

                        aux_p(bound(iboun,ipbou)) = nscbc_p_inf
                        aux_rho(bound(iboun,ipbou)) = nscbc_rho_inf
                        aux_E(bound(iboun,ipbou)) = nscbc_rho_inf*0.5d0*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0d0)
                     end do
                  else if (bcode == 2) then
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),3) = 0.0d0
                        aux_u(bound(iboun,ipbou),3) = 0.0d0
                     end do
                  else if (bcode == 3) then ! non_slip wall
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),1) = 0.0d0
                        aux_q(bound(iboun,ipbou),2) = 0.0d0
                        aux_q(bound(iboun,ipbou),3) = 0.0d0

                        aux_u(bound(iboun,ipbou),1) = 0.0d0
                        aux_u(bound(iboun,ipbou),2) = 0.0d0
                        aux_u(bound(iboun,ipbou),3) = 0.0d0
                     end do
                  end if
               end do
               !$acc end parallel loop
            end if
         end subroutine temporary_bc_routine_dirichlet_prim

      end module mod_bc_routines
