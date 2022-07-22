module mod_bc_routines

   use mod_constants

   implicit none

      contains

         subroutine temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bound,nbnodes,lbnodes,lnbn,aux_rho,aux_q,aux_u,aux_p,aux_E)

            implicit none

            integer(4), intent(in)    :: npoin, nboun, bou_codes(nboun,2), bound(nboun,npbou)
            integer(4), intent(in)    :: nbnodes, lbnodes(nbnodes),lnbn(nboun,npbou)
            real(rp),    intent(inout) :: aux_rho(npoin),aux_q(npoin,ndime),aux_u(npoin,ndime),aux_p(npoin),aux_E(npoin)
            integer(4)                :: iboun, bcode, ipbou
            real(rp)                   :: cin,R_plus,R_minus,v_b,c_b,s_b,rho_b,p_b,rl,rr, sl, sr
            real(rp)                   :: q_hll,rho_hll,E_hll, E_inf

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
                        aux_rho(bound(iboun,ipbou)) = 1.0_rp
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
               aux_q(lbnodes,2) = 0.0_rp
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

                       sl = min(nscbc_u_inf-nscbc_c_inf, aux_u(lnbn(iboun,ipbou),1) - sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))))
                       sr =  max(aux_u(lnbn(iboun,ipbou),1) + sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))), nscbc_u_inf+nscbc_c_inf)

                       E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))

                       rho_hll = (sr*aux_rho(lnbn(iboun,ipbou))-sl*nscbc_rho_inf+nscbc_rho_inf*nscbc_u_inf-aux_q(lnbn(iboun,ipbou),1))/(sr-sl)
                       q_hll   = (sr*aux_q(lnbn(iboun,ipbou),1)-sl*nscbc_rho_inf*nscbc_u_inf+nscbc_rho_inf*nscbc_u_inf**2-aux_u(lnbn(iboun,ipbou),1)*aux_q(lnbn(iboun,ipbou),1))/(sr-sl)
                       E_hll   = (sr*aux_E(lnbn(iboun,ipbou))-sl*E_inf+nscbc_u_inf*(E_inf+nscbc_p_inf)-aux_u(lnbn(iboun,ipbou),1)*(aux_E(lnbn(iboun,ipbou))+aux_p(lnbn(iboun,ipbou))))/(sr-sl)

                       sl = min(nscbc_u_inf-nscbc_c_inf, aux_u(lnbn(iboun,ipbou),1) - sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))))
                       sr =  max(aux_u(lnbn(iboun,ipbou),1) + sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))), nscbc_u_inf+nscbc_c_inf)

                       E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))

                       rho_hll = (sr*aux_rho(lnbn(iboun,ipbou))-sl*nscbc_rho_inf+nscbc_rho_inf*nscbc_u_inf-aux_q(lnbn(iboun,ipbou),1))/(sr-sl)
                       q_hll   = (sr*aux_q(lnbn(iboun,ipbou),1)-sl*nscbc_rho_inf*nscbc_u_inf+nscbc_rho_inf*nscbc_u_inf**2-aux_u(lnbn(iboun,ipbou),1)*aux_q(lnbn(iboun,ipbou),1))/(sr-sl)
                       E_hll   = (sr*aux_E(lnbn(iboun,ipbou))-sl*E_inf+nscbc_u_inf*(E_inf+nscbc_p_inf)-aux_u(lnbn(iboun,ipbou),1)*(aux_E(lnbn(iboun,ipbou))+aux_p(lnbn(iboun,ipbou))))/(sr-sl)


                       aux_q(bound(iboun,ipbou),1) = q_hll

                       aux_rho(bound(iboun,ipbou)) = rho_hll
                       aux_E(bound(iboun,ipbou)) = E_hll

                       sl = min(-nscbc_c_inf, aux_u(lnbn(iboun,ipbou),2) - sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))))
                       sr =  max(aux_u(lnbn(iboun,ipbou),2) + sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))), nscbc_c_inf)
                       q_hll   = (sr*aux_q(lnbn(iboun,ipbou),2)-aux_u(lnbn(iboun,ipbou),2)*aux_q(lnbn(iboun,ipbou),2))/(sr-sl)

                       aux_q(bound(iboun,ipbou),2) = q_hll

                       sl = min(-nscbc_c_inf, aux_u(lnbn(iboun,ipbou),3) - sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))))
                       sr =  max(aux_u(lnbn(iboun,ipbou),3) + sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))), nscbc_c_inf)
                       q_hll   = (sr*aux_q(lnbn(iboun,ipbou),3)-aux_u(lnbn(iboun,ipbou),2)*aux_q(lnbn(iboun,ipbou),3))/(sr-sl)

                       aux_q(bound(iboun,ipbou),3) = q_hll

                       aux_u(bound(iboun,ipbou),1) = aux_q(bound(iboun,ipbou),1)/aux_rho(bound(iboun,ipbou))
                       aux_u(bound(iboun,ipbou),2) = aux_q(bound(iboun,ipbou),2)/aux_rho(bound(iboun,ipbou))
                       aux_u(bound(iboun,ipbou),3) = aux_q(bound(iboun,ipbou),3)/aux_rho(bound(iboun,ipbou))

                       aux_p(bound(iboun,ipbou)) = aux_rho(bound(iboun,ipbou))*(nscbc_gamma_inf-1.0_rp)*((aux_E(bound(iboun,ipbou))/aux_rho(bound(iboun,ipbou)))- &
                          0.5_rp*dot_product(aux_u(bound(iboun,ipbou),:),aux_u(bound(iboun,ipbou),:)))

                     end do
                  else if (bcode == 5) then ! inlet just for aligened inlets with x
                     !$acc loop vector
                     do ipbou = 1,npbou

                       sl = min(nscbc_u_inf-nscbc_c_inf, aux_u(lnbn(iboun,ipbou),1) - sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))))
                       sr =  max(aux_u(lnbn(iboun,ipbou),1) + sqrt(nscbc_gamma_inf*aux_p(lnbn(iboun,ipbou))/aux_rho(lnbn(iboun,ipbou))), nscbc_u_inf+nscbc_c_inf)

                       E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))

                       rho_hll = (sr*aux_rho(lnbn(iboun,ipbou))-sl*nscbc_rho_inf+nscbc_rho_inf*nscbc_u_inf-aux_q(lnbn(iboun,ipbou),1))/(sr-sl)
                       q_hll   = (sr*aux_q(lnbn(iboun,ipbou),1)-sl*nscbc_rho_inf*nscbc_u_inf+nscbc_rho_inf*nscbc_u_inf**2-aux_u(lnbn(iboun,ipbou),1)*aux_q(lnbn(iboun,ipbou),1))/(sr-sl)
                       E_hll   = (sr*aux_E(lnbn(iboun,ipbou))-sl*E_inf+nscbc_u_inf*(E_inf+nscbc_p_inf)-aux_u(lnbn(iboun,ipbou),1)*(aux_E(lnbn(iboun,ipbou))+aux_p(lnbn(iboun,ipbou))))/(sr-sl)


                       aux_q(bound(iboun,ipbou),1) = q_hll
                       aux_q(bound(iboun,ipbou),2) = 0.0_rp
                       aux_q(bound(iboun,ipbou),3) = 0.0_rp

                       aux_rho(bound(iboun,ipbou)) = rho_hll
                       aux_E(bound(iboun,ipbou)) = E_hll

                       aux_u(bound(iboun,ipbou),1) = aux_q(bound(iboun,ipbou),1)/aux_rho(bound(iboun,ipbou))
                       aux_u(bound(iboun,ipbou),2) = aux_q(bound(iboun,ipbou),2)/aux_rho(bound(iboun,ipbou))
                       aux_u(bound(iboun,ipbou),3) = aux_q(bound(iboun,ipbou),3)/aux_rho(bound(iboun,ipbou))

                       aux_p(bound(iboun,ipbou)) = aux_rho(bound(iboun,ipbou))*(nscbc_gamma_inf-1.0_rp)*((aux_E(bound(iboun,ipbou))/aux_rho(bound(iboun,ipbou)))- &
                          0.5_rp*dot_product(aux_u(bound(iboun,ipbou),:),aux_u(bound(iboun,ipbou),:)))

                     end do
                  else if (bcode == 1) then
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),1) = nscbc_rho_inf*nscbc_u_inf
                        aux_u(bound(iboun,ipbou),1) = nscbc_u_inf
                        aux_q(bound(iboun,ipbou),2) = 0.0_rp
                        aux_u(bound(iboun,ipbou),2) = 0.0_rp
                        aux_q(bound(iboun,ipbou),3) = 0.0_rp
                        aux_u(bound(iboun,ipbou),3) = 0.0_rp

                        aux_p(bound(iboun,ipbou)) = nscbc_p_inf
                        aux_rho(bound(iboun,ipbou)) = nscbc_rho_inf
                        aux_E(bound(iboun,ipbou)) = nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
                     end do
                  else if (bcode == 2) then
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),3) = 0.0_rp
                        aux_u(bound(iboun,ipbou),3) = 0.0_rp
                     end do
                  else if (bcode == 3) then ! non_slip wall
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),1) = 0.0_rp
                        aux_q(bound(iboun,ipbou),2) = 0.0_rp
                        aux_q(bound(iboun,ipbou),3) = 0.0_rp

                        aux_u(bound(iboun,ipbou),1) = 0.0_rp
                        aux_u(bound(iboun,ipbou),2) = 0.0_rp
                        aux_u(bound(iboun,ipbou),3) = 0.0_rp
                     end do
                  else if (bcode == 6) then ! non_slip wall hot
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),1) = 0.0_rp
                        aux_q(bound(iboun,ipbou),2) = 0.0_rp
                        aux_q(bound(iboun,ipbou),3) = 0.0_rp

                        aux_u(bound(iboun,ipbou),1) = 0.0_rp
                        aux_u(bound(iboun,ipbou),2) = 0.0_rp
                        aux_u(bound(iboun,ipbou),3) = 0.0_rp

                        aux_rho(bound(iboun,ipbou)) = nscbc_p_inf/nscbc_Rgas_inf/586.0_rp
                        aux_p(bound(iboun,ipbou)) = nscbc_p_inf
                        aux_E(bound(iboun,ipbou)) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
                     end do
                  else if (bcode == 7) then ! non_slip wall cold
                     !$acc loop vector
                     do ipbou = 1,npbou
                        aux_q(bound(iboun,ipbou),1) = 0.0_rp
                        aux_q(bound(iboun,ipbou),2) = 0.0_rp
                        aux_q(bound(iboun,ipbou),3) = 0.0_rp

                        aux_u(bound(iboun,ipbou),1) = 0.0_rp
                        aux_u(bound(iboun,ipbou),2) = 0.0_rp
                        aux_u(bound(iboun,ipbou),3) = 0.0_rp

                        aux_rho(bound(iboun,ipbou)) = nscbc_p_inf/nscbc_Rgas_inf/293.0_rp
                        aux_p(bound(iboun,ipbou)) = nscbc_p_inf
                        aux_E(bound(iboun,ipbou)) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
                    end do
                  end if
               end do
               !$acc end parallel loop
            end if
         end subroutine temporary_bc_routine_dirichlet_prim

      end module mod_bc_routines
