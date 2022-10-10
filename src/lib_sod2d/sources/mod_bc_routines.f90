module mod_bc_routines

   use mod_constants
   use mod_comms
   use mod_nvtx
   use mod_mpi

   implicit none

      contains

         subroutine temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,aux_rho,aux_q,aux_u,aux_p,aux_E)

            implicit none

            integer(4), intent(in)     :: npoin, nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
            integer(4), intent(in)     :: nbnodes, lbnodes(nbnodes),lnbn(nboun,npbou),lnbn_nodes(npoin)
            real(rp),    intent(inout) :: aux_rho(npoin),aux_q(npoin,ndime),aux_u(npoin,ndime),aux_p(npoin),aux_E(npoin)
            real(rp)                   :: aux_rho2(npoin),aux_q2(npoin,ndime),aux_u2(npoin,ndime),aux_p2(npoin),aux_E2(npoin)
            integer(4)                 :: iboun, bcode, ipbou,inode
            real(rp)                   :: cin,R_plus,R_minus,v_b,c_b,s_b,rho_b,p_b,rl,rr, sl, sr
            real(rp)                   :: q_hll,rho_hll,E_hll, E_inf

#if 0
            if (ndime == 3) then
               !
               ! Janky wall BC for 2 codes (1=y, 2=z) in 3D
               ! Nodes belonging to both codes will be zeroed on both directions.
               ! Like this, there's no need to fnd intersections.
               !
               !$acc parallel loop gang
               do iboun = 1,nboun
                  bcode = bou_codes(iboun) ! Boundary element code
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
                  bcode = bou_codes(iboun) ! Boundary element code
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
                  else if (bcode == 5) then ! outlet just for aligened with x
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
                  else if (bcode == 8) then ! outlet just for aligened with y
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

                        aux_rho(bound(iboun,ipbou)) = nscbc_p_inf/nscbc_Rgas_inf/nscbc_T_H
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

                        aux_rho(bound(iboun,ipbou)) = nscbc_p_inf/nscbc_Rgas_inf/nscbc_T_C
                        aux_p(bound(iboun,ipbou)) = nscbc_p_inf
                        aux_E(bound(iboun,ipbou)) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
                    end do
                  end if
               end do
               !$acc end parallel loop
            end if
#endif
            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  aux_q2(inode,1) = aux_q(lnbn_nodes(inode),1)
                  aux_q2(inode,2) = aux_q(lnbn_nodes(inode),2)
                  aux_q2(inode,3) = aux_q(lnbn_nodes(inode),3)
                  aux_rho2(inode) = aux_rho(lnbn_nodes(inode))
                  aux_E2(inode) = aux_E(lnbn_nodes(inode))
               end if
            end do
            !$acc end parallel loop

            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_max_update_float_sendRcv(aux_q2(:,1))
               call mpi_halo_max_update_float_sendRcv(aux_q2(:,2))
               call mpi_halo_max_update_float_sendRcv(aux_q2(:,3))
               call mpi_halo_max_update_float_sendRcv(aux_rho2(:))
               call mpi_halo_max_update_float_sendRcv(aux_E2(:))
               call nvtxEndRange
            end if

            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  aux_u2(inode,1) = aux_q2(inode,1)/aux_rho2(inode)
                  aux_u2(inode,2) = aux_q2(inode,2)/aux_rho2(inode)
                  aux_u2(inode,3) = aux_q2(inode,3)/aux_rho2(inode)

                  aux_p2(inode) = aux_rho2(inode)*(nscbc_gamma_inf-1.0_rp)*((aux_E2(inode)/aux_rho2(inode))- &
                     0.5_rp*dot_product(aux_u2(inode,:),aux_u2(inode,:)))
               end if
            end do
            !$acc end parallel loop

            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  bcode = bou_codes_nodes(inode) ! Boundary element code
                  if (bcode == bc_type_inlet) then ! inlet just for aligened inlets with x

                   sl = min(nscbc_u_inf-nscbc_c_inf, aux_u2(lnbn_nodes(inode),1) - sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))))
                   sr =  max(aux_u2(lnbn_nodes(inode),1) + sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))), nscbc_u_inf+nscbc_c_inf)

                   E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))

                   rho_hll = (sr*aux_rho2(lnbn_nodes(inode))-sl*nscbc_rho_inf+nscbc_rho_inf*nscbc_u_inf-aux_q2(lnbn_nodes(inode),1))/(sr-sl)
                   q_hll   = (sr*aux_q2(lnbn_nodes(inode),1)-sl*nscbc_rho_inf*nscbc_u_inf+nscbc_rho_inf*nscbc_u_inf**2-aux_u2(lnbn_nodes(inode),1)*aux_q2(lnbn_nodes(inode),1))/(sr-sl)
                   E_hll   = (sr*aux_E2(lnbn_nodes(inode))-sl*E_inf+nscbc_u_inf*(E_inf+nscbc_p_inf)-aux_u2(lnbn_nodes(inode),1)*(aux_E2(lnbn_nodes(inode))+aux_p2(lnbn_nodes(inode))))/(sr-sl)

                   aux_q(inode,1) = q_hll

                   aux_rho(inode) = rho_hll
                   aux_E(inode) = E_hll

                   sl = min(-nscbc_c_inf, aux_u2(lnbn_nodes(inode),2) - sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))))
                   sr =  max(aux_u2(lnbn_nodes(inode),2) + sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))), nscbc_c_inf)
                   q_hll   = (sr*aux_q2(lnbn_nodes(inode),2)-aux_u2(lnbn_nodes(inode),2)*aux_q2(lnbn_nodes(inode),2))/(sr-sl)

                   aux_q(inode,2) = q_hll

                   sl = min(-nscbc_c_inf, aux_u2(lnbn_nodes(inode),3) - sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))))
                   sr =  max(aux_u2(lnbn_nodes(inode),3) + sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))), nscbc_c_inf)
                   q_hll   = (sr*aux_q2(lnbn_nodes(inode),3)-aux_u2(lnbn_nodes(inode),2)*aux_q2(lnbn_nodes(inode),3))/(sr-sl)

                     aux_q(inode,3) = q_hll

                    aux_u(inode,1) = aux_q2(inode,1)/aux_rho2(inode)
                    aux_u(inode,2) = aux_q2(inode,2)/aux_rho2(inode)
                    aux_u(inode,3) = aux_q2(inode,3)/aux_rho2(inode)

                    aux_p(inode) = aux_rho2(inode)*(nscbc_gamma_inf-1.0_rp)*((aux_E2(inode)/aux_rho2(inode))- &
                       0.5_rp*dot_product(aux_u2(inode,:),aux_u2(inode,:)))

                  else if (bcode == bc_type_non_slip_adiabatic) then ! non_slip wall adiabatic
                     
                     aux_q(inode,1) = 0.0_rp
                     aux_q(inode,2) = 0.0_rp
                     aux_q(inode,3) = 0.0_rp

                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp

                     aux_rho(inode) = nscbc_rho_inf
                     aux_E(inode) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
                  else if (bcode == bc_type_non_slip_hot) then ! non_slip wall hot

                     aux_q(inode,1) = 0.0_rp
                     aux_q(inode,2) = 0.0_rp
                     aux_q(inode,3) = 0.0_rp

                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp

                     aux_rho(inode) = nscbc_p_inf/nscbc_Rgas_inf/nscbc_T_H
                     !aux_p(inode) = nscbc_p_inf
                     aux_E(inode) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)

                  else if (bcode == bc_type_non_slip_cold) then ! non_slip wall cold

                     aux_q(inode,1) = 0.0_rp
                     aux_q(inode,2) = 0.0_rp
                     aux_q(inode,3) = 0.0_rp

                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp

                     aux_rho(inode) = nscbc_p_inf/nscbc_Rgas_inf/nscbc_T_C
                     !aux_p(inode) = nscbc_p_inf
                     aux_E(inode) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)

                  else if (bcode == bc_type_slip_adiabatic) then !slip wall in x

                     aux_q(inode,1) = nscbc_rho_inf*nscbc_u_inf
                     aux_u(inode,1) = nscbc_u_inf
                     aux_q(inode,2) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_q(inode,3) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp

                     aux_p(inode) = nscbc_p_inf
                     aux_rho(inode) = nscbc_rho_inf
                     aux_E(inode) = nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)

                  else if (bcode == bc_type_outlet) then ! outlet just for aligened with x

                     sl = min(nscbc_u_inf-nscbc_c_inf, aux_u2(lnbn_nodes(inode),1) - sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))))
                     sr =  max(aux_u2(lnbn_nodes(inode),1) + sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))), nscbc_u_inf+nscbc_c_inf)

                     E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))

                     rho_hll = (sr*aux_rho2(lnbn_nodes(inode))-sl*nscbc_rho_inf+nscbc_rho_inf*nscbc_u_inf-aux_q2(lnbn_nodes(inode),1))/(sr-sl)
                     q_hll   = (sr*aux_q2(lnbn_nodes(inode),1)-sl*nscbc_rho_inf*nscbc_u_inf+nscbc_rho_inf*nscbc_u_inf**2-aux_u2(lnbn_nodes(inode),1)*aux_q2(lnbn_nodes(inode),1))/(sr-sl)
                     E_hll   = (sr*aux_E2(lnbn_nodes(inode))-sl*E_inf+nscbc_u_inf*(E_inf+nscbc_p_inf)-aux_u2(lnbn_nodes(inode),1)*(aux_E2(lnbn_nodes(inode))+aux_p2(lnbn_nodes(inode))))/(sr-sl)

                     aux_q(inode,1) = q_hll

                     aux_rho(inode) = rho_hll
                     aux_E(inode) = E_hll

                     sl = min(-nscbc_c_inf, aux_u2(lnbn_nodes(inode),2) - sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))))
                     sr =  max(aux_u2(lnbn_nodes(inode),2) + sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))), nscbc_c_inf)
                     q_hll   = (sr*aux_q2(lnbn_nodes(inode),2)-aux_u2(lnbn_nodes(inode),2)*aux_q2(lnbn_nodes(inode),2))/(sr-sl)

                     aux_q(inode,2) = q_hll

                     sl = min(-nscbc_c_inf, aux_u2(lnbn_nodes(inode),3) - sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))))
                     sr =  max(aux_u2(lnbn_nodes(inode),3) + sqrt(nscbc_gamma_inf*aux_p2(lnbn_nodes(inode))/aux_rho2(lnbn_nodes(inode))), nscbc_c_inf)
                     q_hll   = (sr*aux_q2(lnbn_nodes(inode),3)-aux_u2(lnbn_nodes(inode),2)*aux_q2(lnbn_nodes(inode),3))/(sr-sl)

                     aux_q(inode,3) = q_hll

                     aux_u(inode,1) = aux_q2(inode,1)/aux_rho2(inode)
                     aux_u(inode,2) = aux_q2(inode,2)/aux_rho2(inode)
                     aux_u(inode,3) = aux_q2(inode,3)/aux_rho2(inode)

                     aux_p(inode) = aux_rho2(inode)*(nscbc_gamma_inf-1.0_rp)*((aux_E2(inode)/aux_rho2(inode))- &
                        0.5_rp*dot_product(aux_u2(inode,:),aux_u2(inode,:)))
                  else if (bcode == bc_type_slip_wall_model) then ! slip wall model
                     aux_q(inode,2) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_rho(inode) = nscbc_rho_inf
                     !aux_E(inode) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
                  end if
               end if
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_prim

      end module mod_bc_routines
