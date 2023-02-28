module mod_bc_routines

   use mod_mpi
   use mod_constants
   use mod_comms
   use mod_comms_boundaries
   use mod_nvtx

   implicit none

      contains

         subroutine temporary_bc_routine_dirichlet_prim(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn,lnbn_nodes,normalsAtNodes,aux_rho,aux_q,aux_u,aux_p,aux_E,u_buffer)

            implicit none

            integer(4), intent(in)     :: npoin, nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
            integer(4), intent(in)     :: nbnodes, lbnodes(nbnodes),lnbn(nboun,npbou),lnbn_nodes(npoin)
            real(rp), intent(in)     :: normalsAtNodes(npoin,ndime),u_buffer(npoin,ndime)
            real(rp),    intent(inout) :: aux_rho(npoin),aux_q(npoin,ndime),aux_u(npoin,ndime),aux_p(npoin),aux_E(npoin)
            real(rp)                   :: aux_rho2(npoin),aux_q2(npoin,ndime),aux_u2(npoin,ndime),aux_p2(npoin),aux_E2(npoin)
            integer(4)                 :: iboun,bcode,ipbou,inode,idime,iBoundNode
            real(rp)                   :: cin,R_plus,R_minus,v_b,c_b,s_b,rho_b,p_b,rl,rr, sl, sr
            real(rp)                   :: q_hll,rho_hll,E_hll,E_inf,norm

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
               call mpi_halo_max_boundary_update_real_iSendiRcv(aux_q2(:,1))
               call mpi_halo_max_boundary_update_real_iSendiRcv(aux_q2(:,2))
               call mpi_halo_max_boundary_update_real_iSendiRcv(aux_q2(:,3))
               call mpi_halo_max_boundary_update_real_iSendiRcv(aux_rho2(:))
               call mpi_halo_max_boundary_update_real_iSendiRcv(aux_E2(:))
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
                  if (bcode == bc_type_far_field) then ! inlet just for aligened inlets with x

                     E_inf = (nscbc_rho_inf*0.5_rp*nscbc_u_inf**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp))

                     sl = min(nscbc_u_inf-nscbc_c_inf, aux_u2(inode,1) - sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)))
                     sr =  max(aux_u2(inode,1) + sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)), nscbc_u_inf+nscbc_c_inf)

                     rho_hll = (sr*aux_rho2(inode)-sl*nscbc_rho_inf+nscbc_rho_inf*nscbc_u_inf-aux_q2(inode,1))/(sr-sl)
                     E_hll   = (sr*aux_E2(inode)-sl*E_inf+nscbc_u_inf*(E_inf+nscbc_p_inf)-aux_u2(inode,1)*(aux_E2(inode)+aux_p2(inode)))/(sr-sl)

                     sl = min(u_buffer(inode,1)-nscbc_c_inf, aux_u2(inode,1) - sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)))
                     sr =  max(aux_u2(inode,1) + sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)), u_buffer(inode,1)+nscbc_c_inf)

                     q_hll   = (sr*aux_q2(inode,1)-sl*nscbc_rho_inf*u_buffer(inode,1)+nscbc_rho_inf*u_buffer(inode,1)**2-aux_u2(inode,1)*aux_q2(inode,1))/(sr-sl)

                     aux_rho(inode) = rho_hll
                     aux_E(inode) = E_hll
                     aux_q(inode,1) = q_hll

                     sl = min(u_buffer(inode,2)-nscbc_c_inf, aux_u2(inode,2) - sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)))
                     sr =  max(aux_u2(inode,2) + sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)), u_buffer(inode,2)+nscbc_c_inf)
                     q_hll   = (sr*aux_q2(inode,2)-sl*nscbc_rho_inf*u_buffer(inode,2)+nscbc_rho_inf*u_buffer(inode,2)**2-aux_u2(inode,2)*aux_q2(inode,2))/(sr-sl)


                     aux_q(inode,2) = q_hll

                     sl = min(u_buffer(inode,3)-nscbc_c_inf, aux_u2(inode,3) - sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)))
                     sr =  max(aux_u2(inode,3) + sqrt(nscbc_gamma_inf*aux_p2(inode)/aux_rho2(inode)), u_buffer(inode,3)+nscbc_c_inf)
                     q_hll   = (sr*aux_q2(inode,3)-sl*nscbc_rho_inf*u_buffer(inode,3)+nscbc_rho_inf*u_buffer(inode,3)**2-aux_u2(inode,3)*aux_q2(inode,3))/(sr-sl)

                     aux_q(inode,3) = q_hll

                     aux_u(inode,1) = aux_q(inode,1)/aux_rho(inode)
                     aux_u(inode,2) = aux_q(inode,2)/aux_rho(inode)
                     aux_u(inode,3) = aux_q(inode,3)/aux_rho(inode)

                     aux_p(inode) = aux_rho(inode)*(nscbc_gamma_inf-1.0_rp)*((aux_E(inode)/aux_rho(inode))- &
                        0.5_rp*dot_product(aux_u(inode,:),aux_u(inode,:)))
                  else if (bcode == bc_type_non_slip_adiabatic) then ! non_slip wall adiabatic
                     
                     aux_q(inode,1) = 0.0_rp
                     aux_q(inode,2) = 0.0_rp
                     aux_q(inode,3) = 0.0_rp

                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp

                     aux_rho(inode) = nscbc_rho_inf
                     !aux_E(inode) = nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)
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

                     aux_q(inode,1) = nscbc_rho_inf*u_buffer(inode,1)
                     aux_u(inode,1) = u_buffer(inode,1)
                     aux_q(inode,2) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_q(inode,3) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp

                     aux_p(inode) = nscbc_p_inf
                     aux_rho(inode) = nscbc_rho_inf
                     aux_E(inode) = nscbc_rho_inf*0.5_rp*u_buffer(inode,1)**2 + nscbc_p_inf/(nscbc_gamma_inf-1.0_rp)

                  else if (bcode == bc_type_slip_wall_model) then ! slip wall model
                     norm = dot_product(normalsAtNodes(inode,:),aux_q(inode,:))
                     !$acc loop seq
                     do idime = 1,ndime     
                        aux_q(inode,idime) = aux_q(inode,idime) - norm*normalsAtNodes(inode,idime)
                     end do
                     aux_rho(inode) = nscbc_rho_inf

                     aux_u(inode,1) = aux_q(inode,1)/aux_rho(inode)
                     aux_u(inode,2) = aux_q(inode,2)/aux_rho(inode)
                     aux_u(inode,3) = aux_q(inode,3)/aux_rho(inode)

                     aux_E(inode) = aux_p(inode)/(nscbc_gamma_inf-1.0_rp)+ &
                                    aux_rho(inode)*0.5_rp*dot_product(aux_u(inode,:),aux_u(inode,:))

                  end if
               end if
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_prim

      end module mod_bc_routines
