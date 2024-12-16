module mod_bc_routines_meshElasticity

   use mod_mpi
   use mod_numerical_params
   use mod_comms
   use mod_comms_boundaries
   use mod_nvtx

   implicit none

      contains
            
         subroutine temporary_bc_routine_dirichlet_prim_residual_meshElasticity(npoin,nboun,bou_codes_nodes,normalsAtNodes,aux_u,u_buffer)

            implicit none

            integer(4), intent(in)     :: npoin, nboun,  bou_codes_nodes(npoin)
            real(rp), intent(in)     :: normalsAtNodes(npoin,ndime),u_buffer(npoin,ndime)
            real(rp),    intent(inout) :: aux_u(npoin,ndime)
            integer(4)                 :: iboun,bcode,ipbou,inode,idime,iBoundNode
            real(rp)                   :: cin,R_plus,R_minus,v_b,c_b,s_b,rho_b,p_b,rl,rr, sl, sr
            real(rp)                   :: q_hll,rho_hll,E_hll,E_inf,norm

            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  bcode = bou_codes_nodes(inode) ! Boundary element code
                  if (bcode == bc_type_non_slip_adiabatic ) then ! non_slip wall adiabatic
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  end if
               end if ! This guy
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_prim_residual_meshElasticity

         subroutine temporary_bc_routine_dirichlet_prim_meshElasticity(npoin,nboun,bou_codes_nodes,lnbn_nodes,normalsAtNodes,aux_u,u_buffer)

            implicit none

            integer(4), intent(in)  :: npoin,nboun,bou_codes_nodes(npoin)
            integer(4), intent(in)  :: lnbn_nodes(npoin)
            real(rp), intent(in)    :: normalsAtNodes(npoin,ndime),u_buffer(npoin,ndime)
            real(rp), intent(inout) :: aux_u(npoin,ndime)
            integer(4)              :: iboun,bcode,ipbou,inode,idime,iBoundNode
            real(rp)                :: cin,R_plus,R_minus,v_b,c_b,s_b,rho_b,p_b,rl,rr, sl, sr
            real(rp)                :: q_hll,rho_hll,E_hll,E_inf,norm

            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  bcode = bou_codes_nodes(inode) ! Boundary element code
                  if (bcode == bc_type_non_slip_adiabatic ) then 
                     aux_u(inode,1) = u_buffer(inode,1)
                     aux_u(inode,2) = u_buffer(inode,2)
                     aux_u(inode,3) = u_buffer(inode,3)
                  end if
               end if
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_prim_meshElasticity
          
      end module mod_bc_routines_meshElasticity
