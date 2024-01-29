module mod_bc_routines_incomp

   use mod_mpi
   use mod_numerical_params
   use mod_comms
   use mod_comms_boundaries
   use mod_nvtx

   implicit none

      contains
               subroutine temporary_bc_routine_dirichlet_prim_residual_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,aux_u,u_buffer)

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
                  if (bcode == bc_type_far_field) then ! inlet just for aligened inlets with x
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  else if (bcode == bc_type_non_slip_adiabatic) then ! non_slip wall adiabatic
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  else if ((bcode == bc_type_slip_wall_model) .or. (bcode == bc_type_slip_adiabatic)) then ! slip
                     norm = (normalsAtNodes(inode,1)*aux_u(inode,1)) + (normalsAtNodes(inode,2)*aux_u(inode,2)) + (normalsAtNodes(inode,3)*aux_u(inode,3))
                     !$acc loop seq
                     do idime = 1,ndime     
                        aux_u(inode,idime) = aux_u(inode,idime) - norm*normalsAtNodes(inode,idime)
                     end do
                  end if
               end if ! This guy
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_prim_residual_incomp

         subroutine temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,aux_u,u_buffer)

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
                  if (bcode == bc_type_far_field) then ! inlet just for aligened inlets with x
                     aux_u(inode,1) = u_buffer(inode,1)
                     aux_u(inode,2) = u_buffer(inode,2)
                     aux_u(inode,3) = u_buffer(inode,3)
                  else if (bcode == bc_type_non_slip_adiabatic) then ! non_slip wall adiabatic
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  else if ((bcode == bc_type_slip_wall_model) .or. (bcode == bc_type_slip_adiabatic)) then ! slip
                     norm = (normalsAtNodes(inode,1)*aux_u(inode,1)) + (normalsAtNodes(inode,2)*aux_u(inode,2)) + (normalsAtNodes(inode,3)*aux_u(inode,3))
                     !$acc loop seq
                     do idime = 1,ndime     
                        aux_u(inode,idime) = aux_u(inode,idime) - norm*normalsAtNodes(inode,idime)
                     end do
                  end if
               end if
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_prim_incomp

          subroutine temporary_bc_routine_dirichlet_pressure_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,aux_p)

            implicit none

            integer(4), intent(in)     :: npoin, nboun,bou_codes_nodes(npoin)
            real(rp), intent(in)     :: normalsAtNodes(npoin,ndime)
            real(rp),    intent(inout) :: aux_p(npoin)
            integer(4)                 :: iboun,bcode,ipbou,inode

            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  bcode = bou_codes_nodes(inode) ! Boundary element code
                  if (bcode == bc_type_outlet_incomp) then 
                     aux_p(inode) = 0.0_rp
                  end if
               end if
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_pressure_incomp

         subroutine bc_routine_pressure_flux(nelem,npoin,nboun,connec,bound,point2elem,bou_code,bou_codes_nodes, &
                                             bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,omega,bpress)

            implicit none

            integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun),bou_codes_nodes(npoin)
            integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
            real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
            real(rp),   intent(in)  :: rho(npoin),omega(npoin,ndime),mu_fluid(npoin)
            real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
            real(rp),   intent(inout) :: bpress(npoin)
            integer(4)              :: idime,igaus,bcode, inode, ielem,ibound
            real(rp)                :: bnorm(npbou*ndime),nmag
            real(rp)                :: aux(ndime)
            real(rp)                 :: auxmag, sig            


            !$acc parallel loop gang private(bnorm) 
            do ibound = 1,nboun
               bcode = bou_code(ibound)
               ! Boundary element code
               if (bcode == bc_type_non_slip_adiabatic) then 
                  bnorm(:) = bounorm(ibound,:)
                   !$acc loop vector private(aux)
                  do igaus = 1,npbou
                     aux(1) = bnorm((igaus-1)*ndime+1)
                     aux(2) = bnorm((igaus-1)*ndime+2)
                     aux(3) = bnorm((igaus-1)*ndime+3)
                     auxmag = sqrt(aux(1)*aux(1) + aux(2)*aux(2)  +  aux(3)*aux(3))
                     aux(:) = aux(:)/auxmag
                     inode = bound(ibound,igaus)
                     ielem = point2elem(inode)
                     nmag = dot_product(aux(:), mu_fluid(inode)*omega(inode,:))
                     if(dot_product(coord(connec(ielem,nnode),:)-coord(inode,:), aux(:)) .lt. 0.0_rp ) then
					         sig=-1.0_rp
					      end if
                     !$acc atomic update
                     bpress(inode) = bpress(inode)+auxmag*wgp_b(igaus)*nmag*sig
                     !$acc end atomic
                     !write(111,*) "correcion  ",auxmag*wgp_b(igaus)*nmag*sig
                  end do
               end if
            end do
            !$acc end parallel loop
           

         end subroutine bc_routine_pressure_flux

      end module mod_bc_routines_incomp
