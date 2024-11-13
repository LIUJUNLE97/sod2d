module mod_bc_routines_species

    use mod_mpi
    use mod_numerical_params
    use mod_comms
    use mod_comms_boundaries
    use mod_nvtx
 
    implicit none

    contains

    subroutine temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,aux_Yk,Yk_buffer)

        implicit none

        integer(4), intent(in)     :: npoin, nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
        integer(4), intent(in)     :: nbnodes, lbnodes(nbnodes),lnbn_nodes(npoin)
        real(rp), intent(in)     :: normalsAtNodes(npoin,ndime),Yk_buffer(npoin)
        real(rp),    intent(inout) :: aux_Yk(npoin)
        integer(4)                 :: iboun,bcode,ipbou,inode,idime,iBoundNode

        !$acc parallel loop  
        do inode = 1,npoin
           if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
              bcode = bou_codes_nodes(inode) ! Boundary element code
              if ((bcode == bc_type_far_field) .or. (bcode == bc_type_non_slip_isothermal) .or. (bcode == bc_type_slip_adiabatic).or. (bcode == bc_type_slip_wall_model)) then
                 aux_Yk(inode) = Yk_buffer(inode)
              end if
            end if
        end do
        !$acc end parallel loop

    end subroutine temporary_bc_routine_dirichlet_species

    subroutine temporary_bc_routine_dirichlet_residual_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,aux_Yk,Yk_buffer)

        implicit none

        integer(4), intent(in)     :: npoin, nboun, bou_codes(nboun), bou_codes_nodes(npoin), bound(nboun,npbou)
        integer(4), intent(in)     :: nbnodes, lbnodes(nbnodes),lnbn_nodes(npoin)
        real(rp), intent(in)     :: normalsAtNodes(npoin,ndime),Yk_buffer(npoin)
        real(rp),    intent(inout) :: aux_Yk(npoin)
        integer(4)                 :: iboun,bcode,ipbou,inode,idime,iBoundNode

        !$acc parallel loop  
        do inode = 1,npoin
           if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
              bcode = bou_codes_nodes(inode) ! Boundary element code
              if ((bcode == bc_type_far_field) .or. (bcode == bc_type_non_slip_isothermal)  .or. (bcode == bc_type_slip_adiabatic)) then
                 aux_Yk(inode) = 0.0_rp
              end if
            end if
        end do
        !$acc end parallel loop

    end subroutine temporary_bc_routine_dirichlet_residual_species


end module mod_bc_routines_species