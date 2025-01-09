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
                  if (bcode == bc_type_far_field .or. bcode == bc_type_far_field_SB .or. bcode == bc_type_unsteady_inlet) then ! inlet just for aligened inlets with x
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  else if (bcode == bc_type_non_slip_adiabatic .or. bcode == bc_type_non_slip_isothermal) then ! non_slip wall adiabatic
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  else if (bcode == bc_type_recirculation_inlet) then ! non_slip wall adiabatic
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp                     
                  else if ((bcode == bc_type_slip_wall_model) .or. (bcode == bc_type_slip_adiabatic).or. (bcode == bc_type_symmetry)) then ! slip
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

         subroutine temporary_bc_routine_dirichlet_prim_incomp(npoin,nboun,bou_codes_nodes,lnbn_nodes,normalsAtNodes,aux_u,u_buffer)

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
                  if (bcode == bc_type_far_field .or. bcode == bc_type_far_field_SB .or. bcode == bc_type_unsteady_inlet) then ! inlet just for aligened inlets with x
                     aux_u(inode,1) = u_buffer(inode,1)
                     aux_u(inode,2) = u_buffer(inode,2)
                     aux_u(inode,3) = u_buffer(inode,3)
                  else if (bcode == bc_type_non_slip_adiabatic .or. bcode == bc_type_non_slip_isothermal) then ! non_slip wall adiabatic
                     aux_u(inode,1) = 0.0_rp
                     aux_u(inode,2) = 0.0_rp
                     aux_u(inode,3) = 0.0_rp
                  else if (bcode == bc_type_recirculation_inlet) then ! recirculation inlet
                     aux_u(inode,1) = aux_u(lnbn_nodes(inode),1)
                     aux_u(inode,2) = aux_u(lnbn_nodes(inode),2)
                     aux_u(inode,3) = aux_u(lnbn_nodes(inode),3)
                  else if ((bcode == bc_type_slip_wall_model) .or. (bcode == bc_type_slip_adiabatic).or. (bcode == bc_type_symmetry)) then ! slip
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

          subroutine temporary_bc_routine_dirichlet_pressure_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,aux_p,p_buffer)

            implicit none

            integer(4), intent(in)     :: npoin, nboun,bou_codes_nodes(npoin)
            real(rp), intent(in)     :: normalsAtNodes(npoin,ndime), p_buffer(npoin)
            real(rp),    intent(inout) :: aux_p(npoin)
            integer(4)                 :: iboun,bcode,ipbou,inode

            !$acc parallel loop  
            do inode = 1,npoin
               if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
                  bcode = bou_codes_nodes(inode) ! Boundary element code
                  if (bcode == bc_type_outlet_incomp) then 
                     aux_p(inode) = p_buffer(inode)
                  end if
               end if
            end do
            !$acc end parallel loop

         end subroutine temporary_bc_routine_dirichlet_pressure_incomp

         subroutine temporary_bc_routine_dirichlet_pressure_residual_incomp(npoin,nboun,bou_codes_nodes,normalsAtNodes,aux_p)

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

         end subroutine temporary_bc_routine_dirichlet_pressure_residual_incomp

         subroutine bc_routine_pressure_flux(nelem,npoin,nboun,connec,bound,point2elem,bou_code,bou_codes_nodes,numBoundCodes,bouCodes2BCType, &
                                             bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,mu_sgs,rho,omega,bpress)

            implicit none

            integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun),bou_codes_nodes(npoin),numBoundCodes
            integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin)
            real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode),bouCodes2BCType(numBoundCodes)
            real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
            real(rp),   intent(in)  :: rho(npoin),omega(npoin,ndime),mu_fluid(npoin),mu_sgs(ielem,ngaus)
            real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
            real(rp),   intent(inout) :: bpress(npoin)
            integer(4)              :: idime,igaus,bcode, inode, ielem,ibound,igausV
            real(rp)                :: bnorm(npbou*ndime),nmag
            real(rp)                :: aux(ndime)
            real(rp)                 :: auxmag, sig            

            !$acc kernels
            bpress(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(bnorm) 
            do ibound = 1,nboun
               bcode = bouCodes2BCType(bou_code(ibound))
               ! Boundary element code
               if ((bcode == bc_type_non_slip_adiabatic) .or. (bcode == bc_type_slip_wall_model) .or. (bcode == bc_type_slip_adiabatic).or. (bcode == bc_type_symmetry)) then 
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
                     nmag = dot_product(aux(:), omega(inode,:))
                     igausV = minloc(abs(connec(ielem,:)-inode),1)
                     if(dot_product(coord(connec(ielem,nnode),:)-coord(inode,:), aux(:)) .lt. 0.0_rp ) then
                        sig=1.0_rp
                     else
                        sig=-1.0_rp
                     end if
                     !$acc atomic update
                     bpress(inode) = bpress(inode)+auxmag*wgp_b(igaus)*nmag*sig*(mu_fluid(inode)+mu_sgs(ielem,igausV))
                     !$acc end atomic
                  end do
               end if
            end do
            !$acc end parallel loop
            if(mpi_size.ge.2) then
               call nvtxStartRange("MPI_comms_tI")
               call mpi_halo_atomic_update_real(bpress(:))
               call nvtxEndRange
            end if
            call nvtxEndRange
         end subroutine bc_routine_pressure_flux                   

         subroutine bc_routine_momentum_flux(nelem,npoin,nboun,connec,bound,point2elem,bou_code,bou_codes_nodes,numBoundCodes,bouCodes2BCType, &
                                             bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,rho,u_flux_buffer,R)

            implicit none

            integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun),bou_codes_nodes(npoin),numBoundCodes
            integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin),bouCodes2BCType(numBoundCodes)
            real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
            real(rp),   intent(in)  :: rho(npoin),u_flux_buffer(npoin,ndime),mu_fluid(npoin)
            real(rp),   intent(in)  :: coord(npoin,ndime), gpvol(1,ngaus,nelem)
            real(rp),   intent(inout) :: R(npoin,ndime)
            integer(4)              :: idime,igaus,bcode, inode, ielem,ibound
            real(rp)                :: bnorm(npbou*ndime),nmag
            real(rp)                :: aux(ndime)
            real(rp)                 :: auxmag, sig , fact         

            !$acc kernels
            R(:,:) = 0.0_rp
            !$acc end kernels
            !$acc parallel loop gang private(bnorm) 
            do ibound = 1,nboun
               bcode = bouCodes2BCType(bou_code(ibound))
               ! Boundary element code
               if (bcode == bc_type_outlet_incomp) then 
                  bnorm(:) = bounorm(ibound,:)
                  !$acc loop vector private(aux)
                  do igaus = 1,npbou
                     aux(1) = bnorm((igaus-1)*ndime+1)
                     aux(2) = bnorm((igaus-1)*ndime+2)
                     aux(3) = bnorm((igaus-1)*ndime+3)
                     auxmag = sqrt(aux(1)*aux(1) + aux(2)*aux(2)  +  aux(3)*aux(3))
                     inode = bound(ibound,igaus)
                     ielem = point2elem(inode)
                     if(dot_product(coord(connec(ielem,nnode),:)-coord(inode,:), aux(:)) .lt. 0.0_rp ) then
                        sig=1.0_rp
                     else
                        sig=-1.0_rp
                     end if
                     !$acc loop seq
                     do idime = 1,ndime
                        !$acc atomic update
                        R(inode,idime) = R(inode,idime)+auxmag*wgp_b(igaus)*sig*u_flux_buffer(bound(ibound,igaus),idime)
                        !$acc end atomic
                     end do
                  end do
               end if
            end do
            !$acc end parallel loop
         end subroutine bc_routine_momentum_flux         

         subroutine copy_periodicNodes_for_mappedInlet_incomp(uField)
            implicit none
            real(rp),intent(inout) :: uField(numNodesRankPar,ndime)
            integer(4) :: iPer

            !$acc parallel loop
            do iPer = 1,nPerRankPar
               uField(masSlaRankPar(iPer,2),1) = uField(masSlaRankPar(iPer,1),1)
               uField(masSlaRankPar(iPer,2),2) = uField(masSlaRankPar(iPer,1),2)
               uField(masSlaRankPar(iPer,2),3) = uField(masSlaRankPar(iPer,1),3)
            end do
            !$acc end parallel loop
         end subroutine copy_periodicNodes_for_mappedInlet_incomp

         subroutine evalPAtOutlet(nelem,npoin,npoin_w,nboun,connec,bound,point2elem,bou_code,bou_codes_nodes,lpoin_w, &
            bounorm,normalsAtNodes,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,wgp_b,coord,dlxigp_ip,He,gpvol,mu_fluid,mu_e,mu_sgs,rho,u,p_buffer,u_flux_buffer)
   
         implicit none
   
         integer(4), intent(in)  :: npoin,nboun,bound(nboun,npbou),bou_code(nboun),bou_codes_nodes(npoin),npoin_w
         integer(4), intent(in)  :: nelem,connec(nelem,nnode),point2elem(npoin),lpoin_w(npoin_w)
         real(rp),   intent(in)  :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
         integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
         real(rp),   intent(in)  :: dlxigp_ip(ngaus,ndime,porder+1), He(ndime,ndime,ngaus,nelem)
         real(rp),   intent(in)  :: rho(npoin),u(npoin,ndime),mu_fluid(npoin)
         real(rp),   intent(in)  :: coord(npoin,ndime),gpvol(1,ngaus,nelem),mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
         real(rp),   intent(inout) :: p_buffer(npoin)
         real(rp),   intent(inout) :: u_flux_buffer(npoin,ndime)
         real(rp)                :: gradIsoU(ndime,ndime), gradU(ndime,ndime),tau(ndime,ndime)
         integer(4)              :: iBound,iElem,idime,igaus,iAux,inode,bcode
         integer(4)              :: jdime, isoI, isoJ, isoK,kdime,ii,ipoin
         real(rp)                :: normal(ndime),aux_p, aux(ndime),aux_u,aux_s,aux_ufb,auxmag,mu_fgp,sig
   
         !$acc parallel loop private(normal,gradIsoU,gradU,aux,tau)
         do ipoin = 1,npoin_w
            inode = lpoin_w(ipoin)
            if(bou_codes_nodes(inode) .lt. max_num_bou_codes) then
               bcode = bou_codes_nodes(inode) ! Boundary element code
               if (bcode == bc_type_outlet_incomp) then 
                  iElem = point2elem(inode)
                  normal(1:ndime) = normalsAtNodes(inode,1:ndime)
                  auxmag = sqrt(normal(1)*normal(1) + normal(2)*normal(2)  +  normal(3)*normal(3))
                  normal(:) = normal(:)/auxmag
                  if(dot_product(coord(connec(ielem,nnode),:)-coord(inode,:), normal(:)) .lt. 0.0_rp ) then
                     sig=1.0_rp
                  else
                     sig=-1.0_rp
                  end if
                  normal(:) = sig*normal(:)

                  igaus = minloc(abs(connec(iElem,:)-inode),1)
                  mu_fgp = mu_fluid(inode)+ mu_sgs(iElem,igaus)+mu_e(iElem,igaus)


                  isoI = gmshAtoI(igaus) 
                  isoJ = gmshAtoJ(igaus) 
                  isoK = gmshAtoK(igaus) 
   
                  gradIsoU(:,:) = 0.0_rp
                  !$acc loop seq
                  do ii=1,porder+1
                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*u(connec(iElem,invAtoIJK(ii,isoJ,isoK)),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*u(connec(iElem,invAtoIJK(isoI,ii,isoK)),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*u(connec(iElem,invAtoIJK(isoI,isoJ,ii)),idime)
                     end do
                  end do
                     
                  gradU(:,:) = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                        !$acc loop seq
                        do kdime=1,ndime
                           gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                        end do
                     end do
                  end do

                  !$acc loop seq
                  do idime = 1,ndime
                     !$acc loop seq
                     do jdime = 1,ndime
                          tau(idime,jdime) = (gradU(idime,jdime)+gradU(jdime,idime))
                     end do
                  end do

                  aux(:) = 0.0_rp                    
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                        aux(idime) = aux(idime) + normal(jdime)*tau(idime,jdime)
                     end do
                  end do

                  aux_u = 0.0_rp
                  aux_s = 0.0_rp
                  aux_p = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     aux_p = aux_p + aux(idime)*normal(idime)
                     aux_u = aux_u + u(inode,idime)*u(inode,idime)
                     aux_s = aux_s + normal(idime)*u(inode,idime)
                  end do 
                  aux_p = (mu_fgp*aux_p - 0.5_rp*aux_u*(0.5_rp*(1.0_rp-tanh(aux_s/(nscbc_u_inf*nscbc_delta)))))
                  p_buffer(inode) = aux_p
                  !aux_ufb = gradU(1,1)+gradU(2,2)+gradU(3,3)
                  !!$acc loop seq
                  !do idime=1, ndime
                  !   u_flux_buffer(inode,idime) =  (aux_p*normal(idime) +  0.5_rp*aux_u*(0.5_rp*(1.0_rp-tanh(aux_s/(nscbc_u_inf*nscbc_delta))))*normal(idime) - mu_fgp*aux_ufb*normal(idime))/mu_fgp
                  !   !u_flux_buffer(inode,idime) =  (aux_p*normal(idime) +  0.5_rp*aux_u*(0.5_rp*(1.0_rp-tanh(aux_s/(nscbc_u_inf*nscbc_delta))))*normal(idime))/mu_fgp
                  !end do 
               end if
            end if
         end do
         if(mpi_size.ge.2) then
            call mpi_halo_max_boundary_update_real_iSendiRcv(p_buffer)
         end if
         !!$acc end parallel loop   
         !if(mpi_size.ge.2) then
         !   call nvtxStartRange("MPI_comms_tI")
         !   call mpi_halo_max_boundary_update_real_arrays_iSendiRcv(ndime,u_flux_buffer(:,:))
         !   !call mpi_halo_max_boundary_update_real_iSendiRcv(u_flux_buffer(:,1))
         !   !call mpi_halo_max_boundary_update_real_iSendiRcv(u_flux_buffer(:,2))
         !   !call mpi_halo_max_boundary_update_real_iSendiRcv(u_flux_buffer(:,3))
         !   call nvtxEndRange
         !end if
      end subroutine evalPAtOutlet

      end module mod_bc_routines_incomp
