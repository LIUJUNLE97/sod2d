
module time_integ_species_rk4

   use mod_nvtx
   use elem_convec_species
   use elem_convec , only : generic_scalar_convec_ijk
   use elem_diffu_species
   use mod_solver
   use mod_entropy_viscosity_species
   use mod_numerical_params
   use mod_bc_routines_species

   implicit none

   real(rp), allocatable, dimension(:)   :: RYk,Reta
   real(rp), allocatable, dimension(:)   :: aux_Yk,aux_eta
   real(rp), allocatable, dimension(:)   :: RYk_sum,Rdiff_Yk,Reta_sum
   real(rp), allocatable, dimension(:,:) :: f_eta

   real(rp), allocatable, dimension(:)   :: a_i, b_i, c_i,b_i2
   real(rp), allocatable, dimension(:,:) :: a_ij
   logical :: firstTimeStep = .true.

   contains

   subroutine init_species_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin
      integer(4) :: numSteps

      call nvtxStartRange("Init RK4 species solver")

      allocate(RYk(npoin),Reta(npoin))
      !$acc enter data create(RYk(:))
      !$acc enter data create(Reta(:))

      allocate(aux_Yk(npoin),aux_eta(npoin))
      !$acc enter data create(aux_Yk(:))
      !$acc enter data create(aux_eta(:))

      allocate(RYk_sum(npoin),Reta_sum(npoin),Rdiff_Yk(npoin))
      !$acc enter data create(RYk_sum(:))
      !$acc enter data create(Reta_sum(:))
      !$acc enter data create(Rdiff_Yk(:))
      allocate(f_eta(npoin,ndime))
      !$acc enter data create(f_eta(:,:))

      allocate(a_i(4),b_i(4),c_i(4))
      !$acc enter data create(a_i(:))
      !$acc enter data create(b_i(:))
      !$acc enter data create(c_i(:))

      if (flag_rk_order == 1) then
         a_i = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
         c_i = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
         b_i = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
      else if (flag_rk_order == 2) then
         a_i = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
         c_i = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
         b_i = [0.5_rp, 0.5_rp, 0.0_rp, 0.0_rp]
      else if (flag_rk_order == 3) then
         write(1,*) "--| NOT CODED FOR RK3 YET!"
         stop 1
      else if (flag_rk_order == 4) then
         a_i = [0.0_rp, 0.5_rp, 0.5_rp, 1.0_rp]
         c_i = [0.0_rp, 0.5_rp, 0.5_rp, 1.0_rp]
         b_i = [1.0_rp/6.0_rp, 1.0_rp/3.0_rp, 1.0_rp/3.0_rp, 1.0_rp/6.0_rp]
      else
         write(1,*) "--| NOT CODED FOR RK > 4 YET!"
         stop 1
      end if
      !$acc update device(a_i(:))
      !$acc update device(b_i(:))
      !$acc update device(c_i(:))

      call nvtxEndRange

   end subroutine init_species_solver

   subroutine end_species_solver()
      implicit none

      ! to be done

      !$acc exit data delete(a_i(:))
      !$acc exit data delete(b_i(:))
      !$acc exit data delete(c_i(:))
      deallocate(a_i,b_i,c_i)

   end subroutine end_species_solver

         subroutine rk_4_species_main(ispc,noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
            ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
            rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
            ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
            listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer)  ! Optional args

            implicit none

            integer(4),           intent(in)    :: ispc
            logical,              intent(in)   :: noBoundaries,isWallModelOn
            integer(4),           intent(in)    :: nelem, nboun, npoin
            integer(4),           intent(in)    :: connec(nelem,nnode), npoin_w, lpoin_w(npoin_w),point2elem(npoin),lnbn_nodes(npoin)
            integer(4),           intent(in)    :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            integer(4),           intent(in)    :: ppow
            real(rp),             intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),             intent(in)    :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime)
            real(rp),             intent(in)    :: gpvol(1,ngaus,nelem)
            real(rp),             intent(in)    :: dt, helem(nelem)
            real(rp),             intent(in)    :: helem_l(nelem,nnode)
            real(rp),             intent(in)    :: Ml(npoin)
            real(rp),             intent(in)    :: Cp, Prt
            real(rp),             intent(inout) :: rho(npoin,4)
            real(rp),             intent(inout) :: u(npoin,ndime,4)
            real(rp),             intent(inout) :: Yk(npoin,nspecies,4)
            real(rp),             intent(inout) :: eta_Yk(npoin,nspecies,4)
            real(rp),             intent(inout) :: mu_fluid(npoin)
            real(rp),             intent(inout) :: mu_e_Yk(nelem,ngaus,nspecies)
            real(rp),             intent(inout) :: mu_sgs(nelem,ngaus)
            real(rp),             intent(inout) :: mue_l(nelem,nnode)
            real(rp),             intent(in)    :: coord(npoin,ndime)
            real(rp),             intent(in)  ::  wgp(ngaus)
            integer(4),            intent(in)    :: numBoundsWM
            integer(4), optional, intent(in)    :: ndof, nbnodes, ldof(*), lbnodes(*)
            integer(4), optional, intent(in)    :: bound(nboun,npbou), bou_codes(nboun), bou_codes_nodes(npoin)
            integer(4), optional, intent(in)    :: listBoundsWM(*)
            real(rp), optional, intent(in)      :: wgp_b(npbou), bounorm(nboun,ndime*npbou),normalsAtNodes(npoin,ndime)
            real(rp), optional,   intent(in)    :: Yk_buffer(npoin,nspecies)
            integer(4)                          :: pos
            integer(4)                          :: istep, ipoin, idime,icode


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value
            if(firstTimeStep .eqv. .true.) then
               firstTimeStep = .false.
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*eta_Yk(lpoin_w(ipoin),ispc,1)
                  end do
               end do
               !$acc end parallel loop
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,eta_Yk(:,ispc,1),u(:,:,1),Reta(:))

               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:))
               end if
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:))   

               call nvtxStartRange("Entropy viscosity evaluation")
               call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta(:),Ngp,coord,dNgp,gpvol,wgp, &
                  rho(:,1),u(:,:,1),eta_Yk(:,ispc,1),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
               call nvtxEndRange
            end if


            !
            ! Initialize variables to zero
            !
            call nvtxStartRange("Initialize variables")
            !$acc kernels
            aux_Yk(1:npoin) = 0.0_rp
            aux_eta(1:npoin) = 0.0_rp
            Rdiff_Yk(1:npoin) = 0.0_rp
            RYk(1:npoin) = 0.0_rp
            RYk_sum(1:npoin) = 0.0_rp
            !$acc end kernels

            call nvtxEndRange
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")

            do istep = 1,flag_rk_order
               !
               ! Compute variable at substep (y_i = y_n+dt*A_ij*R_j)
               !
               call nvtxStartRange("Update aux_*")
               !$acc parallel loop
               do ipoin = 1,npoin
                  aux_Yk(ipoin) = Yk(ipoin,ispc,pos) - dt*a_i(istep)*RYk(ipoin)
               end do
               !$acc end parallel loop
               call nvtxEndRange

               !
               ! Apply bcs after update
               !
               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,aux_Yk(:),Yk_buffer(:,ispc))
                  call nvtxEndRange
               end if

               !
               ! Update velocity and equations of state
               !
               call nvtxStartRange("Update u and EOS")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  aux_eta(lpoin_w(ipoin)) = 0.5_rp*aux_Yk(lpoin_w(ipoin))*aux_Yk(lpoin_w(ipoin))
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime,1)*aux_eta(lpoin_w(ipoin))
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               call nvtxStartRange("Entropy convection")
               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta,aux_eta(:),u(:,:,1),Reta(:))
               call nvtxEndRange

               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(Reta(:))
                  call nvtxEndRange
               end if

               if (noBoundaries .eqv. .false.) then
                  call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Reta)
               end if

               call nvtxStartRange("Lumped solver")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta)
               call nvtxEndRange
               
               call nvtxStartRange("DIFFUSIONS")
               call species_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),aux_Yk(:),mu_fluid,mu_e_Yk(:,:,ispc),mu_sgs,Ml,Rdiff_Yk)
               call nvtxEndRange
               !
               !
               ! Compute convective terms
               !
               call nvtxStartRange("CONVECTIONS")
               call species_convec_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),aux_Yk(:),u(:,:,pos),RYk(:))               
               call nvtxEndRange

               call nvtxStartRange("Add convection and diffusion")
               !$acc kernels
               RYk(:) = RYk(:) + Rdiff_Yk(:)
               !$acc end kernels
               call nvtxEndRange

               !TESTING NEW LOCATION FOR MPICOMMS
               if(mpi_size.ge.2) then
                  call nvtxStartRange("MPI_comms_tI")
                  call mpi_halo_atomic_update_real(RYk(:))
                  call nvtxEndRange
               end if

               !
               ! Call lumped mass matrix solver
               !
               call nvtxStartRange("Call solver")
               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,RYk(:))
               call nvtxEndRange
               !
               ! Accumulate the residuals
               !
               call nvtxStartRange("Accumulate residuals")
               !$acc parallel loop
               do ipoin = 1,npoin
                  Reta_sum(ipoin)  = Reta_sum(ipoin)  + b_i(istep)*Reta(ipoin)
                  RYk_sum(ipoin) = RYk_sum(ipoin) + b_i(istep)*RYk(ipoin)
               end do
               !$acc end parallel loop
               call nvtxEndRange
            end do
            call nvtxEndRange
            !
            ! RK update to variables
            !
            call nvtxStartRange("RK_UPDATE")
            !$acc parallel loop
            do ipoin = 1,npoin
               Yk(ipoin,ispc,pos) = Yk(ipoin,ispc,pos)-dt*RYk_sum(ipoin)
            end do
            !$acc end parallel loop
            call nvtxEndRange

            !
            ! Apply bcs after update
            !
            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Yk(:,ispc,pos),Yk_buffer(:,ispc))
               call nvtxEndRange
            end if

            !
            ! Update velocity and equations of state
            !
            call nvtxStartRange("Update u and EOS")
            

            !$acc parallel loop
            do ipoin = 1,npoin_w
               eta_Yk(lpoin_w(ipoin),ispc,pos) = 0.5_rp*Yk(lpoin_w(ipoin),ispc,pos)*Yk(lpoin_w(ipoin),ispc,pos)
            end do
            !$acc end parallel loop
            call nvtxEndRange

            call nvtxStartRange("Entropy residual")
            !$acc parallel loop
            do ipoin = 1,npoin_w
               Reta(lpoin_w(ipoin)) =  -Reta_sum(lpoin_w(ipoin))-(eta_Yk(lpoin_w(ipoin),ispc,2)-eta_Yk(lpoin_w(ipoin),ispc,1))/dt
            end do
            !$acc end parallel loop
            call nvtxEndRange

            if (noBoundaries .eqv. .false.) then
               call nvtxStartRange("BCS_AFTER_UPDATE")
               call bc_fix_dirichlet_residual_entropy(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Reta)
               call nvtxEndRange
            end if
            

            call nvtxStartRange("Entropy viscosity evaluation")
            !
            ! Compute entropy viscosity
            !
            call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Ngp,coord,dNgp,gpvol,wgp, &
               rho(:,2),u(:,:,2),eta_Yk(:,ispc,2),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
            call nvtxEndRange

         end subroutine rk_4_species_main
               
      end module time_integ_species_rk4
