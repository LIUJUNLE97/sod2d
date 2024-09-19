
module time_integ_species

   use mod_nvtx
   use elem_convec_species
   use elem_convec , only : generic_scalar_convec_ijk
   use elem_diffu_species
   use mod_solver
   use mod_entropy_viscosity_species
   use mod_numerical_params
   use mod_bc_routines_species

   implicit none

   real(rp), allocatable, dimension(:,:,:) :: f_eta, Reta
   real(rp), allocatable, dimension(:,:) :: K2spc
   real(rp), allocatable, dimension(:,:) :: Rspc
   real(rp), allocatable, dimension(:) :: auxReta

   real(rp), allocatable, dimension(:)   :: a_i, b_i
   logical :: firstTimeStep = .true.

   contains

   subroutine init_rk4_species_solver(npoin)
      implicit none
      integer(4),intent(in) :: npoin

      write(*,*) "--| Allocating species"

      allocate(auxReta(npoin))
      !$acc enter data create(auxReta(:))

      allocate(K2spc(npoin,ndime),f_eta(npoin,ndime,nspecies),Reta(npoin,nspecies,2),Rspc(npoin,nspecies))
      !$acc enter data create(K2spc(:,:))
      !$acc enter data create(Rspc(:,:))
      !$acc enter data create(f_eta(:,:,:))
      !$acc enter data create(Reta(:,:,:))

      allocate(a_i(flag_rk_ls_stages),b_i(flag_rk_ls_stages))
      !$acc enter data create(a_i(:))
      !$acc enter data create(b_i(:))

      if (flag_rk_ls_stages == 5) then
         a_i = [0.0_rp, real( -0.4178904745d0,rp) , real( -1.192151694643d0,rp), real( -1.697784692471d0,rp), real( -1.514183444257d0,rp)]
         b_i = [real( 0.1496590219993d0,rp), real( 0.3792103129999d0,rp), real( 0.8229550293869d0,rp), real( 0.6994504559488d0,rp), real( 0.1530572479681d0,rp)]
      else if (flag_rk_ls_stages == 14) then

         a_i(1)   = real( 0.0000000000000000d0  ,rp) 
         a_i(2)   = real(-0.7188012108672410d0  ,rp) 
         a_i(3)   = real(-0.7785331173421570d0  ,rp) 
         a_i(4)   = real(-0.0053282796654044d0  ,rp) 
         a_i(5)   = real(-0.8552979934029281d0  ,rp) 
         a_i(6)   = real(-3.9564138245774565d0  ,rp) 
         a_i(7)   = real(-1.5780575380587385d0  ,rp) 
         a_i(8)   = real(-2.0837094552574054d0  ,rp) 
         a_i(9)   = real(-0.7483334182761610d0  ,rp) 
         a_i(10)  = real(-0.7032861106563359d0  ,rp) 
         a_i(11)  = real( 0.0013917096117681d0  ,rp) 
         a_i(12)  = real(-0.0932075369637460d0  ,rp) 
         a_i(13)  = real(-0.9514200470875948d0  ,rp) 
         a_i(14)  = real(-7.1151571693922548d0  ,rp) 

         b_i(1)   = real( 0.0367762454319673d0 ,rp) 
         b_i(2)   = real( 0.3136296607553959d0 ,rp) 
         b_i(3)   = real( 0.1531848691869027d0 ,rp) 
         b_i(4)   = real( 0.0030097086818182d0 ,rp) 
         b_i(5)   = real( 0.3326293790646110d0 ,rp) 
         b_i(6)   = real( 0.2440251405350864d0 ,rp) 
         b_i(7)   = real( 0.3718879239592277d0 ,rp) 
         b_i(8)   = real( 0.6204126221582444d0 ,rp) 
         b_i(9)   = real( 0.1524043173028741d0 ,rp) 
         b_i(10)  = real( 0.0760894927419266d0 ,rp) 
         b_i(11)  = real( 0.0077604214040978d0 ,rp) 
         b_i(12)  = real( 0.0024647284755382d0 ,rp) 
         b_i(13)  = real( 0.0780348340049386d0 ,rp) 
         b_i(14)  = real( 5.5059777270269628d0 ,rp) 
      else
         write(1,*) "--| NOT CODED FOR RK_LS 4 stages > 5 YET!"
         stop 1
      end if
      !$acc update device(a_i(:))
      !$acc update device(b_i(:))

      call nvtxEndRange

   end subroutine init_rk4_species_solver

   subroutine end_rk4_species_solver()
      implicit none

      !$acc exit data delete(auxReta(:))
      deallocate(auxReta)

      !$acc exit data delete(K2spc(:,:))
      !$acc exit data delete(Rspc(:,:))
      !$acc exit data delete(f_eta(:,:,:))
      !$acc exit data delete(Reta(:,:,:))
      deallocate(K2spc,f_eta,Reta,Rspc)

      !$acc exit data delete(a_i(:))
      !$acc exit data delete(b_i(:))
      deallocate(a_i,b_i)

   end subroutine end_rk4_species_solver

         subroutine rk_4_species_main(noBoundaries,isWallModelOn,nelem,nboun,npoin,npoin_w,numBoundsWM,point2elem,lnbn_nodes,dlxigp_ip,xgp,atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,&
                         ppow,connec,Ngp,dNgp,coord,wgp,He,Ml,gpvol,dt,helem,helem_l,Cp,Prt, &
                         rho,u,Yk,eta_Yk,mu_e_Yk,mu_sgs,lpoin_w,mu_fluid,mue_l, &
                         ndof,nbnodes,ldof,lbnodes,bound,bou_codes,bou_codes_nodes,&               ! Optional args
                         listBoundsWM,wgp_b,bounorm,normalsAtNodes,Yk_buffer)  ! Optional args

            implicit none

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
            integer(4)                          :: istep, ipoin, idime,icode,ispc
            real(rp)                            :: umag, rho_min, rho_avg

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! New version of RK4 using loops                 !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            pos = 2 ! Set correction as default value

            if(firstTimeStep .eqv. .true.) then
               firstTimeStep = .false.
               
               !$acc parallel loop collapse(2)
               do ipoin = 1,npoin_w
                  do ispc = 1,nspecies
                     f_eta(lpoin_w(ipoin),1,ispc) = u(lpoin_w(ipoin),1,1)*eta_Yk(lpoin_w(ipoin),ispc,1)
                     f_eta(lpoin_w(ipoin),2,ispc) = u(lpoin_w(ipoin),2,1)*eta_Yk(lpoin_w(ipoin),ispc,1)
                     f_eta(lpoin_w(ipoin),3,ispc) = u(lpoin_w(ipoin),3,1)*eta_Yk(lpoin_w(ipoin),ispc,1)
                  end do
               end do
               !$acc end parallel loop
               
               !$acc loop seq
               do ispc = 1,nspecies
                  call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                     gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:,ispc),eta_Yk(:,ispc,1),u(:,:,1),Reta(:,ispc,1))

                  if(mpi_size.ge.2) then
                     call mpi_halo_atomic_update_real(Reta(:,ispc,1))
                  end if
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,ispc,1))   

                  call nvtxStartRange("Entropy viscosity evaluation")
                  call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta(:,ispc,1),Ngp,coord,dNgp,gpvol,wgp, &
                     rho(:,1),u(:,:,1),eta_Yk(:,ispc,1),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
                  call nvtxEndRange
               end do
            end if

            !$acc parallel loop collapse(2)
            do ipoin = 1,npoin
               do ispc = 1,nspecies
                  K2spc(ipoin,ispc) = 0.0_rp
               end do
            end do
            !$acc end parallel loop
            !
            ! Loop over all RK steps
            !
            call nvtxStartRange("Loop over RK steps")
            do istep = 1,flag_rk_ls_stages
               do ispc = 1,nspecies
                  !
                  ! Apply bcs after update
                  !
                  if (noBoundaries .eqv. .false.) then
                     call nvtxStartRange("BCS_AFTER_UPDATE")
                     call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Yk(:,ispc,pos),Yk_buffer(:,ispc))
                     call nvtxEndRange
                  end if

                  ! Compute diffusion terms with values at current substep
                  !
                  call nvtxStartRange("CONVDIFFS")

                  call species_diffusion_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Prt,rho(:,pos),Yk(:,ispc,pos),mu_fluid,mu_e_Yk(:,:,ispc),mu_sgs,Ml,Rspc(:,ispc),.true.,-1.0_rp)
                  call species_convec_ijk(nelem,npoin,connec,Ngp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,rho(:,pos),Yk(:,ispc,pos),u(:,:,pos),Rspc(:,ispc),.false.,-1.0_rp)               
                  
                  call nvtxEndRange


                  !TESTING NEW LOCATION FOR MPICOMMS
                  if(mpi_size.ge.2) then
                     call nvtxStartRange("MPI_comms_tI")
                     call mpi_halo_atomic_update_real(Rspc(:,ispc))
                     call nvtxEndRange
                  end if

                  !
                  ! Call lumped mass matrix solver
                  !
                  call nvtxStartRange("Call solver")
                  call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Rspc(:,ispc))
                  call nvtxEndRange
                  !
                  ! Accumulate the residuals
                  !
                  call nvtxStartRange("Accumulate residuals")
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     K2spc(lpoin_w(ipoin),ispc) = a_i(istep)*K2spc(lpoin_w(ipoin),ispc) + Rspc(lpoin_w(ipoin),ispc)*dt
                     Yk(lpoin_w(ipoin),ispc,pos) = Yk(lpoin_w(ipoin),ispc,pos) + b_i(istep)*K2spc(lpoin_w(ipoin),ispc)
                  end do
                  !$acc end parallel loop
                  call nvtxEndRange
               end do
            end do
            call nvtxEndRange
            !
            ! Apply bcs after update
            !
            do ispc = 1,nspecies
               if (noBoundaries .eqv. .false.) then
                  call nvtxStartRange("BCS_AFTER_UPDATE")
                  call temporary_bc_routine_dirichlet_species(npoin,nboun,bou_codes,bou_codes_nodes,bound,nbnodes,lbnodes,lnbn_nodes,normalsAtNodes,Yk(:,ispc,pos),Yk_buffer(:,ispc))
                  call nvtxEndRange
               end if
            end do

            !
            ! Update velocity and equations of state
            !
            call nvtxStartRange("Update u and EOS")

            do ispc = 1,nspecies
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  eta_Yk(lpoin_w(ipoin),ispc,1) = eta_Yk(lpoin_w(ipoin),ispc,2)
                  eta_Yk(lpoin_w(ipoin),ispc,2) = 0.5_rp*Yk(lpoin_w(ipoin),ispc,pos)*Yk(lpoin_w(ipoin),ispc,pos)
                  !$acc loop seq
                  do idime = 1,ndime
                     f_eta(lpoin_w(ipoin),idime,ispc) = u(lpoin_w(ipoin),idime,1)*eta_Yk(lpoin_w(ipoin),ispc,1)
                  end do
               end do
               !$acc end parallel loop
               call nvtxEndRange

               call generic_scalar_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He, &
                  gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,f_eta(:,:,ispc),eta_Yk(:,ispc,1),u(:,:,1),Reta(:,ispc,2))

               if(mpi_size.ge.2) then
                  call mpi_halo_atomic_update_real(Reta(:,ispc,2))
               end if

               call lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,Reta(:,ispc,2))

               call nvtxStartRange("Entropy residual")
               !$acc parallel loop
               do ipoin = 1,npoin_w
                  auxReta(lpoin_w(ipoin)) = (1.5_rp*Reta(lpoin_w(ipoin),ispc,2)-0.5_rp*Reta(lpoin_w(ipoin),ispc,1)) + &
                                             (eta_Yk(lpoin_w(ipoin),ispc,2)-eta_Yk(lpoin_w(ipoin),ispc,1))/dt
                  Reta(lpoin_w(ipoin),ispc,1) = Reta(lpoin_w(ipoin),ispc,2)            
               end do
               !$acc end parallel loop
               call nvtxEndRange

               !
               ! Compute entropy viscosity
               !
               call nvtxStartRange("Entropy viscosity evaluation")
               call species_smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,auxReta,Ngp,coord,dNgp,gpvol,wgp, &
                  rho(:,2),u(:,:,2),eta_Yk(:,ispc,2),helem_l,helem,Ml,mu_e_Yk(:,:,ispc),invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,mue_l)
               call nvtxEndRange
            end do

         end subroutine rk_4_species_main

      end module time_integ_species
