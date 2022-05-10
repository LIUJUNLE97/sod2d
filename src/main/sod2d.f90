program sod2d

        !**************************************************************************!
        ! Computes a 2D or 3D compressible flow problem by the Galerkin FE method. !
        ! Stabilized through Entropy viscosity method.                             !
        !**************************************************************************!

        use mod_nvtx
#ifdef GPU
        use mod_gpu_vars
#endif
        use cudafor
        use mod_veclen

        use elem_qua
        use elem_hex
        use jacobian_oper
        use quadrature_rules
        use mesh_reader
        use inicond_reader
        use mass_matrix
        !use mod_graph ! Only useful if CSR model is being used
        use mod_geom
        use mod_output
        use mod_period
        use time_integ
        use mod_analysis
        use mod_constants
        use mod_time_ops
        use mod_fluid_viscosity

        implicit none

        integer(4)                 :: nstep, nper!, nzdom
        integer(4)                 :: idime, inode, igaus, istep, iper!, izdom
        integer(4)                 :: nelem, npoin, nboun, nbcodes
        integer(4)                 :: ielem, ipoin, iboun, ipbou
        integer(4)                 :: idof, ndof, nbnodes, ibnodes
        integer(4)                 :: ppow
        integer(4)                 :: flag_predic, flag_emac
        integer(4)                 :: nsave, nleap
        integer(4)                 :: nsave2, nleap2
        integer(4)                 :: counter
        integer(4)                 :: isPeriodic, npoin_w
        !integer(4), allocatable    :: rdom(:), cdom(:), aux_cdom(:) ! Use with CSR matrices
        integer(4), allocatable    :: connec(:,:), bound(:,:), ldof(:), lbnodes(:), bou_codes(:,:)
        integer(4), allocatable    :: masSla(:,:), connec_orig(:,:), aux1(:), bound_orig(:,:)
        integer(4), allocatable    :: lpoin_w(:), atoIJK(:), listHEX08(:,:), connecLINEAR(:,:)
        real(8),    allocatable    :: coord(:,:), coord_old(:,:), helem(:),helem_l(:)
        real(8),    allocatable    :: xgp(:,:), wgp(:)
        real(8),    allocatable    :: Ngp(:,:), dNgp(:,:,:)
        real(8),    allocatable    :: Ngp_l(:,:), dNgp_l(:,:,:)
        real(8),    allocatable    :: Je(:,:), He(:,:,:,:)
        real(8),    allocatable    :: gpvol(:,:,:)
        real(8),    allocatable    :: u(:,:,:), q(:,:,:), rho(:,:), pr(:,:), E(:,:), Tem(:,:), e_int(:,:), csound(:)
        real(8),    allocatable    :: Ml(:)!, Mc(:)
        real(8),    allocatable    :: mu_e(:,:), mu_fluid(:),mu_sgs(:,:)
        real(8),    allocatable    :: source_term(:)
        real(8),    allocatable    :: aux_1(:,:), aux_2(:)
        real(8)                    :: s, t, z, detJe
        real(8)                    :: dt, he_aux, time, P0, T0, EK, VolTot
        real(8)                    :: cfl_conv, cfl_diff
        character(500)             :: file_path
        character(500)             :: file_name, dumpfile
        character(4)               :: timeStep

        real(8)                    :: Rgas, gamma_gas, Cv
#ifdef CHANNEL
        !channel flow setup
        real(8)  :: vo = 1.0d0
        real(8)  :: M  = 0.1d0
        real(8)  :: delta  = 1.0d0
        real(8)  :: U0     = 1.0d0
        real(8)  :: rho0   = 1.0d0
        real(8)  :: Retau  = 950.0d0
        real(8)  ::  cp = 1004.0d0
        real(8)  ::  gamma_g = 1.40d0
        real(8)  :: yp=0.0d0, ti(3)
        real(8)  :: velo = 0.0d0, vol = 0.0d0
        real(8)  :: Re,mul,utau,Rg,to,po,mur
#else
        real(8)                    :: Cp, rho0, Re, mul,mur,to
#endif

        open(unit=1,file="sod2d.log",status="replace")
#ifdef CHANNEL
        Re     = exp((1.0d0/0.88d0)*log(Retau/0.09d0))
        mul    = (rho0*2.0d0*delta*vo)/Re
        utau   = (Retau*mul)/(delta*rho0)
        Rg = cp*(gamma_g-1.0d0)/gamma_g
        to = vo*vo/(gamma_g*Rg*M*M)
        po = rho0*Rg*to
        mur = 0.000001458d0*(to**1.50d0)/(to+110.40d0)
        flag_mu_factor = mul/mur

        write(1,*) " Gp ", utau*utau*rho0/delta
#else
        Re = 1600.0d0
        mul = 1.0d0*1.0d0*1.0d0/Re
        to = 1.0d0*1.0d0/(1.4d0*287.0d0*0.1d0*0.1d0)
        mur = 0.000001458d0*(to**1.50d0)/(to+110.40d0)
        flag_mu_factor = mul/mur

        !flag_mu_factor = 0.0d0 ! shock
#endif

        !*********************************************************************!
        ! Basic data, hardcoded for now                                       !
        !*********************************************************************!

        write(1,*) "--| ENTER PROBLEM DIMENSION (2 OR 3) :"
        !read(*,*) ndime
        !ndime = 3 ! Nsys
        !nnode = 27 ! TODO: need to allow for mixed elements...
        !porder = 1 ! TODO: make it input
        !npbou = 9 ! TODO: Need to get his from somewhere...
        nstep = 900000 ! TODO: Needs to be input...
#ifdef CHANNEL
        Rgas = Rg
#else
        Rgas = 287.00d0 ! TODO: Make it input
#endif
        Cp = 1004.00d0 ! TODO: Make it input
        gamma_gas = 1.40d0 ! TODO: Make it innput
        Cv = Cp/gamma_gas
        cfl_conv = 0.85d0
        cfl_diff = 2.0d0
        nsave  = 1   ! First step to save, TODO: input
        nsave2 = 1   ! First step to save, TODO: input
        nleap = 800 ! Saving interval, TODO: input
        nleap2 = 10  ! Saving interval, TODO: input
#ifdef CHANNEL
        isPeriodic = 1 ! TODO: make it a read parameter (0 if not periodic, 1 if periodic)
#else
        isPeriodic = 1 ! TODO: make it a read parameter (0 if not periodic, 1 if periodic)
#endif        
        if (isPeriodic == 1) then
#ifdef CHANNEL
           nper = 6499 ! TODO: if periodic, request number of periodic nodes
           !nper = 2145 ! TODO: if periodic, request number of periodic nodes
#else
           !nper = 1387 ! TODO: if periodic, request number of periodic nodes
           nper = 97741  ! TODO: if periodic, request number of periodic nodes
#endif
        else if (isPeriodic == 0) then
           nper = 0 ! Set periodic nodes to zero if case is not periodic
        end if
        flag_emac = 1
        if (flag_emac == 1) then
           write(1,*) "--| RUNNING WITH EMAC CONVECTION"
        else if (flag_emac == 0) then
           write(1,*) "--| RUNNING WITH CONSERV CONVECTION"
        else
           write(1,*) "--| FLAG_EMAC MUST BE EITHER 0 OR 1!"
           STOP(1)
        end if

#ifdef CHANNEL
        allocate(source_term(ndime))
        !set the source term
        source_term(1) = (utau*utau*rho0/delta)/36.0d0
        source_term(2) = 0.00d0
        source_term(3) = 0.00d0
#endif

        !*********************************************************************!
        ! Define vector length to be used                                     !
        !*********************************************************************!
        
        call define_veclen()

        !*********************************************************************!
        ! Read mesh in Alya format                                            !
        !*********************************************************************!

        write(file_path,*) "./mesh/"
        write(1,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(file_path))," !"
        write(1,*) "--| ENTER NAME OF MESH RELATED FILES :"
        call nvtxStartRange("Read mesh")
        !read(*,*) file_name
#ifdef CHANNEL
        write(file_name,*) "channel" ! Nsys
#else
        !write(file_name,*) "shock_tube" ! Nsys
        write(file_name,*) "cube" ! Nsys
#endif
        call read_dims(file_path,file_name,npoin,nelem,nboun)
        allocate(connec(nelem,nnode))
        if (nboun .ne. 0) then
           allocate(bound(nboun,npbou))
           allocate(bou_codes(nboun,2))
           call read_fixbou(file_path,file_name,nboun,nbcodes,bou_codes)
        end if
        allocate(coord(npoin,ndime))
        call read_geo_dat(file_path,file_name,npoin,nelem,nboun,connec,bound,coord)
        if (isPeriodic == 1) then
           allocate(masSla(nper,2))
           allocate(connec_orig(nelem,nnode))
           if (nboun .ne. 0) then
              allocate(bound_orig(nboun,npbou))
           end if
           call read_periodic(file_path,file_name,nper,masSla)
        end if
        call nvtxEndRange

        !*********************************************************************!
        ! Compute characteristic size of elements                             !
        !*********************************************************************!

        call nvtxStartRange("Elem size compute")
        allocate(helem(nelem))
        do ielem = 1,nelem
           call char_length(ielem,nelem,npoin,connec,coord,he_aux)
           helem(ielem) = he_aux
        end do
        call nvtxEndRange

        !*********************************************************************!
        ! Create mesh graph for CSR matrices                                  !
        !*********************************************************************!

        !!write(*,*) "--| PERFORMING GRAPH OPERATIONS..."
        !!allocate(rdom(npoin+1))                                          ! Implicit row indexing
        !!allocate(aux_cdom(nelem*nnode*nnode))                            ! Preliminary cdom for subroutine
        !!call compute_nzdom(npoin,nnode,nelem,connec,nzdom,rdom,aux_cdom) ! Computes nzdom, rdom and aux_cdom
        !!allocate(cdom(nzdom))                                            ! Row indexes with proper array size
        !!do izdom = 1,nzdom
        !!   cdom(izdom) = aux_cdom(izdom)
        !!end do
        !!deallocate(aux_cdom)
        !!write(*,*) "--| END OF GRAPH OPERATIONS!"

        !*********************************************************************!
        ! Generate list of "free" nodes                                       !
        !*********************************************************************!

        if (nboun .ne. 0) then
            write(1,*) "--| SPLITTING BOUNDARY NODES FROM DOFs..."
            call nvtxStartRange("Bnodes split")
            allocate(aux1(npoin))

            !
            ! Fill aux1 with all nodes in order
            !
            !$acc parallel loop
            do ipoin = 1,npoin
               aux1(ipoin) = ipoin
            end do
            !$acc end parallel loop

            !
            ! If node is on boundary, zero corresponding aux1 entry
            !
            !$acc parallel loop gang 
            do iboun = 1,nboun
               !$acc loop vector
               do ipbou = 1,npbou
                  aux1(bound(iboun,ipbou)) = 0
               end do
            end do
            !$acc end parallel loop

            !
            ! Determine how many nodes are boundary nodes
            !
            ndof = 0
            do ipoin = 1,npoin
               if (aux1(ipoin) == 0) then
                  ndof = ndof+1
               end if
            end do

            nbnodes = ndof    ! Nodes on boundaries
            ndof = npoin-ndof ! Free nodes
            write(1,*) '--| TOTAL FREE NODES := ',ndof
            write(1,*) '--| TOTAL BOUNDARY NODES := ',nbnodes

            allocate(ldof(ndof))
            allocate(lbnodes(nbnodes))

            !
            ! Split aux1 into the 2 lists
            !
            idof = 0    ! Counter for free nodes
            ibnodes = 0 ! Counter for boundary nodes
            !$acc parallel loop reduction(+:idof,ibnodes)
            do ipoin = 1,npoin
               if (aux1(ipoin) == 0) then
                  ibnodes = ibnodes+1
                  lbnodes(ibnodes) = ipoin
               else
                  idof = idof+1
                  ldof(idof) = aux1(ipoin)
               end if
            end do
            !$acc end parallel loop
            call nvtxEndRange
        end if

        !*********************************************************************!
        ! Allocate variables                                                  !
        !*********************************************************************!

        WRITE(1,*) "--| ALLOCATING MAIN VARIABLES"
        call nvtxStartRange("Allocate main vars")
        !
        ! Last rank is for prediction-advance related to entropy viscosity,
        ! where 1 is prediction, 2 is final value
        !
        allocate(u(npoin,ndime,2))  ! Velocity
        allocate(q(npoin,ndime,2))  ! momentum
        allocate(rho(npoin,2))      ! Density
        allocate(pr(npoin,2))       ! Pressure
        allocate(E(npoin,2))        ! Total Energy
        allocate(Tem(npoin,2))      ! Temperature
        allocate(e_int(npoin,2))    ! Internal Energy
        allocate(csound(npoin))     ! Speed of sound
        allocate(mu_fluid(npoin))   ! Fluid viscosity
        allocate(mu_e(nelem,ngaus))  ! Elemental viscosity
        allocate(mu_sgs(nelem,ngaus))! SGS viscosity
        call nvtxEndRange

        !*********************************************************************!
        ! Read initial conditions                                             !
        !*********************************************************************!

#ifndef CHANNEL
        call nvtxStartRange("Read ICs")
        call read_veloc(npoin,file_path,u(:,:,2))
        call read_densi(npoin,file_path,rho(:,2))
        call read_press(npoin,file_path,pr(:,2)) ! Can be switched for TEMPE
        call nvtxEndRange
#endif

        !*********************************************************************!
        ! Generate complementary info                                         !
        !*********************************************************************!

        ! Assuming u, rho, p as IC:

        call nvtxStartRange("Additional data")
#ifdef CHANNEL
        if(1) then
           do ipoin = 1,npoin

              if(coord(ipoin,2)<delta) then
                 yp = coord(ipoin,2)*utau/mul
              else
                 yp = abs(coord(ipoin,2)-2.0d0*delta)*utau/mul
              end if

              velo = utau*((1.0d0/0.41d0)*log(1.0d0+0.41d0*yp)+7.8d0*(1.0d0-exp(-yp/11.0d0)-(yp/11.0d0)*exp(-yp/3.0d0))) 

              call random_number(ti)
              u(ipoin,1,2) = velo*(1.0d0 + 0.1d0*(ti(1) -0.5d0))
              u(ipoin,2,2) = velo*(0.1d0*(ti(2) -0.5d0))
              u(ipoin,3,2) = velo*(0.1d0*(ti(3) -0.5d0))

              pr(ipoin,2) = po

              rho(ipoin,2) = po/Rg/to

              e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(gamma_gas-1.0d0))
              Tem(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*Rgas)
              E(ipoin,2) = rho(ipoin,2)*(0.5d0*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
              q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
              csound(ipoin) = sqrt(gamma_gas*pr(ipoin,2)/rho(ipoin,2))
           end do
        else        
           call read_vtk_binary(isPeriodic,230001,npoin,nelem,coord,connec, &
              rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
           do ipoin = 1,npoin
              pr(ipoin,2) = po
              rho(ipoin,2) = po/Rg/to
              e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(gamma_gas-1.0d0))
              Tem(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*Rgas)
              E(ipoin,2) = rho(ipoin,2)*(0.5d0*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
              q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
              csound(ipoin) = sqrt(gamma_gas*pr(ipoin,2)/rho(ipoin,2))
           end do
        endif        
#else
        !$acc parallel loop
        do ipoin = 1,npoin
           e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(gamma_gas-1.0d0))
           Tem(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*Rgas)
           E(ipoin,2) = rho(ipoin,2)*(0.5d0*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
           q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
           csound(ipoin) = sqrt(gamma_gas*pr(ipoin,2)/rho(ipoin,2)) 
        end do
        !$acc end parallel loop
#endif
        !$acc kernels
        mu_e(:,:) = 0.0d0 ! Element syabilization viscosity
        mu_sgs(:,:) = 0.0d0
        !$acc end kernels
        call nvtxEndRange

        if (flag_real_diff == 1) then
           if (flag_diff_suth == 0) then
              call constant_viscosity(npoin,0.000055d0,mu_fluid)
           else
              call sutherland_viscosity(npoin,Tem(:,2),mu_fluid)
           end if
        else if (flag_real_diff == 0) then
           !$acc kernels
           mu_fluid(:) = 0.0d0
           !$acc end kernels
        else
           write(1,*) "--| DIFFUSION FLAG MUST BE EITHER 0 OR 1, NOT: ",flag_real_diff
           STOP(1)
        end if

        !*********************************************************************!
        ! Compute initial time-step size                                      !
        !*********************************************************************!
        
        if (flag_real_diff == 1) then
           call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid,rho(:,2))
           write(1,*) "--| TIME STEP SIZE dt := ",dt,"s"
        else
           call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
           write(1,*) "--| TIME STEP SIZE dt := ",dt,"s"
        end if

        if (((porder+1)**ndime) .le. (3**ndime) .and. flag_spectralElem == 0) then
           !
           ! Call VTK output (0th step)
           !
           write(1,*) "--| GENERATING 1st OUTPUT..."
           call nvtxStartRange("1st write")
           if (isPeriodic == 0) then
              call write_vtk_binary(isPeriodic,0,npoin,nelem,coord,connec, &
                                   rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper)
           else
              print*, 'sub call: ok!'
              call write_vtk_binary(isPeriodic,0,npoin,nelem,coord,connec, &
                                   rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
           end if
        end if
        call nvtxEndRange

        !*********************************************************************!
        ! Generate GLL table                                                  !
        !*********************************************************************!

        write(1,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."

        call nvtxStartRange("Gaussian Quadrature")

        allocate(xgp(ngaus,ndime))
        allocate(wgp(ngaus))

        if (ndime == 2) then
           if (nnode == (porder+1)**2) then ! QUA_XX of order porder
              call gll_qua(xgp,wgp)
           else if (nnode == 3 .or. nnode == 6 .or. nnode == 10) then ! TRI_XX
              write(1,*) '--| NOT CODED YET!'
              STOP 1
           end if
        else if (ndime == 3) then
           if (nnode == (porder+1)**3) then ! HEX_XX
              if (flag_SpectralElem == 0) then
                 write(1,*) "  --| GENERATING GLL TABLE..."
                 call gll_hex(xgp,wgp)
              else if (flag_SpectralElem == 1) then
                 write(1,*) "  --| GENERATING CHEBYSHEV TABLE..."
                 allocate(atoIJK(64))
                 call hex64(1.0d0,1.0d0,1.0d0,atoIJK)
                 call chebyshev_hex(atoIJK,xgp,wgp)
              end if
           else if (nnode == 4 .or. nnode == 10 .or. nnode == 20) then ! TET_XX
              write(1,*) '--| NOT CODED YET!'
              STOP 1
           end if
        end if
        call nvtxEndRange

        !*********************************************************************!
        ! Generate N and dN for all GP                                        !
        !*********************************************************************!

        ! TODO: Allow for more element types

        write(1,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."
        call nvtxStartRange("N and dN")

        allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))
        if (flag_SpectralElem == 1) then
           allocate(Ngp_l(ngaus,nnode),dNgp_l(ndime,nnode,ngaus))
        end if

        do igaus = 1,ngaus
           s = xgp(igaus,1)
           t = xgp(igaus,2)
           if (ndime == 2) then
              if (nnode == 4) then
                 call qua04(s,t,Ngp(igaus,:),dNgp(:,:,igaus))
              else if (nnode == 9) then
                 !call qua09(s,t,Ngp(igaus,:),dNgp(:,:,igaus))
                 write(1,*) '--| NOT CODED YET!'
                 STOP 1
              end if
           else if (ndime == 3) then
              z = xgp(igaus,3)
              if (nnode == 8) then
                 call hex08(s,t,z,Ngp(igaus,:),dNgp(:,:,igaus))
              else if (nnode == 27) then
                 call hex27(s,t,z,Ngp(igaus,:),dNgp(:,:,igaus))
              else if (nnode == 64) then
                 if (flag_spectralElem == 0) then
                    allocate(listHEX08((porder**ndime),2**ndime))
                    call hex64(s,t,z,atoIJK,listHEX08,Ngp(igaus,:),dNgp(:,:,igaus))
                 else if (flag_spectralElem == 1) then
                    allocate(listHEX08((porder**ndime),2**ndime))
                    call hex64(s,t,z,atoIJK,listHEX08,Ngp(igaus,:),dNgp(:,:,igaus),Ngp_l(igaus,:),dNgp_l(:,:,igaus))
                 end if
              else
                 write(1,*) '--| NOT CODED YET!'
                 STOP(1)
              end if
           end if
        end do
        call nvtxEndRange

        !*********************************************************************!
        ! Adjust element nodes and variables if spectral 
        ! element type being used                                             !
        !*********************************************************************!

        if (flag_spectralElem == 1) then
           allocate(aux_1(npoin,ndime))
           aux_1(:,:) = coord(:,:)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 do idime = 1,ndime
                    call var_interpolate(aux_1(connec(ielem,:),idime),Ngp_l(inode,:),coord(connec(ielem,inode),idime))
                 end do
              end do
           end do
           aux_1(:,:) = u(:,:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 do idime = 1,ndime
                    call var_interpolate(aux_1(connec(ielem,:),idime),Ngp_l(inode,:),u(connec(ielem,inode),idime,2))
                 end do
              end do
           end do
           aux_1(:,:) = q(:,:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 do idime = 1,ndime
                    call var_interpolate(aux_1(connec(ielem,:),idime),Ngp_l(inode,:),q(connec(ielem,inode),idime,2))
                 end do
              end do
           end do
           deallocate(aux_1)
           allocate(aux_2(npoin))
           aux_2(:) = rho(:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),rho(connec(ielem,inode),2))
              end do
           end do
           aux_2(:) = pr(:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),pr(connec(ielem,inode),2))
              end do
           end do
           aux_2(:) = E(:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),E(connec(ielem,inode),2))
              end do
           end do
           aux_2(:) = Tem(:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),Tem(connec(ielem,inode),2))
              end do
           end do
           aux_2(:) = e_int(:,2)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),e_int(connec(ielem,inode),2))
              end do
           end do
           aux_2(:) = csound(:)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),csound(connec(ielem,inode)))
              end do
           end do
           aux_2(:) = mu_fluid(:)
           do ielem = 1,nelem
              do inode = (2**ndime)+1,nnode
                 call var_interpolate(aux_2(connec(ielem,:)),Ngp_l(inode,:),mu_fluid(connec(ielem,inode)))
              end do
           end do
           deallocate(aux_2)
        end if
        !charecteristic length for spectral elements for the entropy
        !stablisation
        allocate(helem_l(npoin))
        if (flag_SpectralElem == 1) then
           helem_l(:) = 1000000000000000000000000000000000000000000000000000000000000000000000.0d0
           do ielem = 1,nelem
              call char_length_spectral(ielem,nelem,npoin,connec,coord,helem_l)
           end do
        end if

        !*********************************************************************!
        ! Generate linear mesh and output for spectral case                   !
        !*********************************************************************!
        if (flag_spectralElem == 1) then
           allocate(connecLINEAR(nelem*(porder**ndime),2**ndime))
           call linearMeshOutput(nelem,connec,listHEX08,connecLINEAR)
           !
           ! Call VTK output (0th step)
           !
           write(1,*) "--| GENERATING 1st OUTPUT..."
           call nvtxStartRange("1st write")
           if (isPeriodic == 0) then
              call write_vtk_binary_linearized(isPeriodic,0,npoin,nelem,coord,connecLINEAR,connec, &
                                   rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper)
           else
              call write_vtk_binary_linearized(isPeriodic,0,npoin,nelem,coord,connecLINEAR,connec, &
                                   rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
           end if
        end if

        !*********************************************************************!
        ! Generate Jacobian related information                               !
        !*********************************************************************!

        write(1,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

        call nvtxStartRange("Jacobian info")
        allocate(He(ndime,ndime,ngaus,nelem))
        allocate(gpvol(1,ngaus,nelem))
        call elem_jacobian(nelem,npoin,connec,coord,dNgp,wgp,gpvol,He)
        call  nvtxEndRange
        VolTot = 0.0d0
        do ielem = 1,nelem
           do igaus = 1,ngaus
              VolTot = VolTot+gpvol(1,igaus,ielem)
           end do
        end do
        write(1,*) '--| DOMAIN VOLUME := ',VolTot
        !STOP(1)

        !*********************************************************************!
        ! Treat periodicity                                                   !
        !*********************************************************************!
        
        ! TODO: Verify if lpoin_W is really necessary. Complexity may not be worth a few
        ! less ops.

        if (isPeriodic == 1) then
           if (nboun .eq. 0) then
              call periodic_ops(nelem,npoin,nboun,npoin_w,nper, &
                                lpoin_w,connec,connec_orig,masSla)
           else
              call periodic_ops(nelem,npoin,nboun,npoin_w,nper, &
                                lpoin_w,connec,connec_orig,masSla,bound,bound_orig)
           end if
        else if (isPeriodic == 0) then
           npoin_w = npoin
           allocate(lpoin_w(npoin_w)) ! All nodes are working nodes
           !$acc parallel loop
           do ipoin = 1,npoin_w
              lpoin_w(ipoin) = ipoin
           end do
           !$acc end parallel loop
        end if

        !*********************************************************************!
        ! Compute mass matrix (Lumped and Consistent) and set solver type     !
        !*********************************************************************!

        write(1,*) '--| COMPUTING LUMPED MASS MATRIX...'
        call nvtxStartRange("Lumped mass compute")
        allocate(Ml(npoin))
        if (flag_spectralElem == 0) then
           call lumped_mass(nelem,npoin,connec,gpvol,Ngp,Ml)
        else if (flag_spectralElem == 1) then
           call lumped_mass_spectral(nelem,npoin,connec,gpvol,Ml)
        else
           write(1,*) '--| SPECTRAL ELEMENT FLAG MUST BE 1 OR 0!'
        end if
        call nvtxEndRange

        !
        ! Consisten mass: activate with CSR related operations
        !
        !write(*,*) '--| COMPUTING CONSISTENT MASS MATRIX...'
        !allocate(Mc(nzdom))
        !call consistent_mass(nelem,nnode,npoin,ngaus,connec,nzdom,rdom,cdom,gpvol,Ngp,Mc)

        if (flag_solver_type == 2 .and. flag_spectralElem == 0) then
            write(1,*) '--| ENTER NUMBER OF ITERATIONS FOR APINV SOLVER:'
            !read(*,*) ppow
            ppow = 2 ! Nsys
        end if
        write(1,*) '--| USING SOLVER ',flag_solver_type,' FOR MASS MATRIX'

        !*********************************************************************!
        ! Create variables on GPU                                             !
        !*********************************************************************!

#ifdef GPU

        !! Range with standard color
        !call nvtxStartRange("Memory Management")

        !! Mesh info

        !allocate(connec_d(nelem,nnode))
        !allocate(lbnodes_d(nbnodes))

        !connec_d = connec
        !lbnodes_d = lbnodes

        !! Primary vars.

        !allocate(rho_d(npoin,2))      ! Density
        !allocate(u_d(npoin,ndime,2))  ! Velocity
        !allocate(q_d(npoin,ndime,2))  ! momentum
        !allocate(pr_d(npoin,2))       ! Pressure
        !allocate(E_d(npoin,2))        ! Total Energy
        !allocate(Tem_d(npoin,2))      ! Temperature
        !allocate(e_int_d(npoin,2))    ! Internal Energy
        !allocate(mu_e_d(nelem))       ! Elemental viscosity

        !rho_d = rho
        !u_d = u
        !q_d = q
        !pr_d = pr
        !E_d = E
        !Tem_d = Tem
        !e_int_d = e_int

        !! Mass matrices

        !allocate(Ml_d(npoin))
        !!allocate(Mc_d(nzdom))

        !Ml_d = Ml
        !!Mc_d = Mc

        !! Elemental info

        !allocate(Ngp_d(ngaus,npoin))
        !allocate(gpvol_d(1,ngaus,nelem))
        !allocate(gpcar_d(ndime,nnode,ngaus,nelem))

        !Ngp_d = Ngp
        !gpvol_d = gpvol
        !gpcar_d = gpcar

        !! End nvtx range
        !call nvtxEndRange

#endif

        !*********************************************************************!
        ! Start of time stepping                                              !
        !*********************************************************************!

        !
        ! Write EK to file
        !
        open(unit=666,file="analysis.dat",status="replace")
        time = 0.0d0
        !P0 = 1.0d0
        !T0 = 1.0d0/(1.4d0*Rgas*(0.1**2))
        rho0 = 1.0d0
        call volAvg_EK(nelem,npoin,connec,gpvol,Ngp,rho0,rho(:,2),u(:,:,2),EK)
        call write_EK(time,EK)
        write(1,*) "--| time   ,   EK"
        write(1,*) "--| ",time,"  |  ",EK

        counter = 1

        call nvtxStartRange("Start RK4")
        if (isPeriodic .eq. 0) then ! Case is not periodic
           if (nboun .eq. 0) then ! Case has no boundaries
              write(1,*) '--| ERROR: CASE MUST HAVE BOUNDARIES!'
              STOP(1)
           else ! Case has boundaries
              !TODO: Add call for write_vtk_linearized when mesh is spectral
              write(1,*) '--| NON-PERIODIC CASE WITH BOUNDARIES'
              do istep = 1,nstep

                 write(1,*) '   --| STEP: ', istep

                 !
                 ! Prediction
                 !
                 flag_predic = 1
                 call nvtxStartRange("Init pred "//timeStep,istep)
                 !$acc kernels
                 rho(:,1) = rho(:,2)
                 u(:,:,1) = u(:,:,2)
                 q(:,:,1) = q(:,:,2)
                 pr(:,1) = pr(:,2)
                 E(:,1) = E(:,2)
                 Tem(:,1) = Tem(:,2)
                 e_int(:,1) = e_int(:,2)
                 !$acc end kernels
                 call nvtxEndRange

                 ! nvtx range for full RK
                 write(timeStep,'(i4)') istep
                 call nvtxStartRange("RK step "//timeStep,istep)

#ifndef NOPRED
                 if(flag_rk_order .eq. 3) then
                    ! TODO: Remove rk_3_main
                    call rk_3_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                       ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                       rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                       ndof,nbnodes,ldof,lbnodes,bound,bou_codes) ! Optional args
                 else
                    ! TODO: Rename this to something more descriptive
                    call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                       ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                       rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                       ndof,nbnodes,ldof,lbnodes,bound,bou_codes) ! Optional args
                 end if
#endif
                 !
                 ! Advance with entropy viscosity
                 !
                 flag_predic = 0
                 if(flag_rk_order .eq. 3) then
                    call rk_3_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                       ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                       rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                       ndof,nbnodes,ldof,lbnodes,bound,bou_codes) ! Optional args
                 else
                    call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                       ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                       rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                       ndof,nbnodes,ldof,lbnodes,bound,bou_codes) ! Optional args
                 end if

                 if (flag_real_diff == 1) then
                    call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid,rho(:,2))
                    write(1,*) "DT := ",dt,"s"
                 else
                    call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
                    write(1,*) "DT := ",dt,"s"
                 end if

                 call nvtxEndRange

                 !
                 ! Call VTK output
                 !
                 if (istep == nsave) then
                    call nvtxStartRange("Output "//timeStep,istep)
                    call write_vtk_binary(isPeriodic,counter,npoin,nelem,coord,connec, &
                                         rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper)
                    nsave = nsave+nleap
                    call nvtxEndRange
                 end if

                 counter = counter+1

              end do
           end if
        else if (isPeriodic .eq. 1) then ! Case is periodic
           if (nboun .eq. 0) then ! Case has no boundaries
               write(1,*) '--| PERIODIC CASE WITH NO BOUNDARIES'
               do istep = 1,nstep

                  if (istep == nsave) write(1,*) '   --| STEP: ', istep

                  !
                  ! Prediction
                  !
                  flag_predic = 1
                  call nvtxStartRange("Init pred "//timeStep,istep)
                  !$acc kernels
                  rho(:,1) = rho(:,2)
                  u(:,:,1) = u(:,:,2)
                  q(:,:,1) = q(:,:,2)
                  pr(:,1) = pr(:,2)
                  E(:,1) = E(:,2)
                  Tem(:,1) = Tem(:,2)
                  e_int(:,1) = e_int(:,2)
                  !$acc end kernels
                  call nvtxEndRange

                  ! nvtx range for full RK
                  !write(timeStep,'(i4)') istep
                  call nvtxStartRange("RK4 step "//timeStep,istep)
#ifndef NOPRED
                 if(flag_rk_order .eq. 3) then
                    call rk_3_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                       ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                       rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid)
                 else
                    call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                       ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                       rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid)
                 end if
#endif
                  !
                  ! Advance with entropy viscosity
                  !
                  flag_predic = 0
                  if(flag_rk_order .eq. 3) then
                     call rk_3_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                        ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                        rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid)
                  else
                     call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                        ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                        rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid)
                  end if

                  time = time+dt

                  if (istep == nsave2) then
                     call volAvg_EK(nelem,npoin,connec,gpvol,Ngp,rho0,rho(:,2),u(:,:,2),EK)
                     call write_EK(time,EK)
                     write(1,*) "--| time   ,   EK"
                     write(1,*) "--| ",time,"  |  ",EK
                  end if

                  if (flag_real_diff == 1) then
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid,rho(:,2))
                     if (istep == nsave2) write(1,*) "DT := ",dt,"s"
                  else
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
                     if (istep == nsave2) write(1,*) "DT := ",dt,"s"
                  end if

                  call nvtxEndRange

                  !
                  ! Call VTK output
                  !
                  if (istep == nsave) then
                     call nvtxStartRange("Output "//timeStep,istep)
                     if (flag_spectralElem == 1) then
                        call write_vtk_binary_linearized(isPeriodic,counter,npoin,nelem,coord,connecLINEAR,connec, &
                           rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
                     else 
                        call write_vtk_binary(isPeriodic,counter,npoin,nelem,coord,connec_orig, &
                           rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
                     end if
                     nsave = nsave+nleap
                     call nvtxEndRange
                  end if

                  if(istep==nsave2) then
                     nsave2 = nsave2+nleap2
                     call flush(666)
                     call flush(1)
                  end if

                  counter = counter+1

               end do
               close(666)
            else
               write(1,*) '--| PERIODIC CASE WITH BOUNDARIES'
               do istep = 1,nstep

                  if (istep == nsave)write(1,*) '   --| STEP: ', istep

                  !
                  ! Prediction
                  !
                  flag_predic = 1
                  call nvtxStartRange("Init pred "//timeStep,istep)
                  !$acc kernels
                  rho(:,1) = rho(:,2)
                  u(:,:,1) = u(:,:,2)
                  q(:,:,1) = q(:,:,2)
                  pr(:,1) = pr(:,2)
                  E(:,1) = E(:,2)
                  Tem(:,1) = Tem(:,2)
                  e_int(:,1) = e_int(:,2)
                  !$acc end kernels
                  call nvtxEndRange

                  ! nvtx range for full RK
                 ! write(timeStep,'(i4)') istep
                  call nvtxStartRange("RK4 step "//timeStep,istep)
#ifndef NOPRED
                if(flag_rk_order .eq. 3) then
                   call rk_3_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                      ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                      rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                      ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args
                else
                   call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                      ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                      rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                      ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args
                end if
#endif
                  !
                  ! Advance with entropy viscosity
                  !
                  flag_predic = 0
                  if(flag_rk_order .eq. 3) then
                     call rk_3_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                        ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                        rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                        ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args
                  else
                     call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                        ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,helem_l,Rgas,gamma_gas, &
                        rho,u,q,pr,E,Tem,e_int,mu_e,mu_sgs,lpoin_w,mu_fluid, &
                        ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args
                  end if

                  time = time+dt
                  if (flag_real_diff == 1) then
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid,rho(:,2))
                     if (istep == nsave2)write(1,*) "DT := ",dt,"s time := ",time,"s"
                  else
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
                     if (istep == nsave2)write(1,*) "DT := ",dt,"s time := ",time,"s"
                  end if

                  call nvtxEndRange

                  !
                  ! Call VTK output
                  !
                  if (istep == nsave) then
                     call nvtxStartRange("Output "//timeStep,istep)
                     if (flag_spectralElem == 1) then
                        call write_vtk_binary_linearized(isPeriodic,counter,npoin,nelem,coord,connecLINEAR,connec, &
                           rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
                     else 
                        call write_vtk_binary(isPeriodic,counter,npoin,nelem,coord,connec_orig, &
                           rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_fluid,mu_e,mu_sgs,nper,masSla)
                     end if
                     nsave = nsave+nleap
                     call nvtxEndRange
                  end if

                  if(istep==nsave2) then
                     nsave2 = nsave2+nleap2
                     call flush(1)
                  end if

                  counter = counter+1

               end do
            end if
        end if
        call nvtxEndRange

end program sod2d
