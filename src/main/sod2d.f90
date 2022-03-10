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
        integer(4)                 :: counter
        integer(4)                 :: isPeriodic, npoin_w
        !integer(4), allocatable    :: rdom(:), cdom(:), aux_cdom(:) ! Use with CSR matrices
        integer(4), allocatable    :: connec(:,:), bound(:,:), ldof(:), lbnodes(:), bou_codes(:,:)
        integer(4), allocatable    :: masSla(:,:), connec_orig(:,:), aux1(:), bound_orig(:,:)
        integer(4), allocatable    :: lpoin_w(:)
        real(8),    allocatable    :: coord(:,:), helem(:)
        real(8),    allocatable    :: xgp(:,:), wgp(:)
        real(8),    allocatable    :: Ngp(:,:), dNgp(:,:,:)
        real(8),    allocatable    :: Je(:,:), He(:,:,:,:)
        real(8),    allocatable    :: gpvol(:,:,:)
        real(8),    allocatable    :: u(:,:,:), q(:,:,:), rho(:,:), pr(:,:), E(:,:), Tem(:,:), e_int(:,:), csound(:)
        real(8),    allocatable    :: Ml(:)!, Mc(:)
        real(8),    allocatable    :: mu_e(:), mu_fluid(:)
        real(8),    allocatable    :: source_term(:)
        real(8)                    :: s, t, z, detJe
        real(8)                    :: dt, he_aux, time, P0, T0, EK, VolTot
        real(8)                    :: cfl_conv, cfl_diff
        character(500)             :: file_path
        character(500)             :: file_name, dumpfile
        character(5)               :: matrix_type, solver_type
        character(4)               :: timeStep

        real(8)                    :: Rgas, gamma_gas, Cv
#ifdef CHANNEL
        !channel flow setup
        real(8)  :: vo = 27.263d0
        real(8)  :: M  = 0.2d0
        real(8)  :: delta  = 1.0d0
        real(8)  :: U0     = 1.0d0
        real(8)  :: rho0   = 1.0d0
        real(8)  :: Retau  = 180.0d0
        real(8)  ::  cp = 1004.0d0
        real(8)  ::  gamma_g = 1.40d0
        real(8)  :: yp=0.0d0, ti(3)
        real(8)  :: velo = 0.0d0, vol = 0.0d0
        real(8)  :: Re,mul,utau,Rg,to,po
#else
        real(8)                    :: Cp, rho0
#endif

#ifdef CHANNEL
        Re     = exp((1.0d0/0.88d0)*log(Retau/0.09d0))
        mul    = (rho0*2.0d0*delta*vo)/Re
        utau   = (Retau*mul)/(delta*rho0)
        Rg = cp*(gamma_g-1.0d0)/gamma_g
        to = vo*vo/(gamma_g*Rg*M*M)
        po = rho0*Rg*to

        write(*,*) " Gp ", utau*utau*rho0/delta
#endif

        !*********************************************************************!
        ! Basic data, hardcoded for now                                       !
        !*********************************************************************!

        write(*,*) "--| ENTER PROBLEM DIMENSION (2 OR 3) :"
        !read(*,*) ndime
        !ndime = 3 ! Nsys
        !nnode = 27 ! TODO: need to allow for mixed elements...
        !porder = 2 ! TODO: make it input
        !npbou = 9 ! TODO: Need to get his from somewhere...
        nstep = 10 ! TODO: Needs to be input...
#ifdef CHANNEL
        Rgas = Rg
#else
        Rgas = 287.00d0 ! TODO: Make it input
#endif
        Cp = 1004.00d0 ! TODO: Make it input
        gamma_gas = 1.40d0 ! TODO: Make it innput
        Cv = Cp/gamma_gas
        cfl_conv = 0.80d0
        cfl_diff = 0.20d0
        nsave = 1 ! First step to save, TODO: input
        nleap = 1 ! Saving interval, TODO: input
        isPeriodic = 1 ! TODO: make it a read parameter (0 if not periodic, 1 if periodic)
        if (isPeriodic == 1) then
           nper = 49537 ! TODO: if periodic, request number of periodic nodes
        else if (isPeriodic == 0) then
           nper = 0 ! Set periodic nodes to zero if case is not periodic
        end if
        flag_emac = 0
        if (flag_emac == 1) then
           write(*,*) "--| RUNNING WITH EMAC CONVECTION"
        else if (flag_emac == 0) then
           write(*,*) "--| RUNNING WITH CONSERV CONVECTION"
        else
           write(*,*) "--| FLAG_EMAC MUST BE EITHER 0 OR 1!"
           STOP(1)
        end if

#ifdef CHANNEL
        allocate(source_term(ndime))
        !set the source term
        source_term(1) = utau*utau*rho0/delta !this is for Retau 180
        source_term(2) = 0.00d0
        source_term(3) = 0.00d0
#endif

        !*********************************************************************!
        ! Read mesh in Alya format                                            !
        !*********************************************************************!

        write(file_path,*) "./mesh/"
        write(*,*) "--| ALL MESH FILES MUST BE IN ",trim(adjustl(file_path))," !"
        write(*,*) "--| ENTER NAME OF MESH RELATED FILES :"
        call nvtxStartRange("Read mesh")
        !read(*,*) file_name
        write(file_name,*) "cube" ! Nsys
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
            write(*,*) "--| SPLITTING BOUNDARY NODES FROM DOFs..."
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
            write(*,*) '--| TOTAL FREE NODES := ',ndof
            write(*,*) '--| TOTAL BOUNDARY NODES := ',nbnodes

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

        WRITE(*,*) "--| ALLOCATING MAIN VARIABLES"
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
        allocate(mu_e(nelem))       ! Elemental viscosity
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
        do ipoin = 1,npoin

           if(coord(ipoin,2)<delta) then
               yp = coord(ipoin,2)*utau/mul
           else
              yp = abs(coord(ipoin,2)-2.0d0*delta)*utau/mul
           end if

          velo = utau*((1.0d0/0.41d0)*log(1.0d0+0.41d0*yp)+7.8d0*(1.0d0-exp(-yp/11.0d0)-(yp/11.0d0)*exp(-yp/3.0d0))) 

          call random_number(ti)
          !u(ipoin,1,2) = velo*(1.0d0 + 0.05d0*(sin(5.0d0*coord(ipoin,1)) -0.5d0))
          !u(ipoin,2,2) = velo*(0.05d0*(sin(5.0d0*coord(ipoin,2)) -0.5d0))
          !u(ipoin,3,2) = velo*(0.05d0*(sin(5.0d0*coord(ipoin,3)) -0.5d0))
          u(ipoin,1,2) = velo*(1.0d0 + 0.05d0*(ti(1) -0.5d0))
          u(ipoin,2,2) = velo*(0.05d0*(ti(2) -0.5d0))
          u(ipoin,3,2) = velo*(0.05d0*(ti(3) -0.5d0))

          pr(ipoin,2) = po

          rho(ipoin,2) = po/Rg/to
          
          e_int(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*(gamma_gas-1.0d0))
          Tem(ipoin,2) = pr(ipoin,2)/(rho(ipoin,2)*Rgas)
          E(ipoin,2) = rho(ipoin,2)*(0.5d0*dot_product(u(ipoin,:,2),u(ipoin,:,2))+e_int(ipoin,2))
          q(ipoin,1:ndime,2) = rho(ipoin,2)*u(ipoin,1:ndime,2)
        end do
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
        mu_e(:) = 0.0d0 ! Element syabilization viscosity
        !$acc end kernels
        call nvtxEndRange

        if (flag_real_diff == 1) then
           call constant_viscosity(npoin,0.000055d0,mu_fluid)
        else if (flag_real_diff == 0) then
           !$acc kernels
           mu_fluid(:) = 0.0d0
           !$acc end kernels
        else
           write(*,*) "--| DIFFUSION FLAG MUST BE EITHER 0 OR 1, NOT: ",flag_real_diff
           STOP(1)
        end if

        !*********************************************************************!
        ! Compute initial time-step size                                      !
        !*********************************************************************!
        
        if (flag_real_diff == 1) then
           call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid)
           write(*,*) "--| TIME STEP SIZE dt := ",dt,"s"
        else
           call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
           write(*,*) "--| TIME STEP SIZE dt := ",dt,"s"
        end if

        !
        ! Call VTK output (0th step)
        !
        write(*,*) "--| GENERATING 1st OUTPUT..."
        call nvtxStartRange("1st write")
        if (isPeriodic == 0) then
           call write_vtk_binary(isPeriodic,0,npoin,nelem,coord,connec, &
                                rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,nper)
        else
           print*, 'sub call: ok!'
           call write_vtk_binary(isPeriodic,0,npoin,nelem,coord,connec, &
                                rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,nper,masSla)
        end if
        call nvtxEndRange

        !*********************************************************************!
        ! Generate GLL table                                                  !
        !*********************************************************************!

        write(*,*) "--| GENERATING GAUSSIAN QUADRATURE TABLE..."

        call nvtxStartRange("Gaussian Quadrature")
        ! TODO: allow for more element types...

     !  if (ndime == 2) then ! 2D elements
     !     if (nnode == 3) then ! TRI03
     !         ngaus = 3
     !         write(*,*) 'NOT CODED YET!'
     !         STOP 1
     !     else if (nnode == 6) then ! TRI06
     !         ngaus = 7
     !         write(*,*) 'NOT CODED YET!'
     !         STOP 1
     !     else if (nnode == 4) then ! QUA04
     !         ngaus = 4
     !     else if (nnode == 9) then ! QUA09
     !         ngaus = 9
     !     end if
     !  else if (ndime == 3) then ! 3D elements
     !     if (nnode == 4) then ! TET04
     !        ngaus = 4
     !        write(*,*) 'NOT CODED YET!'
     !        stop 1
     !     else if (nnode == 8) then ! HEX08
     !        ngaus = 8
     !     else if (nnode == 27) then ! HEX27
     !        ngaus = 27
     !        write(*,*) '--| USING 27 GAUSS NODES PER ELEMENT!'
     !     else if (nnode == 64) then ! HEX64
     !        ngaus = 64
     !     else
     !        write(*,*) 'ELEMENT DOES NOT EXIST, OR IS NOT CODED YET!!'
     !        stop 1
     !     end if
     !  else
     !     write(*,*) 'ONLY 2D AND 3D ELEMMENTS SUPPORTED!'
     !     stop 1
     !  end if

        ! Option to run with 1 Gauss point per element, for testing and debugging. Will override the previous section.

#ifdef TGAUS
        ngaus = 1 ! Test value
#endif

        allocate(xgp(ngaus,ndime))
        allocate(wgp(ngaus))

        if (ndime == 2) then
           if (nnode == (porder+1)**2) then ! QUA_XX of order porder
              call gll_qua(xgp,wgp)
           else if (nnode == 3 .or. nnode == 6 .or. nnode == 10) then ! TRI_XX
              write(*,*) '--| NOT CODED YET!'
              STOP 1
           end if
        else if (ndime == 3) then
           if (nnode == (porder+1)**3) then ! HEX_XX
              call gll_hex(xgp,wgp)
           else if (nnode == 4 .or. nnode == 10 .or. nnode == 20) then ! TET_XX
              write(*,*) '--| NOT CODED YET!'
              STOP 1
           end if
        end if
        call nvtxEndRange

        !*********************************************************************!
        ! Generate N and dN for all GP                                        !
        !*********************************************************************!

        ! TODO: Allow for more element types

        write(*,*) "--| GENERATING SHAPE FUNCTIONS AND ISOPAR. DERIVATIVES..."
        call nvtxStartRange("N and dN")

        allocate(Ngp(ngaus,nnode),dNgp(ndime,nnode,ngaus))

        do igaus = 1,ngaus
           s = xgp(igaus,1)
           t = xgp(igaus,2)
           if (ndime == 2) then
              if (nnode == 4) then
                 call qua04(s,t,Ngp(igaus,:),dNgp(:,:,igaus))
              else if (nnode == 9) then
                 !call qua09(s,t,Ngp(igaus,:),dNgp(:,:,igaus))
                 write(*,*) '--| NOT CODED YET!'
                 STOP 1
              end if
           else if (ndime == 3) then
              z = xgp(igaus,3)
              if (nnode == 8) then
                 call hex08(s,t,z,Ngp(igaus,:),dNgp(:,:,igaus))
              else if (nnode == 27) then
                 call hex27(s,t,z,Ngp(igaus,:),dNgp(:,:,igaus))
              else if (nnode == 64) then
                 !call hex64(s,t,z,Ngp(igaus,:),dNgp(:,:,igaus))
              else
              end if
           end if
        end do
        call nvtxEndRange

        !*********************************************************************!
        ! Generate Jacobian related information                               !
        !*********************************************************************!

        ! TODO: create as a subroutine
        
        write(*,*) "--| GENERATING JACOBIAN RELATED INFORMATION..."

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
        write(*,*) '--| DOMAIN VOLUME := ',VolTot

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

        write(*,*) '--| COMPUTING LUMPED MASS MATRIX...'
        call nvtxStartRange("Lumped mass compute")
        allocate(Ml(npoin))
        call lumped_mass(nelem,npoin,connec,gpvol,Ngp,Ml)
        solver_type = 'LUMSO'
        call nvtxEndRange

        !
        ! Consisten mass: activate with CSR related operations
        !
        !write(*,*) '--| COMPUTING CONSISTENT MASS MATRIX...'
        !allocate(Mc(nzdom))
        !call consistent_mass(nelem,nnode,npoin,ngaus,connec,nzdom,rdom,cdom,gpvol,Ngp,Mc)

        !
        ! Set solver type
        !
        write(*,*) '--| ENTER REQUIRED SOLVER FOR MASS MATRIX:'
        write(*,*) '--| AVAILABLE SOLVERS ARE: LUMSO, APINV:'
        !read(*,*) solver_type
        write(solver_type,'(a)') "APINV" ! Nsys
        if (solver_type == 'APINV') then
            write(*,*) '--| ENTER NUMBER OF ITERATIONS FOR APINV SOLVER:'
            !read(*,*) ppow
            ppow = 2 ! Nsys
        end if
        write(*,*) '--| USING SOLVER ',solver_type,' FOR MASS MATRIX'

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
        open(unit=666,file="analysis.dat",status="new")
        time = 0.0d0
        P0 = 1.0d0
        T0 = 1.0d0/(1.4d0*Rgas*(0.1**2))
        rho0 = P0/(Rgas*T0)
        call volAvg_EK(nelem,npoin,connec,gpvol,Ngp,rho0,rho(:,2),u(:,:,2),EK)
        call write_EK(time,EK)
        write(*,*) "--| time   ,   EK"
        write(*,*) "--| ",time,"  |  ",EK

        counter = 1

        call nvtxStartRange("Start RK4")
        if (isPeriodic .eq. 0) then ! Case is not periodic
           if (nboun .eq. 0) then ! Case has no boundaries
              write(*,*) '--| ERROR: CASE MUST HAVE BOUNDARIES!'
              STOP(1)
           else ! Case has boundaries
              write(*,*) '--| NON-PERIODIC CASE WITH BOUNDARIES'
              do istep = 1,nstep

                 write(*,*) '   --| STEP: ', istep

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
                 call nvtxStartRange("RK4 step "//timeStep,istep)

                 call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                           ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                           rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w,mu_fluid, &
                           ndof,nbnodes,ldof,lbnodes,bound,bou_codes) ! Optional args

                 !
                 ! Advance with entropy viscosity
                 !
                 flag_predic = 0
                 call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                           ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                           rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w,mu_fluid, &
                           ndof,nbnodes,ldof,lbnodes,bound,bou_codes) ! Optional args

                 if (flag_real_diff == 1) then
                    call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid)
                    write(*,*) "DT := ",dt,"s"
                 else
                    call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
                    write(*,*) "DT := ",dt,"s"
                 end if

                 call nvtxEndRange

                 !
                 ! Call VTK output
                 !
                 if (istep == nsave) then
                    call nvtxStartRange("Output "//timeStep,istep)
                    call write_vtk_binary(isPeriodic,counter,npoin,nelem,coord,connec, &
                                         rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,nper)
                    nsave = nsave+nleap
                    call nvtxEndRange
                 end if

                 counter = counter+1

              end do
           end if
        else if (isPeriodic .eq. 1) then ! Case is periodic
           if (nboun .eq. 0) then ! Case has no boundaries
               write(*,*) '--| PERIODIC CASE WITH NO BOUNDARIES'
               do istep = 1,nstep

                  write(*,*) '   --| STEP: ', istep

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
                  call nvtxStartRange("RK4 step "//timeStep,istep)

                  call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                           ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                           rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w,mu_fluid)

                  !
                  ! Advance with entropy viscosity
                  !
                  flag_predic = 0
                  call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                           ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                           rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w,mu_fluid)

                  time = time+dt
                  call volAvg_EK(nelem,npoin,connec,gpvol,Ngp,rho0,rho(:,2),u(:,:,2),EK)
                  call write_EK(time,EK)
                  write(*,*) "--| time   ,   EK"
                  write(*,*) "--| ",time,"  |  ",EK

                  if (flag_real_diff == 1) then
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid)
                     write(*,*) "DT := ",dt,"s"
                  else
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
                     write(*,*) "DT := ",dt,"s"
                  end if

                  call nvtxEndRange

                  !
                  ! Call VTK output
                  !
                  if (istep == nsave) then
                     call nvtxStartRange("Output "//timeStep,istep)
                     call write_vtk_binary(isPeriodic,counter,npoin,nelem,coord,connec_orig, &
                                          rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,nper,masSla)
                     nsave = nsave+nleap
                     call nvtxEndRange
                  end if

                  counter = counter+1

               end do
               close(666)
            else
               write(*,*) '--| PERIODIC CASE WITH BOUNDARIES'
               do istep = 1,nstep

                  if (istep == nsave) write(*,*) '   --| STEP: ', istep

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
                  if (istep == nsave) write(timeStep,'(i4)') istep
                  call nvtxStartRange("RK4 step "//timeStep,istep)

                  call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                           ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                           rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w,mu_fluid, &
                           ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args

                  !
                  ! Advance with entropy viscosity
                  !
                  flag_predic = 0
                  call rk_4_main(flag_predic,flag_emac,nelem,nboun,npoin,npoin_w, &
                           ppow,connec,Ngp,dNgp,He,Ml,gpvol,dt,helem,Rgas,gamma_gas, &
                           rho,u,q,pr,E,Tem,e_int,mu_e,lpoin_w,mu_fluid, &
                           ndof,nbnodes,ldof,lbnodes,bound,bou_codes,source_term) ! Optional args

                  if (flag_real_diff == 1) then
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt,cfl_diff,mu_fluid)
                     write(*,*) "DT := ",dt,"s"
                  else
                     call adapt_dt_cfl(nelem,npoin,connec,helem,u(:,:,2),csound,cfl_conv,dt)
                     write(*,*) "DT := ",dt,"s"
                  end if

                  call nvtxEndRange

                  !
                  ! Call VTK output
                  !
                  if (istep == nsave) then
                     call nvtxStartRange("Output "//timeStep,istep)
                     call write_vtk_binary(isPeriodic,counter,npoin,nelem,coord,connec_orig, &
                                          rho(:,2),u(:,:,2),pr(:,2),E(:,2),mu_e,nper,masSla)
                     nsave = nsave+nleap
                     call nvtxEndRange
                  end if

                  counter = counter+1

               end do
            end if
        end if
        call nvtxEndRange

end program sod2d
