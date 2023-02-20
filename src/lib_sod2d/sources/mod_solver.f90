module mod_solver

      use mod_constants
      use mod_nvtx

      implicit none
      !-----------------------------------------------------------------------
      ! Creating variables for the GMRES solver
      !-----------------------------------------------------------------------
      integer(4), parameter                     :: maxSave=5
      real(rp)  , allocatable, dimension(:)     :: Jy_mass, Jy_ener, ymass, yener
      real(rp)  , allocatable, dimension(:)     :: Rmass_fix, Rener_fix, Dmass, Dener, Rmass, Rener
      real(rp)  , allocatable, dimension(:,:)   :: Jy_mom, ymom, Rmom_fix, Dmom, Rmom
      real(rp)  , allocatable, dimension(:,:)   :: Q_Mass, Q_Ener
      real(rp)  , allocatable, dimension(:,:,:) :: Q_Mom

      contains

              subroutine lumped_solver_scal(npoin,npoin_w,lpoin_w,Ml,R)

                      implicit none

                      integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
                      real(rp),    intent(in)    :: Ml(npoin)
                      real(rp),    intent(inout) :: R(npoin)
                      integer(4)                :: ipoin

                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         R(lpoin_w(ipoin)) = R(lpoin_w(ipoin))/Ml(lpoin_w(ipoin))
                      end do
                      !$acc end parallel

              end subroutine lumped_solver_scal

              subroutine lumped_solver_vect(npoin,npoin_w,lpoin_w,Ml,R)

                      implicit none

                      integer(4), intent(in)    :: npoin, npoin_w, lpoin_w(npoin_w)
                      real(rp),    intent(in)    :: Ml(npoin)
                      real(rp),    intent(inout) :: R(npoin,ndime)
                      integer(4)                :: idime, ipoin

                      !$acc parallel loop collapse(2)
                      do ipoin = 1,npoin_w
                         do idime = 1,ndime
                            R(lpoin_w(ipoin),idime) = R(lpoin_w(ipoin),idime)/Ml(lpoin_w(ipoin))
                         end do
                      end do
                      !$acc end  parallel loop

              end subroutine lumped_solver_vect

              subroutine approx_inverse_scalar(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,R)

                      use mass_matrix

                      implicit none

                      integer(4), intent(in)       :: nelem,npoin,ppow,npoin_w
                      integer(4), intent(in)       :: connec(nelem,nnode),lpoin_w(npoin_w)
                      real(rp),    intent(in)       :: gpvol(1,ngaus,nelem),Ngp(ngaus,nnode),Ml(npoin)
                      real(rp),    intent(inout)    :: R(npoin)
                      integer(4)                   :: ipoin, ipow
                      real(rp),    dimension(npoin) :: b, v, x

                      call nvtxStartRange("Scalar APINV")

                      !
                      ! Initialize series at k=0
                      !
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         b(lpoin_w(ipoin)) = R(lpoin_w(ipoin))
                         v(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                         x(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop

                      !
                      ! Step over sucessive powers
                      !
                      do ipow = 1,ppow
                         call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,v,b)
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            v(lpoin_w(ipoin)) = v(lpoin_w(ipoin))- &
                               (b(lpoin_w(ipoin))/Ml(lpoin_w(ipoin)))
                            x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+v(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end do
                      !$acc parallel loop
                      do ipoin = 1,npoin_w
                         R(lpoin_w(ipoin)) = x(lpoin_w(ipoin))
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine approx_inverse_scalar

              subroutine approx_inverse_vect(nelem,npoin,npoin_w,lpoin_w,connec,gpvol,Ngp,ppow,Ml,R)

                      use mass_matrix

                      implicit none

                      integer(4), intent(in)       :: nelem,npoin,ppow,npoin_w
                      integer(4), intent(in)       :: connec(nelem,nnode),lpoin_w(npoin_w)
                      real(rp),    intent(in)       :: gpvol(1,ngaus,nelem),Ngp(ngaus,nnode),Ml(npoin)
                      real(rp),    intent(inout)    :: R(npoin,ndime)
                      integer(4)                   :: ipoin, idime, ipow
                      real(rp),    dimension(npoin) :: b, v, x

                      call nvtxStartRange("Vector APINV")

                      do idime = 1,ndime
                         !
                         ! Initialize series at k=0
                         !
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            b(lpoin_w(ipoin)) = R(lpoin_w(ipoin),idime)
                            v(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                            x(lpoin_w(ipoin)) = b(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop

                         !
                         ! Step over sucessive powers
                         !
                         do ipow = 1,ppow
                            call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,v,b)
                            !$acc parallel loop
                            do ipoin = 1,npoin_w
                               v(lpoin_w(ipoin)) = v(lpoin_w(ipoin))- &
                                  (b(lpoin_w(ipoin))/Ml(lpoin_w(ipoin)))
                               x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+v(lpoin_w(ipoin))
                            end do
                            !$acc end parallel loop
                         end do
                         !$acc parallel loop
                         do ipoin = 1,npoin_w
                            R(lpoin_w(ipoin),idime) = x(lpoin_w(ipoin))
                         end do
                         !$acc end parallel loop
                      end do
                      call nvtxEndRange

              end subroutine approx_inverse_vect

              subroutine CSR_SpMV_scal(npoin,nzdom,rdom,cdom,Mc,v,u)

                      implicit none

                      integer(4), intent(in)  :: npoin, nzdom
                      integer(4), intent(in)  :: rdom(npoin+1), cdom(nzdom)
                      real(rp),    intent(in)  :: Mc(nzdom), v(npoin)
                      real(rp),    intent(out) :: u(npoin)
                      integer(4)              :: ipoin, izdom, jpoin, rowb, rowe

                      call nvtxStartRange("SPMV")
                      !$acc kernels
                      u(:) = 0.0_rp
                      !$acc end kernels
                      !
                      !$acc parallel loop
                      do ipoin = 1,npoin
                         !
                         ! Get CSR section for row ipoin
                         !
                         rowb = rdom(ipoin)+1
                         rowe = rdom(ipoin+1)
                         !
                         ! Loop inside CSR section
                         !
                         !!$acc loop seq
                         do izdom = rowb,rowe
                            jpoin = cdom(izdom) ! Col. index
                            u(ipoin) = u(ipoin)+Mc(izdom)*v(jpoin) ! Dot product
                         end do
                      end do
                      !$acc end parallel loop
                      call nvtxEndRange

              end subroutine CSR_SpMV_scal

              subroutine conjGrad_scalar(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,R)

                 use mass_matrix

                 implicit none

                 integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
                 real(rp)   , intent(in)    :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                 real(rp)   , intent(inout) :: R(npoin)
                 integer(4)                :: ipoin, iter
                 real(rp), dimension(npoin) :: x, r0, p0, q, v, b
                 real(rp)                   :: Q1, Q2, T1, alpha, beta


                 call nvtxStartRange("CG solver scalar")
                 !$acc kernels
                 x(:) = 0.0_rp
                 r0(:) = 0.0_rp
                 p0(:) = 0.0_rp
                 q(:) = 0.0_rp
                 v(:) = 0.0_rp
                 b(:) = 0.0_rp
                 !$acc end kernels
                 !
                 ! Initialize solver
                 !
                 !$acc parallel loop
                 do ipoin = 1,npoin_w
                    x(lpoin_w(ipoin)) = R(lpoin_w(ipoin))
                    b(lpoin_w(ipoin)) = R(lpoin_w(ipoin))
                 end do
                 !$acc end parallel loop
                 call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,x,q) ! A*x0
                 !$acc parallel loop
                 do ipoin = 1,npoin_w
                    r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-q(lpoin_w(ipoin)) ! b-A*x0
                    p0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin)) ! s0 = r0
                 end do

                 !
                 ! Start iterations
                 !
                 do iter = 1,maxSave
                    call nvtxStartRange("Iteration")
                    call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,p0,q) ! A*s_k-1
                    Q1 = 0.0_rp
                    Q2 = 0.0_rp
                    !$acc parallel loop reduction(+:Q1,Q2)
                    do ipoin = 1,npoin_w
                       Q1 = Q1+p0(lpoin_w(ipoin))*r0(lpoin_w(ipoin)) ! <s_k-1,r_k-1>
                       Q2 = Q2+p0(lpoin_w(ipoin))*q(lpoin_w(ipoin)) ! <s_k-1,A*s_k-1>
                    end do
                    !$acc end parallel loop
                    alpha = Q1/Q2
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+alpha*p0(lpoin_w(ipoin)) ! x_k = x_k-1 + alpha*s_k-1
                    end do
                    !$acc end parallel loop
                    call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,x,v) ! A*x_k
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-v(lpoin_w(ipoin)) ! b-A*x_k
                    end do
                    !$acc end parallel loop
                    T1 = 0.0_rp
                    !$acc parallel loop reduction(+:T1)
                    do ipoin = 1,npoin
                       T1 = T1+r0(ipoin)*r0(ipoin)
                    end do
                    !$acc end parallel loop
                    !
                    ! Stop cond
                    !
                    if (sqrt(T1) .lt. tol) then
                       call nvtxEndRange
                       exit
                    end if
                    !
                    ! Update p
                    !
                    T1 = 0.0_rp
                    !$acc parallel loop reduction(+:T1)
                    do ipoin = 1,npoin
                       T1 = T1+r0(ipoin)*q(ipoin) ! <r_k,A*s_k-1>
                    end do
                    !$acc end parallel loop
                    beta = T1/Q1
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       p0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))-beta*p0(lpoin_w(ipoin)) ! s_k = r_k+beta*s_k-1
                    end do
                    !$acc end parallel loop
                    call nvtxEndRange
                 end do
                 if (iter == maxSave) then
                    write(1,*) "--| TOO MANY ITERATIONS!"
                    call nvtxEndRange
                    stop 1
                 end if
                 !$acc kernels
                 R(:) = x(:)
                 !$acc end kernels
                 call nvtxEndRange

              end subroutine conjGrad_scalar

              subroutine conjGrad_vector(nelem,npoin,npoin_w,connec,lpoin_w,gpvol,Ngp,R)

                 use mass_matrix

                 implicit none

                 integer(4), intent(in)    :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
                 real(rp)   , intent(in)    :: gpvol(1,ngaus,nelem), Ngp(ngaus,nnode)
                 real(rp)   , intent(inout) :: R(npoin,ndime)
                 integer(4)                :: ipoin, idime, iter
                 real(rp), dimension(npoin) :: x, r0, p0, q, v, b
                 real(rp)                   :: Q1, Q2, T1, alpha, beta


                 call nvtxStartRange("CG solver vector")
                 do idime = 1,ndime
                    call nvtxStartRange("Dimension")
                    !$acc kernels
                    x(:) = 0.0_rp
                    r0(:) = 0.0_rp
                    p0(:) = 0.0_rp
                    q(:) = 0.0_rp
                    v(:) = 0.0_rp
                    b(:) = 0.0_rp
                    !$acc end kernels
                    !
                    ! Initialize solver
                    !
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       x(lpoin_w(ipoin)) = R(lpoin_w(ipoin),idime)
                       b(lpoin_w(ipoin)) = R(lpoin_w(ipoin),idime)
                    end do
                    !$acc end parallel loop
                    call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,x,q) ! A*x0
                    !$acc parallel loop
                    do ipoin = 1,npoin_w
                       r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-q(lpoin_w(ipoin)) ! b-A*x0
                       p0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin)) ! s0 = r0
                    end do
                    !$acc end parallel loop
                    !
                    ! Start iterations
                    !
                    do iter = 1,maxSave
                       call nvtxStartRange("Iteration")
                       call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,p0,q) ! A*s_k-1
                       Q1 = 0.0_rp
                       Q2 = 0.0_rp
                       !$acc parallel loop reduction(+:Q1,Q2)
                       do ipoin = 1,npoin
                          Q1 = Q1+p0(ipoin)*r0(ipoin) ! <s_k-1,r_k-1>
                          Q2 = Q2+p0(ipoin)*q(ipoin) ! <s_k-1,A*s_k-1>
                       end do
                       !$acc end parallel loop
                       alpha = Q1/Q2
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          x(lpoin_w(ipoin)) = x(lpoin_w(ipoin))+alpha*p0(lpoin_w(ipoin)) ! x_k = x_k-1 + alpha*s_k-1
                       end do
                       !$acc end parallel loop
                       call cmass_times_vector(nelem,npoin,connec,gpvol,Ngp,x,v) ! A*x_k
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          r0(lpoin_w(ipoin)) = b(lpoin_w(ipoin))-v(lpoin_w(ipoin)) ! b-A*x_k
                       end do
                       !$acc end parallel loop
                       T1 = 0.0_rp
                       !$acc parallel loop reduction(+:T1)
                       do ipoin = 1,npoin
                          T1 = T1+r0(ipoin)*r0(ipoin)
                       end do
                       !$acc end parallel loop
                       !
                       ! Stop cond
                       !
                       if (sqrt(T1) .lt. tol) then
                          call nvtxEndRange
                          exit
                       end if
                       !
                       ! Update p
                       !
                       T1 = 0.0_rp
                       !$acc parallel loop reduction(+:T1)
                       do ipoin = 1,npoin
                          T1 = T1+r0(ipoin)*q(ipoin) ! <r_k,A*s_k-1>
                       end do
                       !$acc end parallel loop
                       beta = T1/Q1
                       !$acc parallel loop
                       do ipoin = 1,npoin_w
                          p0(lpoin_w(ipoin)) = r0(lpoin_w(ipoin))-beta*p0(lpoin_w(ipoin)) ! s_k = r_k+beta*s_k-1
                       end do
                       !$acc end parallel loop
                       call nvtxEndRange
                    end do
                    if (iter == maxSave) then
                       write(1,*) "--| TOO MANY ITERATIONS!"
                       call nvtxEndRange
                       stop 1
                    end if
                    !$acc parallel loop
                    do ipoin = 1,npoin
                       R(ipoin,idime) = x(ipoin)
                    end do
                    !$acc end parallel loop
                    call nvtxEndRange
                 end do
                 call nvtxEndRange

              end subroutine conjGrad_vector

              subroutine gmres_full(npoin, flag_gmres_mem_alloc, flag_gmres_mem_free)
                  implicit none
                  logical   , intent(in)    :: flag_gmres_mem_alloc, flag_gmres_mem_free
                  integer(4), intent(in) :: npoin

                  ! Allocate the memory for the gmres solver if not yet allocated
                  if (flag_gmres_mem_alloc .eqv. .true.) then
                     allocate(Jy_mass(npoin), Jy_mom(npoin,ndime), Jy_ener(npoin))
                     allocate(ymass(npoin), ymom(npoin,ndime), yener(npoin))
                     allocate(Rmass_fix(npoin), Rmom_fix(npoin,ndime), Rener_fix(npoin))
                     allocate(Rmass(npoin), Rmom(npoin,ndime), Rener(npoin))
                     allocate(Dmass(npoin), Dmom(npoin,ndime), Dener(npoin))
                     allocate(Q_Mass(npoin,maxSave+1), Q_Mom(npoin,ndime,maxSave+1), Q_Ener(npoin,maxSave+1))
                  end if

                  ! Form the approximate inv(Ml)*J*y for all equations
                  !call form_approx_Jy()

                  ! Initialize the solver
                  !call init_gmres()

                  ! If memory not needed anymore, deallocate arrays
                  if (flag_gmres_mem_free .eqv. .true.) then
                     deallocate(Jy_mass, Jy_mom, Jy_ener)
                     deallocate(ymass, ymom, yener)
                     deallocate(Rmass_fix, Rmom_fix, Rener_fix)
                     deallocate(Rmass, Rmom, Rener)
                     deallocate(Dmass, Dmom, Dener)
                     deallocate(Q_Mass, Q_Mom, Q_Ener)
                  end if

              end subroutine gmres_full

              !-----------------------------------------------------------------------
              ! Auxiliary routines for the gmres solver
              !-----------------------------------------------------------------------

              subroutine arnoldi_iter()
                  implicit none
              end subroutine arnoldi_iter

              subroutine form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                        rho,u,q,pr,E,Tem,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                        zmass,zmom,zener,flag_gmres_form_fix)
                  use mod_comms
                  use elem_convec
                  use elem_diffu
                  implicit none
                  logical   , intent(in) :: flag_gmres_form_fix
                  integer(4), intent(in) :: nelem, npoin, npoin_w, connec(nelem,nnode), lpoin_w(npoin_w)
                  integer(4), intent(in) :: atoIJK(nnode),invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                  real(rp)  , intent(in) :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus), gpvol(1,ngaus,nelem), xgp(ngaus,ndime)
                  real(rp)  , intent(in) :: He(ndime,ndime,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), Ml(npoin)
                  real(rp)  , intent(in) :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin), Cp, Prt
                  real(rp)  , intent(in) :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
                  real(rp)  , intent(in) :: zmass(npoin), zener(npoin), zmom(npoin,ndime)
                  integer(4)             :: idime
                  real(rp)  , parameter  :: eps = 1.0e-7

                  ! Form the R(u^n) arrays if not formed already
                  if (flag_gmres_form_fix .eqv. .true.) then
                     call full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                                          gmshAtoI, gmshAtoJ, gmshAtoK, u, q, rho, pr, E, Rmass, Rmom, Rener)
                     call full_diffusion_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                                             gmshAtoI, gmshAtoJ, gmshAtoK, Cp, Prt, rho, u, Tem, &
                                             mu_fluid, mu_e, mu_sgs, Ml, Dmass, Dmom, Dener)
                     Rmass_fix(:) = Rmass(:) + Dmass(:)
                     Rener_fix(:) = Rener(:) + Dener(:)
                     Rmom_fix(:,:) = Rmom(:,:) + Dmom(:,:)
                  end if

                  ! Form the R(u^n + eps*y0) arrays
                  call full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                                       gmshAtoI, gmshAtoJ, gmshAtoK, u, q+eps*zmom, rho+eps*zmass, pr, E+eps*zener, Rmass, Rmom, Rener)
                  call full_diffusion_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                                          gmshAtoI, gmshAtoJ, gmshAtoK, Cp, Prt, rho+eps*zmass, u, Tem, &
                                          mu_fluid, mu_e, mu_sgs, Ml, Dmass, Dmom, Dener)
                  Rmass(:) = Rmass(:) + Dmass(:)
                  Rener(:) = Rener(:) + Dener(:)
                  Rmom(:,:) = Rmom(:,:) + Dmom(:,:)

                  ! Form the J*y arrays
                  Jy_mass(:) = (Rmass(:) - Rmass_fix(:))/eps
                  Jy_ener(:) = (Rener(:) - Rener_fix(:))/eps
                  Jy_mom(:,:) = (Rmom(:,:) - Rmom_fix(:,:))/eps

                  ! Communicate before applying ML
                  if(mpi_size.ge.2) then
                     call nvtxStartRange("MPI_comms_tI")
                     call mpi_halo_atomic_update_float(Jy_mass)
                     call mpi_halo_atomic_update_float(Jy_ener)
                     do idime = 1,ndime
                        call mpi_halo_atomic_update_float(Jy_mom(:,idime))
                     end do
                     call nvtxEndRange
                  end if

                  ! Multiply residuals by inverse Ml
                  call lumped_solver_scal(npoin, npoin_w, lpoin_w, Ml, Jy_mass)
                  call lumped_solver_scal(npoin, npoin_w, lpoin_w, Ml, Jy_ener)
                  call lumped_solver_vect(npoin, npoin_w, lpoin_w, Ml, Jy_mom)

              end subroutine form_approx_Jy

              subroutine init_gmres(npoin, bmass, bmom, bener, dt, gammaRK)
                  implicit none
                  integer(4), intent(in) :: npoin
                  real(rp)  , intent(in) :: bmass(npoin), bmom(npoin,ndime), bener(npoin), dt, gammaRK
                  integer(4)             :: ipoin, idime
                  real(rp)               :: aux

                  ! Zero Q_* arrays
                  Q_Mass(:,:) = 0.0_rp
                  Q_Ener(:,:) = 0.0_rp
                  Q_Mom(:,:,:) = 0.0_rp

                  ! Add the remaining terms
                  ymass(:) = ymass(:)/gammaRK*dt + Jy_mass(:)
                  yener(:) = yener(:)/gammaRK*dt + Jy_ener(:)
                  ymom(:,:) = ymom(:,:)/gammaRK*dt + Jy_mom(:,:)

                  ! Compute each r = b-Ax
                  Q_Mass(:,1) = bmass(:) - ymass(:)
                  Q_Ener(:,1) = bener(:) - yener(:)
                  Q_Mom(:,:,1) = bmom(:,:) - ymom(:,:)

                  ! Normalize each residual
                  aux = 0.0
                  do ipoin = 1,npoin
                     aux = aux + Q_Mass(ipoin,1)**2
                  end do
                  Q_Mass(:,1) = Q_Mass(:,1)/sqrt(aux)
                  aux = 0.0
                  do ipoin = 1,npoin
                     aux = aux + Q_Ener(ipoin,1)**2
                  end do
                  Q_Ener(:,1) = Q_Ener(:,1)/sqrt(aux)
                  do idime = 1,ndime
                     aux = 0.0
                     do ipoin = 1,npoin
                        aux = aux + Q_Mom(ipoin,idime,1)**2
                     end do
                     Q_Mom(:,idime,1) = Q_Mom(:,idime,1)/sqrt(aux)
                  end do

              end subroutine init_gmres

end module mod_solver
