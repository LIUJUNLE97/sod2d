module mod_solver

      use mod_constants
      use mod_comms
      use mod_mpi
      use mod_nvtx

      implicit none
      !-----------------------------------------------------------------------
      ! Creating variables for the GMRES solver
      !-----------------------------------------------------------------------
      real(rp)                                  :: norm_bmass, norm_bener, norm_bmom(ndime),epsQ(3),epsR,epsE
      real(rp)                                  :: err_mass, err_ener, err_mom(ndime)
      real(rp)                                  :: e1_mass(maxIter+1), e1_ener(maxIter+1), e1_mom(maxIter+1,ndime)
      real(rp)                                  :: beta_mass(maxIter+1), beta_ener(maxIter+1), beta_mom(maxIter+1,ndime)
      real(rp)  , dimension(maxIter)            :: cs_mass, cs_ener, sn_mass, sn_ener
      real(rp)  , dimension(maxIter,ndime)      :: cs_mom, sn_mom
      real(rp)  , dimension(maxIter)            :: updMass, updEner
      real(rp)  , dimension(maxIter,ndime)      :: updMom
      real(rp)  , allocatable, dimension(:)     :: Jy_mass, Jy_ener, ymass, yener
      real(rp)  , allocatable, dimension(:)     :: Rmass_fix, Rener_fix, Dmass, Dener, Rmass, Rener,pEner,pMass
      real(rp)  , allocatable, dimension(:,:)   :: Jy_mom, ymom, Rmom_fix, Dmom, Rmom,pMom
      real(rp)  , allocatable, dimension(:,:)   :: Q_Mass, Q_Ener, H_mass, H_ener
      real(rp)  , allocatable, dimension(:,:,:) :: Q_Mom, H_mom
      logical                                   :: flag_gmres_mem_alloc=.true.
      logical                                   :: flag_gmres_mem_free=.false.
      real(rp)                                  :: eps=1e-4

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
                 do iter = 1,maxIter
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
                 if (iter == maxIter) then
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
                    do iter = 1,maxIter
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
                    if (iter == maxIter) then
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

              subroutine jacobi_full(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                    atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                    rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                    gammaRK,dt,bmass,bmom,bener,mass_sol,mom_sol,ener_sol,istep)
                  implicit none
                  integer(4), intent(in)    :: nelem, npoin, npoin_w, lpoin_w(npoin_w), connec(nelem,nnode)
                  integer(4), intent(in)    :: atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                  real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,ngaus,nnode), He(ndime,ndime,ngaus,nelem)
                  real(rp)  , intent(in)    :: gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
                  real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin)
                  real(rp)  , intent(in)    :: Rgas, gamma_gas, Cp, Prt, gammaRK, dt, Ml(npoin)
                  real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,nnode), mu_sgs(nelem,nnode)
                  real(rp)  , intent(in)    :: bmass(npoin), bmom(npoin,ndime), bener(npoin)
                  real(rp)  , intent(inout) :: mass_sol(npoin), mom_sol(npoin,ndime), ener_sol(npoin)
                  integer(4), intent(in)    :: istep
                  integer(4)                :: ik, jk, kk, ipoin, idime
                  real(rp)                  :: errMax, errTol, relax=1.0_rp
                  real(8)                  :: auxN(5),auxN2(5),aux

                  errTol = tol

                  ! Allocate the memory for the gmres solver if not yet allocated
                  if (flag_gmres_mem_alloc .eqv. .true.) then
                     allocate(Jy_mass(npoin), Jy_mom(npoin,ndime), Jy_ener(npoin))
                     allocate(ymass(npoin), ymom(npoin,ndime), yener(npoin))
                     allocate(Rmass_fix(npoin), Rmom_fix(npoin,ndime), Rener_fix(npoin))
                     allocate(Rmass(npoin), Rmom(npoin,ndime), Rener(npoin))
                     allocate(Dmass(npoin), Dmom(npoin,ndime),Dener(npoin),pEner(npoin),pMass(npoin),pMom(npoin,ndime))
                     allocate(Q_Mass(npoin,maxIter+1), Q_Mom(npoin,ndime,maxIter+1), Q_Ener(npoin,maxIter+1))
                     allocate(H_mass(maxIter+1,maxIter), H_mom(maxIter+1,maxIter,ndime), H_ener(maxIter+1,maxIter))
                     flag_gmres_mem_alloc = .false.
                  end if

                  epsR = 1e-3
                  epsE = 1e-3
                  epsQ(1) = 1e-3
                  epsQ(2) = 1e-3
                  epsQ(3) = 1e-3
                  if(istep == 1) then 
                     call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                        rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                        mass_sol,mom_sol,ener_sol,.true.)
                  end if

                  do ik = 1,maxIter
                     if(ik .lt. 3) then 
                        auxN(:) = 0.0_rp
                        !$acc parallel loop
                        do ipoin = 1,npoin_w
                           auxN(1) = auxN(1) + real(mass_sol(lpoin_w(ipoin))**2,8)
                           auxN(2) = auxN(2) + real(ener_sol(lpoin_w(ipoin))**2,8)
                           !$acc lloop seq
                           do idime = 1,ndime
                              auxN(idime+2) = auxN(idime+2) + real(mom_sol(lpoin_w(ipoin),idime)**2,8)
                           end do
                        end do
                        !$acc end parallel loop

                        call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                        if(auxN2(1)<epsilon(errTol)) auxN2(1) = 1.0_rp
                        if(auxN2(2)<epsilon(errTol)) auxN2(2) = 1.0_rp
                        if(auxN2(3)<epsilon(errTol)) auxN2(3) = 1.0_rp
                        if(auxN2(4)<epsilon(errTol)) auxN2(4) = 1.0_rp
                        if(auxN2(5)<epsilon(errTol)) auxN2(5) = 1.0_rp

                        epsR = eps/sqrt(real(auxN2(1),rp))
                        epsE = eps/sqrt(real(auxN2(2),rp))
                        epsQ(1) = eps/sqrt(real(auxN2(3),rp))
                        epsQ(2) = eps/sqrt(real(auxN2(4),rp))
                        epsQ(3) = eps/sqrt(real(auxN2(5),rp))
                     end if
                     call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                        rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                        mass_sol,mom_sol,ener_sol,.false.)

                     auxN(:) = 0.0_rp
                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        aux = mass_sol(lpoin_w(ipoin))
                        mass_sol(lpoin_w(ipoin)) = mass_sol(lpoin_w(ipoin))*(1.0_rp-relax) + relax*gammaRK*dt*(bmass(lpoin_w(ipoin)) - Jy_mass(lpoin_w(ipoin)))
                        auxN(1) = auxN(1) + real((aux-mass_sol(lpoin_w(ipoin)))**2 ,8)
                        aux = ener_sol(lpoin_w(ipoin))
                        ener_sol(lpoin_w(ipoin)) = ener_sol(lpoin_w(ipoin))*(1.0_rp-relax) + relax*gammaRK*dt*(bener(lpoin_w(ipoin)) - Jy_ener(lpoin_w(ipoin)))
                        auxN(2) = auxN(2) + real((aux-ener_sol(lpoin_w(ipoin)))**2 ,8)
                        !$acc loop seq
                        do idime = 1,ndime
                           aux = mom_sol(lpoin_w(ipoin),idime)
                           mom_sol(lpoin_w(ipoin),idime) = mom_sol(lpoin_w(ipoin),idime)*(1.0_rp-relax) + relax*gammaRK*dt*(bmom(lpoin_w(ipoin),idime) - Jy_mom(lpoin_w(ipoin),idime))
                           auxN(2+idime) = auxN(2+idime) + real((aux-mom_sol(lpoin_w(ipoin),idime))**2 ,8)
                        end do
                     end do

                     call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                     errMax = real(sqrt(maxval(auxN2(:))),rp)

                     if(mpi_rank.eq.0)print*, " err ",errMax," it ",ik,' emass ',sqrt(auxN2(1))," eener ",sqrt(auxN2(2))," emom ", sqrt(auxN2(3))," ", sqrt(auxN2(4))," ",sqrt(auxN2(5))

                     if(errMax .lt. tol) exit

                  end do


                  ! If memory not needed anymore, deallocate arrays
                  if (flag_gmres_mem_free .eqv. .true.) then
                     deallocate(Jy_mass, Jy_mom, Jy_ener)
                     deallocate(ymass, ymom, yener)
                     deallocate(Rmass_fix, Rmom_fix, Rener_fix)
                     deallocate(Rmass, Rmom, Rener)
                     deallocate(Dmass, Dmom, Dener)
                     deallocate(Q_Mass, Q_Mom, Q_Ener)
                  end if

              end subroutine jacobi_full

              subroutine gmres_full(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                    atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                    rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                    gammaRK,dt,bmass,bmom,bener,mass_sol,mom_sol,ener_sol,istep)
                  implicit none
                  integer(4), intent(in)    :: nelem, npoin, npoin_w, lpoin_w(npoin_w), connec(nelem,nnode)
                  integer(4), intent(in)    :: atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                  real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,ngaus,nnode), He(ndime,ndime,ngaus,nelem)
                  real(rp)  , intent(in)    :: gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
                  real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin)
                  real(rp)  , intent(in)    :: Rgas, gamma_gas, Cp, Prt, gammaRK, dt, Ml(npoin)
                  real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,nnode), mu_sgs(nelem,nnode)
                  real(rp)  , intent(in)    :: bmass(npoin), bmom(npoin,ndime), bener(npoin)
                  real(rp)  , intent(inout) :: mass_sol(npoin), mom_sol(npoin,ndime), ener_sol(npoin)
                  integer(4), intent(in)    :: istep
                  integer(4)                :: ik, jk, kk, ipoin, idime
                  real(rp)                  :: errMax, errTol, aux
                  real(8)                  :: auxN(5),auxN2(5)

                  errTol = tol

                  !$acc kernels
                  ymass(:) = mass_sol(:)
                  yener(:) = ener_sol(:)
                  ymom(:,:) = mom_sol(:,:)
                  !$acc end kernels

                  ! Allocate the memory for the gmres solver if not yet allocated
                  if (flag_gmres_mem_alloc .eqv. .true.) then
                     allocate(Jy_mass(npoin), Jy_mom(npoin,ndime), Jy_ener(npoin))
                     allocate(ymass(npoin), ymom(npoin,ndime), yener(npoin))
                     allocate(Rmass_fix(npoin), Rmom_fix(npoin,ndime), Rener_fix(npoin))
                     allocate(Rmass(npoin), Rmom(npoin,ndime), Rener(npoin))
                     allocate(Dmass(npoin), Dmom(npoin,ndime),Dener(npoin),pEner(npoin),pMass(npoin),pMom(npoin,ndime))
                     allocate(Q_Mass(npoin,maxIter+1), Q_Mom(npoin,ndime,maxIter+1), Q_Ener(npoin,maxIter+1))
                     allocate(H_mass(maxIter+1,maxIter), H_mom(maxIter+1,maxIter,ndime), H_ener(maxIter+1,maxIter))
                     flag_gmres_mem_alloc = .false.
                  end if

                  auxN(:) = 0.0_rp
                  aux = 0.0_rp
                  !$acc parallel loop reduction(+:aux)
                  do ipoin = 1,npoin_w
                     aux = aux + real(bmass(lpoin_w(ipoin))**2,8)
                  end do
                  !$acc end parallel loop
                  auxN(1) = aux

                  aux = 0.0_rp
                  !$acc parallel loop reduction(+:aux)
                  do ipoin = 1,npoin_w
                     aux = aux + real(bener(lpoin_w(ipoin))**2,8)
                  end do
                  !$acc end parallel loop
                  auxN(2) = aux

                  do idime = 1,ndime
                     aux = 0.0_rp
                     !$acc parallel loop reduction(+:aux)
                     do ipoin = 1,npoin_w
                        aux = aux + real(bmom(lpoin_w(ipoin),idime)**2,8)
                     end do
                     !$acc end parallel loop
                     auxN(idime+2) = aux
                  end do

                  call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                  if(auxN2(1)<epsilon(errTol)) auxN2(1) = 1.0_rp
                  if(auxN2(2)<epsilon(errTol)) auxN2(2) = 1.0_rp
                  if(auxN2(3)<epsilon(errTol)) auxN2(3) = 1.0_rp
                  if(auxN2(4)<epsilon(errTol)) auxN2(4) = 1.0_rp
                  if(auxN2(5)<epsilon(errTol)) auxN2(5) = 1.0_rp

                  norm_bmass   = sqrt(real(auxN2(1),rp))
                  norm_bener   = sqrt(real(auxN2(2),rp))
                  norm_bmom(1) = sqrt(real(auxN2(3),rp))
                  norm_bmom(2) = sqrt(real(auxN2(4),rp))
                  norm_bmom(3) = sqrt(real(auxN2(5),rp))

                  epsR = 1e-3
                  epsE = 1e-3
                  epsQ(1) = 1e-3
                  epsQ(2) = 1e-3
                  epsQ(3) = 1e-3

                  if(istep == 1) then 
                     call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                        rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                        mass_sol,mom_sol,ener_sol,.true.)
                  end if

                  call init_gmres(npoin,npoin_w,lpoin_w,bmass,bmom,bener,dt,gammaRK)

                  outer:do ik = 1,maxIter
                     ! TODO: put this inside Jy and make gpu friendly
                     ! Verify: if this is properly done, ||Q|| = 1, so step should be
                     auxN(:) = 0.0_rp
                     if(ik .lt. 3) then 
                        !$acc parallel loop
                        do ipoin = 1,npoin_w
                           auxN(1) = auxN(1) + real(Q_mass(lpoin_w(ipoin),ik)**2,8)
                           auxN(2) = auxN(2) + real(Q_ener(lpoin_w(ipoin),ik)**2,8)
                           !$acc lloop seq
                           do idime = 1,ndime
                              auxN(idime+2) = auxN(idime+2) + real(Q_mom(lpoin_w(ipoin),idime,ik)**2,8)
                           end do
                        end do
                        !$acc end parallel loop

                        call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                        if(auxN2(1)<epsilon(errTol)) auxN2(1) = 1.0_rp
                        if(auxN2(2)<epsilon(errTol)) auxN2(2) = 1.0_rp
                        if(auxN2(3)<epsilon(errTol)) auxN2(3) = 1.0_rp
                        if(auxN2(4)<epsilon(errTol)) auxN2(4) = 1.0_rp
                        if(auxN2(5)<epsilon(errTol)) auxN2(5) = 1.0_rp

                        epsR = eps/sqrt(real(auxN2(1),rp))
                        epsE = eps/sqrt(real(auxN2(2),rp))
                        epsQ(1) = eps/sqrt(real(auxN2(3),rp))
                        epsQ(2) = eps/sqrt(real(auxN2(4),rp))
                        epsQ(3) = eps/sqrt(real(auxN2(5),rp))

                     end if

                     !$acc kernels
                     pMass(:) = Q_Mass(:,ik)
                     pEner(:) = Q_Ener(:,ik)
                     pMom(:,:) = Q_Mom(:,:,ik)
                     !$acc end kernels

                     call  smooth_gmres(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                      atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                      rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                      gammaRK,dt,3)

                     call arnoldi_iter(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                        rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                        gammaRK,dt)

                     call apply_givens_rotation(ik)

                     aux = cs_mass(ik)*beta_mass(ik) + sn_mass(ik)*beta_mass(ik+1)
                     beta_mass(ik+1) = -sn_mass(ik)*beta_mass(ik) + cs_mass(ik)*beta_mass(ik+1)
                     beta_mass(ik) = aux
                     err_mass = abs(beta_mass(ik+1))/norm_bmass
                     aux = cs_ener(ik)*beta_ener(ik) + sn_ener(ik)*beta_ener(ik+1)
                     beta_ener(ik+1) = -sn_ener(ik)*beta_ener(ik) + cs_ener(ik)*beta_ener(ik+1)
                     beta_ener(ik) = aux
                     err_ener = abs(beta_ener(ik+1))/norm_bener
                     do idime = 1,ndime
                        aux = cs_mom(ik,idime)*beta_mom(ik,idime) + sn_mom(ik,idime)*beta_mom(ik+1,idime)
                        beta_mom(ik+1,idime) = -sn_mom(ik,idime)*beta_mom(ik,idime) + cs_mom(ik,idime)*beta_mom(ik+1,idime)
                        beta_mom(ik,idime) = aux
                        err_mom(idime) = abs(beta_mom(ik+1,idime))/norm_bmom(idime)
                     end do

                     errmax = 0.0_rp
                     errmax = max(errmax,abs(err_mass))
                     errmax = max(errmax,abs(err_ener))
                     do idime = 1,ndime
                        errmax = max(errmax,abs(err_mom(idime)))
                     end do

                     if(mpi_rank.eq.0)print*, " err ",errMax," it ",ik,' emass ',err_mass," eener ",err_ener," emom ",err_mom

                     if (errMax .lt. errTol) then
                        ! Set upd* to beta_*
                        !$acc kernels
                        updMass(1:maxIter) = beta_mass(1:maxIter)
                        updEner(1:maxIter) = beta_ener(1:maxIter)
                        updMom(1:maxIter,1:ndime) = beta_mom(1:maxIter,1:ndime)
                        !$acc end kernels
                        do jk=ik,1,-1
                           updMass(jk) = updMass(jk)/H_mass(jk,jk)
                           updEner(jk) = updEner(jk)/H_ener(jk,jk)
                           do kk = jk-1,1,-1
                              updMass(kk) = updMass(kk) - H_mass(kk,jk)*updMass(jk)
                              updEner(kk) = updEner(kk) - H_ener(kk,jk)*updEner(jk)
                           end do
                        end do
                        do idime = 1,ndime
                           do jk=ik,1,-1
                              updMom(jk,idime) = updMom(jk,idime)/H_mom(jk,jk,idime)
                              do kk = jk-1,1,-1
                                 updMom(kk,idime) = updMom(kk,idime) - H_mom(kk,jk,idime)*updMom(jk,idime)
                              end do
                           end do
                        end do

                        !$acc parallel loop
                        do ipoin = 1,npoin_w
                           aux = 0.0_rp
                           do kk = 1,ik
                              aux = aux + Q_Mass(lpoin_w(ipoin),kk)*updMass(kk)
                           end do
                           mass_sol(lpoin_w(ipoin)) = mass_sol(lpoin_w(ipoin)) + aux
                           aux = 0.0_rp
                           do kk = 1,ik
                              aux = aux + Q_Ener(lpoin_w(ipoin),kk)*updEner(kk)
                           end do
                           ener_sol(lpoin_w(ipoin)) = ener_sol(lpoin_w(ipoin)) + aux
                           do idime = 1,ndime
                              aux = 0.0_rp
                              do kk = 1,ik
                                 aux = aux + Q_Mom(lpoin_w(ipoin),idime,kk)*updMom(kk,idime)
                              end do
                              mom_sol(lpoin_w(ipoin),idime) = mom_sol(lpoin_w(ipoin),idime) + aux
                           end do
                        end do
                        !$acc end parallel loop
                        exit outer
                     end if

                  end do outer

                  if(ik == (maxIter+1)) then
                     !$acc kernels
                     updMass(1:maxIter) = beta_mass(1:maxIter)
                     updEner(1:maxIter) = beta_ener(1:maxIter)
                     updMom(1:maxIter,1:ndime) = beta_mom(1:maxIter,1:ndime)
                     !$acc end kernels
                     do jk=ik,1,-1
                        updMass(jk) = updMass(jk)/H_mass(jk,jk)
                        updEner(jk) = updEner(jk)/H_ener(jk,jk)
                        do kk = jk-1,1,-1
                           updMass(kk) = updMass(kk) - H_mass(kk,jk)*updMass(jk)
                           updEner(kk) = updEner(kk) - H_ener(kk,jk)*updEner(jk)
                        end do
                     end do
                     do idime = 1,ndime
                        do jk=ik,1,-1
                           updMom(jk,idime) = updMom(jk,idime)/H_mom(jk,jk,idime)
                           do kk = jk-1,1,-1
                              updMom(kk,idime) = updMom(kk,idime) - H_mom(kk,jk,idime)*updMom(jk,idime)
                           end do
                        end do
                     end do

                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        aux = 0.0_rp
                        do kk = 1,maxIter
                           aux = aux + Q_Mass(lpoin_w(ipoin),kk)*updMass(kk)
                        end do
                        mass_sol(lpoin_w(ipoin)) = mass_sol(lpoin_w(ipoin)) + aux
                        aux = 0.0_rp
                        do kk = 1,maxIter
                           aux = aux + Q_Ener(lpoin_w(ipoin),kk)*updEner(kk)
                        end do
                        ener_sol(lpoin_w(ipoin)) = ener_sol(lpoin_w(ipoin)) + aux
                        do idime = 1,ndime
                           aux = 0.0_rp
                           do kk = 1,maxIter
                              aux = aux + Q_Mom(lpoin_w(ipoin),idime,kk)*updMom(kk,idime)
                           end do
                           mom_sol(lpoin_w(ipoin),idime) = mom_sol(lpoin_w(ipoin),idime) + aux
                        end do
                     end do
                     !$acc end parallel loop
                  end if

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

              subroutine arnoldi_iter(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                      atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                      rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                      gammaRK,dt)
                  implicit none
                  integer(4), intent(in)    :: ik
                  integer(4), intent(in)    :: nelem, npoin, npoin_w, atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1)
                  integer(4), intent(in)    :: connec(nelem,nnode), lpoin_w(npoin_w), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                  real(rp)  , intent(in)    :: gammaRK, dt, gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
                  real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus), He(ndime,ndime,ngaus,nelem), Ml(npoin)
                  real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin),Rgas,gamma_gas, Cp, Prt
                  real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
                  integer(4)                :: idime, jk, ipoin
                  real(rp)                  :: zmass(npoin), zmom(npoin,ndime), zener(npoin)
                  real(8)                  :: aux(5),aux2(5)
                  
                  ! Compute the new J*Q(:,ik) arrays
                  !$acc kernels
                  zmass(:) = pMass(:)
                  zmom (:,:) = pMom(:,:)
                  zener(:) = pEner(:)
                  !$acc end kernels
                  call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                      atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                      rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                      zmass,zmom,zener,.false.)

                  ! Update Jy with the complement
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     Q_Mass(lpoin_w(ipoin),ik+1) = zmass(lpoin_w(ipoin))/(gammaRK*dt) + Jy_mass(lpoin_w(ipoin))
                     Q_Ener(lpoin_w(ipoin),ik+1) = zener(lpoin_w(ipoin))/(gammaRK*dt) + Jy_ener(lpoin_w(ipoin))
                     !$acc loop seq
                     do idime = 1,ndime
                        Q_Mom(lpoin_w(ipoin),idime,ik+1) = zmom(lpoin_w(ipoin),idime)/(gammaRK*dt) + Jy_mom(lpoin_w(ipoin),idime)
                     end do
                  end do
                  !$acc end parallel loop

                  ! Compute the new H matrix
                  ! TODO: make it gpu friendly
                  do jk = 1,ik
                     aux(:) = 0.0_rp
                     do ipoin=1,npoin_w
                        aux(1) = aux(1) + real((Q_mass(lpoin_w(ipoin),jk)*Q_mass(lpoin_w(ipoin),ik+1)),8)
                        aux(2) = aux(2) + real((Q_ener(lpoin_w(ipoin),jk)*Q_ener(lpoin_w(ipoin),ik+1)),8)
                        do idime = 1,ndime
                           aux(idime+2) = aux(idime+2) + real((Q_mom(lpoin_w(ipoin),idime,jk)*Q_mom(lpoin_w(ipoin),idime,ik+1)),8)
                        end do
                     end do

                     call MPI_Allreduce(aux,aux2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                     H_mass(jk,ik) = real(aux2(1),rp)
                     H_ener(jk,ik) = real(aux2(2),rp)
                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        Q_Mass(lpoin_w(ipoin),ik+1) = Q_Mass(lpoin_w(ipoin),ik+1) - H_mass(jk,ik)*Q_Mass(lpoin_w(ipoin),jk)
                        Q_Ener(lpoin_w(ipoin),ik+1) = Q_Ener(lpoin_w(ipoin),ik+1) - H_Ener(jk,ik)*Q_Ener(lpoin_w(ipoin),jk)
                     end do
                     !$acc end parallel loop
                     do idime = 1,ndime
                        H_mom(jk,ik,idime) = real(aux2(idime+2),rp)
                        !$acc parallel loop
                        do ipoin = 1,npoin_w
                           Q_Mom(lpoin_w(ipoin),idime,ik+1) = Q_Mom(lpoin_w(ipoin),idime,ik+1) - H_mom(jk,ik,idime)*Q_Mom(lpoin_w(ipoin),idime,jk)
                        end do
                        !$acc end parallel loop
                     end do
                  end do

                  ! Fill H(ik+1,ik) with the norms of Q(:,ik+1)
                  aux(:) = 0.0_rp
                  do ipoin=1,npoin_w
                     aux(1) = aux(1) + real((Q_mass(lpoin_w(ipoin),ik+1)*Q_mass(lpoin_w(ipoin),ik+1)),8)
                     aux(2) = aux(2) + real((Q_ener(lpoin_w(ipoin),ik+1)*Q_ener(lpoin_w(ipoin),ik+1)),8)
                     do idime = 1,ndime
                        aux(idime+2) = aux(idime+2) + real((Q_mom(lpoin_w(ipoin),idime,ik+1)*Q_mom(lpoin_w(ipoin),idime,ik+1)),8)
                     end do
                  end do

                  call MPI_Allreduce(aux,aux2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                  H_mass(ik+1,ik) = sqrt(real(aux2(1),rp))
                  H_ener(ik+1,ik) = sqrt(real(aux2(2),rp))
                  do idime = 1,ndime
                     H_mom(ik+1,ik,idime) = sqrt(real(aux2(idime+2),rp))
                  end do

                  ! Normalize every Q(:,ik+1)
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     Q_mass(lpoin_w(ipoin),ik+1) = Q_mass(lpoin_w(ipoin),ik+1)/H_mass(ik+1,ik)
                     Q_ener(lpoin_w(ipoin),ik+1) = Q_ener(lpoin_w(ipoin),ik+1)/H_ener(ik+1,ik)
                     !$acc loop seq
                     do idime = 1,ndime
                        Q_mom(lpoin_w(ipoin),idime,ik+1) = Q_mom(lpoin_w(ipoin),idime,ik+1)/H_mom(ik+1,ik,idime)
                     end do
                  end do
                  !$acc end parallel loop
                  
              end subroutine arnoldi_iter

              subroutine apply_givens_rotation(ik)
                  implicit none
                  integer(4), intent(in) :: ik
                  integer(4)             :: jk, idime
                  real(rp)                :: aux

                  ! Modify the ik_th column of H_* matrices
                  do jk = 1,ik-1
                     aux = cs_mass(jk)*H_mass(jk,ik) + sn_mass(jk)*H_mass(jk+1,ik)
                     H_mass(jk+1,ik) = -sn_mass(jk)*H_mass(jk,ik) + cs_mass(jk)*H_mass(jk+1,ik)
                     H_mass(jk,ik) = aux
                     aux = cs_ener(jk)*H_ener(jk,ik) + sn_ener(jk)*H_ener(jk+1,ik)
                     H_ener(jk+1,ik) = -sn_ener(jk)*H_ener(jk,ik) + cs_ener(jk)*H_ener(jk+1,ik)
                     H_ener(jk,ik) = aux
                     do idime = 1,ndime
                        aux = cs_mom(jk,idime)*H_mom(jk,ik,idime) + sn_mom(jk,idime)*H_mom(jk+1,ik,idime)
                        H_mom(jk+1,ik,idime) = -sn_mom(jk,idime)*H_mom(jk,ik,idime) + cs_mom(jk,idime)*H_mom(jk+1,ik,idime)
                        H_mom(jk,ik,idime) = aux
                     end do
                  end do

                  ! Apply the Givens rotation to the ik_th column of H_* matrices
                  call givens_rotation_full(H_mass(ik,ik),H_mass(ik+1,ik), &
                                            H_ener(ik,ik),H_ener(ik+1,ik), &
                                            H_mom(ik,ik,:),H_mom(ik+1,ik,:),ik)

                  ! Eliminate the ik+1_th row of H_* matrices
                  aux = cs_mass(ik)*H_mass(ik,ik) + sn_mass(ik)*H_mass(ik+1,ik)
                  H_mass(ik+1,ik) = -sn_mass(ik)*H_mass(ik,ik) + cs_mass(ik)*H_mass(ik+1,ik)
                  H_mass(ik,ik) = aux
                  aux = cs_ener(ik)*H_ener(ik,ik) + sn_ener(ik)*H_ener(ik+1,ik)
                  H_ener(ik+1,ik) = -sn_ener(ik)*H_ener(ik,ik) + cs_ener(ik)*H_ener(ik+1,ik)
                  H_ener(ik,ik) = aux
                  do idime = 1,ndime
                     aux = cs_mom(ik,idime)*H_mom(ik,ik,idime) + sn_mom(ik,idime)*H_mom(ik+1,ik,idime)
                     H_mom(ik+1,ik,idime) = -sn_mom(ik,idime)*H_mom(ik,ik,idime) + cs_mom(ik,idime)*H_mom(ik+1,ik,idime)
                     H_mom(ik,ik,idime) = aux
                  end do

              end subroutine apply_givens_rotation

              subroutine givens_rotation_full(v1mass,v2mass,v1ener,v2ener,v1mom,v2mom,ik)
                  implicit none
                  integer(4), intent(in) :: ik
                  real(rp)   , intent(in) :: v1mass,v2mass,v1ener,v2ener,v1mom(ndime),v2mom(ndime)
                  integer(4)             :: idime
                  real(rp)                :: tmass,tener,tmom

                  ! Mass ops.
                  if (v2mass .lt. 1e-7) then
                     cs_mass(ik) = 0.0
                     sn_mass(ik) = 1.0
                  else if (abs(v2mass) .gt. abs(v1mass)) then
                     tmass = v1mass/v2mass
                     sn_mass(ik) = 1.0/sqrt(1.0 + tmass**2)
                     cs_mass(ik) = tmass*sn_mass(ik)
                  else
                     tmass = v2mass/v1mass
                     cs_mass(ik) = 1.0/sqrt(1.0 + tmass**2)
                     sn_mass(ik) = tmass*cs_mass(ik)
                  end if

                  ! Ener ops.
                  if (v2ener .lt. 1e-7) then
                     cs_ener(ik) = 0.0
                     sn_ener(ik) = 1.0
                  else if (abs(v2ener) .gt. abs(v1ener)) then
                     tener = v1ener/v2ener
                     sn_ener(ik) = 1.0/sqrt(1.0 + tener**2)
                     cs_ener(ik) = tener*sn_ener(ik)
                  else
                     tener = v2ener/v1ener
                     cs_ener(ik) = 1.0/sqrt(1.0 + tener**2)
                     sn_ener(ik) = tmass*cs_ener(ik)
                  end if

                  ! Mom ops.
                  do idime = 1,ndime
                     if (v2mom(idime) .lt. 1e-7) then
                        cs_mom(ik,idime) = 0.0
                        sn_mom(ik,idime) = 1.0
                     else if (abs(v2mom(idime)) .gt. abs(v1mom(idime))) then
                        tmom = v1mom(idime)/v2mom(idime)
                        sn_mom(ik,idime) = 1.0/sqrt(1.0 + tmom**2)
                        cs_mom(ik,idime) = tmom*sn_mom(ik,idime)
                     else
                        tmom = v2mom(idime)/v1mom(idime)
                        cs_mom(ik,idime) = 1.0/sqrt(1.0 + tmom**2)
                        sn_mom(ik,idime) = tmom*cs_mom(ik,idime)
                     end if
                  end do

              end subroutine givens_rotation_full

              subroutine form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                        rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
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
                  real(rp)  , intent(in) :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin), Rgas,gamma_gas,Cp, Prt
                  real(rp)  , intent(in) :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
                  real(rp)  , intent(in) :: zmass(npoin), zener(npoin), zmom(npoin,ndime)
                  real(rp)               :: zpres(npoin),zu(npoin,ndime),ztemp(npoin),zeint(npoin)
                  real(rp)               :: auxRmass(npoin),auxRmom(npoin,ndime),auxRener(npoin),auxQ(npoin,ndime),auxU(npoin,ndime)
                  integer(4)             :: idime,ipoin
                  real(rp)               :: aux

                  ! Form the R(u^n) arrays if not formed already
                  if (flag_gmres_form_fix .eqv. .true.) then
                     call full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                                          gmshAtoI, gmshAtoJ, gmshAtoK, u, q, rho, pr, E, Rmass, Rmom, Rener)
                      call full_diffusion_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                                              gmshAtoI, gmshAtoJ, gmshAtoK, Cp, Prt, rho, u, Tem, &
                                              mu_fluid, mu_e, mu_sgs, Ml, Dmass, Dmom, Dener)
                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        Rmass_fix(lpoin_w(ipoin)) = Rmass(lpoin_w(ipoin))  + Dmass(lpoin_w(ipoin))
                        Rener_fix(lpoin_w(ipoin)) = Rener(lpoin_w(ipoin))  + Dener(lpoin_w(ipoin))
                        !$acc loop seq
                        do idime = 1,ndime
                           Rmom_fix(lpoin_w(ipoin),idime) = Rmom(lpoin_w(ipoin),idime) + Dmom(lpoin_w(ipoin),idime)
                        end do
                     end do
                     !$acc end parallel loop
                  end if

                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     !$acc loop seq
                     do idime = 1,ndime
                        zu(lpoin_w(ipoin),idime) = zmom(lpoin_w(ipoin),idime)/rho(lpoin_w(ipoin))
                        auxU(lpoin_w(ipoin),idime) = u(lpoin_w(ipoin),idime)+epsQ(idime)*zu(lpoin_w(ipoin),idime)
                        auxQ(lpoin_w(ipoin),idime) = q(lpoin_w(ipoin),idime)+epsQ(idime)*zmom(lpoin_w(ipoin),idime)
                     end do 
                     zeint(lpoin_w(ipoin)) = (zener(lpoin_w(ipoin))/rho(lpoin_w(ipoin)))- &
                        0.5_rp*dot_product(zu(lpoin_w(ipoin),:),zu(lpoin_w(ipoin),:))
                     zpres(lpoin_w(ipoin)) = zmass(lpoin_w(ipoin))*(gamma_gas-1.0_rp)*zeint(lpoin_w(ipoin))
                     ztemp(lpoin_w(ipoin)) = zpres(lpoin_w(ipoin))/(rho(lpoin_w(ipoin))*Rgas)
                  end do
                  !$acc end parallel loop

                  call full_convec_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                     gmshAtoI, gmshAtoJ, gmshAtoK, auxU, auxQ, rho+epsR*zmass, pr+epsE*zpres, E+epsE*zener, Rmass, Rmom, Rener)

                   call full_diffusion_ijk(nelem, npoin, connec, Ngp, dNgp, He, gpvol, dlxigp_ip, xgp, atoIJK, invAtoIJK, &
                      gmshAtoI, gmshAtoJ, gmshAtoK, Cp, Prt, rho+epsR*zmass, auxU, Tem+epsE*ztemp, &
                      mu_fluid, mu_e, mu_sgs, Ml, Dmass, Dmom, Dener)

                  !$acc end parallel loop

                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     Rmass(lpoin_w(ipoin)) = Rmass(lpoin_w(ipoin)) + Dmass(lpoin_w(ipoin))
                     Rener(lpoin_w(ipoin)) = Rener(lpoin_w(ipoin)) + Dener(lpoin_w(ipoin))
                     !$acc loop seq
                     do idime = 1,ndime
                        Rmom(lpoin_w(ipoin),idime) = Rmom(lpoin_w(ipoin),idime) + Dmom(lpoin_w(ipoin),idime)
                     end do
                  end do
                  !$acc end parallel loop

                  ! Form the J*y arrays
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     Jy_mass(lpoin_w(ipoin)) = (Rmass(lpoin_w(ipoin)) - Rmass_fix(lpoin_w(ipoin)))/epsR
                     Jy_ener(lpoin_w(ipoin)) = (Rener(lpoin_w(ipoin)) - Rener_fix(lpoin_w(ipoin)))/epsE
                     !$acc loop seq
                     do idime = 1,ndime
                        Jy_mom(lpoin_w(ipoin),idime) = (Rmom(lpoin_w(ipoin),idime) - Rmom_fix(lpoin_w(ipoin),idime))/epsQ(idime)
                     end do
                  end do
                  !$acc end parallel loop
                  ! Communicate before applying ML
                  if(mpi_size.ge.2) then
                     call nvtxStartRange("MPI_comms_tI")
                     call mpi_halo_atomic_update_real(Jy_mass)
                     call mpi_halo_atomic_update_real(Jy_ener)
                     do idime = 1,ndime
                        call mpi_halo_atomic_update_real(Jy_mom(:,idime)) ! TODO: need to pass the fulll aray?
                     end do
                     call nvtxEndRange
                  end if

                  ! Multiply residuals by inverse Ml
                  call lumped_solver_scal(npoin, npoin_w, lpoin_w, Ml, Jy_mass)
                  call lumped_solver_scal(npoin, npoin_w, lpoin_w, Ml, Jy_ener)
                  call lumped_solver_vect(npoin, npoin_w, lpoin_w, Ml, Jy_mom)

              end subroutine form_approx_Jy

              subroutine init_gmres(npoin,npoin_w,lpoin_w,bmass,bmom,bener,dt,gammaRK)
                  implicit none
                  integer(4), intent(in) :: npoin, npoin_w, lpoin_w(npoin_w)
                  real(rp)  , intent(in) :: bmass(npoin), bmom(npoin,ndime), bener(npoin), dt, gammaRK
                  integer(4)             :: ipoin, idime
                  real(rp)               :: amass(npoin), aener(npoin), amom(npoin,ndime)
                  real(8)                :: aux(5), aux2(5), aux3

                  ! Zero Q_* arrays
                  !$acc kernels
                  Q_Mass(:,:) = 0.0_rp
                  Q_Ener(:,:) = 0.0_rp
                  Q_Mom(:,:,:) = 0.0_rp
                  !$acc end kernels

                  ! Zero cs and sn arrays
                  !$acc kernels
                  cs_mass(:) = 0.0_rp
                  cs_ener(:) = 0.0_rp
                  cs_mom(:,:) = 0.0_rp
                  sn_mass(:) = 0.0_rp
                  sn_ener(:) = 0.0_rp
                  sn_mom(:,:) = 0.0_rp
                  !$acc end kernels

                  ! Zero H_* mmatrices
                  !$acc kernels
                  H_Mass(:,:) = 0.0_rp
                  H_Ener(:,:) = 0.0_rp
                  H_Mom(:,:,:) = 0.0_rp
                  !$acc end kernels

                  ! Initialize e1_* arrays
                  !$acc kernels
                  e1_mass(:) = 0.0_rp
                  e1_ener(:) = 0.0_rp
                  e1_mom(:,:) = 0.0_rp
                  !$acc end kernels
                  e1_mass(1) = 1.0_rp
                  e1_ener(1) = 1.0_rp
                  e1_mom(1,:) = 1.0_rp

                  ! Add the remaining terms too form thhe L*y arrays
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     amass(lpoin_w(ipoin)) = ymass(lpoin_w(ipoin))/(gammaRK*dt) + Jy_mass(lpoin_w(ipoin))
                     aener(lpoin_w(ipoin)) = yener(lpoin_w(ipoin))/(gammaRK*dt) + Jy_ener(lpoin_w(ipoin))
                     !$acc loop seq
                     do idime = 1,ndime
                        amom(lpoin_w(ipoin),idime) = ymom(lpoin_w(ipoin),idime)/(gammaRK*dt) + Jy_mom(lpoin_w(ipoin),idime)
                     end do
                  end do
                  !$acc end parallel loop

                  ! Compute each r = b-Ax
                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     Q_Mass(lpoin_w(ipoin),1) = bmass(lpoin_w(ipoin)) - amass(lpoin_w(ipoin))
                     Q_Ener(lpoin_w(ipoin),1) = bener(lpoin_w(ipoin)) - aener(lpoin_w(ipoin))
                     !$acc loop seq
                     do idime = 1,ndime
                        Q_Mom(lpoin_w(ipoin),idime,1) = bmom(lpoin_w(ipoin),idime) - amom(lpoin_w(ipoin),idime)
                     end do
                  end do
                  !$acc end parallel loop

                  ! Normalize each residual
                  ! TODO: make it GPU compatible without using atomic update
                  !$acc kernels
                  aux(:) = 0.0_rp
                  aux2(:) = 0.0_rp
                  !$acc end kernels
                  aux3 = 0.0_rp
                  !$acc parallel loop reduction(+:aux3)
                  do ipoin = 1,npoin_w
                     aux3 = aux3 + real(Q_Mass(lpoin_w(ipoin),1)**2,8)
                  end do
                  !$acc end parallel loop
                  aux(1) = aux3

                  aux3 = 0.0_rp
                  !$acc parallel loop reduction(+:aux3)
                  do ipoin = 1,npoin_w
                     aux3 = aux3 + real(Q_Ener(lpoin_w(ipoin),1)**2,8)
                  end do
                  aux(2) = aux3
                  !$acc end parallel loop

                  do idime = 1,ndime
                     aux3 = 0.0_rp
                     !$acc parallel loop reduction(+:aux3)
                     do ipoin = 1,npoin_w
                        aux3 = aux3 + real(Q_Mom(lpoin_w(ipoin),idime,1)**2,8)
                     end do
                     !$acc end parallel loop
                     aux(idime+2) = aux3
                  end do

                  call MPI_Allreduce(aux,aux2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                  if(aux2(1)<1e-10) aux2(1) = 1.0_rp
                  if(aux2(2)<1e-10) aux2(2) = 1.0_rp
                  if(aux2(3)<1e-10) aux2(3) = 1.0_rp
                  if(aux2(4)<1e-10) aux2(4) = 1.0_rp
                  if(aux2(5)<1e-10) aux2(5) = 1.0_rp
                  
                  !$acc kernels
                  beta_mass(:) = sqrt(real(aux2(1),rp))*e1_mass(:)
                  beta_ener(:) = sqrt(real(aux2(2),rp))*e1_ener(:)
                  !$acc end kernels

                  !$acc parallel loop
                  do ipoin = 1,npoin_w
                     Q_Mass(lpoin_w(ipoin),1) = Q_Mass(lpoin_w(ipoin),1)/sqrt(real(aux2(1),rp))
                     Q_Ener(lpoin_w(ipoin),1) = Q_Ener(lpoin_w(ipoin),1)/sqrt(real(aux2(2),rp))
                  end do
                  !$acc end parallel loop

                  do idime = 1,ndime
                     !$acc kernels
                     beta_mom(:,idime) = sqrt(real(aux2(idime+2),rp))*e1_mom(:,idime)
                     !$acc end kernels
                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        Q_Mom(lpoin_w(ipoin),idime,1) = Q_Mom(lpoin_w(ipoin),idime,1)/sqrt(real(aux2(idime+2),rp))
                     end do
                     !$acc end parallel loop
                  end do
                  
                  !$acc end parallel loop
                  !print*," beta mass ",beta_mass
                  !print*," beta ener ",beta_ener                  !print*," beta mom ",beta_mom
                  !print*," beta mom ",beta_mom
              end subroutine init_gmres

              subroutine smooth_gmres(ik,nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                                      atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                                      rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                                      gammaRK,dt,iters)
                  implicit none
                  integer(4), intent(in)    :: ik,iters
                  integer(4), intent(in)    :: nelem, npoin, npoin_w, atoIJK(nnode), invAtoIJK(porder+1,porder+1,porder+1)
                  integer(4), intent(in)    :: connec(nelem,nnode), lpoin_w(npoin_w), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                  real(rp)  , intent(in)    :: gammaRK, dt, gpvol(1,ngaus,nelem), dlxigp_ip(ngaus,ndime,porder+1), xgp(ngaus,ndime)
                  real(rp)  , intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus), He(ndime,ndime,ngaus,nelem), Ml(npoin)
                  real(rp)  , intent(in)    :: rho(npoin), u(npoin,ndime), q(npoin,ndime), pr(npoin), E(npoin), Tem(npoin),Rgas,gamma_gas, Cp, Prt
                  real(rp)  , intent(in)    :: mu_fluid(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
                  integer(4)                :: idime, jk, ipoin,i
                  real(rp)                  :: bmass(npoin), bmom(npoin,ndime), bener(npoin),relax=0.2_rp
                  real(8)                  :: auxN(5),auxN2(5)

                  !$acc kernels
                  bmass(:) = pMass(:)
                  bener(:) = pEner(:)
                  bmom(:,:) = pMom(:,:)
                  !$acc end kernels

                  do jk = 1,iters
                     auxN(:) = 0.0_rp
                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        auxN(1) = auxN(1) + real(pMass(lpoin_w(ipoin))**2,8)
                        auxN(2) = auxN(2) + real(pEner(lpoin_w(ipoin))**2,8)
                        !$acc lloop seq
                        do idime = 1,ndime
                           auxN(idime+2) = auxN(idime+2) + real(pMom(lpoin_w(ipoin),idime)**2,8)
                        end do
                     end do
                     !$acc end parallel loop

                     call MPI_Allreduce(auxN,auxN2,5,mpi_datatype_real8,MPI_SUM,MPI_COMM_WORLD,mpi_err)

                     if(auxN2(1)<epsilon(eps)) auxN2(1) = 1.0_rp
                     if(auxN2(2)<epsilon(eps)) auxN2(2) = 1.0_rp
                     if(auxN2(3)<epsilon(eps)) auxN2(3) = 1.0_rp
                     if(auxN2(4)<epsilon(eps)) auxN2(4) = 1.0_rp
                     if(auxN2(5)<epsilon(eps)) auxN2(5) = 1.0_rp

                     epsR = eps/sqrt(real(auxN2(1),rp))
                     epsE = eps/sqrt(real(auxN2(2),rp))
                     epsQ(1) = eps/sqrt(real(auxN2(3),rp))
                     epsQ(2) = eps/sqrt(real(auxN2(4),rp))
                     epsQ(3) = eps/sqrt(real(auxN2(5),rp))
                     call form_approx_Jy(nelem,npoin,npoin_w,lpoin_w,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp, &
                        atoIJK,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK, &
                        rho,u,q,pr,E,Tem,Rgas,gamma_gas,Cp,Prt,mu_fluid,mu_e,mu_sgs,Ml, &
                        pMass,pMom,pEner,.false.)

                     !$acc parallel loop
                     do ipoin = 1,npoin_w
                        pMass(lpoin_w(ipoin)) = pMass(lpoin_w(ipoin))*(1.0_rp-relax) +  relax*gammaRK*dt*(bmass(lpoin_w(ipoin)) - Jy_mass(lpoin_w(ipoin)))
                        pEner(lpoin_w(ipoin)) = pEner(lpoin_w(ipoin))*(1.0_rp-relax) +  relax*gammaRK*dt*(bener(lpoin_w(ipoin)) - Jy_ener(lpoin_w(ipoin)))
                        !$acc loop seq
                        do idime = 1,ndime
                           pMom(lpoin_w(ipoin),idime) = pMom(lpoin_w(ipoin),idime)*(1.0_rp-relax) + relax*gammaRK*dt*(bmom(lpoin_w(ipoin),idime) - Jy_mom(lpoin_w(ipoin),idime))
                        end do
                     end do

                  end do
              end subroutine smooth_gmres
end module mod_solver
