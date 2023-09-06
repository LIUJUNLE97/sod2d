module mod_entropy_viscosity

   use mod_numerical_params
   use mod_nvtx
   
   use mod_mpi
   use mod_mpi_mesh
   use mod_comms

   implicit none

   real(rp),  allocatable :: mue_l(:,:),al(:),am(:),an(:)
   integer(4),allocatable :: convertIJK(:)
   logical  :: allocate_memory_entropyvisc = .true.

   contains
      subroutine smart_visc_spectral(nelem,npoin,npoin_w,connec,lpoin_w,Reta,Rrho,Ngp,coord,dNgp,gpvol,wgp, &
                            gamma_gas,rho,u,csound,Tem,eta,helem,helem_k,Ml,mu_e,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK)

              ! TODO: Compute element size h

              implicit none

              integer(4), intent(in)   :: nelem, npoin,npoin_w, connec(nelem,nnode),lpoin_w(npoin_w)
              real(rp),    intent(in)  :: Reta(npoin), Rrho(npoin), Ngp(ngaus,nnode),gamma_gas
              real(rp),    intent(in)  :: rho(npoin), u(npoin,ndime),csound(npoin), Tem(npoin), eta(npoin),helem(nelem,nnode),helem_k(nelem),Ml(npoin)
              real(rp),    intent(out) :: mu_e(nelem,ngaus)
              real(rp),   intent(in)  :: coord(npoin,ndime), dNgp(ndime,nnode,ngaus), wgp(ngaus)
              real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
              integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1), gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
              integer(4)               :: ielem, inode,igaus,ipoin,npoin_w_g,idime,jdime
              real(rp)                 :: R1, R2, Ve
              real(rp)                 :: betae,mu,vol,vol2
              real(rp)                 :: L3, aux1, aux2, aux3
              real(rp)                 :: maxEta_r,maxEta, maxRho, norm_r,norm, Rgas, maxV, maxC
              real(rp)                :: Je(ndime,ndime), maxJe, minJe,ced,magJe, M, ceM
              integer(4)              :: ii,jj,kk,mm,nn,ll

            if(allocate_memory_entropyvisc) then
               allocate_memory_entropyvisc = .false.
               
               allocate(mue_l(nelem,nnode))
               allocate(al(-1:1))
               allocate(am(-1:1))
               allocate(an(-1:1))
               allocate(convertIJK(0:porder+2))
               !$acc enter data create(mue_l(:,:))
               !$acc enter data create(al(:))
               !$acc enter data create(am(:))
               !$acc enter data create(an(:))
               !$acc enter data create(convertIJK(:))

               do ii=3,porder+1
                  convertIJK(ii-1) = ii
               end do
               convertIJK(0) = 3
               convertIJK(1) = 1
               convertIJK(porder+1) = 2
               convertIJK(porder+2) = porder
               !$acc update device(convertIJK(:))

               al(-1) = 1.0_rp/4.0_rp
               al(0)  = 2.0_rp/4.0_rp
               al(1)  = 1.0_rp/4.0_rp
               !$acc update device(al(:))

               am(-1) = 1.0_rp/4.0_rp
               am(0)  = 2.0_rp/4.0_rp
               am(1)  = 1.0_rp/4.0_rp
               !$acc update device(am(:))

               an(-1) = 1.0_rp/4.0_rp
               an(0)  = 2.0_rp/4.0_rp
               an(1)  = 1.0_rp/4.0_rp
               !$acc update device(an(:))
            end if

             Rgas = nscbc_Rgas_inf

             !$acc kernels
             mue_l(:,:) = mu_e(:,:)
             !$acc end kernels

              if(flag_normalise_entropy .eq. 1) then
                 maxEta_r = 0.0_rp
                 !$acc parallel loop reduction(+:maxEta_r)
                 do ipoin = 1,npoin_w
                    maxEta_r = maxEta_r + eta(lpoin_w(ipoin))
                 end do
                 !$acc end parallel loop

                 call MPI_Allreduce(maxEta_r,maxEta,1,mpi_datatype_real,MPI_SUM,app_comm,mpi_err)
                 call MPI_Allreduce(npoin_w,npoin_w_g,1,mpi_datatype_int,MPI_SUM,app_comm,mpi_err)

                 maxEta = maxEta/real(npoin_w_g,rp)

                 norm_r = 0.0_rp
                 !$acc parallel loop reduction(max:norm_r)
                 do ipoin = 1,npoin_w
                    norm_r = max(norm_r, abs(eta(lpoin_w(ipoin))-maxEta))
                 end do
                 !$acc end parallel loop

                 call MPI_Allreduce(norm_r,norm,1,mpi_datatype_real,MPI_MAX,app_comm,mpi_err)
              else
                 norm = 1.0_rp
              end if

              !$acc parallel loop gang
              do ielem = 1,nelem
                maxJe=0.0_rp
                minJe=1000000.0_rp
                maxV = 0.0_rp
                maxC = 0.0_rp
                !$acc loop seq
                do igaus = 1,ngaus
                   minJe = min(minJe,gpvol(1,igaus,ielem)/wgp(igaus))
                   maxJe = max(maxJe,gpvol(1,igaus,ielem)/wgp(igaus))
                   maxV = max(maxV,sqrt(dot_product(u(connec(ielem,igaus),:),u(connec(ielem,igaus),:))))
                   maxC = max(maxC,csound(connec(ielem,igaus)))
                end do
                M = maxV/maxC
                ced = max(tanh((M**15)*v_pi),ce)
                ceM = max(tanh((M**15)*v_pi),ce)
                ced = max(1.0_rp-(minJe/maxJe)**2,ce)
                ced = max(ced,ceM) 

                mu = 0.0_rp
                betae = 0.0_rp
                !$acc loop vector reduction(max:betae)
                do inode = 1,nnode
                   aux2 = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:))) ! Velocity mag. at element node
                   aux3 = sqrt(Rgas*gamma_gas*Tem(connec(ielem,inode)))     ! Speed of sound at node
                   aux1 = aux2+aux3
                   betae = max(betae,(rho(connec(ielem,inode))*helem_k(ielem))*(cmax/real(porder,rp))*aux1)
                end do
                !$acc loop vector
                do inode = 1,nnode
                   R1 = rho(connec(ielem,inode))*abs(Reta(connec(ielem,inode)))/norm
                   Ve = ced*R1*(helem(ielem,inode))**2
                   mue_l(ielem,inode) = cglob*min(Ve,betae)
                end do
                !$acc loop vector collapse(3)
                do ii=1,porder+1
                   do jj=1,porder+1
                      do kk=1,porder+1
                         aux1 = 0.00_rp
                         !$acc loop seq
                         do ll=-1,1
                            !$acc loop seq
                            do mm=-1,1
                               !$acc loop seq
                               do nn=-1,1
                                  aux1 =   aux1 +  al(ll)*am(mm)*an(nn)*mue_l(ielem,invAtoIJK(convertIJK(ii+ll),convertIJK(jj+mm),convertIJK(kk+nn)))
                               end do
                            end do
                         end do
                         mu_e(ielem,invAtoIJK(convertIJK(ii),convertIJK(jj),convertIJK(kk))) = aux1
                      end do
                   end do
                end do
              end do
              !$acc end parallel loop
      end subroutine smart_visc_spectral
end module mod_entropy_viscosity
