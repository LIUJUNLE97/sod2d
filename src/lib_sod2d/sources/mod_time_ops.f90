module mod_time_ops

   use mod_numerical_params
   
   use mod_nvtx
   use mod_mpi
   use mod_comms

contains

   subroutine adapt_dt_cfl(nelem,npoin,connec,helem,u,csound,cfl_conv,dt,cfl_diff,mu_fluid,mu_sgs,mu_e,rho,Cp,Pr)

      implicit none

      integer(4), intent(in)            :: nelem, npoin, connec(nelem,nnode)
      real(rp)   , intent(in)           :: helem(nelem)
      real(rp)   , intent(in)           :: u(npoin,ndime), csound(npoin)
      real(rp)   , intent(in)           :: cfl_conv
      real(rp)   , intent(in), optional :: cfl_diff,mu_fluid(npoin),mu_sgs(nelem,nnode),mu_e(nelem,nnode),rho(npoin),Cp,Pr
      real(rp)  , intent(inout)         :: dt

      integer(4)                        :: inode,iElem
      real(rp)                          :: dt_l,dt_conv,dt_diff
      real(rp)                          :: umag, aux1, aux2, aux4, L3, max_MU

      call nvtxStartRange("Adapting dt")
      dt_conv = 100000000000000.0_rp
      dt_diff = 100000000000000.0_rp
      dt_l    = 100000000000000.0_rp

      if(flag_use_constant_dt == 0) then
         !$acc parallel loop gang  reduction(min:dt_conv,dt_diff,dt_l)
         do ielem = 1,nelem
            L3 = 0.0_rp
            !$acc loop vector reduction(max:L3)
            do inode = 1,nnode
               umag = abs(u(connec(ielem,inode),1))
               umag = max(umag,abs(u(connec(ielem,inode),2)))
               umag = max(umag,abs(u(connec(ielem,inode),3)))
               aux1 = umag+csound(connec(ielem,inode))
               L3 = max(L3,aux1)
            end do
            aux2 = cfl_conv*(helem(ielem)/real(2.0_rp*porder+1,rp))/L3
            dt_conv = min(dt_conv,aux2)
            if(present(cfl_diff) .and. present(mu_fluid) .and. present(mu_sgs)  .and.  present(rho)) then
               max_MU = 0.0_rp
               !$acc loop vector reduction(max:max_MU)
               do inode = 1,nnode
                  max_MU = max( max_MU, (Cp*rho(connec(ielem,inode))*mu_sgs(ielem,inode)/0.9_rp+mu_fluid(connec(ielem,inode))*(Cp/Pr) + mu_e(ielem,inode)*(Cp/Pr))/rho(connec(ielem,inode)))
               end do
               aux4 = cfl_diff*((helem(ielem)/real(2.0_rp*porder+1,rp))**2)/max_MU
               dt_diff = min(dt_diff,aux4)
            end if
            dt_l = min(dt_conv,dt_diff)
         end do
         !$acc end parallel loop
         call MPI_Allreduce(dt_l,dt,1,mpi_datatype_real,MPI_MIN,app_comm,mpi_err)
      end if
      call nvtxEndRange

      end subroutine adapt_dt_cfl

   subroutine adapt_local_dt_cfl(nelem,npoin,connec,helem,u,csound,cfl_conv,dt,cfl_diff,mu_fluid,mu_sgs,rho)

      implicit none

      integer(4), intent(in)            :: nelem, npoin, connec(nelem,nnode)
      real(rp)   , intent(in)           :: helem(nelem)
      real(rp)   , intent(in)           :: u(npoin,ndime), csound(npoin)
      real(rp)   , intent(in)           :: cfl_conv
      real(rp)   , intent(in), optional :: cfl_diff,mu_fluid(npoin),mu_sgs(nelem,nnode),rho(npoin)
      real(rp)  , intent(out)           :: dt(npoin)

      integer(4)                        :: inode,iElem
      real(rp)                          :: umag, aux1, aux2, aux4, L3, max_MU

      call nvtxStartRange("Adapting dt")

      !$acc kernels
      dt(:) = 100000000000000.0_rp
      !$acc end kernels

      !$acc parallel loop gang
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            umag = abs(u(connec(ielem,inode),1))
            umag = max(umag,abs(u(connec(ielem,inode),2)))
            umag = max(umag,abs(u(connec(ielem,inode),3)))
            aux1 = umag+csound(connec(ielem,inode))
            !$acc atomic update
            dt(connec(ielem,inode)) = min(dt(connec(ielem,inode)) , cfl_conv*(helem(ielem)/real(2.0_rp*porder+1,rp))/aux1)
            !$acc end atomic
         end do
         if(present(cfl_diff) .and. present(mu_fluid) .and. present(mu_sgs)  .and.  present(rho)) then
            !$acc loop vector
            do inode = 1,nnode
               max_MU =(rho(connec(ielem,inode))*mu_sgs(ielem,inode)+mu_fluid(connec(ielem,inode)))/rho(connec(ielem,inode))
               !$acc atomic update
               dt(connec(ielem,inode)) = min(dt(connec(ielem,inode)) , cfl_diff*((helem(ielem)/real(2.0_rp*porder+1,rp))**2)/max_MU)
               !$acc end atomic
            end do
         end if
      end do
      !$acc end parallel loop
      call nvtxEndRange
      if(mpi_size.ge.2) then
         call nvtxStartRange("MPI_comms_tI")
            call mpi_halo_atomic_min_update_real_sendRcv(dt)
         call nvtxEndRange
      end if
      end subroutine adapt_local_dt_cfl

      subroutine adapt_local_dt_cfl_incomp(nelem,npoin,connec,helem,u,cfl_conv,dt,cfl_diff,mu_fluid,mu_sgs,rho)

      implicit none

      integer(4), intent(in)            :: nelem, npoin, connec(nelem,nnode)
      real(rp)   , intent(in)           :: helem(nelem)
      real(rp)   , intent(in)           :: u(npoin,ndime)
      real(rp)   , intent(in)           :: cfl_conv
      real(rp)   , intent(in), optional :: cfl_diff,mu_fluid(npoin),mu_sgs(nelem,nnode),rho(npoin)
      real(rp)  , intent(out)           :: dt(npoin)

      integer(4)                        :: inode,iElem
      real(rp)                          :: umag, aux1, aux2, aux4, L3, max_MU

      call nvtxStartRange("Adapting dt")

      !$acc kernels
      dt(:) = 100000000000000.0_rp
      !$acc end kernels

      !$acc parallel loop gang
      do ielem = 1,nelem
         !$acc loop vector
         do inode = 1,nnode
            umag = dot_product(u(connec(ielem,inode),1:ndime),u(connec(ielem,inode),1:ndime))
            if(umag .gt. epsilon(umag)) then
               !$acc atomic update
               dt(connec(ielem,inode)) = min(dt(connec(ielem,inode)) , cfl_conv*(helem(ielem)/real(2.0_rp*porder+1,rp))/umag)
               !$acc end atomic
            end if
         end do
         if(present(cfl_diff) .and. present(mu_fluid) .and. present(mu_sgs)  .and.  present(rho)) then
            !$acc loop vector
            do inode = 1,nnode
               max_MU =(rho(connec(ielem,inode))*mu_sgs(ielem,inode)+mu_fluid(connec(ielem,inode)))/rho(connec(ielem,inode))
               !$acc atomic update
               dt(connec(ielem,inode)) = min(dt(connec(ielem,inode)) , cfl_diff*((helem(ielem)/real(2.0_rp*porder+1,rp))**2)/max_MU)
               !$acc end atomic
            end do
         end if
      end do
      !$acc end parallel loop
      call nvtxEndRange
      if(mpi_size.ge.2) then
         call nvtxStartRange("MPI_comms_tI")
            call mpi_halo_atomic_min_update_real_sendRcv(dt)
         call nvtxEndRange
      end if
      end subroutine adapt_local_dt_cfl_incomp

subroutine expl_adapt_dt_cfl(nelem,npoin,connec,helem,u,csound,cfl_conv,dt,cfl_diff,mu_fluid,mu_sgs,rho)

      implicit none

      integer(4), intent(in)            :: nelem, npoin, connec(nelem,nnode)
      real(rp)   , intent(in)           :: helem(nelem)
      real(rp)   , intent(in)           :: u(npoin,ndime), csound(npoin)
      real(rp)   , intent(in)           :: cfl_conv
      real(rp)   , intent(in), optional :: cfl_diff,mu_fluid(npoin),mu_sgs(nelem,nnode),rho(npoin)
      real(rp)  , intent(out)           :: dt

      integer(4)                        :: inode,iElem
      real(rp)                          :: dt_l,dt_conv,dt_diff
      real(rp)                          :: umag, aux1, aux2, aux4, L3, max_MU

      call nvtxStartRange("Adapting dt")
      dt_conv = 100000000000000.0_rp
      dt_diff = 100000000000000.0_rp
      dt_l    = 100000000000000.0_rp


      !$acc parallel loop gang  reduction(min:dt_conv,dt_diff,dt_l)
      do ielem = 1,nelem
         L3 = 0.0_rp
         !$acc loop vector reduction(max:L3)
         do inode = 1,nnode
            umag = abs(u(connec(ielem,inode),1))
            umag = max(umag,abs(u(connec(ielem,inode),2)))
            umag = max(umag,abs(u(connec(ielem,inode),3)))
            aux1 = umag+csound(connec(ielem,inode))
            L3 = max(L3,aux1)
         end do
         aux2 = cfl_conv*(helem(ielem)/real(2.0_rp*porder+1,rp))/L3
         dt_conv = min(dt_conv,aux2)
         if(present(cfl_diff) .and. present(mu_fluid) .and. present(mu_sgs)  .and.  present(rho)) then
            max_MU = 0.0_rp
            !$acc loop vector reduction(max:max_MU)
            do inode = 1,nnode
               max_MU = max( max_MU, (rho(connec(ielem,inode))*mu_sgs(ielem,inode)+mu_fluid(connec(ielem,inode)))/rho(connec(ielem,inode)))
            end do
            aux4 = cfl_diff*((helem(ielem)/real(2.0_rp*porder+1,rp))**2)/max_MU
            dt_diff = min(dt_diff,aux4)
         end if
         dt_l = min(dt_conv,dt_diff)
      end do
      !$acc end parallel loop
      call MPI_Allreduce(dt_l,dt,1,mpi_datatype_real,MPI_MIN,app_comm,mpi_err)
      call nvtxEndRange

   end subroutine expl_adapt_dt_cfl

end module mod_time_ops
