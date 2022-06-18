module mod_time_ops

   use mod_constants
   use mod_veclen
   use mod_nvtx

   contains

      subroutine adapt_dt_cfl(nelem,npoin,connec,helem,u,csound,cfl_conv,dt,cfl_diff,mu_fluid,mu_sgs,rho)

         implicit none

         integer(4), intent(in)           :: nelem, npoin, connec(nelem,nnode)
         real(8)   , intent(in)           :: helem(nelem)
         real(8)   , intent(in)           :: u(npoin,ndime), csound(npoin)
         real(8)   , intent(in)           :: cfl_conv
         real(8)   , intent(in), optional :: cfl_diff,mu_fluid(npoin),mu_sgs(nelem,nnode),rho(npoin)
         real(8)  , intent(out)           :: dt

         integer(4)                       :: inode, ielem
         real(8)                          :: dt_conv, dt_diff
         real(8)                          :: umag, aux1, aux2, aux4, L3, max_MU

         call nvtxStartRange("Adapting dt")
         dt_conv = 100000000000000.0d0
         dt_diff = 100000000000000.0d0
         dt      = 100000000000000.0d0
         !$acc parallel loop gang  reduction(min:dt_conv,dt_diff,dt) vector_length(vecLength)
         do ielem = 1,nelem
            L3 = 0.0d0
            !$acc loop vector reduction(max:L3)
            do inode = 1,nnode
               umag = abs(u(connec(ielem,inode),1))
               umag = max(umag,abs(u(connec(ielem,inode),2)))
               umag = max(umag,abs(u(connec(ielem,inode),3)))
               aux1 = umag+csound(connec(ielem,inode))
               L3 = max(L3,aux1)
            end do
            aux2 = cfl_conv*(helem(ielem)/dble(2*porder+1))/L3
            dt_conv = min(dt_conv,aux2)
            if(present(cfl_diff) .and. present(mu_fluid) .and. present(mu_sgs)  .and.  present(rho)) then
               max_MU = 0.0d0
               !$acc loop vector reduction(max:max_MU)
               do inode = 1,nnode
                  max_MU = max(max_MU,mu_sgs(ielem,inode)+mu_fluid(connec(ielem,inode))/rho(connec(ielem,inode)))
               end do
               aux4 = cfl_diff*((helem(ielem)/dble(2*porder+1))**2)/max_MU
               dt_diff = min(dt_diff,aux4)
            end if
            dt = min(dt_conv,dt_diff)
         end do
         !$acc end parallel loop
         call nvtxEndRange

      end subroutine adapt_dt_cfl

end module mod_time_ops
