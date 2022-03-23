module mod_time_ops

   use mod_constants
   use mod_nvtx

   contains

      subroutine adapt_dt_cfl(nelem,npoin,connec,helem,u,csound,cfl_conv,dt,cfl_diff,mu_fluid,rho)

         implicit none

         integer(4), intent(in)           :: nelem, npoin, connec(nelem,nnode)
         real(8)   , intent(in)           :: helem(nelem)
         real(8)   , intent(in)           :: u(npoin,ndime), csound(npoin)
         real(8)   , intent(in)           :: cfl_conv
         real(8)   , intent(in), optional :: cfl_diff,mu_fluid(npoin),rho(npoin)
         real(8)  , intent(out)           :: dt

         integer(4)                       :: inode, ielem
         real(8)                          :: dt_conv, dt_diff
         real(8)                          :: umag, aux1, aux2, aux3, aux4, L3, max_MU

         call nvtxStartRange("Adapting dt")
         dt_conv = 100000000000000.0d0
         dt_diff = 100000000000000.0d0
         dt      = 100000000000000.0d0
         do ielem = 1,nelem
            L3 = 0.0d0
            do inode = 1,nnode
               umag = sqrt(dot_product(u(connec(ielem,inode),:),u(connec(ielem,inode),:)))
               aux1 = umag+csound(connec(ielem,inode))
               L3 = max(L3,aux1)
            end do
            aux2 = cfl_conv*(helem(ielem)/dble(porder))/L3
            dt_conv = min(dt_conv,aux2)
            if(present(cfl_diff) .and. present(mu_fluid) .and. present(rho)) then
               aux3 = 0.0d0
            do inode = 1,nnode
               max_MU = max(aux3,mu_fluid(connec(ielem,inode))/rho(connec(ielem,inode)))
               end do
               dt_diff = min(dt_diff,aux4)
            end if
            dt = min(dt_conv,dt_diff)
         end do
         call nvtxEndRange

      end subroutine adapt_dt_cfl

end module mod_time_ops
