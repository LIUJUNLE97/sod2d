module mod_time_ops

   use mod_constants
   use mod_veclen
   use mod_nvtx
   use mod_mpi

   contains

      subroutine adapt_dt_cfl(nelem,npoin,connec,helem,u,csound,cfl_conv,dt,cfl_diff,mu_fluid,mu_sgs,rho)

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
            !aux2 = cfl_conv*(helem(ielem))/L3
            aux2 = cfl_conv*(helem(ielem)/real(2.0_rp*porder+1,rp))/L3
           ! aux2 = cfl_conv*(helem(ielem)/real(porder**2,rp))/L3
            dt_conv = min(dt_conv,aux2)
            if(present(cfl_diff) .and. present(mu_fluid) .and. present(mu_sgs)  .and.  present(rho)) then
               max_MU = 0.0_rp
               !$acc loop vector reduction(max:max_MU)
               do inode = 1,nnode
                  max_MU = max( max_MU, (rho(connec(ielem,inode))*mu_sgs(ielem,inode)+mu_fluid(connec(ielem,inode)))/rho(connec(ielem,inode)))
               end do
               !aux4 = cfl_diff*((helem(ielem))**2)/max_MU
               !aux4 = cfl_diff*((helem(ielem)/real(porder**2,rp))**2)/max_MU
               aux4 = cfl_diff*((helem(ielem)/real(2.0_rp*porder+1,rp))**2)/max_MU
               dt_diff = min(dt_diff,aux4)
            end if
            dt_l = min(dt_conv,dt_diff)
         end do
         !$acc end parallel loop

         call MPI_Allreduce(dt_l,dt,1,mpi_datatype_real,MPI_MIN,MPI_COMM_WORLD,mpi_err)

         call nvtxEndRange

      end subroutine adapt_dt_cfl

end module mod_time_ops
