module elem_source

   use mod_numerical_params
   use mod_nvtx
   
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5
   use mod_comms

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Computes source term integration                                           !
   ! Added to the rhs                                                           !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

   ! integrates a constant source term (s[ndime]) for each cartessian
   ! direction in the momentum equations 
   subroutine mom_source_const_vect(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,s,Rmom,fact)


      implicit none

      integer(4), intent(in)    :: nelem, npoin
      integer(4), intent(in)    :: connec(nelem,nnode)
      real(rp),    intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)    :: He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)    :: gpvol(1,ngaus,nelem)
      real(rp),    intent(in)    :: s(npoin,ndime), u(npoin,ndime)
      real(rp),    intent(inout) :: Rmom(npoin,ndime)
      real(rp), optional, intent(in)  :: fact
      integer(4)                :: ielem, igaus, idime, inode
      real(rp)                   :: Re(nnode,ndime)
      real(rp)  :: aux_fact = 1.0_rp

      call nvtxStartRange("Momentum source term")

      !oriol: I will assue that you will call
      !this subroutine at least having convection so Rmom is
      !already initialized

      

      if(present(fact)) then
         aux_fact = fact
      end if

      !$acc parallel loop gang private(Re) 
      do ielem = 1,nelem
         !$acc loop vector collapse(2)
         do idime = 1,ndime
            do inode = 1,nnode
               !$acc atomic update
               Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)-aux_fact*gpvol(1,inode,ielem)*s(connec(ielem,inode),idime)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

   end subroutine mom_source_const_vect

   subroutine ener_source(nelem,npoin,connec,Ngp,dNgp,He,gpvol,q,Rener,fact)


      implicit none

      integer(4), intent(in)    :: nelem, npoin
      integer(4), intent(in)    :: connec(nelem,nnode)
      real(rp),    intent(in)    :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
      real(rp),    intent(in)    :: He(ndime,ndime,ngaus,nelem)
      real(rp),    intent(in)    :: gpvol(1,ngaus,nelem)
      real(rp),    intent(in)    :: q(npoin,ndime)
      real(rp),    intent(inout) :: Rener(npoin)
      real(rp), optional, intent(in)  :: fact
      integer(4)                :: ielem, igaus, idime, inode
      real(rp)  :: aux_fact = 1.0_rp

      call nvtxStartRange("Momentum source term")

      !oriol: I will assue that you will call
      !this subroutine at least having convection so Rmom is
      !already initialized

      if(present(fact)) then
         aux_fact = fact
      end if

      !$acc parallel loop gang 
      do ielem = 1,nelem
         !$acc loop vector 
        do inode = 1,nnode
            !$acc atomic update
            Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+aux_fact*gpvol(1,inode,ielem)*q(connec(ielem,inode),3)*nscbc_g
            !$acc end atomic
         end do
      end do
      !$acc end parallel loop
      call nvtxEndRange

   end subroutine ener_source

end module elem_source
