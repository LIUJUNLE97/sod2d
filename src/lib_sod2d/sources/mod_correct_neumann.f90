module mod_correct_neumann

   use mod_numerical_params
   use mod_veclen
   use mod_nvtx
   use mod_mpi

contains

   subroutine evalCorrectNeumann(npoin,nboun,bound,wgp_b,bounorm,u_buffer_flux,Rdiff)

      implicit none

      integer(4), intent(in)     :: npoin,nboun,bound(nboun,npbou)
      real(rp),   intent(in)     :: wgp_b(npbou), bounorm(nboun,ndime*npbou)
      real(rp),   intent(in)     :: u_buffer_flux(npoin,ndime)
      real(rp),   intent(inout)  :: Rdiff(npoin,ndime)
      integer(4)              :: iBound,igaus,idime
      real(rp)                :: aux(ndime),auxmag
      real(rp)                :: bnorm(npbou*ndime)

      !$acc parallel loop gang  private(bnorm)
      do iBound = 1,nboun
         bnorm(1:npbou*ndime) = bounorm(iBound,1:npbou*ndime)
         !$acc loop vector
         do igaus = 1,npbou
            auxmag = sqrt(bnorm((igaus-1)*ndime+1)**2 + bnorm((igaus-1)*ndime+2)**2 + bnorm((igaus-1)*ndime+3)**2)
            !$acc loop seq
            do idime = 1,ndime
               !$acc atomic update
               Rdiff(bound(iBound,igaus),idime) = Rdiff(bound(iBound,igaus),idime)-auxmag*u_buffer_flux(bound(iBound,igaus),idime)*wgp_b(igaus)
               !$acc end atomic
            end do
         end do
      end do
      !$acc end parallel loop

   end subroutine evalCorrectNeumann

end module mod_correct_neumann
