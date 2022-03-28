module mod_fluid_viscosity

   use mod_constants

   contains

   subroutine constant_viscosity(npoin,val,mu_fluid)

      implicit none

      integer(4), intent(in)  :: npoin
      real(8)   , intent(in)  :: val
      real(8)   , intent(inout) :: mu_fluid(npoin)

      !$acc kernels
      mu_fluid(:) = val
      !$acc end kernels

   end subroutine constant_viscosity

   subroutine sutherland_viscosity(npoin,Tem,mu_fluid)

      implicit none
      
      integer(4), intent(in)    :: npoin
      real(8),    intent(in)    :: Tem(npoin)
      real(8),    intent(inout) :: mu_fluid(npoin)
      integer(4)                :: ipoin
      real(8)                   :: C1, S

      C1 = 0.000001458d0
      S = 110.4d0

      !$acc parallel loop
      do ipoin = 1,npoin
         mu_fluid(ipoin) = flag_mu_factor*C1*(Tem(ipoin)**1.5d0)/(Tem(ipoin)+S)
      end do
      !$acc end parallel loop

   end subroutine sutherland_viscosity

end module mod_fluid_viscosity
