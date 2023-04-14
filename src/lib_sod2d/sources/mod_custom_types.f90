module mod_custom_types
   use mod_constants

   implicit none

   type ptr_array1d_rp
      real(rp),dimension(:),pointer :: ptr
   end type ptr_array1d_rp

   type ptr_array2d_rp
      real(rp),dimension(:,:),pointer :: ptr
   end type ptr_array2d_rp

   !type domainptr
      !  type(ptr_array), pointer :: p
   !end type domainptr

end module mod_custom_types