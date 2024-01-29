module mod_custom_types
   use mod_constants

   implicit none

   type ptr_array1d_rp
      real(rp),dimension(:),pointer :: ptr
   end type ptr_array1d_rp

   type ptr_array2d_rp
      real(rp),dimension(:,:),pointer :: ptr
   end type ptr_array2d_rp

   type ptr_array1d_rpavg
      real(rp_avg),dimension(:),pointer :: ptr
   end type ptr_array1d_rpavg

   type ptr_array2d_rpavg
      real(rp_avg),dimension(:,:),pointer :: ptr
   end type ptr_array2d_rpavg

   type ptr_array1d_rp_save
      real(rp),dimension(:),pointer :: ptr_rp
      real(rp_avg),dimension(:),pointer :: ptr_avg
      character(128),public :: nameField
   end type ptr_array1d_rp_save

   type ptr_array2d_rp_save
      real(rp),dimension(:,:),pointer :: ptr_rp
      real(rp_avg),dimension(:,:),pointer :: ptr_avg
      character(128),public :: nameField
   end type ptr_array2d_rp_save
!------------------------------------------------------------------------------------------------------

   type vector_int4
      integer(4),dimension(:), allocatable :: elems
   end type vector_int4

   type vector_int8
      integer(8),dimension(:), allocatable :: elems
   end type vector_int8

   type matrix_int4
      integer(4),dimension(:,:), allocatable :: elems
   end type matrix_int4

   type matrix_int8
      integer(8),dimension(:,:), allocatable :: elems
   end type matrix_int8

   type vector_rp
      real(rp),dimension(:), allocatable :: elems
   end type vector_rp

   type matrix_rp
      real(rp),dimension(:,:), allocatable :: elems
   end type matrix_rp

   type vector_real8
      real(8),dimension(:), allocatable :: elems
   end type vector_real8

   type matrix_real8
      real(8),dimension(:,:), allocatable :: elems
   end type matrix_real8

   type jagged_vector_int4
      type(vector_int4),dimension(:), allocatable :: vector
   end type jagged_vector_int4

   type jagged_matrix_int4
      type(matrix_int4),dimension(:), allocatable :: matrix
   end type jagged_matrix_int4

   type jagged_vector_int8
      type(vector_int8),dimension(:), allocatable :: vector
   end type jagged_vector_int8

   type jagged_matrix_int8
      type(matrix_int8),dimension(:), allocatable :: matrix
   end type jagged_matrix_int8

   type jagged_vector_rp
      type(vector_rp),dimension(:), allocatable :: vector
   end type jagged_vector_rp

   type jagged_matrix_rp
      type(matrix_rp),dimension(:), allocatable :: matrix
   end type jagged_matrix_rp

   type jagged_vector_real8
      type(vector_real8),dimension(:), allocatable :: vector
   end type jagged_vector_real8

   type jagged_matrix_real8
      type(matrix_real8),dimension(:), allocatable :: matrix
   end type jagged_matrix_real8

!------------------------------------------------------------------------------------------------------

end module mod_custom_types
