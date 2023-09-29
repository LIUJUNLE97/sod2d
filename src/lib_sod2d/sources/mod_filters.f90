module mod_filters
   use mod_constants
   use mod_nvtx

   implicit none

   ! integer ---------------------------------------------------
   integer(4), allocatable :: convertIJK(:)

   ! real ------------------------------------------------------
   real(rp), allocatable :: al_weights(:),am_weights(:),an_weights(:)

contains

   subroutine init_filters()
      implicit none
      integer(4) :: ii

      call nvtxStartRange("Init filters vars")

      allocate(al_weights(-1:1))
      !$acc enter data create(al_weights(:))
      al_weights(-1) = 1.0_rp/4.0_rp
      al_weights(0)  = 2.0_rp/4.0_rp
      al_weights(1)  = 1.0_rp/4.0_rp
      !$acc update device(al_weights(:))

      allocate(am_weights(-1:1))
      !$acc enter data create(am_weights(:))
      am_weights(-1) = 1.0_rp/4.0_rp
      am_weights(0)  = 2.0_rp/4.0_rp
      am_weights(1)  = 1.0_rp/4.0_rp
      !$acc update device(am_weights(:))

      allocate(an_weights(-1:1))
      !$acc enter data create(an_weights(:))
      an_weights(-1) = 1.0_rp/4.0_rp
      an_weights(0)  = 2.0_rp/4.0_rp
      an_weights(1)  = 1.0_rp/4.0_rp
      !$acc update device(an_weights(:))
     
      allocate(convertIJK(0:porder+2))
      !$acc enter data create(convertIJK(:))
      do ii=3,porder+1
         convertIJK(ii-1) = ii
      end do
      convertIJK(0) = 3
      convertIJK(1) = 1
      convertIJK(porder+1) = 2
      convertIJK(porder+2) = porder
      !$acc update device(convertIJK(:))

      call nvtxEndRange

   end subroutine init_filters

   subroutine deallocate_filters()

      call nvtxStartRange("Deallocating filters vars")

      !$acc exit data delete(al_weights(:))
      deallocate(al_weights)
      !$acc exit data delete(am_weights(:))
      deallocate(am_weights)
      !$acc exit data delete(an_weights(:))
      deallocate(an_weights)
      !$acc exit data delete(convertIJK(:))    
      deallocate(convertIJK)

      call nvtxEndRange

   end subroutine deallocate_filters



end module mod_filters