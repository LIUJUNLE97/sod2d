module mod_utils
    implicit none

contains

!########################################################################
subroutine quicksort_matrix_int(int_matrix,sort_col,firstrow,lastrow)
   integer, intent(inout)::int_matrix(:,:)
   integer, intent(in) :: sort_col
   integer, intent(in),optional::firstrow,lastrow
   integer, dimension(:), allocatable :: temp
   integer :: i,j,left,right,low,high,pivot
   ! If your compiler lacks storage_size(), replace
   ! storage_size(i) by 64
   integer :: stack(2,storage_size(i)),stack_ptr,num_rows

   num_rows=size(int_matrix(1,:))
   allocate(temp(num_rows))

   if(present(firstrow)) then
      low=firstrow
   else
      low=1
   endif
   if(present(lastrow)) then
      high =lastrow
   else
      high=size(int_matrix(:,sort_col))
   end if
   stack_ptr=1

   do

      if (high-low.lt.50) then ! use insertion sort on small int_arrays
         do i=low+1,high
            temp=int_matrix(i,:)
            do j=i-1,low,-1
               if (int_matrix(j,sort_col).le.temp(sort_col)) exit
               int_matrix(j+1,:)=int_matrix(j,:)
            enddo
            int_matrix(j+1,:)=temp(:)
         enddo
         ! now pop from stack
         if (stack_ptr.eq.1) return
         stack_ptr=stack_ptr-1
         low=stack(1,stack_ptr)
         high=stack(2,stack_ptr)
         cycle
      endif

      ! find median of three pivot
      ! and place sentinels at first and last elements
      temp(:)=int_matrix((low+high)/2,:)
      int_matrix((low+high)/2,:)=int_matrix(low+1,:)
      if (temp(sort_col).gt.int_matrix(high,sort_col)) then
         int_matrix(low+1,:)=int_matrix(high,:)
         int_matrix(high,:)=temp(:)
      else
         int_matrix(low+1,:)=temp(:)
      endif
      if (int_matrix(low,sort_col).gt.int_matrix(high,sort_col)) then
         temp(:)=int_matrix(low,:)
         int_matrix(low,:)=int_matrix(high,:)
         int_matrix(high,:)=temp(:)
      endif
      if (int_matrix(low,sort_col).gt.int_matrix(low+1,sort_col)) then
         temp(:)=int_matrix(low,:)
         int_matrix(low,:)=int_matrix(low+1,:)
         int_matrix(low+1,:)=temp(:)
      endif
      pivot=int_matrix(low+1,sort_col)

      left=low+2
      right=high-1
      do
         do while(int_matrix(left,sort_col).lt.pivot)
            left=left+1
         enddo
         do while(int_matrix(right,sort_col).gt.pivot)
            right=right-1
         enddo
         if (left.ge.right) exit
         temp(:)=int_matrix(left,:)
         int_matrix(left,:)=int_matrix(right,:)
         int_matrix(right,:)=temp(:)
         left=left+1
         right=right-1
      enddo
      if (left.eq.right) left=left+1
      !          call quicksort(int_array(1:left-1))
      !          call quicksort(int_array(left:))
      if (left.lt.(low+high)/2) then
         stack(1,stack_ptr)=left
         stack(2,stack_ptr)=high
         stack_ptr=stack_ptr+1
         high=left-1
      else
         stack(1,stack_ptr)=low
         stack(2,stack_ptr)=left-1
         stack_ptr=stack_ptr+1
         low=left
      endif

   enddo
end subroutine quicksort_matrix_int
!########################################################################



end module mod_utils