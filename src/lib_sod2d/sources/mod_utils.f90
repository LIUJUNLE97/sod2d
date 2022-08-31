module mod_parTimer
    use mpi
    implicit none
    private
    public :: parTimer

    type parTimer
        private
        real(8) :: total_time, wc_start, wc_end, current_time
    contains
        procedure, public :: init_timer
        procedure, public :: start_timer, stop_timer
        procedure, public :: get_totalTime
    end type parTimer

contains

    subroutine init_timer(this)
        class(parTimer), intent(inout) :: this
        this%total_time = 0.d0
        this%wc_start = 0.d0
        this%wc_end = 0.d0
        this%current_time = 0.d0
    end subroutine init_timer

    subroutine start_timer(this)
        class(parTimer), intent(inout) :: this

        this%wc_start = MPI_Wtime()

    end subroutine start_timer

    subroutine stop_timer(this)
        class(parTimer), intent(inout) :: this
        
        this%wc_end = MPI_Wtime()
        this%current_time = this%wc_end - this%wc_start
        this%total_time = this%total_time + this%current_time

    end subroutine stop_timer

    real(8) function get_totalTime(this) result(val)
        class(parTimer), intent(inout) :: this
        val = this%total_time
    end function get_totalTime


end module mod_parTimer

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

subroutine quicksort_array_int(int_array,firstrow,lastrow)
   integer, intent(inout)::int_array(:)
   !integer, intent(in) :: sort_col
   integer, intent(in),optional::firstrow,lastrow
   integer :: temp
   integer :: i,j,left,right,low,high,pivot
   ! If your compiler lacks storage_size(), replace
   ! storage_size(i) by 64
   integer :: stack(2,storage_size(i)),stack_ptr

   if(present(firstrow)) then
      low=firstrow
   else
      low=1
   endif
   if(present(lastrow)) then
      high =lastrow
   else
      high=size(int_array(:))
   end if
   stack_ptr=1

   do
      if (high-low.lt.50) then ! use insertion sort on small int_arrays
         do i=low+1,high
            temp=int_array(i)
            do j=i-1,low,-1
               if (int_array(j).le.temp) exit
               int_array(j+1)=int_array(j)
            enddo
            int_array(j+1)=temp
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
      temp=int_array((low+high)/2)
      int_array((low+high)/2)=int_array(low+1)
      if (temp.gt.int_array(high)) then
         int_array(low+1)=int_array(high)
         int_array(high)=temp
      else
         int_array(low+1)=temp
      endif
      if (int_array(low).gt.int_array(high)) then
         temp=int_array(low)
         int_array(low)=int_array(high)
         int_array(high)=temp
      endif
      if (int_array(low).gt.int_array(low+1)) then
         temp=int_array(low)
         int_array(low)=int_array(low+1)
         int_array(low+1)=temp
      endif
      pivot=int_array(low+1)

      left=low+2
      right=high-1
      do
         do while(int_array(left).lt.pivot)
            left=left+1
         enddo
         do while(int_array(right).gt.pivot)
            right=right-1
         enddo
         if (left.ge.right) exit
         temp=int_array(left)
         int_array(left)=int_array(right)
         int_array(right)=temp
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
end subroutine quicksort_array_int


!########################################################################



end module mod_utils