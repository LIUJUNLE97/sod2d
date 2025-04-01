!-----------------------------------------------------------------------------------------
! # MODULE FOR I/O UTILS
!-----------------------------------------------------------------------------------------
module mod_ioutils
   use mod_mpi
   use json_module
   implicit none

contains

   !-------------------------------------------------------------
   !    JSON
   !-------------------------------------------------------------
   subroutine trim_extension_json_file(json_filename)
      implicit none
      character(len=*),intent(inout) :: json_filename
     
      ! Append the extension
      json_filename = trim(json_filename)
      if (json_filename == "") then 
         if(mpi_rank.eq.0) write(*,*) "Error! No configuration JSON file given on command line!"
         stop 1
      end if

      if(len_trim(json_filename) < 4 .or. json_filename(len_trim(json_filename) - 4:) /= ".json") then
         json_filename = trim(json_filename) // ".json"
      end if
      !if(mpi_rank.eq.0) write(*,*), "Reading the JSON file : ", json_filename
   end subroutine

   subroutine open_json_file(json_filename,json_f)
      implicit none
      character(len=*),intent(inout) :: json_filename
      type(json_file),intent(inout) :: json_f
  
      call trim_extension_json_file(json_filename)
  
      call json_f%initialize()
      call json_f%load_file(json_filename)
      if(json_f%failed()) then 
        if(mpi_rank.eq.0) then
            write(*,*) "Error loading json file: ",trim(adjustl(json_filename))
            write(*,*) "Aborting!"
            call MPI_Abort(app_comm,-1,mpi_err)
        end if
      end if
   end subroutine

   subroutine close_json_file(json_f)
      implicit none
      type(json_file),intent(inout) :: json_f
  
       call json_f%destroy()
   end subroutine

   !-------------------------------------------------------------
   !    TXT FILE
   !-------------------------------------------------------------

   elemental subroutine str2int(str,int,stat)
      implicit none
      character(len=*),intent(in) :: str
      integer,intent(out)         :: int
      integer,intent(out)         :: stat

      read(str,*,iostat=stat)  int
   end subroutine str2int

   subroutine open_inputFile(path_input_file)
      implicit none
      character(len=*),intent(in) :: path_input_file

      if(mpi_rank.eq.0) write(*,*) 'Reading input file: ',trim(adjustl(path_input_file))
      open(unit=99,file=path_input_file,status="old")
   end subroutine

   subroutine close_inputFile()
      implicit none

      close(unit=99)
   end subroutine

   subroutine read_inputFile_integer(numLine,parameter2read,int2read)
      implicit none
      integer,intent(inout) :: numLine
      character(256),intent(in) :: parameter2read
      integer,intent(inout) :: int2read
      character(256) :: read_sec,read_val
      integer :: aux_int,stat

      read(99,*) read_sec,read_val! Section header
      if(read_sec.eq.trim(adjustl(parameter2read))) then
          call str2int(read_val,aux_int,stat)
          if(stat.eq.0) then
              int2read = aux_int
          else
              if(mpi_rank.eq.0) then
                  write(*,*) 'ERROR!',trim(adjustl(parameter2read)),' must be an integer',&
                             ' | wrong value ',int2read
                  call MPI_Abort(app_comm,-1,mpi_err)
              end if
          end if
          if(mpi_rank.eq.0) write(*,*) trim(adjustl(parameter2read)),': ',int2read
      else
          if(mpi_rank.eq.0) write(*,*) 'Error! Line ',numLine,' must be',trim(adjustl(parameter2read)),' value'
      endif

      numLine=numLine+1

   end subroutine read_inputFile_integer

   subroutine read_inputFile_string(numLine,parameter2read,string2read)
      implicit none
      integer,intent(inout) :: numLine
      character(256),intent(in) :: parameter2read
      character(512),intent(inout) :: string2read
      character(256) :: read_sec,read_val
      integer :: aux_int,stat

        read(99,*) read_sec,read_val! Section header
        if(read_sec.eq.trim(adjustl(parameter2read))) then
            string2read = trim(adjustl(read_val))
            if(mpi_rank.eq.0) write(*,*) trim(adjustl(parameter2read)),': ',trim(adjustl(read_val))
        else
            if(mpi_rank.eq.0) write(*,*) 'Error! Line ',numLine,' must be ',trim(adjustl(parameter2read)),' value'
            call MPI_Abort(app_comm,-1,mpi_err)
        endif

      numLine=numLine+1

   end subroutine read_inputFile_string

   subroutine read_inputFile_logical(numLine,parameter2read,logical2read)
      implicit none
      integer,intent(inout) :: numLine
      character(256),intent(in) :: parameter2read
      logical,intent(inout) :: logical2read
      character(256) :: read_sec,read_val
      integer :: aux_int,stat

      read(99,*) read_sec,read_val! Section header
      if(read_sec.eq.trim(adjustl(parameter2read))) then
          call str2int(read_val,aux_int,stat)
          if(stat.eq.0) then
              if(aux_int.eq.0) then
                  logical2read = .false.
              elseif(aux_int .eq. 1) then
                  logical2read = .true.
              else
                  if(mpi_rank.eq.0) then
                      write(*,*) 'ERROR!',trim(adjustl(parameter2read)),'must be 0(false)/1(true)',&
                                 ' | wrong value ',trim(adjustl(read_val))
                      call MPI_Abort(app_comm,-1,mpi_err)
                  end if
              end if
          else
              if(mpi_rank.eq.0) then
                  write(*,*) 'ERROR!',trim(adjustl(parameter2read)),'must be 0(false)/1(true)',&
                             ' | wrong value ',trim(adjustl(read_val))
                  call MPI_Abort(app_comm,-1,mpi_err)
              end if
          end if
          if(mpi_rank.eq.0) write(*,*) trim(adjustl(parameter2read)),': ',logical2read
      else
          if(mpi_rank.eq.0) write(*,*) 'Error! Input File Line ',numLine,' must be ',trim(adjustl(parameter2read)),' 0/1'
      endif

      numLine=numLine+1

   end subroutine read_inputFile_logical

end module mod_ioutils