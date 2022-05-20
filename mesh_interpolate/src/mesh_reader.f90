module mesh_reader

      use mod_constants
      contains
      
         subroutine read_dims(file_path,file_name,npoin,nelem,nboun)
         
            implicit none
            
            character(500), intent(in)  :: file_path, file_name
            integer(4)    , intent(out) :: npoin, nelem, nboun
            character(500)              :: file_type, line
            
            write(file_type,*) ".dims.dat"
            
            write(*,*) "--| READING DIMS FILE..."
            open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
            read(99,*) line, npoin
            write(*,*) "--| NODES ON MESH : ",npoin
            read(99,*) line, nelem
            write(*,*) "--| ELEMENTS ON MESH : ",nelem
            read(99,*) line, nboun
            write(*,*) "--| BOUNDARY ELEMENTS ON MESH : ",nboun
            write(*,*) "--| END OF DIMS FILE!"
            close(99)
         
         end subroutine read_dims
         
         subroutine read_geo_dat(file_path,file_name,npoin,nelem,nboun,nnodes,npbous,connec,bound,coord)
         
            implicit none
            
            character(500), intent(in)  :: file_path, file_name
            integer(4)    , intent(in)  :: npoin, nelem, nboun, nnodes, npbous
            integer(4)    , intent(out) :: connec(nelem,nnodes), bound(nboun,npbous)
            real(8)       , intent(out) :: coord(npoin,ndime)
            integer(4)                  :: iline, int1, inode, idime, aux(nnodes+1), bou_aux(npbous+1)
            character(2000)             :: file_type, line
            
            write(file_type,*) ".geo.dat"
            
            write(*,*) "--| READING GEO.DAT FILE..."
            open(99,file=trim(adjustl(file_path))//trim(adjustl(file_name))//trim(adjustl(file_type)),status="old")
            !
            ! Nodes/element section, blank for now
            !
            read(99,*) ! Section header
            do iline = 1,nelem
               read(99,*)
            end do
            read(99,*) ! Section ender
            !
            ! Connectivity table section
            !
            write(*,*) "--| READING ELEMENT TABLE..."
            read(99,*) ! Section header
            do iline = 1,nelem
               read(99,'(a)') line
               read(line,*) (aux(inode), inode=1,nnodes+1)
               connec(iline,1:nnodes) = aux(2:nnodes+1)
            end do
            read(99,*) ! Section ender
            !
            ! Nodal coordinates section
            !
            write(*,*) "--| READING COORDINATES..."
            read(99,*) ! Section header
            do iline = 1,npoin
               if (ndime == 2) then
                  read(99,*) int1, coord(iline,1), coord(iline,2)
               else if (ndime == 3) then
                  read(99,*) int1, coord(iline,1), coord(iline,2), coord(iline,3)
               end if
            end do
            read(99,*) ! Section ender
            !
            ! Boundary nodes section
            !
            write(*,*) "--| READING BOUNDARIES..."
            read(99,*) line! Section header
            do iline = 1,nboun
               read(99,'(a)')
            end do
            close(99)
            write(*,*) "--| END OF GEO.DAT FILE!"
         
         end subroutine read_geo_dat

end module mesh_reader
