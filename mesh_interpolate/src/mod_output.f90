module mod_output

   use mod_constants

   contains

      subroutine read_vtk_binary(npoin,nelem,coord,connec, &
                                 rho,u,pr,E,mu_fluid,mu_e,mu_sgs)
         implicit none
      
         integer(4), intent(in)                            :: npoin, nelem
         integer(4), intent(in)                            :: connec(nelem,nnode)
         real(8)   , intent(in)                            :: coord(npoin,ndime)
         real(8)   , intent(inout), dimension(npoin)       :: rho, pr, E, mu_fluid
         real(8)   , intent(inout), dimension(npoin,ndime) :: u
         real(8)   , intent(inout), dimension(nelem,ngaus) :: mu_e
         real(8)   , intent(inout), dimension(nelem,ngaus) :: mu_sgs
         integer(4)                                        :: i, iper, ivtk=99
         integer(4)            , dimension(nelem,nnode+1)  :: cells
         integer(4)            , dimension(nelem)          :: cellTypes
         real(8)               , dimension(npoin,3)        :: points
         character(500)                                    :: filename

         write(*,*) " Begining the reading of the vtk binary "
         call flush(1)

         !
         ! Open file with ascii input
         !
         write(filename,'("vtkTstep_",i0,".vtk")') istep
         open(unit=ivtk,file=filename,status='old',access='stream',convert='BIG_ENDIAN',action='read') ! Binary file access with stream
         
         !
         ! Read header in ascii format
         !
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         
         !
         ! Read points
         !
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) points(i,:)
         end do
         
         !
         ! Read cells
         !
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,nelem
            read(ivtk) cells(i,:)
         end do
         
         !
         ! Read cell types
         !
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,nelem
            read(ivtk) cellTypes(i)
         end do
         
         !
         ! Read point scalar data
         !
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) rho(i)
         end do
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) pr(i)
         end do
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) E(i)
         end do
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) mu_fluid(i)
         end do
         !
         ! Write point vector data
         !
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) u(i,:)
         end do

         write(*,*) " end the reading of the vtk binary "
         call flush(1)
         
         close(ivtk)
      
      end subroutine read_vtk_binary

      subroutine skip_line(ivtk)
         implicit none
      
         integer(4), intent(in)                            :: ivtk
         character(1)                                      :: lf, aux

         lf = achar(10)
         do while(.true.)
            read(ivtk) aux
            if(aux .eq. lf) exit
         end do
      end subroutine skip_line
end module mod_output
