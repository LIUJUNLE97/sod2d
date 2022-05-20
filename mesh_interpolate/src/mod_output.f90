module mod_output

   use mod_constants

   contains

      subroutine write_ascii_fields(npoin,rho,u,pr)

         implicit none

         integer(4), intent(in)  :: npoin
         real(8),    intent(out) :: rho(npoin), u(npoin,ndime), pr(npoin)
         integer(4)              :: ipoin
         integer(4), parameter   :: fileID=667
         character(500)          :: filename

         ! Write velocity vector u to VELOC.alya file
         write(filename,'("VELOC.alya")')
         open(unit=fileID,file=filename,status='unknown')
         do ipoin = 1, npoin
            write(fileID,'(I7,2X,F16.8,2X,F16.8,2X,F16.8)') ipoin, u(ipoin,1), u(ipoin,2), u(ipoin,3)
         end do
         close(fileID)

         ! Write density rho to DENSITY.alya file
         write(filename,'("DENSI.alya")')
         open(unit=fileID,file=filename,status='unknown')
         do ipoin = 1, npoin
            write(fileID,*) ipoin, rho(ipoin)
         end do
         close(fileID)

         ! Write pressure pr to PRESSURE.alya file
         write(filename,'("PRESS.alya")')
         open(unit=fileID,file=filename,status='unknown')
         do ipoin = 1, npoin
            write(fileID,*) ipoin, pr(ipoin)
         end do
         close(fileID)

      end subroutine write_ascii_fields

      subroutine write_vtk_binary_linearized(npoin,nelem,coord,connecLINEAR,connec,rho,u,pr)
         implicit none
      
         integer(4), intent(in)                                              :: npoin, nelem
         integer(4), intent(in)                                              :: connecLINEAR(nelem*(porder**ndime),2**ndime)
         integer(4), intent(in)                                              :: connec(nelem,nnode)
         real(8)   , intent(in)                                              :: coord(npoin,ndime)
         real(8)   , intent(inout), dimension(npoin)                         :: rho, pr
         real(8)   , intent(inout), dimension(npoin,ndime)                   :: u
         integer(4)                                                          :: i,j, iper, ivtk=9, nelem_l, nnode_l
         integer(4)            , dimension(nelem*(porder**ndime),2**ndime+1) :: cells
         integer(4)            , dimension(nelem*(porder**ndime))            :: cellTypes
         real(8)               , dimension(npoin,3)                          :: points, u3d
         character(500)                                                      :: filename
         character(80)                                                       :: buffer
         character(8)                                                        :: str1, str2
         character(1)                                                        :: lf

         lf = achar(10)

         nelem_l = nelem*(porder**ndime)
         nnode_l = 2**ndime 
      
         !
         ! Pass coordinates to a suitable 3D generic format
         !
         !$acc kernels
         points(:,:) = 0.0d0
         points(:,1:ndime) = coord(:,1:ndime)
         !$acc end kernels
      
         !
         ! Pass vector data to a suitable 3D generic format
         !
         !$acc kernels
         u3d(:,:) = 0.0d0
         u3d(:,1:ndime) = u(:,1:ndime)
         !$acc end kernels
      
         !
         ! Pass cell list to VTK format
         !
         !$acc kernels
         cells(:,1) = nnode_l
         cells(:,2:nnode_l+1) = connecLINEAR(:,1:nnode_l)-1
         !cells(:,2:nnode+1) = cells(:,2:nnode+1)-1 ! maybe can be removed?
         !$acc end kernels
      
         !
         ! Define cell types
         !
         if (ndime .eq. 2) then
            cellTypes = 9
         else if (ndime .eq. 3) then
            cellTypes = 12
         end if
      
         !
         ! Open file with ascii input
         !
         write(filename,'("vtkInterp_",i0,".vtk")') istep
         open(unit=ivtk,file=filename,status='replace',access='stream',convert='BIG_ENDIAN') ! Binary file access with stream
         
         !
         ! Write header in ascii format
         !
         write(ivtk) '# vtk DataFile Version 3.0'//lf
         write(ivtk) 'unstr_grid'//lf
         write(ivtk) 'BINARY'//lf
         write(ivtk) 'DATASET UNSTRUCTURED_GRID'//lf//lf
         
         !
         ! Write points
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) 'POINTS '//str1//'  double'//lf
         do i = 1,npoin
            write(ivtk) points(i,:)
         end do
         
         !
         ! Write cells
         !
         write(str1(1:8),'(i8)') nelem_l
         write(str2(1:8),'(i8)') nelem_l*(nnode_l+1)
         write(ivtk) lf//lf//'CELLS '//str1//' '//str2//lf
         do i = 1,nelem_l
            write(ivtk) cells(i,:)
         end do
         
         !
         ! Write cell types
         !
         write(str1(1:8),'(i8)') nelem_l
         write(ivtk) lf//lf//'CELL_TYPES '//str1//lf
         do i = 1,nelem_l
            write(ivtk) cellTypes(i)
         end do
         
         !
         ! Write point scalar data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'POINT_DATA '//str1//lf
         write(ivtk) 'SCALARS DENSI double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) rho(i)
         end do
         write(ivtk) lf//lf//'SCALARS PRESS double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) pr(i)
         end do
         
         !
         ! Write point vector data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS VELOC double'//lf
         do i = 1,npoin
            write(ivtk) u3d(i,:)
         end do
         
         close(ivtk)
      
      end subroutine write_vtk_binary_linearized

      subroutine read_vtk_binary(npoin,nelem,coord,connec,rho,u,pr,E,mu_fluid)
      
         implicit none
      
         integer(4), intent(in)                            :: npoin, nelem
         integer(4), intent(in)                            :: connec(nelem,nnode)
         real(8)   , intent(in)                            :: coord(npoin,ndime)
         real(8)   , intent(inout), dimension(npoin)       :: rho, pr, E, mu_fluid
         real(8)   , intent(inout), dimension(npoin,ndime) :: u
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
