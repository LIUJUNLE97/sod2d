module mod_output

   use mod_constants
   use mod_mpi
   use mod_mpi_mesh
   use mod_hdf5

   contains

      subroutine write_vtk_ascii(istep,npoin,nelem,coord,connec, &
                                 rho,u,pr,E,mu_e, mu_sgs)
         implicit none
      
         integer(4), intent(in)                           :: istep, npoin, nelem
         integer(4), intent(in)                           :: connec(nelem,nnode)
         real(rp)   , intent(in)                           :: coord(npoin,ndime)
         real(rp)   , intent(in), dimension(npoin)         :: rho, pr, E
         real(rp)   , intent(in), dimension(npoin,ndime)   :: u
         real(rp)   , intent(in), dimension(nelem,ngaus)   :: mu_e
         real(rp)   , intent(in), dimension(nelem,ngaus)   :: mu_sgs
         integer(4)                                       :: i, ivtk=9
         integer(4)            , dimension(nelem,nnode+1) :: cells
         integer(4)            , dimension(nelem)         :: cellTypes
         real(rp)               , dimension(npoin,3)       :: points, u3d
         character(500)                                   :: filename
         character(16)                                   :: str1, str2
      
         !
         ! Pass coordinates to a suitable 3D generic format
         !
         points = 0.0_rp
         points(:,1:ndime) = coord(:,1:ndime)
      
         !
         ! Pass vector data to a suitable 3D generic format
         !
         u3d = 0.0_rp
         u3d(:,1:ndime) = u(:,1:ndime)
      
         !
         ! Pass cell list to VTK format
         !
         cells(:,1) = nnode
         cells(:,2:nnode+1) = connec(:,1:nnode)
         cells(:,2:nnode+1) = cells(:,2:nnode+1)-1
      
         !
         ! Define cell types
         !
         if (ndime .eq. 2) then
            if (nnode .eq. 4) then ! QUA04
               cellTypes = 9
            end if
         else if (ndime .eq. 3) then
            if (nnode .eq. 8) then ! HEX08
               cellTypes = 12
            else if (nnode .eq. 27) then ! HEX27
               cellTypes = 29
            end if
         end if
      
         !
         ! Open file with ascii input
         !
         call set_vtk_filename(filename,'results',istep)
         !write(filename,'("vtkTstep_",i0,".vtk")') istep
         open(unit=ivtk,file=filename,status='replace') ! Binary file access with stream
         
         !
         ! Write header in ascii format
         !
         write(ivtk,'(a)') '# vtk DataFile Version 3.0'
         write(ivtk,'(a)') 'unstr_grid'
         write(ivtk,'(a)') 'ASCII' 
         write(ivtk,'(a)') 'DATASET UNSTRUCTURED_GRID'
         
         !
         ! Write points
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk,'(a)') 'POINTS '//trim(str1)//'  double'
         do i = 1,npoin
            write(ivtk,*) points(i,:)
         end do
         
         !
         ! Write cells
         !
         write(str1(1:8),'(i8)') nelem
         write(str2(1:8),'(i8)') nelem*(nnode+1)
         write(ivtk,'(a)') 'CELLS '//trim(str1)//trim(str2)
         do i = 1,nelem
            write(ivtk,*) cells(i,:)
         end do
         
         
         ! Write cell types
         !
         write(str1(1:8),'(i8)') nelem
         write(ivtk,'(a)') 'CELL_TYPES '//trim(str1)
         do i = 1,nelem
            write(ivtk,*) cellTypes(i)
         end do
         
         !
         ! Write point scalar data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk,'(a)') 'POINT_DATA '//trim(str1)
         write(ivtk,'(a)') 'SCALARS '//' DENSI '//' double '//' 1'
         write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
         do i = 1,npoin
            write(ivtk,*) rho(i)
         end do
         write(ivtk,'(a)') 'SCALARS '//' PRESS '//' double '//' 1'
         write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
         do i = 1,npoin
            write(ivtk,*) pr(i)
         end do
         write(ivtk,'(a)') 'SCALARS '//' TENER '//' double '//' 1'
         write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
         do i = 1,npoin
            write(ivtk,*) E(i)
         end do
         
         !
         ! Write point vector data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk,'(a)') 'VECTORS '//' VELOC '//' double'
         do i = 1,npoin
            write(ivtk,*) u3d(i,1:3)
         end do
         
         !!
         !! Write cell scalar data
         !!
         write(str1(1:8),'(i8)') nelem
         write(ivtk,'(a)') 'CELL_DATA '//trim(str1)
         write(ivtk,'(a)') 'SCALARS '//' ENVIT '//' double'
         write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
         do i = 1,nelem
            write(ivtk,*) mu_e(i,1)
         end do
         write(ivtk,'(a)') 'SCALARS '//' SGSVI '//' double'
         write(ivtk,'(a)') 'LOOKUP_TABLE '//' default'
         do i = 1,nelem
            write(ivtk,*) mu_sgs(i,1)
         end do
         
         close(ivtk)
      
      end subroutine

      subroutine write_vtk_binary(isPeriodic,istep,npoin,nelem,coord,connec, &
                                 rho,u,pr,E,csound,machno, &
                                 gradRho,curlU,divU,Qcrit,mu_fluid,mu_e,mu_sgs,nper,masSla)
         implicit none
      
         logical, intent(in)                               :: isPeriodic
         integer(4), intent(in)                            :: nper, istep, npoin, nelem
         integer(4), intent(in)                            :: connec(nelem,nnode)
         integer(4), intent(in), optional                  :: masSla(nper,2)
         real(rp)   , intent(in)                            :: coord(npoin,ndime)
         real(rp)   , intent(inout), dimension(npoin)       :: rho, pr, E, csound, machno, mu_fluid, divU, Qcrit
         real(rp)   , intent(inout), dimension(npoin,ndime) :: u
         real(rp)   , intent(inout), dimension(nelem,ngaus) :: mu_e
         real(rp)   , intent(inout), dimension(nelem,ngaus) :: mu_sgs
         real(rp)   , intent(inout), dimension(npoin,ndime) :: gradRho, curlU
         integer(4)                                        :: i, j, iper, ivtk=9
         integer(4)            , dimension(nelem,nnode+1)  :: cells
         integer(4)            , dimension(nelem)          :: cellTypes
         real(rp)               , dimension(npoin)          :: envit, mut
         real(rp)               , dimension(npoin,3)        :: points, u3d
         character(500)                                    :: filename
         character(80)                                     :: buffer
         character(8)                                      :: str1, str2
         character(1)                                      :: lf

         lf = achar(10)
      
         !
         ! Pass coordinates to a suitable 3D generic format
         !
         !$acc kernels
         points(:,:) = 0.0_rp
         points(:,1:ndime) = coord(:,1:ndime)
         !$acc end kernels

         !
         ! Adjust mu_sgs and envit to be nodal
         !
         do i = 1,nelem
            do j = 1, nnode
               envit(connec(i,j)) =  mu_e(i,j)
               mut(connec(i,j))   =  mu_sgs(i,j)
            end do
         end do
         !
         ! If case is periodic, adjust slave nodes
         !
         if (isPeriodic .and. present(masSla)) then
            !$acc parallel loop
            do iper = 1,nper
               u(masSla(iper,2),1) = u(masSla(iper,1),1)
               u(masSla(iper,2),2) = u(masSla(iper,1),2)
               u(masSla(iper,2),3) = u(masSla(iper,1),3)
               gradRho(masSla(iper,2),1) = gradRho(masSla(iper,1),1)
               gradRho(masSla(iper,2),2) = gradRho(masSla(iper,1),2)
               gradRho(masSla(iper,2),3) = gradRho(masSla(iper,1),3)
               curlU(masSla(iper,2),1) = curlU(masSla(iper,1),1)
               curlU(masSla(iper,2),2) = curlU(masSla(iper,1),2)
               curlU(masSla(iper,2),3) = curlU(masSla(iper,1),3)
               rho(masSla(iper,2)) = rho(masSla(iper,1))
               pr(masSla(iper,2)) = pr(masSla(iper,1))
               E(masSla(iper,2)) = E(masSla(iper,1))
               E(masSla(iper,2)) = E(masSla(iper,1))
               csound(masSla(iper,2)) = csound(masSla(iper,1))
               machno(masSla(iper,2)) = machno(masSla(iper,1))
               mu_fluid(masSla(iper,2)) = mu_fluid(masSla(iper,1))
               mut(masSla(iper,2)) = mut(masSla(iper,1))
               envit(masSla(iper,2)) = envit(masSla(iper,1))
               divU(masSla(iper,2)) = divU(masSla(iper,1))
               Qcrit(masSla(iper,2)) = Qcrit(masSla(iper,1))
            end do
            !$acc end parallel loop
         end if

         !
         ! Pass vector data to a suitable 3D generic format
         !
         !$acc kernels
         u3d(:,:) = 0.0_rp
         u3d(:,1:ndime) = u(:,1:ndime)
         !$acc end kernels
      
         !
         ! Pass cell list to VTK format
         !
         !$acc kernels
         cells(:,1) = nnode
         cells(:,2:nnode+1) = connec(:,1:nnode)-1
         !$acc end kernels
      
         !
         ! Define cell types
         !
         if (ndime .eq. 2) then
            if (nnode .eq. 4) then ! QUA04
               cellTypes = 9
            end if
         else if (ndime .eq. 3) then
            if (nnode .eq. 8) then ! HEX08
               cellTypes = 12
            else if (nnode .eq. 27) then ! HEX27
               cellTypes = 29
            else if (nnode .eq. 64) then ! HEX64
               cellTypes = 72
            end if
         end if
      
         !
         ! Open file with ascii input
         !
         call set_vtk_filename(filename,'results',istep)
         !call set_vtk_filename(filename,istep)
         !write(filename,'("vtkTstep_",i0,"-",i0,".vtk")'),istep,mpi_rank
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
            write(ivtk) dble(points(i,:))
         end do
         
         !
         ! Write cells
         !
         write(str1(1:8),'(i8)') nelem
         write(str2(1:8),'(i8)') nelem*(nnode+1)
         write(ivtk) lf//lf//'CELLS '//str1//' '//str2//lf
         do i = 1,nelem
            write(ivtk) cells(i,:)
         end do
         
         !
         ! Write cell types
         !
         write(str1(1:8),'(i8)') nelem
         write(ivtk) lf//lf//'CELL_TYPES '//str1//lf
         do i = 1,nelem
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
            write(ivtk) dble(rho(i))
         end do
         write(ivtk) lf//lf//'SCALARS PRESS double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(pr(i))
         end do
         write(ivtk) lf//lf//'SCALARS TENER double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(E(i))
         end do
         write(ivtk) lf//lf//'SCALARS CSOUND double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(csound(i))
         end do
         write(ivtk) lf//lf//'SCALARS DIVEU double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(divU(i))
         end do
         write(ivtk) lf//lf//'SCALARS QCRIT double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(Qcrit(i))
         end do
         write(ivtk) lf//lf//'SCALARS MACHNO double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(machno(i))
         end do
         write(ivtk) lf//lf//'SCALARS VISCO double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(mu_fluid(i))
         end do
         write(ivtk) lf//lf//'SCALARS SGSVI double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(mut(i))
         end do
         write(ivtk) lf//lf//'SCALARS ENVIT double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(envit(i))
         end do
         !
         ! Write point vector data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS VELOC double'//lf
         do i = 1,npoin
            write(ivtk) dble(u3d(i,:))
         end do
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS GRDEN double'//lf
         do i = 1,npoin
            write(ivtk) dble(gradRho(i,:))
         end do
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS VORTI double'//lf
         do i = 1,npoin
            write(ivtk) dble(curlU(i,:))
         end do

         !
         ! Write cell scalar data
         !
         !write(str1(1:8),'(i8)') nelem
         !write(ivtk) lf//lf//'CELL_DATA '//str1//lf
         !write(ivtk) 'SCALARS ENVIT double'//lf
         !write(ivtk) 'LOOKUP_TABLE default'//lf
         !do i = 1,nelem
         !   write(ivtk) mu_e(i,1)
         !end do
         !write(ivtk) 'SCALARS SGSVI double'//lf
         !write(ivtk) 'LOOKUP_TABLE default'//lf
         !do i = 1,nelem
         !   write(ivtk) mu_sgs(i,1)
         !end do
         
         close(ivtk)

      end subroutine write_vtk_binary

      subroutine write_vtkAVG_binary(isPeriodic,istep,npoin,nelem,coord,connec, &
                                 acuvel,acuve2,acurho,acupre,acumueff,acutim,nper,masSla)

         implicit none
         logical, intent(in)                                 :: isPeriodic
         integer(4), intent(in)                              :: nper, istep, npoin, nelem
         integer(4), intent(in)                              :: connec(nelem,nnode)
         integer(4), intent(in), optional                    :: masSla(nper,2)
         real(rp)   , intent(in)                              :: coord(npoin,ndime)
         real(rp)   , intent(inout)                           :: acutim
         real(rp)   , intent(inout), dimension(npoin)         :: acurho, acupre, acumueff
         real(rp)   , intent(inout), dimension(npoin,ndime)   :: acuvel, acuve2
         integer(4)                                          :: i, j, iper, ivtk=9, idime, ipoin
         integer(4)               , dimension(nelem,nnode+1) :: cells
         integer(4)               , dimension(nelem)         :: cellTypes
         real(rp)                  , dimension(npoin)         :: avrho, avpre, avmueff
         real(rp)                  , dimension(npoin,ndime)   :: avvel, avve2
         real(rp)                  , dimension(npoin,3)       :: points, u3d
         character(500)                                      :: filename
         character(80)                                       :: buffer
         character(8)                                        :: str1, str2
         character(1)                                        :: lf

         lf = achar(10)

         !
         ! Pass coordinates to a suitable 3D generic format
         !
         !$acc kernels
         points(:,:) = 0.0_rp
         points(:,1:ndime) = coord(:,1:ndime)
         !$acc end kernels

         !
         ! If case is periodic, adjust slave nodes
         !
         if (isPeriodic .and. present(masSla)) then
            !$acc parallel loop
            do iper = 1,nper
               acuvel(masSla(iper,2),1) = acuvel(masSla(iper,1),1)
               acuvel(masSla(iper,2),2) = acuvel(masSla(iper,1),2)
               acuvel(masSla(iper,2),3) = acuvel(masSla(iper,1),3)
               acuve2(masSla(iper,2),1) = acuve2(masSla(iper,1),1)
               acuve2(masSla(iper,2),2) = acuve2(masSla(iper,1),2)
               acuve2(masSla(iper,2),3) = acuve2(masSla(iper,1),3)
               acurho(masSla(iper,2)) = acurho(masSla(iper,1))
               acupre(masSla(iper,2)) = acupre(masSla(iper,1))
               acumueff(masSla(iper,2)) = acumueff(masSla(iper,1))
            end do
            !$acc end parallel loop
         end if

         !
         ! Divide accumulated vars by the accumulated time
         !
         !$acc kernels
         avrho(:) = acurho(:) / acutim
         avpre(:) = acupre(:) / acutim
         avmueff(:) = acumueff(:) / acutim
         avvel(:,:) = acuvel(:,:) / acutim
         avve2(:,:) = acuve2(:,:) / acutim
         !$acc end kernels

         !
         ! Favre average the rho*phi reynolds-averagedd variables
         !
         !$acc parallel loop
         do ipoin = 1,npoin
            !$acc loop seq
            do idime = 1,ndime
               avvel(ipoin,idime) = avvel(ipoin,idime)!/avrho(ipoin)
               avve2(ipoin,idime) = avve2(ipoin,idime)!/avrho(ipoin)
            end do
         end do
         !$acc end parallel loop

         !
         ! Reset the accumulated variables
         !
         !$acc kernels
         acurho(:) = 0.0_rp
         acupre(:) = 0.0_rp
         acumueff(:) = 0.0_rp
         acuvel(:,:) = 0.0_rp
         acuve2(:,:) = 0.0_rp
         !$acc end kernels
         acutim = 0.0_rp

         !
         ! Pass cell list to VTK format
         !
         !$acc kernels
         cells(:,1) = nnode
         cells(:,2:nnode+1) = connec(:,1:nnode)-1
         !$acc end kernels

         !
         ! Define cell types
         !
         if (ndime .eq. 2) then
            if (nnode .eq. 4) then ! QUA04
               cellTypes = 9
            end if
         else if (ndime .eq. 3) then
            if (nnode .eq. 8) then ! HEX08
               cellTypes = 12
            else if (nnode .eq. 27) then ! HEX27
               cellTypes = 29
            else if (nnode .eq. 64) then ! HEX64
               cellTypes = 72
            end if
         end if

         !
         ! Open file with ascii input
         !
         call set_vtk_filename(filename,'resAVG',istep)
         !write(filename,'("vtkAVG_step_",i0,".vtk")') istep
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
            write(ivtk) dble(points(i,:))
         end do

         !
         ! Write cells
         !
         write(str1(1:8),'(i8)') nelem
         write(str2(1:8),'(i8)') nelem*(nnode+1)
         write(ivtk) lf//lf//'CELLS '//str1//' '//str2//lf
         do i = 1,nelem
            write(ivtk) cells(i,:)
         end do

         !
         ! Write cell types
         !
         write(str1(1:8),'(i8)') nelem
         write(ivtk) lf//lf//'CELL_TYPES '//str1//lf
         do i = 1,nelem
            write(ivtk) cellTypes(i)
         end do

         !
         ! Write point scalar data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'POINT_DATA '//str1//lf
         write(ivtk) 'SCALARS AVDEN double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(avrho(i))
         end do
         write(ivtk) lf//lf//'SCALARS AVPRE double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(avpre(i))
         end do
         write(ivtk) lf//lf//'SCALARS AVMUEFF double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) dble(avmueff(i))
         end do

         !
         ! Write point vector data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS FAVVEL double'//lf
         do i = 1,npoin
            write(ivtk) dble(avvel(i,:))
         end do
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS FAVVE2 double'//lf
         do i = 1,npoin
            write(ivtk) dble(avve2(i,:))
         end do

         !
         ! Write cell scalar data
         !
         !write(str1(1:8),'(i8)') nelem
         !write(ivtk) lf//lf//'CELL_DATA '//str1//lf
         !write(ivtk) 'SCALARS ENVIT double'//lf
         !write(ivtk) 'LOOKUP_TABLE default'//lf
         !do i = 1,nelem
         !   write(ivtk) mu_e(i,1)
         !end do
         !write(ivtk) 'SCALARS SGSVI double'//lf
         !write(ivtk) 'LOOKUP_TABLE default'//lf
         !do i = 1,nelem
         !   write(ivtk) mu_sgs(i,1)
         !end do

         close(ivtk)

      end subroutine write_vtkAVG_binary

      subroutine write_vtk_binary_linearized(isPeriodic,istep,npoin,nelem,coord,connecLINEAR,connec, &
                                 rho,u,pr,E,csound,machno,mu_fluid,mu_e,mu_sgs,nper,masSla)
         implicit none
      
         logical, intent(in)                               :: isPeriodic
         integer(4), intent(in)                            :: nper,istep,npoin,nelem
         integer(4), intent(in)                            :: connecLINEAR(nelem*(porder**ndime),2**ndime)
         integer(4), intent(in)                            :: connec(nelem,nnode)
         integer(4), intent(in), optional                  :: masSla(nper,2)
         real(rp)   , intent(in)                            :: coord(npoin,ndime)
         real(rp)   , intent(inout), dimension(npoin)       :: rho, pr, E, mu_fluid, csound, machno
         real(rp)   , intent(inout), dimension(npoin,ndime) :: u
         real(rp)   , intent(inout), dimension(nelem,ngaus) :: mu_e
         real(rp)   , intent(inout), dimension(nelem,ngaus) :: mu_sgs
         integer(4)                                        :: i,j, iper, ivtk=9, nelem_l, nnode_l
         integer(4)            , dimension(nelem*(porder**ndime),2**ndime+1)  :: cells
         integer(4)            , dimension(nelem*(porder**ndime))          :: cellTypes
         real(rp)               , dimension(npoin,3)        :: points, u3d
         real(rp)               , dimension(npoin)          :: envit, mut
         character(500)                                    :: filename
         character(80)                                     :: buffer
         character(8)                                      :: str1
         character(16)                                     :: str2
         character(1)                                      :: lf

         lf = achar(10)

         nelem_l = nelem*(porder**ndime)
         nnode_l = 2**ndime 
      
         !
         ! Pass coordinates to a suitable 3D generic format
         !
         !$acc kernels
         points(:,:) = 0.0_rp
         points(:,1:ndime) = coord(:,1:ndime)
         !$acc end kernels
      
         !
         ! If case is periodic, adjust slave nodes
         !
         do i = 1,nelem
            do j = 1, nnode
               envit(connec(i,j)) =  mu_e(i,j)
               mut(connec(i,j))   =  mu_sgs(i,j)
            end do
         end do
         if (isPeriodic .and. present(masSla)) then
            !$acc parallel loop
            do iper = 1,nper
               u(masSla(iper,2),1) = u(masSla(iper,1),1)
               u(masSla(iper,2),2) = u(masSla(iper,1),2)
               u(masSla(iper,2),3) = u(masSla(iper,1),3)
               rho(masSla(iper,2)) = rho(masSla(iper,1))
               pr(masSla(iper,2)) = pr(masSla(iper,1))
               E(masSla(iper,2)) = E(masSla(iper,1))
               csound(masSla(iper,2)) = csound(masSla(iper,1))
               machno(masSla(iper,2)) = machno(masSla(iper,1))
               mu_fluid(masSla(iper,2)) = mu_fluid(masSla(iper,1))
               mut(masSla(iper,2)) = mut(masSla(iper,1))
               envit(masSla(iper,2)) = envit(masSla(iper,1))
            end do
            !$acc end parallel loop
         end if

         !
         ! Pass vector data to a suitable 3D generic format
         !
         !$acc kernels
         u3d(:,:) = 0.0_rp
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
         call set_vtk_filename(filename,'results',istep)
         !call set_vtk_filename(filename,istep)
         !write(filename,'("vtkTstep_",i0,".vtk")') istep
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
         write(str2(1:16),'(i16)') nelem_l*(nnode_l+1)
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
         write(ivtk) lf//lf//'SCALARS TENER double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) E(i)
         end do
         write(ivtk) lf//lf//'SCALARS CSOUND double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) csound(i)
         end do
         write(ivtk) lf//lf//'SCALARS MACHNO double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) machno(i)
         end do
         write(ivtk) lf//lf//'SCALARS ENVIT double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) envit(i)
         end do
         write(ivtk) lf//lf//'SCALARS SGSVI double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) mut(i)
         end do
         write(ivtk) lf//lf//'SCALARS VISCO double '//lf
         write(ivtk) 'LOOKUP_TABLE default'//lf
         do i = 1,npoin
            write(ivtk) mu_fluid(i)
         end do
         
         !
         ! Write point vector data
         !
         write(str1(1:8),'(i8)') npoin
         write(ivtk) lf//lf//'VECTORS VELOC double'//lf
         do i = 1,npoin
            write(ivtk) u3d(i,:)
         end do
         
         !!
         !! Write cell scalar data
         !!
         !write(str1(1:8),'(i8)') nelem_l
         !write(ivtk) lf//lf//'CELL_DATA '//str1//lf
         !write(ivtk) 'SCALARS ENVIT double'//lf
         !write(ivtk) 'LOOKUP_TABLE default'//lf
         !do i = 1,nelem_l
         !   write(ivtk) mu_e(i,1)
         !end do
         !write(ivtk) 'SCALARS SGSVI double'//lf
         !write(ivtk) 'LOOKUP_TABLE default'//lf
         !do i = 1,nelem_l
         !   write(ivtk) mu_sgs(i,1)
         !end do
         
         close(ivtk)
      
      end subroutine write_vtk_binary_linearized

      subroutine read_vtk_binary(isPeriodic,istep,npoin,nelem,coord,connec, &
                                 rho,u,pr,E,csound,machno,mu_fluid,mu_e,mu_sgs,nper,masSla)
         implicit none
      
         logical, intent(in)                               :: isPeriodic
         integer(4), intent(in)                            :: nper,istep,npoin,nelem
         integer(4), intent(in)                            :: connec(nelem,nnode)
         integer(4), intent(in), optional                  :: masSla(nper,2)
         real(rp)   , intent(in)                            :: coord(npoin,ndime)
         real(rp)   , intent(inout), dimension(npoin)       :: rho, pr, E, csound, machno, mu_fluid
         real(rp)   , intent(inout), dimension(npoin,ndime) :: u
         real(rp)   , intent(inout), dimension(nelem,ngaus) :: mu_e
         real(rp)   , intent(inout), dimension(nelem,ngaus) :: mu_sgs
         integer(4)                                        :: i, iper, ivtk=99
         integer(4)            , dimension(nelem,nnode+1)  :: cells
         integer(4)            , dimension(nelem)          :: cellTypes
         integer(4)            , dimension(npoin)          :: envit, mut
         real(rp)               , dimension(npoin,3)        :: points
         character(500)                                    :: filename

         write(1,*) " begining the reading of the vtk binary "
         call flush(111)

         !
         ! Open file with ascii input
         !
         call set_vtk_filename(filename,'results',istep)
         !call set_vtk_filename(filename,istep)
         !write(filename,'("vtkTstep_",i0,".vtk")') istep
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
            read(ivtk) csound(i)
         end do
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) machno(i)
         end do
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) envit(i)
         end do
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         call skip_line(ivtk)
         do i = 1,npoin
            read(ivtk) mut(i)
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

         write(1,*) " end the reading of the vtk binary "
         call flush(111)
         
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

      subroutine set_vtk_filename(filename,prefix,istep)
         implicit none
         character(len=*), intent(out) :: filename
         character(len=*), intent(in) :: prefix
         integer, intent(in) :: istep

         write(filename,'(A,"_step_",i0,"_s",i0,"-r",i0,".vtk")'),prefix,istep,mpi_size,mpi_rank
         !write(filename,'("vtkTstep_",i0,".vtk")') istep         !old, for serial output

      end subroutine set_vtk_filename




end module mod_output
