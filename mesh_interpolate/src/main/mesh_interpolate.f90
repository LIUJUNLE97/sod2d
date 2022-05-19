program mesh_interpolate

   use mod_constants
   use mesh_reader
   use quadrature_rules
   use elem_hex
   use inicond_reader
   use mod_output
   use mod_maths

   implicit none

   integer(4)              :: igaus, ielem, inode, idime, jelem
   integer(4)              :: lin_npoin, lin_nelem, lin_nboun
   integer(4)              :: npoin, nelem, nboun
   integer(4)              :: counter
   integer(4), allocatable :: lin_connec(:,:), lin_bound(:,:)
   integer(4), allocatable :: atoIJK(:), listHEX08(:,:)
   integer(4), allocatable :: connec(:,:), bound(:,:)
   integer(4), allocatable :: lmatch(:,:)
	real(8)                 :: xi, eta, zeta
   real(8),    allocatable :: lin_xyz(:,:), xyz(:,:), aux(:,:)
   real(8),    allocatable :: xnodes(:,:), wgp(:)
   real(8),    allocatable :: Ngp(:,:), dNgp(:,:,:)
   real(8),    allocatable :: lin_u(:,:), lin_rho(:), lin_pr(:), lin_E(:), lin_mu_fluid(:)
   real(8),    allocatable :: u(:,:), rho(:), pr(:), E(:), mu_fluid(:)
   character(500)          :: lin_file_path, lin_file_name
   character(500)          :: file_path, file_name

   ! Define case name
   write(lin_file_name,*) "cube"
   file_name = lin_file_name

   ! Read the linear mesh
   write(lin_file_path,*) "./mesh_lin/"
   call read_dims(lin_file_path,lin_file_name,lin_npoin,lin_nelem,lin_nboun)
   allocate(lin_connec(lin_nelem,lin_nnode), lin_bound(lin_nboun,lin_npbou), lin_xyz(lin_npoin,ndime))
   call read_geo_dat(lin_file_path,lin_file_name,lin_npoin,lin_nelem,lin_nboun,lin_connec,lin_bound,lin_xyz)

   ! Read the high-order mesh
   write(file_path,*) "./mesh/"
   call read_dims(file_path,file_name,npoin,nelem,nboun)
   allocate(connec(nelem,nnode), bound(nboun,npbou), xyz(npoin,ndime))
   call read_geo_dat(file_path,file_name,npoin,nelem,nboun,connec,bound,xyz)

   ! Check if the mesh is compatible
   if (nelem .ne. lin_nelem) then
      write(*,*) "The number of elements in the high-order mesh is different from the number of elements in the linear mesh"
      error stop
   end if
   
   ! Create a list of matching pairs of high-order and linear elements
   ! Col. 1: low-order element number
   ! Col. 2: high-order element number
   allocate(lmatch(nelem,2))
   allocate(aux(lin_nnode,ndime))
   linear: do ielem = 1, nelem
      lmatch(ielem,1) = ielem
      aux(1:lin_nnode,1:ndime) = lin_xyz(lin_connec(ielem,1:lin_nnode),1:ndime)
      nonlinear: do jelem = 1,nelem
         counter = 0
         corners: do inode = 1,nncorner
            do idime = 1,ndime
               if (aux(inode,idime) .ne. xyz(connec(jelem,inode),idime)) then
                  exit corners
               else if (aux(inode,idime) .eq. xyz(connec(jelem,inode),idime)) then
                  counter = counter + 1
               end if
            end do
         end do corners
         if (counter .eq. ndime*nncorner) then
            lmatch(ielem,2) = jelem
         else
            write(*,*) "Did not found a matching pair of elements! Check your mesh!"
            error stop
         end if
      end do nonlinear
   end do linear

   ! Allocate data for linear and high order  variables
   allocate(lin_u(lin_npoin,ndime))
   allocate(lin_rho(lin_npoin))
   allocate(lin_pr(lin_npoin))
   allocate(lin_E(lin_npoin))
   allocate(lin_mu_fluid(lin_npoin))
   allocate(u(npoin,ndime))
   allocate(rho(npoin))
   allocate(pr(npoin))
   allocate(E(npoin))
   allocate(mu_fluid(npoin))

   ! Read initial condition from linear mesh
   if (flag_readVTK == 1) then
      call read_vtk_binary(lin_npoin,lin_nelem,lin_xyz,lin_connec,lin_rho,lin_u,lin_pr,lin_E,lin_mu_fluid)
   else
      ! Read ascii file *.alya
      call read_veloc(lin_npoin,lin_file_path,lin_u)
      call read_densi(lin_npoin,lin_file_path,lin_rho)
      call read_press(lin_npoin,lin_file_path,lin_pr)
   end if

   ! Create master element geometries (HEX08 as a subset of HEX high order)
   allocate(atoIJK(nnode))
   allocate(listHEX08(porder**ndime,nncorner))
   allocate(xnodes(ngaus,ndime))
   allocate(wgp(ngaus))
   if (nnode == 64) then
      call hex64(1.0d0,1.0d0,1.0d0,atoIJK)
   end if
   call lagrange_hex(atoIJK,xnodes,wgp)

   ! Generate linear shape functions at each high-order node
	allocate(Ngp(ngaus,lin_nnode),dNgp(ndime,lin_nnode,ngaus))
	do igaus = 1, ngaus
		xi = xnodes(igaus,1)
		eta = xnodes(igaus,2)
		zeta = xnodes(igaus,3)
      call hex08(xi,eta,zeta,Ngp(igaus,:),dNgp(:,:,igaus))
	end do

   ! Loop over all elements and interpolate variables
   !$acc parallel loop gang
   do ielem = 1,nelem
      !$acc loop vector
      do inode = nncorner,nnode
         !$acc loop seq
         do idime = 1,ndime
            call var_interpolate(lin_u(lin_connec(ielem,:),idime),Ngp(inode,:),u(connec(lmatch(ielem,2),inode),idime))
         end do
         call var_interpolate(lin_rho(lin_connec(ielem,:)),Ngp(inode,:),rho(connec(lmatch(ielem,2),inode)))
         call var_interpolate(lin_pr(lin_connec(ielem,:)),Ngp(inode,:),pr(connec(lmatch(ielem,2),inode)))
      end do
   end do
   !$acc end parallel loop

end program mesh_interpolate
