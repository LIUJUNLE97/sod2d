program mesh_interpolate

    use mod_constants
    use mesh_reader
    use quadrature_rules
    use elem_hex

    implicit none

    integer(4)              :: lin_npoin, lin_nelem, lin_nboun
    integer(4)              :: npoin, nelem, nboun
    integer(4), allocatable :: lin_connec(:,:), lin_bound(:,:)
    integer(4), allocatable :: atoIJK(:), listHEX08(:,:)
    integer(4), allocatable :: connec(:,:), bound(:,:)
    real(8),    allocatable :: lin_xyz(:,:), xyz(:,:)
    real(8),    allocatable :: xnodes(:,:), wgp(:)
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

    ! Create master element geometries (HEX08 as a subset of HEX high order)
    allocate(atoIJK(nnode))
    allocate(listHEX08(porder**ndime,nncorner))
    allocate(xnodes(ngaus,ndime))
    allocate(wgp(ngaus))
    if (ndime == 3) then
       if (nnode == 64) then
          call hex64(1.0d0,1.0d0,1.0d0,atoIJK)
       end if
    end if
    call lagrange_hex(atoIJK,xnodes,wgp)

end program mesh_interpolate
