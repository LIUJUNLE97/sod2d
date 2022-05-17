program mesh_interpolate

    use mod_constants
    use mesh_reader

    implicit none

    integer(4)     :: lin_npoin, lin_nelem, lin_nboun
    integer(4)     :: npoin, nelem, nboun
    character(500) :: lin_file_path, lin_file_name
    character(500) :: file_path, file_name

    ! Define case name
    write(lin_file_name,*) "cube"
    file_name = lin_file_name

    ! Read the linear mesh
    write(lin_file_path,*) "./mesh_lin/"
    call read_dims(lin_file_path,lin_file_name,lin_npoin,lin_nelem,lin_nboun)

    ! Read the high-order mesh
    write(file_path,*) "./mesh/"
    call read_dims(file_path,file_name,npoin,nelem,nboun)

end program mesh_interpolate
