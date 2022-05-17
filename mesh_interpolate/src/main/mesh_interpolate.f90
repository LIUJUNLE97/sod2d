program mesh_interpolate

    use mod_constants
    use mesh_reader

    implicit none

    integer(4)     :: lin_npoin, lin_nelem, lin_nboun
    integer(4)     :: npoin, nelem, nboun
    character(500) :: file_path, file_name

    ! Read the linear mesh
    write(file_path,*) "./mesh/"
    write(file_name,*) "cube"
    call read_dims(file_path,file_name,lin_npoin,lin_nelem,nboun)
end program mesh_interpolate
