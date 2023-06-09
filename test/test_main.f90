program main

    implicit none

    write(*,*) 'Test code for libpaw'
    call prepare_libpaw
    call get_vloc_ncoret
    call get_nhat
    call get_dij
end program