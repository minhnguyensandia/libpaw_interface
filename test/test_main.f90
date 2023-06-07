program main

    implicit none

    write(*,*) 'Test code for libpaw'
    call prepare_libpaw
    call get_nhat
end program