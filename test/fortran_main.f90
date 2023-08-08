! The input variables:
! 1. ecut, ecutpaw : kinetic energy cutoff of the planewave basis set
! there will be one coarse grid for density/potential, and a fine grid for PAW
! the unit is in Hartree
! 2. rprimd, gprimd : real and reciprocal space lattice vectors, respectively
! unit for rprimd is in Bohr, and for gprimd is in Bohr^-1
! 3. gmet : reciprocal space metric (bohr^-2)
! 4. ucvol : volume of unit cell (Bohr^3)
! 5. ngfft, ngfftdg : dimension of FFT grids of the corase and fine grids
! 6. natom, ntypat, typat: #. atoms, #. element types
! and typat records the type of each atom
! 7. xred : coordinate of each atom, in terms of rprimd (namely, direct coordinate)
! 8. filename_list : filename of the PAW xml files for each element

subroutine fortran_main(ecut,ecutpaw,gmet,rprimd,gprimd,ucvol, &
    ngfft,ngfftdg,natom,ntypat,typat,xred,filename_list)
    implicit none

    real*8  :: ecut, ecutpaw !a coarse grid, and a fine grid for PAW
    real*8  :: gmet(3,3) !reciprocal space metric
    real*8  :: rprimd(3,3) !lattice vectors
    real*8  :: gprimd(3,3) !reciprocal space lattice vectors
    integer :: ngfft(3),ngfftdg(3)
    real*8  :: ucvol !volume of unit cell
    integer :: ixc, xclevel !functional; will be set externally in practice
    integer :: ntypat, natom
    !integer, allocatable :: typat(:)
    !real*8,  allocatable :: xred(:,:)

    integer :: typat(natom)
    real*8  :: xred(3,natom)
    
    character(len=264) :: filename_list(ntypat)

    write(*,*) '1. Setting up libpaw'
    !open(unit=10,file='pawfiles')
    !open(unit=11,file='input')

    ! Read some input variables
    !call scan_input_double_scalar('ecut',ecut)
    !call scan_input_double_scalar('ecutpaw',ecutpaw)
    !call scan_input_double('gmet',gmet,9)
    !call scan_input_double('rprimd',rprimd,9)
    !call scan_input_double('gprimd',gprimd,9)
    !call scan_input_double_scalar('ucvol',ucvol)
    !call scan_input_int('ngfft',ngfft,3)
    !call scan_input_int('ngfftdg',ngfftdg,3)
    !call scan_input_int_scalar('natom',natom)
    !call scan_input_int_scalar('ntypat',ntypat)

    ! (Temporary) set xc functional type
    ixc = 7 ! corresponds to PW92 LDA functional
    xclevel = 1

    !allocate(typat(natom), xred(3,natom))
    
    !call scan_input_int('typat',typat,natom)
    !call scan_input_double('xred',xred,3*natom)

    write(*,*) 'Test code for libpaw'
    call prepare_libpaw(ecut,ecutpaw,gmet,rprimd,gprimd,ucvol,ngfft,ngfftdg, &
        natom,ntypat,typat,xred,ixc,xclevel,filename_list)

    call get_vloc_ncoret(ngfftdg,ngfft,natom,ntypat,rprimd,gprimd,gmet,ucvol,xred)
    call get_nhat(natom,ntypat,xred,ngfft,ngfftdg,gprimd,rprimd,ucvol)
    call get_dij(natom,ntypat,ixc,xclevel,ngfft,ngfftdg,xred,ucvol,gprimd)

!contains
!    subroutine scan_input_double(name_in,value,n)
!        character(*)      :: name_in
!        character(len=20) :: name
!        integer           :: n
!        real*8            :: value(n)
!
!        read(11,*) name, value
!        if(trim(name)/=trim(name_in)) then
!            write(*,*) 'variable name does not match : ',trim(name),trim(name_in)
!        endif
!    end subroutine
!
!    subroutine scan_input_double_scalar(name_in,value)
!        character(*)      :: name_in
!        character(len=20) :: name
!        real*8            :: value
!
!        read(11,*) name, value
!        if(trim(name)/=trim(name_in)) then
!            write(*,*) 'variable name does not match : ',trim(name),trim(name_in)
!        endif
!    end subroutine
!
!    subroutine scan_input_int(name_in,value,n)
!        character(*)      :: name_in
!        character(len=20) :: name
!        integer           :: n
!        integer           :: value(n)
!
!        read(11,*) name, value
!        if(trim(name)/=trim(name_in)) then
!            write(*,*) 'variable name does not match : ',trim(name),trim(name_in)
!        endif
!    end subroutine
!
!    subroutine scan_input_int_scalar(name_in,value)
!        character(*)      :: name_in
!        character(len=20) :: name
!        integer           :: value
!
!        read(11,*) name, value
!        if(trim(name)/=trim(name_in)) then
!            write(*,*) 'variable name does not match : ',trim(name),trim(name_in)
!        endif
!    end subroutine

end subroutine
