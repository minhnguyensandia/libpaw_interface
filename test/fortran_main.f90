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
    use libpaw_mod, only : cplex, pawrhoij
    implicit none

    real*8  :: ecut, ecutpaw !a coarse grid, and a fine grid for PAW
    real*8  :: gmet(3,3) !reciprocal space metric
    real*8  :: rprimd(3,3) !lattice vectors
    real*8  :: gprimd(3,3) !reciprocal space lattice vectors
    integer :: ngfft(3),ngfftdg(3)
    real*8  :: ucvol !volume of unit cell
    integer :: ixc, xclevel !functional; will be set externally in practice
    integer :: ntypat, natom, iatom
    integer :: nfft
    integer :: nspden,nsppol
    integer :: nrhoijsel,size_rhoij
    integer, allocatable :: rhoijselect(:)
    real*8,  allocatable :: rhoijp(:,:)
    real*8,  allocatable :: vtrial(:,:), vxc(:,:)
    real*8,  allocatable :: vloc(:), ncoret(:)
    real*8,  allocatable :: nhat(:,:), nhatgr(:,:,:)
    !integer, allocatable :: typat(:)
    !real*8,  allocatable :: xred(:,:)

    integer :: typat(natom)
    real*8  :: xred(3,natom)
    
    character(len=264) :: filename_list(ntypat)

    write(*,*) '1. Setting up libpaw'

    ! (Temporary) set xc functional type
    ixc = 7 ! corresponds to PW92 LDA functional
    xclevel = 1
    nspden = 1
    nsppol = 1

    write(*,*) 'Test code for libpaw'
    call prepare_libpaw(ecut,ecutpaw,gmet,rprimd,gprimd,ucvol,ngfft,ngfftdg, &
        natom,ntypat,typat,xred,ixc,xclevel,filename_list,nspden,nsppol)

    nfft = ngfftdg(1) * ngfftdg(2) * ngfftdg(3)
    allocate(ncoret(nfft),vloc(nfft))
    allocate(vtrial(cplex*nfft,nspden),vxc(cplex*nfft,nspden))
    allocate(nhat(nfft,nspden), nhatgr(nfft,nspden,3))

    call get_vloc_ncoret(ngfftdg,nfft,natom,ntypat,rprimd,gprimd,gmet,ucvol,xred,ncoret,vloc)
    
    open(unit=10,file='rhoij')
    do iatom = 1, natom
        read(10,*) nrhoijsel

        size_rhoij = size(pawrhoij(iatom)%rhoijselect)
        allocate(rhoijselect(size_rhoij),rhoijp(size_rhoij,nspden))
        read(10,*) rhoijselect
        read(10,*) rhoijp
        call set_rhoij(iatom,nrhoijsel,size_rhoij,nspden,rhoijselect,rhoijp)
        deallocate(rhoijselect,rhoijp)
    enddo
    close(10)

    call get_nhat(natom,ntypat,xred,ngfft,nfft,nspden,gprimd,rprimd,ucvol,nhat,nhatgr)

    open(unit=10,file='veff')
    read(10,*) vtrial
    read(10,*) vxc
    close(10)

    call calculate_dij(natom,ntypat,ixc,xclevel,nfft,nspden,xred,ucvol,gprimd,vtrial,vxc)

end subroutine
