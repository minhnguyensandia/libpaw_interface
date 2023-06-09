! Following m_prcref.F90

subroutine get_vloc_ncoret
    use libpaw_mod
    use m_atm2fft
    use m_kg

    implicit none

    real*8, allocatable :: ncoret(:), vloc(:), ph1d(:,:)
    real*8, allocatable :: dummy(:),dummy1(:),dummy2(:),dummy3(:),dummy4(:),dummy5(:),dummy8(:),dummy9(:)
    real*8  :: dummy6(6),dummy7(6),rcut,vprtrb(2)
    integer :: nfft, mgrid, qprtrb(3)

    write(*,*) '2. Generating pseudo charge density and vloc on FFT grid'

    nfft = ngfftdg(1) * ngfftdg(2) * ngfftdg(3)
    mgrid = maxval(ngfftdg)

    allocate(ncoret(nfft),vloc(nfft),ph1d(2,3*(2*mgrid+1)*natom))
    ph1d = 0d0
    qprtrb = 0
    rcut = 0d0
    vprtrb = 0d0

    call getph(atindx,natom,ngfftdg(1),ngfftdg(2),ngfftdg(3),ph1d,xred)

    call atm2fft(atindx1,ncoret,vloc,dummy,dummy2,dummy9,dummy1,gmet,gprimd,dummy3,dummy4,gsqcutdg, &
        & mgrid,mqgrid,natom,nattyp,nfft,ngfftdg,ntypat,1,0,0,0,1,1,0,1, &
        & pawtab,ph1d,qgrid_vl,qprtrb,rcut,dummy5,rprimd,dummy6,dummy7,&
        & ucvol,1,dummy8,dummy8,dummy8,vprtrb,vlspl,&
        & ngfftdg(2),fftn2_distrib,ffti2_local,ngfftdg(3),fftn3_distrib,ffti3_local)
    
end subroutine