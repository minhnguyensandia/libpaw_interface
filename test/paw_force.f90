! This subroutine calls pawgrnl, which gives term F2 (Eq. 31)
! of the PAW paper (Torrent 2008)

subroutine paw_force(natom,ntypat,typat,ngfftdg,nfft,nhat,xred,rprimd,gmet,ucvol,vtrial,vxc,nspden, &
    grnl,nlstr)
    use libpaw_mod
    use m_kg
    use m_paw_dfpt
    implicit none

    integer :: ngfftdg(3)
    integer :: nfft,nspden,natom

    real*8  :: nhat(nfft,nspden),vtrial(nfft,nspden),vxc(nfft,nspden)

    real*8  :: grnl(3*natom),nlstr(6)
    integer :: mgrid, ntypat, typat(natom)
    real*8, allocatable :: ph1d(:,:)

    real*8  :: qphon(3)
    real*8  :: rprimd(3,3),xred(3,natom),gmet(3,3),ucvol

    mgrid = maxval(ngfftdg)
    allocate(ph1d(2,3*(2*mgrid+1)*natom))
    ph1d = 0d0
    qphon = 0.0

    call getph(atindx,natom,ngfftdg(1),ngfftdg(2),ngfftdg(3),ph1d,xred)

    call pawgrnl(atindx1,nspden,grnl,gsqcut,mgrid,natom,natom,nattyp,nfft,ngfftdg,nhat, &
        nlstr,nspden,ntypat,1,0,1,0,pawang,pawfgrtab,pawrhoij,pawtab,ph1d, &
        qphon,rprimd,gmet,typat,ucvol,vtrial,vxc,xred)
end subroutine