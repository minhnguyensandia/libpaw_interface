subroutine get_dij(natom,ntypat,ixc,xclevel,nfft,nspden,xred,ucvol,gprimd,vtrial,vxc)
    use libpaw_mod
    use m_paw_denpot
    use m_pawdij

    implicit none

    real*8  :: compch_sph, epaw, epawdc
    real*8  :: nucdipmom(3,natom), qphon(3)
    real*8  :: vtrial(cplex*nfft,nspden),vxc(cplex*nfft,nspden)

    integer :: natom,ntypat,ixc,xclevel
    integer :: nfft, nspden
    real*8  :: xred(3,natom),ucvol,gprimd(3,3)

    write(*,*) '4. Generating dij from on-site rhoij, v_ks and v_xc'

    ! Note : in PAW calculations, v_h is calculated with ncomp + ntilde + ncoretilde
    ! while v_xc is calculated with ntilde + ncoretilde

    nucdipmom = 0.0
    qphon = 0.0

    call pawdenpot(compch_sph, epaw, epawdc, 0, ixc, natom, natom, & !ipert
        & nspden, ntypat, nucdipmom, 0, 0, paw_an, paw_an, paw_ij, & !nzlmopt, option
        & pawang, 0, pawrad, pawrhoij, 0, pawtab, xcdev, 1d0, xclevel, & !pawprtvol, pawspnorb, spnorbscl
        & denpos, ucvol, znucl)
    
    !write(21,*) 'epaw',epaw
    !write(21,*) 'epawdc',epawdc

    call pawdij(cplex, 0, gprimd, 0, natom, natom, nfft, nfft, nspden, ntypat, & !enunit, ipert
        & paw_an, paw_ij, pawang, pawfgrtab, 0, pawrad, pawrhoij, 0, pawtab, & !pawprtvol, pawspnorb
        & xcdev, qphon, 1d0, ucvol, 0d0, vtrial, vxc, xred) !spnorbscl, charge

end subroutine
