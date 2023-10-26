subroutine prepare_libpaw(ecut,ecutpaw,gmet,rprimd,gprimd,ucvol,ngfft,ngfftdg, &
        natom,ntypat,typat,xred,ixc,xclevel,filename_list,nspden,nsppol,epsatm)
    use m_pawpsp
    use m_pawxmlps
    use m_paw_init
    use m_kg
    use m_paw_occupancies
    use m_paw_nhat
    use libpaw_mod

    implicit none

    integer :: stat
    integer :: n3, i3, n2, i2
    integer :: it, ia

    real*8  :: ecut, ecutpaw !a coarse grid, and a fine grid for PAW
    real*8  :: gmet(3,3) !reciprocal space lattice vectors
    real*8  :: rprimd(3,3) !lattice vectors
    real*8  :: gprimd(3,3) !reciprocal space lattice vectors
    integer :: ngfft(3),ngfftdg(3)
    real*8  :: ucvol !volume of unit cell
    integer :: ixc, xclevel !functional; will be set externally in practice
    integer :: ntypat, natom
    integer :: typat(natom)
    real*8  :: xred(3,natom)
    real*8  :: epsatm(ntypat)
    integer :: nspden,nsppol!number of spin components for density and wavefunctions
    integer :: llmax
    !for normal nspin = 1,2 calculations, nsppol is set to be same as nspden

    character(len=264) :: filename_list(ntypat)

    hyb_mixing = 0.0
    hyb_range_fock = 0.0

    allocate(zion(ntypat),znucl(ntypat), nattyp(ntypat), lexexch(ntypat), lpawu(ntypat), l_size_atm(natom))
    allocate(pawrad(ntypat), pawtab(ntypat), pawrhoij(natom), paw_ij(natom), &
        & paw_an(natom), pawfgrtab(natom))
    allocate(atindx(natom),atindx1(natom))

    ! not supporting dftu for now!
    lexexch = 0
    lpawu = 0

    ! Process atomic information
    call map_atom_index(ntypat,natom,typat)

    ! Process energy cutoff
    call getcut(ecut,gmet,gsqcut,iboxcut,ngfft)
    call getcut(ecutpaw,gmet,gsqcutdg,iboxcut,ngfftdg)
    
    allocate(qgrid_ff(mqgrid),qgrid_vl(mqgrid),stat=stat)
    if(stat/=0) then
        write(*,*) 'problem allocating mqgrid'
        call exit(1)
    endif
    
    call generate_qgrid(gsqcut,qgrid_ff,mqgrid)
    call generate_qgrid(gsqcutdg,qgrid_vl,mqgrid)

    allocate(ffspl(mqgrid,2,lnmax), vlspl(mqgrid,2,ntypat))

    llmax = 0

    do it = 1, ntypat
        ! Read paw input files
        call rdpawpsxml(filename_list(it), pawsetup)
        call rdpawpsxml(filename_list(it), paw_setuploc)
        call pawpsp_read_header_xml(lloc, lmax, pspcod, pspxc,&
            & pawsetup, r2well, zion(it), znucl(it))
        call pawpsp_read_pawheader(pawpsp_header%basis_size,&
            &   lmax,pawpsp_header%lmn_size,&
            &   pawpsp_header%l_size, pawpsp_header%mesh_size,&
            &   pawpsp_header%pawver, pawsetup,&
            &   pawpsp_header%rpaw, pawpsp_header%rshp, pawpsp_header%shape_type)

        ! Process onsite information
        call pawtab_nullify(pawtab(it))
        call pawtab_set_flags(pawtab(it),has_tvale=1,has_vhnzc=1,has_vhtnzc=1)
        call pawpsp_17in(epsatm(it), ffspl, icoulomb, ipsp, hyb_mixing, ixc, lmax,&
                &       lnmax, pawpsp_header%mesh_size, mqgrid, mqgrid, pawpsp_header,&
                &       pawrad(it), pawtab(it), xcdev, qgrid_ff, qgrid_vl, usewvl, usexcnhat,&
                &       vlspl(:,:,it), xcccrc, xclevel, denpos, zion(it), znucl(it))
        call paw_setup_free(pawsetup)
        call paw_setup_free(paw_setuploc)

        !write(13,*) 'dij0',pawtab(1)%dij0

        llmax = max(lmax,llmax)
    enddo

    do ia = 1, natom
        it = typat(ia)
        l_size_atm(ia) = pawtab(it)%l_size
    enddo

    mpsang = llmax + 1
    call pawinit(effmass_free, gnt_option,gsqcut_eff,hyb_range_fock,lcutdens,lmix,mpsang,nphi,nsym,ntheta,&
        &     pawang,pawrad,pawspnorb,pawtab,xcdev,xclevel,usepotzero)
    !write(13,*) 'eijkl',pawtab(1)%eijkl

    !See m_scfcv_core.f90, lines 669 and forth
    call pawfgrtab_init(pawfgrtab, cplex, l_size_atm, nspden, typat)

    call paw_an_nullify(paw_an)
    call paw_ij_nullify(paw_ij)

    call paw_an_init(paw_an, natom, ntypat, 0, 0, nspden, cplex, &
        & xcdev, typat, pawang, pawtab, has_vxc = 1, has_vxc_ex = 1)
    call paw_ij_init(paw_ij, cplex, nspinor, nsppol, nspden, pawspnorb, &
        & natom, ntypat, typat, pawtab, &
        & has_dij = 1, has_dijhartree = 1, has_dijso = 1, has_pawu_occ = 1, has_exexch_pot = 1)
    
    call initrhoij(cplex, lexexch, lpawu, natom, natom, nspden, nspinor, &
        & nsppol, ntypat, pawrhoij, pawspnorb, pawtab, 1, spinat, typat)

    n3 = ngfftdg(3)
    allocate(fftn3_distrib(n3), ffti3_local(n3))
    !This is the case when running serially
    fftn3_distrib = 0
    
    do i3 = 1, n3
        ffti3_local(i3) = i3
    enddo

    n2 = ngfftdg(2)
    allocate(fftn2_distrib(n2), ffti2_local(n2))
    !This is the case when running serially
    fftn2_distrib = 0

    do i2 = 1, n2
        ffti2_local(i2) = i2
    enddo

    call nhatgrid(atindx1, gmet, natom, natom, nattyp, ngfftdg, ntypat, &
        & 0, 1, 0, 0, 0, & !optcut, optgr0, optgr1, optgr2, optrad
        & pawfgrtab, pawtab, rprimd, typat, ucvol, xred, &
        & n3, fftn3_distrib, ffti3_local)

    !do ia = 1, natom
    !    write(16,*) pawfgrtab(ia)%gylm
    !    write(16,*)
    !enddo
contains
    subroutine generate_qgrid(gsqcut,qgrid,mqgrid)
        real*8  :: gsqcut
        real*8  :: qmax, dq
        real*8  :: qgrid(mqgrid)

        integer :: mqgrid, iq

        qmax = 1.2d0 * sqrt(gsqcut)
        dq = qmax/(1.0*(mqgrid-1))
        do iq = 1,mqgrid
            qgrid(iq) = (iq-1)*dq
        enddo
    end subroutine
end subroutine
