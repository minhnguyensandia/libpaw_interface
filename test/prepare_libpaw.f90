program main
    use m_pawpsp
    use m_pawxmlps
    use m_paw_init
    use m_kg
    use m_paw_occupancies, only : initrhoij
    use libpaw_mod

    implicit none

    integer :: stat

    write(*,*) 'Test code for setting up libpaw'
    open(unit=10,file='pawfiles')
    open(unit=11,file='input')

    ! Read some input variables
    call scan_input_double('ecut',ecut,1)
    call scan_input_double('ecutpaw',ecutpaw,1)
    call scan_input_double('gmet',gmet,9)
    call scan_input_int('ngfft',ngfft,3)
    call scan_input_int('ngfftdg',ngfftdg,3)
    
    ! (Temporary) set xc functional type
    ixc = 7 ! corresponds to PW92 LDA functional
    xclevel = 1
    hyb_mixing = 0.0
    hyb_range_fock = 0.0

    ! Process energy cutoff
    call getcut(ecut(1),gmet,gsqcut,iboxcut,ngfft)
    call getcut(ecutpaw(1),gmet,gsqcutdg,iboxcut,ngfftdg)
    
    allocate(qgrid_ff(mqgrid),qgrid_vl(mqgrid),stat=stat)
    if(stat/=0) then
        write(*,*) 'problem allocating mqgrid'
        call exit(1)
    endif
    
    call generate_qgrid(gsqcut,qgrid_ff,mqgrid)
    call generate_qgrid(gsqcutdg,qgrid_vl,mqgrid)

    allocate(ffspl(mqgrid,2,lnmax), vlspl(mqgrid,2))

    !do
        read(10,*,iostat=stat) filename
        !if(stat/=0) exit

        ! Read paw input files
        call rdpawpsxml(filename, pawsetup)
        call rdpawpsxml(filename, paw_setuploc)
        call pawpsp_read_header_xml(lloc, lmax, pspcod, pspxc,&
            & pawsetup, r2well, zion, znucl)
        call pawpsp_read_pawheader(pawpsp_header%basis_size,&
            &   lmax,pawpsp_header%lmn_size,&
            &   pawpsp_header%l_size, pawpsp_header%mesh_size,&
            &   pawpsp_header%pawver, pawsetup,&
            &   pawpsp_header%rpaw, pawpsp_header%rshp, pawpsp_header%shape_type)

        ! Process onsite information
        call pawtab_set_flags(pawtab,has_tvale=1,has_vhnzc=1,has_vhtnzc=1)
        call pawpsp_17in(epsatm, ffspl, icoulomb, ipsp, hyb_mixing, ixc, lmax,&
                &       lnmax, pawpsp_header%mesh_size, mqgrid, mqgrid, pawpsp_header,&
                &       pawrad(1), pawtab(1), xcdev, qgrid_ff, qgrid_vl, usewvl, usexcnhat,&
                &       vlspl(:,:), xcccrc, xclevel, denpos, zion, znucl)
        call paw_setup_free(pawsetup)
        call paw_setup_free(paw_setuploc)

        !write(13,*) 'dij0',pawtab(1)%dij0

    !enddo

    mpsang = lmax + 1
    call pawinit(effmass_free, gnt_option,gsqcut_eff,hyb_range_fock,lcutdens,lmix,mpsang,nphi,nsym,ntheta,&
        &     pawang,pawrad,pawspnorb,pawtab,xcdev,xclevel,usepotzero)
    !write(13,*) 'eijkl',pawtab(1)%eijkl

    !See m_scfcv_core.f90, lines 669 and forth

    l_size_atm(1) = pawtab(1)%l_size
    typat(1) = 1 !typat stores the atomic number, I'm using H as test case so set to 1
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

    

    close(10)
    close(11)
contains
    subroutine scan_input_double(name_in,value,n)
        character(*)      :: name_in
        character(len=20) :: name
        integer           :: n
        real*8            :: value(n)

        read(11,*) name, value
        if(trim(name)/=trim(name_in)) then
            write(*,*) 'variable name does not match : ',trim(name),trim(name_in)
        endif
    end subroutine

    subroutine scan_input_int(name_in,value,n)
        character(*)      :: name_in
        character(len=20) :: name
        integer           :: n
        integer           :: value(n)

        read(11,*) name, value
        if(trim(name)/=trim(name_in)) then
            write(*,*) 'variable name does not match : ',trim(name),trim(name_in)
        endif
    end subroutine

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
end program
