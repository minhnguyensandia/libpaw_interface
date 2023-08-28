subroutine set_rhoij(iatom,nrhoijsel,size_rhoij,nspden,rhoijselect,rhoijp)
    use libpaw_mod
    implicit none

    integer :: nrhoijsel,iatom,size_rhoij,nspden
    integer :: rhoijselect(size_rhoij)
    real*8  :: rhoijp(size_rhoij,nspden)

    pawrhoij(iatom)%nrhoijsel = nrhoijsel
    pawrhoij(iatom)%rhoijselect = rhoijselect
    pawrhoij(iatom)%rhoijp = rhoijp
end subroutine