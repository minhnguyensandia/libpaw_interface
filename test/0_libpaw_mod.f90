module libpaw_mod
      use m_pawpsp
      use m_pawtab
      use m_pawrad
      use m_pawxmlps
    
      character(264)           :: filename !name of PAW xml file
      type(paw_setup_t)        :: pawsetup
      type(pawpsp_header_type) :: pawpsp_header
      type(pawtab_type)        :: pawtab
      type(pawrad_type)        :: pawrad

      integer :: lloc, lmax, pspcod, pspxc
      integer :: ixc, xclevel
      real*8  :: hyb_mixing
      real*8  :: r2well, zion, znucl

      integer :: ngfft(3),ngfftdg(3)
      real*8  :: gmet(3,3) !reciprocal space lattice vectors
      
      real*8  :: gsqcut, gsqcutdg
      real*8,allocatable  :: qgrid_ff(:), qgrid_vl(:) !ff:'corase', vl:'fine'
      real*8,allocatable  :: ffspl(:,:,:), vlspl(:,:)

      real*8  :: epsatm
      real*8  :: xcccrc

      ! Some default values in ABINIT
      integer :: iboxcut = 0
      integer :: mqgrid = 3001, lnmax = 4, ipsp = 1
      integer :: usexcnhat = 0, xcdev = 1, usewvl = 0, icoulomb = 0
      real*8  :: denpos = 1d-14

end module
