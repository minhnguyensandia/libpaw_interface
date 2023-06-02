module libpaw_mod
      use m_pawpsp
      use m_pawtab
      use m_pawrad
      use m_pawxmlps
      use m_pawang
    
      !In practice, this part will be obtained from
      !the program calling libpaw_interface
      real*8  :: ecut(1), ecutpaw(1) !a coarse grid, and a fine grid for PAW
      real*8  :: gmet(3,3) !reciprocal space lattice vectors

      character(264)           :: filename !name of PAW xml file
      type(paw_setup_t)        :: pawsetup
      type(pawpsp_header_type) :: pawpsp_header
      type(pawtab_type)        :: pawtab(1) !set to 1 for now, will expand later
      type(pawrad_type)        :: pawrad(1)
      type(pawang_type)        :: pawang

      integer :: lloc, lmax, pspcod, pspxc
      integer :: ixc, xclevel !functional; will be set externally in practice
      real*8  :: hyb_mixing, hyb_range_fock
      real*8  :: r2well, zion, znucl
      integer :: pawspnorb = 0 !not considering soc yet

      integer :: ngfft(3),ngfftdg(3)
      
      real*8  :: gsqcut, gsqcutdg, gsqcut_eff
      real*8,allocatable  :: qgrid_ff(:), qgrid_vl(:) !ff:'corase', vl:'fine'
      real*8,allocatable  :: ffspl(:,:,:), vlspl(:,:)

      real*8  :: epsatm
      real*8  :: xcccrc

      ! Some default values in ABINIT
      integer :: iboxcut = 0
      integer :: mqgrid = 3001, lnmax = 4, ipsp = 1
      integer :: usexcnhat = 0, xcdev = 1, usewvl = 0, icoulomb = 0, usepotzero = 0
      real*8  :: denpos = 1d-14
      integer :: gnt_option = 1, lcutdens = 10, lmix = 10
      integer :: mpsang, nphi = 13, ntheta = 12, nsym = 1 !the spherical grid
      real*8  :: effmass_free = 1.0

end module
