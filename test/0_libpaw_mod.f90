module libpaw_mod
      use m_pawpsp
      use m_pawtab
      use m_pawrad
      use m_pawxmlps
      use m_pawang
      use m_pawfgrtab
      use m_paw_ij
      use m_paw_an
      use m_pawrhoij
    
      !In practice, this part will be obtained from
      !the program calling libpaw_interface
      real*8  :: ecut, ecutpaw !a coarse grid, and a fine grid for PAW
      real*8  :: gmet(3,3) !reciprocal space lattice vectors
      real*8  :: rprimd(3,3) !lattice vectors
      real*8  :: gprimd(3,3) !reciprocal space lattice vectors
      real*8  :: ucvol !volume of unit cell
      integer :: nspden = 1, nsppol = 1 !number of spin components for density and wavefunctions
      !for normal nspin = 1,2 calculations, nsppol is set to be same as nspden
      real*8  :: spinat(3,1) = 0.0 !no atomic magnetization for now

      character(264)           :: filename !name of PAW xml file
      type(paw_setup_t)        :: pawsetup
      type(pawpsp_header_type) :: pawpsp_header
      type(pawang_type)        :: pawang

      type(pawtab_type),    allocatable     :: pawtab(:) !set to 1 for now, will expand later
      type(pawrad_type),    allocatable     :: pawrad(:)
      type(pawfgrtab_type), allocatable     :: pawfgrtab(:)
      type(paw_ij_type),    allocatable     :: paw_ij(:)
      type(paw_an_type),    allocatable     :: paw_an(:)
      type(pawrhoij_type),  allocatable     :: pawrhoij(:)

      integer :: ntypat, natom
      integer, allocatable :: typat(:), l_size_atm(:), nattyp(:)
      integer, allocatable :: lexexch(:), lpawu(:) !no exact exchange, no dft+u

      real*8,  allocatable :: znucl(:), xred(:,:)
      integer, allocatable :: atindx1(:),atindx(:)

      integer :: lloc, lmax, pspcod, pspxc
      integer :: ixc, xclevel !functional; will be set externally in practice
      real*8  :: hyb_mixing, hyb_range_fock
      real*8  :: r2well, zion

      integer :: pawspnorb = 0, nspinor = 1 !not considering soc yet

      integer :: ngfft(3),ngfftdg(3)
      
      real*8  :: gsqcut, gsqcutdg, gsqcut_eff
      real*8,allocatable  :: qgrid_ff(:), qgrid_vl(:) !ff:'corase', vl:'fine'
      real*8,allocatable  :: ffspl(:,:,:), vlspl(:,:,:)

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
      integer :: cplex = 1

end module
