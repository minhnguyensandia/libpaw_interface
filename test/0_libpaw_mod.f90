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
      real*8  :: spinat(3,1) = 0.0 !no atomic magnetization for now

      type(paw_setup_t)        :: pawsetup
      type(pawpsp_header_type) :: pawpsp_header
      type(pawang_type)        :: pawang

      type(pawtab_type),    allocatable     :: pawtab(:) !set to 1 for now, will expand later
      type(pawrad_type),    allocatable     :: pawrad(:)
      type(pawfgrtab_type), allocatable     :: pawfgrtab(:)
      type(paw_ij_type),    allocatable     :: paw_ij(:)
      type(paw_an_type),    allocatable     :: paw_an(:)
      type(pawrhoij_type),  allocatable     :: pawrhoij(:)

      integer, allocatable :: l_size_atm(:), nattyp(:)
      integer, allocatable :: lexexch(:), lpawu(:) !no exact exchange, no dft+u

      real*8,  allocatable :: znucl(:), zion(:)
      integer, allocatable :: atindx1(:),atindx(:)

      integer :: lloc, lmax, pspcod, pspxc
      real*8  :: hyb_mixing, hyb_range_fock
      real*8  :: r2well

      integer :: pawspnorb = 0, nspinor = 1 !not considering soc yet
      
      real*8  :: gsqcut, gsqcutdg, gsqcut_eff
      real*8,allocatable  :: qgrid_ff(:), qgrid_vl(:) !ff:'corase', vl:'fine'
      real*8,allocatable  :: ffspl(:,:,:), vlspl(:,:,:)

      integer, allocatable :: fftn3_distrib(:), ffti3_local(:)
      integer, allocatable :: fftn2_distrib(:), ffti2_local(:) !For FFT parallelization

      real*8  :: xcccrc

      ! Some default values in ABINIT
      integer :: iboxcut = 0
      integer :: mqgrid = 3001, lnmax = 6, ipsp = 1
      integer :: usexcnhat = 0, xcdev = 1, usewvl = 0, icoulomb = 0, usepotzero = 0
      real*8  :: denpos = 1d-14
      integer :: gnt_option = 1, lcutdens = 10, lmix = 10
      integer :: mpsang, nphi = 13, ntheta = 12, nsym = 1 !the spherical grid
      real*8  :: effmass_free = 1.0
      integer :: cplex = 1

end module
