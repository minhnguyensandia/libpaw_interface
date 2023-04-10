!!From 56_recipspace/m_kg.F90, but only the getcut subroutine
!!note that size of ngfft has been set to 3, because only nx, ny, nz are used

#include "libpaw.h"

!!****m* ABINIT/m_fftcore
!! NAME
!!  m_fftcore
!!
!! FUNCTION
!!  Low-level tools for FFT (sequential and MPI parallel version)
!!  It also provides helper functions to set up the list of G vectors
!!  inside a sphere or to count them.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2022 ABINIT group (SG, XG, AR, MG, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  1) Pass distribfft instead of MPI_enreg to simplify the API and facilitate code-reuse.
!!
!!  2) Merge this module with m_distribfft
!!
!!  3) Get rid of paral_kgb and MPI_type! This is a low-level module that may be called by other
!!     code in which paral_kgb is meaningless! FFT tables and a MPI communicator are sufficient.
!!
!! SOURCE


module m_fftcore

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING
 implicit none

 private

 public :: bound               ! Find distance**2 to boundary point of fft box nearest to kpt

contains
!!***

!!****f* m_fftcore/bound
!! NAME
!! bound
!!
!! FUNCTION
!! For given kpt, ngfft, and gmet,
!!  Find distance**2 to boundary point of fft box nearest to kpt
!!  Find distance**2 to boundary point of fft box farthest to kpt
!!
!! INPUTS
!!  kpt(3)=real input k vector (reduced coordinates)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gmet(3,3)=reciprocal space metric (currently in Bohr**-2)
!!
!! OUTPUT
!!  dsqmax=maximum distance**2 from k to boundary in Bohr**-2.
!!  dsqmin=minimum distance**2 from k to boundary in Bohr**-2.
!!  gbound(3)=coords of G on boundary (correspnding to gsqmin)
!!  plane=which plane min occurs in (1,2, or 3 for G1,etc).
!!
!! NOTES
!! Potential trouble: this routine was written assuming kpt lies inside
!! first Brillouin zone.  No measure is taken to fold input kpt back
!! into first zone.  Given arbitrary kpt, this will cause trouble.
!!
!! SOURCE

subroutine bound(dsqmax,dsqmin,gbound,gmet,kpt,ngfft,plane)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: plane
 real(dp),intent(out) :: dsqmax,dsqmin
!arrays
 integer,intent(in) :: ngfft(3)
 integer,intent(out) :: gbound(3)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: i1,i1min,i2,i2min,i3,i3min
 real(dp) :: dsm,dsp
 character(len=500) :: msg

! *************************************************************************

!Set plane to impossible value
 plane=0

!look at +/- g1 planes:
 dsqmax=zero
 dsqmin=dsq(ngfft(1)/2,-ngfft(2)/2,-ngfft(3)/2,gmet,kpt)+0.01_dp
 do i2=-ngfft(2)/2,ngfft(2)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(ngfft(1)/2, i2, i3,gmet,kpt)
     dsm = dsq( - ngfft(1)/2, i2, i3,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min =  - ngfft(1)/2
       i2min = i2
       i3min = i3
       plane=1
     end if
   end do
 end do
!
!+/- g2 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i3=-ngfft(3)/2,ngfft(3)/2
     dsp = dsq(i1,ngfft(2)/2,i3,gmet,kpt)
     dsm = dsq(i1,-ngfft(2)/2,i3,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = ngfft(2)/2
       i3min = i3
       plane=2
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min =  - ngfft(2)/2
       i3min = i3
       plane=2
     end if
   end do
 end do
!
!+/- g3 planes:
 do i1=-ngfft(1)/2,ngfft(1)/2
   do i2=-ngfft(2)/2,ngfft(2)/2
     dsp = dsq(i1,i2,ngfft(3)/2,gmet,kpt)
     dsm = dsq(i1,i2,-ngfft(3)/2,gmet,kpt)
     if (dsp>dsqmax) dsqmax = dsp
     if (dsm>dsqmax) dsqmax = dsm
     if (dsp<dsqmin) then
       dsqmin = dsp
       i1min = i1
       i2min = i2
       i3min = ngfft(3)/2
       plane=3
     end if
     if (dsm<dsqmin) then
       dsqmin = dsm
       i1min = i1
       i2min = i2
       i3min =  - ngfft(3)/2
       plane=3
     end if
   end do
 end do

 if (plane==0) then
!  Trouble: missed boundary somehow
   write(msg, '(a,a,a,3f9.4,a,3(i0,1x),a,a,a,a,a)' )&
   'Trouble finding boundary of G sphere for',ch10,&
   'kpt=',kpt(:),' and ng=',ngfft(1:3),ch10,&
   'Action : check that kpt lies',&
   'reasonably within first Brillouin zone; ',ch10,&
   'else code bug, contact ABINIT group.'
   ABI_BUG(msg)
 end if

 gbound(1)=i1min
 gbound(2)=i2min
 gbound(3)=i3min

 contains

   function dsq(i1,i2,i3,gmet,kpt)

     integer :: i1,i2,i3
     real(dp) :: dsq
     real(dp) :: kpt(3),gmet(3,3)

     dsq=gmet(1,1)*(kpt(1)+dble(i1))**2&
&      +gmet(2,2)*(kpt(2)+dble(i2))**2&
&      +gmet(3,3)*(kpt(3)+dble(i3))**2&
&      +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
&      +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
&      +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))
   end function dsq

end subroutine bound
!!***

END MODULE m_fftcore
!!***
