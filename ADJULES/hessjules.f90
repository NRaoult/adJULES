!======================================================================
!   file:               hessjules.f90
!   purpose:            for optimisation using R
!                       allows compilation of a shared object file
!                       that can be called from within R
!======================================================================


!
! make sure these modules are defined and available for use by other subroutines
!


      module mo_h
      implicit none
      integer, parameter :: mh = 1001
      integer, parameter :: mx = 100
      integer ih
      real hx(mx,mh)
      real hf(mh)
      real hg(mx,mh)
      real lx(mx)
      end module mo_h



