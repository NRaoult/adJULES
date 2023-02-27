!===============================================
! this subroutine provides a wrapper so that 
! it can be called from within R
!===============================================

      subroutine vfunc_wrapper_adm(x1,v)
      use mo_n
      use inout, ONLY : nts

      implicit none
      real   x1(n)
      real   v
      real   vfunc      
      external vfunc
      v = vfunc(x1)
      
      end subroutine


!===============================================
! this subroutine evaluates the function
!===============================================
real function vfunc( x)
  use mo_n

  use inout, only: weight

  implicit none
  real    x(n), y(m)
  real    fc
! call model

  y = 0.
  call func( n, x, m, y )

  fc = y(1)

! avoid stop due to overflow
  if(.not.fc.lt.huge(fc)/1.e50) then
     fc = huge(fc)/1.e50
     print*, 'PRGDFPMIN: avoiding function overflow'
     call flush(6)
  endif

  if(.not.fc.gt.-huge(fc)/1.e50) then
     fc = huge(fc)/1.e50
     print*, 'PRGDFPMIN: avoiding function overflow'
     call flush(6)
  endif

  vfunc = fc

end function vfunc
