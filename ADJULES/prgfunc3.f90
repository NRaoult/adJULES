!**************************************************************
!   Copyright (C) 2000, 2001, 2010
!   FastOpt GmbH, Hamburg, Germany
!
!   Contact: http://www.FastOpt.com
!
!   All rights reserved.
!**************************************************************

      program main
!*************************************************************
      implicit none
!=========================================
! declaration
!=========================================
      integer n, m

!==============================================================
! get the number of the independent and dependent variables
!==============================================================
      call setfunc( n, m )
!-----------------------------------------
! call the subroutine
!-----------------------------------------
      call prgfunction( n, m )

      end program main


      subroutine prgfunction( n, m )
!==============================================================
! declaration
!==============================================================
      implicit none
      
      integer n, m
      
      
      real :: x(n)
      real :: Y(m),y1,y2,y3,y4,eps
      integer :: nmax, mmax, i, my_pe, stdout

!==============================================================
! allocate memory for the control variables
! and the gradient vector
!==============================================================

      nmax = n
      mmax = m
      my_pe = 0
      stdout = 6
      eps = 1.e-6
      
!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================

      call initfunc( n, x )

!==============================================================
! compute function twice
!==============================================================
      call func( n, x, m, y )
      write(stdout,"(/,a,/)") 'after first function call :'
      print 722,sum(y)
      y1= sum(y)
      y=0.
      call func( n, x, m, y )
      write(stdout,"(/,a,/)") 'after second function call :'
      print 722,sum(y)
      y2= sum(y)
      y=0.
      write(stdout,"(/,a,e12.6)") 'perturbing x by :',eps
      x = x + eps
      call func( n, x, m, y )
      write(stdout,"(/,a,/)") 'after function call w perturbed x:'
      print 722,sum(y)
      y3= sum(y)
      y=0.
      call func( n, x, m, y )
      write(stdout,"(/,a,/)")                                            &
     &     'after second function call w perturbed x:'
      print 722,sum(y)
      y4= sum(y)

!===============================================
! compare
!===============================================
      eps = 10*epsilon(eps)

      write(*,*) ' cost2 with eps = ',eps
      if (abs(y1-y2).lt.abs(eps*(1+y1))) then
         write(*,*) ' cost2: o.k.'
      else
         write(*,*) ' cost2: DIFFERENT'
      endif
      if (abs(y3-y4).lt.abs(eps*(1+y3))) then
         write(*,*) ' cost2+eps: o.k.'
      else
         write(*,*) ' cost2+eps: DIFFERENT'
      endif

 722  format ("sum: ",e38.30)
!==============================================================
! postprocessing set
!==============================================================
!      call postfunc( n, x, m, y)

      end subroutine prgfunction

