!======================================================================
!   file:               prgopti.f90
!   purpose:            driver for optimisation using dpfmin
!======================================================================


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

      program main   ! program
!*************************************************************
      use mo_n
      implicit none

!==============================================================
! get the number of control variables
!==============================================================
      call setfunc( n, m )

!-----------------------------------------
! call the subroutine
!-----------------------------------------
      call opti

      end program

      subroutine opti
      use mo_n
      implicit none

      REAL    objf

      REAL   x(n), y(m)

      REAL   adx(n)
      integer ifail

!==============================================================
! initialisize the model
! and set the first guess of the control variables
!==============================================================
      call initfunc( n, x)
      ifail=0
      OPEN (39,file='output_FO/xf.b',status='old',form='unformatted',iostat=ifail)
      if (ifail.gt.0) then
         print*,' OPTI : parameter values set by initmod'
      else
         read(39) x
         print*,' OPTI : parameter values read from x.b'
      endif
      CLOSE(39)
      call flush(6)

!==============================================================
! optimize the costfunction using control variables
!==============================================================
      call optimum( n, x, objf, adx )
      y(1) = objf

!==============================================================
! start of postprocessor
!==============================================================
! save optimim
      OPEN (39,file='xf.b',status='unknown',form='unformatted')
      WRITE(39) x
      CLOSE(39) 

! save optimim
      OPEN (39,file='xf.dat',status='unknown',form='formatted')
      WRITE(39,'(20e26.16)') x
      CLOSE(39)

! save final gradient
      OPEN (39,file='gf.b',status='unknown',form='unformatted')
      WRITE(39) adx
      CLOSE(39) 

      call postfunc( n, x, m , y)

      end

      subroutine optimum( n, x, fc, adx )
      use mo_h
      implicit none

!==============================================================
! declare call parameters
!==============================================================
      integer n
      REAL    x(n), x_true(n)
      REAL    fc
      REAL    adx(n)

!==============================================================
! declare local variables
!==============================================================
      integer    ioptns, iter, i, j, unith
      parameter( ioptns = 87 )
      parameter( unith = 88 )
      REAL epsg
      external  vfunc,dfunc
      logical :: cook = .false.
      namelist/ OPTIM/ epsg

      epsg  = 1.E-2

!======================================================
! read new optimisation switches from file
!======================================================
      open( ioptns, FILE   = 'opti.par'&
                 , ACCESS = 'SEQUENTIAL'&
                 , FORM   = 'FORMATTED'  )
      read( ioptns, OPTIM, end=9000, ERR=8000 )
      print *, ' OPTI : options have been read'
 8000 continue
 9000 continue
      close( ioptns )
      print *, ' OPTI : Running with epsg = ', epsg
      call flush(6)

!==============================================================
! perform optimisation
!==============================================================
      ih = 0 ! our own count of iterations
      hg = 0. 
      hf = 0.
      hx = 0.
      lx = 0.
      call dfpmin(x,n,epsg,iter,fc,vfunc,dfunc)
      print*,' OPTI : exit after ',iter,' iterations.'
      call flush(6)
      fc = hf(ih)
      adx = hg(1:n,ih) ! save final gradient

!==============================================================
! output convergence overview
!==============================================================
      open (unith,file='h.dat',form='formatted')
      do i=1,ih
         write(unith,'(i5,2(3x,e12.6))')i,hf(i),sum(sqrt(hg(:,i)**2))
      enddo
      close (unith)

      end
