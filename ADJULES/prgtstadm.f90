!**************************************************************
!   file: 	        prgtstadm.f90
!   purpose: 	        driver to check gradient from scalar 
!                       valued adjoint against finite 
!                       differences and print timing information
!   creation date:	06 02
!
!   Copyright (C) 2000-2008
!   FastOpt, GmbH, Hamburg, Germany
!   http://FastOpt.com
!   All rights reserved.
!**************************************************************
program prgtstadm
  implicit none
  integer :: n,m

n=0  ! initialise 
m=0

!==============================================================
! get the number of control variables
!==============================================================
  call setfunc( n, m )

!==============================================================
! do the actual computation, see subroutine below
!==============================================================
  call doit( n, m )

end program prgtstadm

module mo_timing
  integer isec0, isec, crate, isecmax
  integer dt
  integer iclock0(8), iclock(8)

contains

  subroutine inittiming
    implicit none
    call system_clock( count=isec0 )
    call system_clock( count=isec )
    dt = isec - isec0
    call system_clock( count_rate=crate, count_max=isecmax )
  end subroutine inittiming

  subroutine btiming
    implicit none
    call date_and_time( values=iclock0 )
    call system_clock( count=isec0 )
  end subroutine btiming
  
  subroutine etiming(deltat)
    implicit none
    real :: deltat
    call system_clock( count=isec )
    call date_and_time( values=iclock )
    if (iclock(3).eq.iclock0(3)) then
       deltat=float((isec-isec0) - dt)/crate
    elseif (iclock(3).gt.iclock0(3)) then
       isec=isec+(iclock(3)-iclock0(3))*isecmax
       deltat=float((isec-isec0) - dt)/crate   
    else
       deltat=0.0
       print*,'mo_timing: timing problem'
       print*,' maybe a new month has begun overnight?'
       print*,' to avoid crash set deltat = ',deltat
       print*,' isec0 = ', isec0
       print*,' isec = ', isec
       print*,' iclock0(3) = ', iclock0(3)
       print*,' iclock(3) = ', iclock(3)
    endif
  end subroutine etiming

end module mo_timing

      subroutine doit( n, m )
      use mo_timing
      implicit none
      integer n, m
      integer, parameter :: nmax=100000

      REAL   ady(m)
      REAL   x(n)
      REAL   adx(n)
      REAL   gfd(nmax,m)

      REAL        epsmach
      parameter ( epsmach = 1.E-20 )
      real ::     epsbeg
      parameter ( epsbeg = 1.E-4 )
      integer i
      REAL    y(m), yh(m), yh2(m), xmemo, rel(nmax,m), deltay(m), absmax
      REAL    eps, acc, tlspeed, adspeed
      integer nbeg, nend, nstep
      logical :: lexist
      real :: neps,epsstep ! dummies for compatibility of tst.par
      namelist /check/ eps, acc, nbeg, nend, nstep, tlspeed, adspeed, neps,epsstep
      character(7) cfile
      parameter ( cfile = 'tst.par')

      real :: tfunc, tgrad, tforw, deltat
      integer istat
      integer niter, ncfunc, ncgrad

      INTEGER nn
      INTEGER ifail, iostat
      LOGICAL cold
      INTEGER ifunc, j
      
      logical :: lsuccess = .true.

  real  y_ad(1), x_ad(n), x_tl(n), x_ad_tl(n), hess(n,n) !WORK 

!==============================================================
! parameters for check
!==============================================================
      eps   = epsbeg
      nbeg  = 1
      nend  = n
      nstep = 1
      acc   = 1.e-3
      adspeed = 10

      open (UNIT=87, FILE=cfile )
      read (UNIT=87, nml=check, iostat=iostat )
      close(UNIT=87)
      if (iostat==0) then
         print '(A)',' parameters for gradient check read:'
      else if (iostat<0) then
         print '(A)',' default parameters for gradient check:'
      else
         print '(A)','problem reading '//cfile
         stop 'problem reading parameter file in doit!'
      end if
      print check

      nbeg = max(1,nbeg)
      nbeg = min(n,nbeg)
      nend = max(1,nend)
      nend = min(n,nend)

      if (nend.gt.nmax) then
         print*, 'increase nmax = ',nmax,' to be bigger than nend =',nend
         stop
      endif
!======================================================
! initialisations for timing
!======================================================
      tfunc = 0.0
      tgrad = 0.0
      niter = 1
      ncgrad = 0
      call inittiming

!==============================================================
! initialise the model and set control variables
!==============================================================
      print*, 'starting initfunc'
! portability     call flush(6)
      call initfunc( n, x)
! get initial state if present
      inquire(file='x.b',exist=lexist)
      if (lexist) then
         OPEN (39,file='x.b')
         read(39,*) x
         CLOSE(39)
         print*,'Point read from "x.b":',x
      end if

!==============================================================
! compare components of gradient against finite differences
! and do timing for admodel and each call of the model
!==============================================================
! run model for initial values of control variables
      y=0.
      print*, 'starting func'
! portability     call flush(6)
      call btiming

      call func( n, x, m, y )
      call etiming(deltat)
      ncfunc = 1
      tfunc = tfunc + deltat

! initialise adjoint
      adx(:) = 0.0
      ady(:) = 0.0
      ady(1) = 1.0

! run adjoint for initial values of control variables
      call btiming
      call func_ad( n, x, adx, m, yh2, ady)
      call etiming(deltat)
      ncgrad = 1
      tgrad = tgrad + deltat
      print *,'checking function value of ADM'
      print '(a6,4(1x,a19))','row','y(row) fu','y(row) ad','deltay'
      do j=1,m
         print '(i6,4(1x,e19.12)  )',j,y(j),yh2(j),yh2(j)-y(j)
      enddo


hess = 0.0 ! initialise

write(*,*) ' calling hessian....'
do i=1,n
  x_tl     = 0.0
  x_tl(i)  = 1.0
  x_ad     = 0.0
  y        = 0.0 
  y_ad     = 1.0
  x_ad_tl  = 0.0
  
  call func_ad_tl( n, x, x_tl, x_ad, x_ad_tl, m, y, y_ad )

  call flush(6)
  
  hess(:,i) = x_ad_tl


 
end do

write(*,*) hess

write(*,*) '....done'








      print *
      print *, '**************************************************'
      print '(A,E11.3)', ' CHECK OF ADJOINT USING eps = ', eps
      print *, '**************************************************'
      print '(2a6,9(1x,a19))', 'col i', 'row j','eps','x(i)','y2(j)','y(j)', 'delty(j)/eps','dy(j)/dx(i) ','rel diff'

! loop over components of the gradient 
      do i = nbeg, nend, nstep

! perturb one control variable
         xmemo = x(i)
         x(i) = x(i) + eps

! run model for perturbed values of control variables
         call btiming
         call func( n, x, m, yh )
         call etiming(deltat)

         ncfunc=ncfunc+1
         tfunc = tfunc + deltat
         x(i) = xmemo
        
! compute finite difference approximation
         deltay = yh - y
         do j=1,m
            if (abs(deltay(j)) .lt. epsmach) then
               gfd(i,j) = 0.0
            else
               gfd(i,j) = deltay(j) / eps
            end if
! initialise adjoint
            adx(:) = 0.0
            ady(:) = 0.0
            ady(j) = 1.0
! run adjoint 
            call btiming
            call func_ad( n, x, adx, m, yh2, ady)
            call etiming(deltat)
            ncgrad = ncgrad + 1
            tgrad = tgrad + deltat
! compute normalised difference
            rel(i,j) = abs(gfd(i,j)-adx(i))
            absmax = max( abs(gfd(i,j)), abs(adx(i)) )
            if (absmax .lt. epsmach) then
               rel(i,j) = 0.
            else
               rel(i,j) = rel(i,j) / absmax
            end if

! print result of test
            if (rel(i,j) .lt. acc.and.rel(i,j).gt.-huge(rel(i,j)).and.rel(i,j).lt.huge(rel(i,j))) then
               print '(2i6,7(1x,e19.12)  )',i,j,eps,x(i),yh(j),y(j),gfd(i,j),adx(i),rel(i,j)
            else
               print '(2i6,7(1x,e19.12),a  )',i,j,eps,x(i),yh(j),y(j),gfd(i,j),adx(i),rel(i,j), &
                    ' DIFFERENT'
               lsuccess = .false.
            endif
         end do
      end do

      tfunc=tfunc/ncfunc
      tgrad=tgrad/ncgrad
      print *, '**************************************************'
      print *

!===============================================
! print timing 
!===============================================
      print *
      print *, '**************************************************'
      print *, '  TIMING OF ADJOINT '
      print *, '**************************************************'
      print '(a,i8)'   , ' based on function calls      : ', ncfunc
      print '(a,f12.3)', ' run time function            : ', tfunc
      print '(a,f12.3)', ' run time function + adjoint  : ', tgrad

      if (tfunc .gt. 0.) then
         print '(2a,f12.3)', ' rel. run time reverse mode'&
           , ' (FUNC+GRAD)/FUNC : ', tgrad/tfunc
         if (tgrad/tfunc.gt.adspeed) print*, 'TOO SLOW'
      endif
      print *, '**************************************************'
      print *

!==============================================================
! initialise the model and set control variables
!==============================================================
      call POSTFUNC(  n, x, m, y )
  if (lsuccess) then
     print *,'tstadm: o.k.'
  else
     print *,'tstadm: differences'
  end if
end subroutine doit
