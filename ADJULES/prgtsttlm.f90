!**************************************************************
!   file: 	        prgtsttlm.f90
!   purpose: 	        driver to check gradient of scalar valued function
!                       from multiple runs of product tangent linear * vector 
!                       against finite differences 
!                       and print timing information
!   creation date:	01 07
!
!   Copyright (C) 2000-2011
!   FastOpt, GmbH, Hamburg, Germany
!   http://FastOpt.com
!   All rights reserved.
!**************************************************************
program main
  implicit none
  integer :: n, m
! get the number of the independent and dependent variables
  call setfunc( n, m )
! trick for dynamic allocation
  call doit( n, m )
end program main

!**************************************************************
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
    real deltat
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

!**************************************************************
subroutine doit( n, m )
  use mo_timing
  implicit none
  integer n, m
  real :: x(n), y(m), x_tl(n), y_tl(m)
  real :: yh(n,m), fdy(m), dyfd(n,m), dyad(n,m), rel(n,m)
  character*(9) :: flag(n,m)
  character*(25) :: lineout, headout
  real, parameter ::epsmach = 1.E-20 
  real, parameter ::epsbeg = 1.E-4 
  integer :: i, j
  integer :: nbeg, nend, nstep
  real :: eps, acc, tlspeed, adspeed
  real :: xmemo, absmax
  real :: neps,epsstep ! dummies for compatibility of tst.par
  namelist /check/ eps, acc, nbeg, nend, nstep, tlspeed, adspeed, neps,epsstep
  integer :: iostat
  character*(*) cfile
  logical :: lexist
  parameter ( cfile = 'tst.par')
  ! for timing
  real ::    tfunc, tgrad, tforw, deltat
  integer :: ncfunc, ncgrad
  ! initialisations
  flag = ''
  lineout = '(2i6,7(1x,e19.12),a)'
  headout = '(2a6,7(1x,a19   ),a)'
  ! set/read parameters for check
  eps   = epsbeg
  nbeg  = 1
  nend  = n
  nstep = 1
  acc   = 1.e-3
  tlspeed = 2.0
  adspeed = 10.0
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
  ! initialise timing
  tfunc = 0.0
  tgrad = 0.0
  ncgrad = 0
  call inittiming
  ! initialise the model and set control variables
  call initfunc( n, x )
! get initial state if present
  inquire(file='x.b',exist=lexist)
  if (lexist) then
     OPEN (39,file='x.b',status='old',form='unformatted')
     read(39) x
     CLOSE(39)
     print*,'Point read from "x.b":',x(1:min(5,n))
  end if
  ! fd test
!  print *
!  print *, '**************************************************'
!  print '(A,E9.3)', ' CHECK OF TANGENT'
!  print *, '**************************************************'
!  write (*,headout) 'col i', 'row j','eps','x(i)','y2(j)','y(j)', 'delty(j)/eps','dy(j)/dx(i) ','rel diff'
  ! run model for initial values of control variables
  y=0.
  call btiming
  call func( n, x, m, y(:))
  call etiming(deltat)
  ncfunc = 1
  tfunc = tfunc + deltat
  ! loop over components of the gradient 
  do i = nbeg, nend, nstep
     ! initialise tangent linear: perturb in i-th direction
     x_tl = 0.0
     x_tl(i) = 1.0
     ! run TLM for initial values of control variables
     yh(i,:)=0.
     call btiming
     call func_tl( n, x, x_tl, m, yh(i,:), y_tl )
     call etiming(deltat)
     dyad(i,:) = y_tl
     ncgrad=ncgrad+1
     tgrad = tgrad + deltat
     ! checking function value of tangent
     if(i.eq.nbeg) then
        print *,'checking function value of tangent'
        print '(a6,4(1x,a19))','row','y(row) fu','y(row) ad','deltay'
        do j=1,m
           print '(i6,4(1x,e19.12)  )',j,y(j),yh(i,j),yh(i,j)-y(j)
        enddo
        call btiming
        call func_tl( n, x, x_tl, m, yh(i,:), y_tl )
        call etiming(deltat)
        dyad(i,:) = y_tl
        ncgrad=ncgrad+1
        tgrad = tgrad + deltat
        print *,'checking function value of tangent again'
        print '(a6,4(1x,a19))','row','y(row) fu','y(row) ad','deltay'
        do j=1,m
           print '(i6,4(1x,e19.12)  )',j,y(j),yh(i,j),yh(i,j)-y(j)
        enddo
        call btiming
        call func( n, x, m, y(:))
        call etiming(deltat)
        ncfunc = 1
        tfunc = tfunc + deltat
        print *,'checking function value of model after call of tangent'
        print '(a6,4(1x,a19))','row','y(row) fu','y(row) ad','deltay'
        do j=1,m
           print '(i6,4(1x,e19.12)  )',j,y(j),yh(i,j),yh(i,j)-y(j)
        enddo
     endif
     ! perturb one control variable
     xmemo = x(i)
     x(i) = x(i) + eps
     ! run model for perturbed values of control variables
     yh(i,:)=0.
     call btiming
     call func( n, x, m, yh(i,:))
     call etiming(deltat)
     ncfunc=ncfunc+1
     tfunc = tfunc + deltat
     x(i) = xmemo
     ! compute finite difference approximation
     fdy(:) = yh(i,:) - y(:)
     ! look over function components
     do j = 1, m
        ! compute y finite difference
        if (abs(fdy(j)) .lt. epsmach) then
           dyfd(i,j) = 0.0
        else
           dyfd(i,j) = fdy(j) / eps
        end if
        ! compute normalised difference
        rel(i,j) = abs(dyfd(i,j)-dyad(i,j))
        absmax = max( abs(dyfd(i,j)), abs(dyad(i,j)) )
        if (absmax .lt. epsmach.and.rel(i,j).gt.-huge(rel(i,j)).and.rel(i,j).lt.huge(rel(i,j))) then
           rel(i,j) = 0.0
        else
           rel(i,j) = rel(i,j) / absmax
           if (abs(rel(i,j)) .gt. acc) then
              flag(i,j) = 'DIFFERENT'
           endif
        end if
        ! print result of test
!        write (*,lineout) i,j,eps,x(i),yh(j),y(j),dyfd(i,j),dyad(i,j),rel(i,j),flag(i,j)
     enddo
  end do
  print *, '**************************************************'
  tfunc=tfunc/ncfunc
  tgrad=tgrad/ncgrad
  ! print timing 
  print *
  print *, '**************************************************'
  print *, '  TIMING OF TANGENT '
  print *, '**************************************************'
  print '(a,i8)'   , ' based on function calls      : ', ncfunc
  print '(a,i8)'   , ' based on tangent linear calls: ', ncgrad
  print '(a,f12.3)', ' run time function            : ', tfunc
  print '(a,f12.3)', ' run time function + tangent  : ', tgrad
  if (tfunc .gt. 0.0) then
     print '(2a,f12.3)', ' rel. run time forw mode'&
          , ' (FUNC+GRAD)/FUNC          : ', tgrad/tfunc
     if (tgrad/tfunc.gt.tlspeed) print*, 'TOO SLOW'
  endif
  print *, '**************************************************'
  ! repeat check output
  print *
  print *, '**************************************************'
  print '(A,E9.3)', ' CHECK OF TANGENT '
  print *, '**************************************************'
  write (*,headout) 'col i', 'row j','eps','x(i)','y2(j)','y(j)', 'delty(j)/eps','dy(j)/dx(i) ','rel diff'
  do i = nbeg, nend, nstep
     do j = 1, m
        write (*,lineout) i,j,eps,x(i),yh(i,j),y(j),dyfd(i,j),dyad(i,j),rel(i,j),flag(i,j)
     enddo
  enddo
  print *, '**************************************************'
  if (all(flag(:,:)=='')) then
     print *,'tsttlm: o.k.'
  else
     print *,'tsttlm: differences'
  end if
end subroutine doit
