program line
  implicit none
  integer :: n, m

! get the number of control variables and dependents
  call setfunc( n, m )

! action
  call linerun( n, m )

end program line

subroutine linerun( n, m )
  implicit none

  integer :: n, m
  integer, parameter  :: uopt = 7        ! unit for options
  integer, parameter  :: usect = 22      ! unit for section
  integer, parameter  :: upar1 = 23      ! unit for parameter values 
  integer, parameter  :: upar2 = 24      ! unit for parameter values 
  real    :: fc, fc0, y(m)
  real    :: x0(n), xend(n), x(n), t, ndx



  integer :: nt                          ! number of steps in search direction
  integer :: it1,it2                     ! start and end point 
  real :: a                              ! scaling factor default: 1.

  integer :: i, it, ix
  namelist/nlist/a,nt,it1,it2

! check for too large n
  if (n.gt.99) print *, 'Line: n = ',n,' too large, change print format'

! default values for options
  a = 1. 
  nt = 10
  it1 = -2
  it2 = nt

! read namelist for options
      open( uopt, file   = 'line.par', access = 'sequential', form = 'formatted'  )
      read( uopt, nlist, end=9000, err=8000 )
      print *, ' Line : options have been read'
 8000 continue
 9000 continue
      close( uopt )

! initialize the model
  call initfunc( n, x0 )

! get test interval
  
  print*, 'Read previous point from line search '
  OPEN (39,file='lx.b',status='unknown',form='unformatted')
  READ(39) xend
  CLOSE(39)

  print*, 'Read final point '
  OPEN (39,file='xf.b',status='unknown',form='unformatted')
  READ(39) x0
  CLOSE(39)

! compute norm of delta x
  ndx =0.
  do ix=1,n
     ndx=ndx+(xend(ix)-x0(ix))**2
  enddo
  ndx = sqrt(ndx)

!  open units
  open (usect,file='line.set',form='formatted')
  open (upar1,file='parset.dat',form='formatted')
  open (upar2,file='parlist.dat',form='formatted')

!  xmgrace header
  write (usect,'(A)') '@WITH G00'
  write (usect,'(A)') '@G00 ON'
  write (usect,'(A)') '@TYPE xy'
  write (usect,'(A,E16.8)') '@WORLD XMAX', ndx
  write (usect,'(A)') '@s0 symbol 1'
  write (usect,'(A)') '@s0 line type 0'
  write (usect,'(A)') '@xaxis  label "norm of difference of control vector to base point"'
  write (usect,'(A)') '@TITLE "cost function"'

  write (upar1,'(A,A19,99I20)') '#', 't', (i,i=1,n)

  fc0 = 0.
! evaluate at end point of optimisation
  IF (it2-it1.GT.1) y = 0.
  IF (it2-it1.GT.1) call func( n, x0, m, y)
  IF (it2-it1.GT.1) fc0 = y(1)

! this is the loop that does the section
  do it=it1,it2
     t=a*float(it)/nt
     x=x0+t*(xend(3)-x0(3))
     y = 0.
     call func( n, x, m, y)
     fc = y(1)
     write (usect,*) t*ndx,fc-fc0
     write (upar1,'(100e20.6)') t*ndx,x(:)
     write (upar2,*) '# parameter vector for it = ',it
     do i=1,n
        write (upar2,*) i, x(i)
     enddo
  enddo

!  close units
  close (usect)
  close (upar1)
  close (upar2)
end subroutine linerun
