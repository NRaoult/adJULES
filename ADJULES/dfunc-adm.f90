! *****************************COPYRIGHT*******************************
! (c) University of Exeter/FastOpt 2011. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237] 
! *****************************COPYRIGHT*******************************

! evaluates the gradient using the adjoint model

  subroutine dfuncadm( x, grady )
  use mo_n
  use mo_h
  use inout, only : echo, nts, weight
  implicit none
!
! Description:
!   calculates gradient grady of cost function y at parameter vector x for the adjoint model
!
! Method:
!   
!  calls func_ad
!
! Current Code Owner: Tim Jupp, t.e.jupp@exeter.ac.uk
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!

! Subroutine arguments
! 
  real    x(n)                        !IN parameter vector
  real    grady(n)                    !OUT gradient of cost 

! local variables
!
  real y(m), yh2(m), adx(n), ady(m)   !WORK
  character*23 outfile,itername       !WORK
  integer i, icount                   !WORK

  save icount
  call setfunc(n,m) !  set m,n    

! call tangent 
  icount=icount+1
      
! initialise adjoint
  adx(:)   = 0.0
  ady(:)   = 0.0
  ady(1)   = 1.0
  y        = 0.0 
  yh2      = 0.0
           
  call func_ad(  n, x, adx, m, yh2, ady)

!      print *, ' OPTI : derivative evaluated', ih + 1
  call flush(6)
  grady(1:n) = adx(1:n) 
      
!      write(*,*) 'grady is', m,n,grady(:)
!      print*, 'done'
! avoid stop due to overflow
  do i = 1,n
     if(.not.grady(i).lt.huge(grady(i))) then
        print*, 'PRGDFPMIN: avoiding derivative overflow for component ',i,grady(i)
        grady(i) = huge(grady(i))/1.e50
        print*, 'reset to ',grady(i)
     endif
     if(.not.-grady(i).lt.huge(grady(i))) then
        print*, 'PRGDFPMIN: avoiding derivative underflow for component ',i,grady(i)
        grady(i) = -huge(grady(i))/1.e50
        print*, 'reset to ',grady(i)
     endif
  enddo
! store current point 
!  ih = ih + 1
!  do i = 1,n
!     hx(i,ih)=x(i)
!     hg(i,ih)=grady(i)
!  enddo
!  hf(ih)=yh2(1)
    if (echo) write(*,*) 'Opti: ',ih,hf(ih),sqrt(sum(hg(:,ih)**2))
    if (echo) write(itername,'(i18)') ih 
    if (echo) outfile='output_FO/opti-'//trim(adjustl(itername))//'.b'
    if (echo) OPEN (39,file=outfile,status='unknown',form='unformatted')
    if (echo)  WRITE(39) x
    if (echo) CLOSE(39)
    if (echo) outfile='output_FO/opti-'//trim(adjustl(itername))//'.dat'
    if (echo) OPEN (39,file=outfile,status='unknown',form='formatted')
    if (echo) WRITE(39,'(10(1x,e12.6))') x
    if (echo)  CLOSE(39)
  end subroutine
