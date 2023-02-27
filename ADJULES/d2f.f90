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

  subroutine d2func_wrapper( x, hess )
  use mo_n
  use mo_h
  use inout, only : echo, nts, weight, xvary
  implicit none
!
! Description:
!   calculates hessian hess of cost function y at parameter vector x for the adjoint model
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
  real    hess(n,n)                   !OUT hessian of cost 
 

! local variables
!
  real y(1), y_ad(1), x_ad(n), x_tl(n), x_ad_tl(n) !WORK  check these

  character*23 outfile,itername       !WORK
  integer i, icount                   !WORK
  save icount

  call setfunc(n,m) !  set m,n    

! call tangent 
  icount=icount+1
      
! initialise adjoint

hess = 0.0 ! initialise

do i=1,n
  x_tl     = 0.0
  x_tl(i)  = 1.0
  x_ad     = 0.0
  y        = 0.0 
  y_ad     = 1.0
  x_ad_tl  = 0.0
  
  if (xvary(i)) call func_ad_tl( n, x, x_tl, x_ad, x_ad_tl, m, y, y_ad )

  call flush(6)
  
  hess(:,i) = x_ad_tl 
end do

  end subroutine
