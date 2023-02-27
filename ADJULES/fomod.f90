!*****************************************************************
!   file: 	        fomod.f90
!   purpose: 	        make main loop counter available globally.
!
!   creation date:	05 09
!
!   Copyright (C) 2009
!   FastOpt, GmbH, Hamburg, Germany
!   http://FastOpt.com
!   All rights reserved.
!*****************************************************************
    module fomod
    implicit none
    integer :: iloopcount ! main loop counter
    integer :: nloopcount ! upper limit of main loop
    integer :: call_fcdch_sea  ! holds local call number on repeated calls
    integer :: call_fcdch_land ! holds local call number on repeated calls
    integer :: call_sf_flux_land ! holds local call number on repeated calls
    integer :: call_sf_flux_land_max ! holds maximum local call number on repeated calls
    integer :: call_data_init_vals ! holds local call number on repeated calls
    integer :: spin_a_step         ! records value of a_step at end of spin-up

!Luke
    real,allocatable :: obsts(:,:) ! globally available observed time-series
    real,allocatable :: modts(:,:) ! globally available modelled time-series
    real :: window                 ! window length for averaging
    logical,allocatable :: monthstarts(:) ! Vector of when months start
!Luke

    logical :: l_smooth   ! if FALSE force use of max/min intrinsics
    logical :: l_mask     ! if TRUE mask nighttime fluxes (sw_down < 20)

    contains
!Luke function for smoothing MAX intrinsic
    function softmax(x,y,kk) result(smx)
       real, intent(in) :: x,y,kk  ! input
       real             :: smx ! output
       real             ::  mx !maximum
       real             :: k
!write(*,*) 'inside softmax',x,y,k
       mx = max(x,y)
       k = 50./ sqrt(sqrt(sqrt(tiny(1.)))+x**2+y**2)                               ! Jupp hack to reset k since need arguments of exp to be less than about 50
       if (l_smooth) then
 !       smx = (log(exp(k*(x-mx))+exp(k*(y-mx))))/k + mx      ! one way of doing it
         smx = (x*exp(k*x)+y*exp(k*y)) / (exp(k*x)+exp(k*y))  ! another way of doing it
       else
         smx = mx
       end if
!write(*,*) ' leaving softmax', smx
    end function softmax

!Luke function for smoothing MIN intrinsic
    function softmin(x,y,k) result(smn)
       real, intent(in) :: x,y,k ! input
       real             :: smn   ! output
    smn = -softmax(-x,-y,k)
    end function softmin
    

    end module fomod

    module mo_n
    implicit none
    integer :: n, m
    end module mo_n
