                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine func( n, x, m, y ) 
                                                                        
      use keep_init 
      use observation 
                                                                        
! edithere begins                                                       
! get access to parameters                                              
      use SURF_PARAM, only  :  n_its, iter, q10_leaf, q10_soil !Luke 
      use pftparm, only     :  alpha, nl0, f0, tlow, tupp, rootd_ft, dcatch_dlai, dqcrit !Luke
      use prognostics, only : lai, canht_ft, cs 
      use p_s_parms, only   : b,satcon,albsoil,smvccl,smvcwt,sathh,smvcst,hcap,hcon 
! edithere ends                                                         
                                                                                                                                         
      use fomod, only      : nloopcount,modts,obsts,window,monthstarts
      use soil_param, only : mmax 
      USE Ancil_info, ONLY : sm_levels,ntiles,land_pts 
      use switches, only   : ilayers 
      use descent , only   : iter_eq 
      USE timeConst, ONLY  : iSecInDay 
      USE veg_io_vars, ONLY : nvegVarMax 
      use spin_mod, only    : nspinvarmax
      use nstypes, only     : ntype
      use switches, only    : ilayers
      use inout, only : ts_bias, ts_scal, nts, echo, PM,xoff ! Luke

      use drive_io_vars, only: driveTimeIndex !Luke

      use time_loc, only: date


                                                                       
      implicit none 
                                                                        
      integer :: n, m, i, j,k ,nwin !Luke added counters i,j,k
      real    :: x(n), y(m), u, v, p(n), parcost ,temp,tcost,temp2,temp3,ks,kc,kc2,tt,ty !Luke added variables for parcost and tcost 

      real    :: smvcst_old(sm_levels), smvccl_old(sm_levels) ,smvcwt_old(sm_levels)!Luke for new definition of soil moisture as frac. of sat.               
!      real, parameter :: window=1.0

      real,allocatable    :: newmodts(:,:),newobsts(:,:)

      integer :: dateYear, dateMonth




!     local                                                            
                                                                        
                                                                        
!$openad INDEPENDENT(X)                                                                        
                                                                        
!#define 2 1                                                            
                                                                                                                                                             
!$TAF INIT tape_n            = static, nloopcount                                  

!$TAF INIT tape_x            = static, 1                                  

!$TAF INIT tape_n_sfe        = static, nloopcount                                  

!$TAF INIT tape_n_ta         = static, nloopcount                                  

!$TAF INIT tape_n_physiol      = static, nloopcount
!$TAF INIT tape_n_npft_physiol = static, nloopcount*npft

!$TAF INIT tape_n_control0   = static, nloopcount
!$TAF INIT tape_n_control1   = static, nloopcount
!$TAF INIT tape_n_control2   = static, nloopcount
!$TAF INIT tape_n_control3   = static, nloopcount
!$TAF INIT tape_n_control4   = static, nloopcount
!$TAF INIT tape_n_control5   = static, nloopcount
!$TAF INIT tape_n_control6   = static, nloopcount
!$TAF INIT tape_n_control7   = static, nloopcount
!$TAF INIT tape_n_control8   = static, nloopcount
!$TAF INIT tape_n_control9   = static, nloopcount
!$TAF INIT tape_n_control10  = static, nloopcount
!$TAF INIT tape_n_control11  = static, nloopcount

!$TAF INIT tape_n_sfstom     = static, nloopcount

!$TAF INIT tape_as           = static, land_pts*ntype*nloopcount

!$TAF INIT tape_albp         = static, ilayers*land_pts*2*npft*nloopcount

!$TAF INIT tape_ta           = static, nloopcount * ntype
!$TAF INIT tape_ta1          = static, 4*land_pts*ntype*nloopcount
!$TAF INIT tape_ta2          = static, land_pts*ntype*nloopcount

!$TAF INIT tape_n_ta0        = static, nloopcount
!$TAF INIT tape_n_ta1        = static, nloopcount
!$TAF INIT tape_n_ta2        = static, nloopcount
!$TAF INIT tape_n_ta3        = static, nloopcount
!$TAF INIT tape_n_ta4        = static, nloopcount

!$TAF INIT tape_veg1a        = static, nloopcount*land_pts*npft

!$TAF INIT tape_n_np_spin    = static, land_pts*nspinvarmax*nloopcount

!$TAF INIT tape_n_nt         = static, ntiles * nloopcount   

!$TAF INIT tape_n_compete    = static, nloopcount                                  
!$TAF INIT tape_n_vegcarb    = static, nloopcount   

!$TAF INIT tape_n_phenol1    = static, nloopcount  
                               
!$TAF INIT tape_n_lotka      = static, npft * nloopcount    
!$TAF INIT tape_n_triffid    = static, npft * nloopcount    
!$TAF INIT tape_n_triffid1   = static, nloopcount    

!$TAF INIT tape_n_m          = static, nloopcount*mmax   
                        
!$TAF INIT tape_n_nshyd      = static, nloopcount*sm_levels                  
!$TAF INIT tape_n_soilcarb   = static, nloopcount*sm_levels       ! Jupp 

!$TAF INIT tape_rf      = static, nloopcount*sm_levels                   
               
!$TAF INIT tape_n_sfstom0    = static, iter * ilayers * nloopcount            
!$TAF INIT tape_n_sfstom1    = static, iter * npft * nloopcount 
!$TAF INIT tape_n_sfstom2    = static, iter * npft * nloopcount 
  
!$TAF INIT tape_veg2a        = static,           nloopcount/(INT(REAL(iSecInDay)*TRIFFID_PERIOD/TIMESTEP))                   
!$TAF INIT tape_veg2b        = static, iter_eq * nloopcount/(INT(REAL(iSecInDay)*TRIFFID_PERIOD/TIMESTEP))       
            
!$TAF INIT tape_n1           = static,  1 + nloopcount     
                        
!$TAF INIT tape_nv           = static, (1 + nloopcount)*nvegVarMax   
             
!$TAF INIT tape_spin         = static,  1 + nspin                             
  
!$TAF INIT tape_calczw       = static,  nloopcount * soil_pts * 3                              
!$TAF INIT tape_calczw1     = static,  nloopcount * soil_pts 		!Luke

!$TAF INIT tape_n_soil     = static,  nloopcount * soil_pts 		!Luke

!$TAF INIT tape_n_gauss      = static, nloopcount*soil_pts*sm_levels       ! Jupp 
!$TAF INIT tape_n_gauss1     = static, nloopcount*sm_levels       ! Jupp

!$TAF INIT tape_soilhyd      = static, nloopcount*sm_levels*land_pts                
!$TAF INIT tape_soilhyd1     = static, nloopcount*sm_levels*land_pts                

!$TAF INIT tape_surfhyd      = static, ntiles*nloopcount*land_pts 

!$TAF INIT tape_calcbfj      = static, nloopcount*soil_pts       ! Luke

!$TAF INIT tape_snow      = static, nloopcount*land_pts       ! Luke
!$TAF INIT tape_snow1     = static, nloopcount*ntiles		!Luke
!$TAF INIT tape_snow2     = static, nloopcount*ntiles*nsmax		!Luke

!$TAF INIT tape_rlsnow1      = static, nloopcount*land_pts*(nsmax+1)       ! Luke
!$TAF INIT tape_rlsnow2      = static, nloopcount*land_pts*(nsmax+1)*nsmax       ! Luke

!$TAF INIT tape_com_snow      = static, nloopcount*land_pts*nsmax       ! Luke
!$TAF INIT tape_snowgrain      = static, nloopcount*land_pts*nsmax       ! Luke
!$TAF INIT tape_snowpack      = static, nloopcount*land_pts*nsmax       ! Luke
!$TAF INIT tape_tridag      = static, nloopcount*nsmax       ! Luke

!$TAF INIT tape_sfimpl      = static, nloopcount*nice       ! Luke
!$TAF INIT tape_sfimpl1      = static, nloopcount*ntiles       ! Luke

!$TAF INIT tape_sfevap      = static, nloopcount*ntiles       ! Luke

!$TAF INIT tape_fcdch       = static, 1		!Luke
!$TAF INIT tape_sfexch       = static, nloopcount*rows*row_length*5			!Luke

!$TAF INIT tape_leaf      = static, nloopcount*land_pts       ! Luke

!$TAF INIT tape_sfstom        = static, nloopcount*ilayers    !Luke
!$TAF INIT tape_sfstom1        = static, 1    !Luke
!$TAF INIT tape_sfstom2    = static, nloopcount*ilayers*iter*land_pts	!Luke            
!$TAF INIT tape_sfstom3    = static, nloopcount*ilayers*iter	!Luke
!$TAF INIT tape_sfstom4    = static, nloopcount*iter*land_pts	!Luke
!$TAF INIT tape_sfstom5    = static, nloopcount*iter	!Luke                    

!$TAF INIT tape_update     = static, nloopcount*driveTimeIndex(2) !Luke

!    !$TAF INIT tape_n_fcdchs     = static, nloopcount*n_its*3              
!    !$TAF INIT tape_n_fcdchl     = static, ntiles*nloopcount*n_its 
!    !$TAF INIT tape_n_phenol     = static, nloopcount 
!    !$TAF INIT tape_n_sfstom     = static, iter * npft * nloopcount         
!    !$TAF INIT tape_n_sfflux     = static, ntiles * nloopcount  

!$TAF INIT tape_iter=static,1 !Luke

!$TAF INIT tape_nts=static,nts !Luke
!$TAF INIT tape_n_nts=static,nts*nloopcount !Luke
!$TAF INIT tape_win_nts=static,nts*12 !Luke
                                                                      
!                                                                       
! begin executable statements                                           
!     
      tcost=0.0            !initialise tcost                                                                  
      parcost=0.0          !initialise parcost                                                    
!                                                                       
!     reset the initial state                                           
!                                                                       

      if (echo) write(*,*) 'entering func at top.f90 line 140'
      call reinit 
      if (echo) write(*,*) 'top.f90 line 142 finished call to reinit'


      dateYear = date / 10000
      dateMonth = ( date - dateYear*10000 ) / 100 ! will be useful for monthly LAI


                                                                        
!                                                                       
!     map control variables                                             
!                                                                       
                                                                        
! edithere begins                                                       
                                                                        
      q10_leaf    = x(1) 
      nl0(1)      = x(2) 
      nl0(2)      = x(3) 
      nl0(3)      = x(4) 
      nl0(4)      = x(5) 
      nl0(5)      = x(6) 
      alpha(1)    = x(7) 
      alpha(2)    = x(8) 
      alpha(3)    = x(9) 
      alpha(4)    = x(10) 
      alpha(5)    = x(11) 
      f0(1)       = x(12) 
      f0(2)       = x(13) 
      f0(3)       = x(14) 
      f0(4)       = x(15) 
      f0(5)       = x(16) 
      tlow(1)     = x(17) 
      tlow(2)     = x(18) 
      tlow(3)     = x(19) 
      tlow(4)     = x(20) 
      tlow(5)     = x(21) 
      tupp(1)     = x(22) 
      tupp(2)     = x(23) 
      tupp(3)     = x(24) 
      tupp(4)     = x(25) 
      tupp(5)     = x(26) 

      lai(1,1)    = x(27)     ! multiply by month factors * x(94+dateMonth)
      lai(1,2)    = x(28) 
      lai(1,3)    = x(29) 
      lai(1,4)    = x(30) 
      lai(1,5)    = x(31) 

      canht_ft(1,1) = x(32) 
      canht_ft(1,2) = x(33) 
      canht_ft(1,3) = x(34) 
      canht_ft(1,4) = x(35) 
      canht_ft(1,5) = x(36) 
      satcon(1,1) = x(37) 
      satcon(1,2) = x(38) 
      satcon(1,3) = x(39) 
      satcon(1,4) = x(40) 

      b(1,1)      = x(41) 
      b(1,2)      = x(42) 
      b(1,3)      = x(43) 
      b(1,4)      = x(44) 

      cs(1,1)     = x(45) 
      rootd_ft(1) = x(46) 
      rootd_ft(2) = x(47) 
      rootd_ft(3) = x(48) 
      rootd_ft(4) = x(49) 
      rootd_ft(5) = x(50) 
      ts_bias(1)  = x(51)
      ts_bias(2)  = x(52)
      ts_bias(3)  = x(53)
      ts_bias(4)  = x(54)
      ts_scal(1)  = x(55)
      ts_scal(2)  = x(56)
      ts_scal(3)  = x(57)
      ts_scal(4)  = x(58)
      albsoil(1)  = x(59)

      smvcst(1,1)  = x(82)        !Jupp  !Luke perturb saturation
      smvcst(1,2)  = x(83)        !Jupp
      smvcst(1,3)  = x(84)        !Jupp
      smvcst(1,4)  = x(85)        !Jupp

      smvccl(1,1) = smvcst(1,1)  *  x(60) 	!Luke perturb critical
      smvccl(1,2) = smvcst(1,2)  *  x(61) 	
      smvccl(1,3) = smvcst(1,3)  *  x(62) 
      smvccl(1,4) = smvcst(1,4)  *  x(63) 

      smvcwt(1,1) = smvccl(1,1)  *  x(64) 	!Luke perturb wilt
      smvcwt(1,2) = smvccl(1,2)  *  x(65) 	
      smvcwt(1,3) = smvccl(1,3)  *  x(66) 
      smvcwt(1,4) = smvccl(1,4)  *  x(67) 
      
      dcatch_dlai(1) = x(68)	!Luke       
      dcatch_dlai(2) = x(69)  	!Luke     
      dcatch_dlai(3) = x(70) 	!Luke   
      dcatch_dlai(4) = x(71)    !Luke    
      dcatch_dlai(5) = x(72)  	!Luke 
        
      dqcrit(1) = x(73)       	!Luke
      dqcrit(2) = x(74)         !Luke
      dqcrit(3) = x(75)	        !Luke 
      dqcrit(4) = x(76)     	!Luke  
      dqcrit(5) = x(77) 	!Luke
      
      sathh(1,1)  = x(78)        !Jupp
      sathh(1,2)  = x(79)        !Jupp
      sathh(1,3)  = x(80)        !Jupp
      sathh(1,4)  = x(81)        !Jupp

      hcap(1,1)  = x(86)        !Jupp
      hcap(1,2)  = x(87)        !Jupp
      hcap(1,3)  = x(88)        !Jupp
      hcap(1,4)  = x(89)        !Jupp      
      
      hcon(1,1)  = x(90)        !Jupp
      hcon(1,2)  = x(91)        !Jupp
      hcon(1,3)  = x(92)        !Jupp
      hcon(1,4)  = x(93)        !Jupp

      q10_soil    = x(94) 	!Luke


!write(*,*) 'top',ts_bias,ts_scal
! edithere ends                                                         
                                                                        
                                                                        
! allocate space for modts and obsts plus monthstarts

      allocate(obsts(nloopcount,nts))
      allocate(modts(nloopcount,nts))
      allocate(monthstarts(nloopcount))                                                                                                                             
                                                               
      call initmodts 
      call initobsts 
      call initavmodts 
      call initavobsts

      if (echo) write(*,*) 'about to call jules top line 219'
! rewind() ! need this to make sure everything rewound?
        call jules
      if (echo) write(*,*) 'have called jules, top lin 221'
! Luke calculate parameter cost
      do i=1,n
        temp=0.0
        do j=1,n
!           write(*,*) 'parcost loop', i,j,PM(i,j),x(j),xoff(j), temp
	  temp = temp + PM(i,j) * (x(j)-xoff(j))
        end do
        parcost=parcost+temp**2
      end do

      nwin=floor(nloopcount/window)
      allocate(newobsts(nwin,nts))
      allocate(newmodts(nwin,nts))

      call getobs !Jupp

      call avgts(newmodts,newobsts,nwin) ! need kahan summation in here?

      call writeavts(nwin,newmodts,newobsts)

!      write(*,*)'averaged'
!      write(*,*) newmodts(1:10,:), newobsts(1:10,:)
!Luke calculate time-series cost

      kc = 0.0 ! Kahan compensation
      do k=1,nwin
        do i=1,nts
            temp2=0.0
            kc2  =0.0
	    do j=1,nts
              if (newobsts(k,j) .ne. -9999.99) then ! more robust than isnan()
              ty = weight(i,j) * (newobsts(k,j)-newmodts(k,j))  - kc2
              tt = temp2 + ty
              kc2 = (tt-temp2)-ty
              temp2 = tt
              end if
	    end do
          temp  = (temp2)**2 - kc        ! Jupp: code rewritten as a Kahan summation
          temp3 = tcost + temp           ! to improve roundoff error - need to rewrite avgts too I wonder? 
          kc    = (temp3-tcost) - temp
          tcost = temp3
!          write(*,*),tcost
        end do
      end do   

!write(*,*) 'parcost ',parcost
!write(*,*) 'tcost   ',tcost
      y(1) = tcost+parcost !Luke 

     if(.not.y(1).lt.huge(y(1))/1.e50) then
     open(300,file='overflow.txt')
     write(300,*) modts,"gap",newmodts
     close(300)
     endif

     if(.not.y(1).gt.-huge(y(1))/1.e50) then
     open(300,file='overflow.txt')
     write(300,*) modts,"gap",newmodts
     close(300)
     endif


!      write(*,*)'cost',y(1)
!         y(1) = cost

      call postmodts 
      call postobsts
      call postavmodts 
      call postavobsts

      deallocate(obsts)
      deallocate(modts)
      deallocate(newobsts)
      deallocate(newmodts)
      deallocate(monthstarts)
!$openad DEPENDENT(Y)                                                                     
      end subroutine func 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine setfunc(nn,mm) 
      use mo_n
      implicit none 
      integer :: nn,mm
! edithere begins                                                       
! user may want to change number of parameters to vary: n               
      n=94      !Luke
! edithere ends                                                         
      m=1 
      nn=n 
      mm=m 
      end subroutine setfunc 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
     subroutine getxoff(xout)
     use inout
     real :: xout(npars)
     xout  = xoff
     end subroutine                                                                 
                                                                        
                                                                        
                                                                        
                                                  
                                                                        
                                                                        
      subroutine initfunc(n,x) 
! initialisation                        
      use inout, only :  stdIn,nts,weight,PM,xoff,npars
      use keep_init 
      use inout, only: formatasc,xoff
      use fomod, only: l_smooth,window
      use file_utils
      implicit none 
      integer :: n, junit, i 
      real    :: x(n), pert(n), xoffout(npars,npars)
      character(8) :: pertfile = 'pert.txt' 
      logical :: ex = .false. 
!                                                                       
! open the run control file on the Standard Input channel (set as "54" e
!                                                                       
                                                                        
      open(stdIn,file='control.jin',action='read',status='unknown') 
                                      ! make sure the standard input is 
      rewind(stdIn) 
                                      ! initialise all variables in the 
      call INIT 
                               ! keep initial state of variables for lat
      call keepinit(n) 

                                                                        
!     set values of control variables                                   
      x(:) = xoff(:) 

!     initialise pert to zero
      pert(:)=0.0                                                                        
                                                                        
                                                                        
                                                                        
!     allow to change the first guess, e.g. for identical twin experimen
      inquire (file=pertfile,exist=ex) 
      if (ex) then 
                                   ! always check whether unit is unconn
         junit=fileunit(formatasc) 
         open (unit=junit,file=pertfile,form='formatted') 
         read (junit,*) pert 
         close (junit) 
         x(:) = x(:) + pert(:) 
         write(*,*)  'adding perturbation of ',pert 
      else 
!         write(*,*)  'no perturbation added' 
      endif

!     Luke: set weight matrix to default identity matrix
      weight(:,:) = 0.0
      do i = 1,nts
        weight(i,i) = 1.0
      enddo
!     Luke initialise smooth logical
      l_smooth=.TRUE.            
!Define matrix for penalisation of parameter departure
      PM(:,:)=0.0

!     Luke set averaging window length to default 1.0
      window=1.0


      end subroutine initfunc 
                                                                        
                                                                        
                                                                        
      subroutine postfunc(n,x,m,y) 
!     postprocessing                                                    
      use inout,       only :  stdIn 
                                        ! FastOpt: force final dump for 
      use inout,       only :  dumpFreq 
      use keep_init 
      implicit none 
      integer      :: n,m 
      real         :: x(n),y(m) 
      integer      :: i, j 
!                                                                       
! close the standard input channel containing the run control file      
!                                                                       
      close(stdIn) 
!      if (dumpFreq==0) then                                ! FastOpt: f
!         write(*,*) 'changing inout::dumpFreq from 0 to 1',! FastOpt: f
!     &        ' in postfunc for debugging. FastOpt'        ! FastOpt: f
!         dumpFreq = 1                                      ! FastOpt: f
!      end if                                                           
                                                                        
!      write(*,*)  
!      write(*,*)  'independent variables :' 
!      print '(a20,2x,a20)', 'variable number', 'value' 

      do i = 1, n 
!         print '(i20,10x,e12.5)', i,x(i) 
      enddo 

!      write(*,*)  
!      write(*,*)  'function values :' 
!      print '(a20,2x,a20)', 'function component', 'value' 

      do j = 1, m 
!         print '(i20,10x,e12.5)', j ,y(j) 
      enddo 

      call postjules 
                         ! deallocate variables recording the initial st
      call postkeep 
                                                                        
      end subroutine postfunc 
                                                                        
    


















                                                                    
      subroutine posttlm(n,x,x_tl,m,y,y_tl) 
!     postprocessing derivative                                         
      implicit none 
      integer      :: n,m 
      real         :: x(n), y(m), x_tl(n), y_tl(m) 
      integer      :: j 
      write(*,*)  
      write(*,*)  'partial derivative value :' 
      print '(a20,2x,a20)', 'function component', 'value' 
      do j = 1, m 
         print '(i20,10x,e12.5)', j ,y_tl(j) 
      enddo 
      end subroutine posttlm 
     
















                                                                   
                                                                        
      subroutine postadm(n,x,x_ad,m,y,y_ad) 
!     postprocessing derivative                                         
      implicit none 
      integer      :: n,m 
      real         :: x(n), y(m), x_ad(n), y_ad(m) 
      integer      :: i 
      write(*,*)  
      write(*,*)  'gradient value :' 
      print '(a20,2x,a20)', 'variable number', 'derivative value' 
      do i = 1, n 
         print '(i20,10x,e12.5)', i, x_ad(i) 
      enddo 
      END











      subroutine setw(w)
      use inout, only: weight,nts
      implicit none
      real :: w(nts,nts)
      weight(:,:)=w(:,:)
      end subroutine setw













      subroutine accessvars(time_step, start_date, start_time)
      use time_loc, only: timestep,dateMainRun,timeRun
      implicit none
      real :: time_step, start_date, start_time
      time_step=timestep
      start_date=dateMainRun(1)
      start_time=timeRun(1)
      end subroutine accessvars










      subroutine setPM(M)
      use inout, only: PM,npars
      implicit none
      real :: M(npars,npars)
      PM(:,:)=M(:,:)
      end subroutine setPM



      subroutine setxvary(xvaryin)
      use inout, only: xvary,npars
      implicit none
      logical :: xvaryin(npars)
      xvary = xvaryin
      end subroutine setxvary











      subroutine setsmooth(lsmooth)
      use fomod, only: l_smooth
      implicit none
      real :: lsmooth
      if (lsmooth==1) then
      l_smooth=.TRUE.
      else if (lsmooth==2) then
      l_smooth=.FALSE.
      end if
      end subroutine setsmooth













      subroutine setwindow(Rwindow)
      use fomod, only: window
      implicit none
      real :: Rwindow
      window=Rwindow
!      write(*,*) window
      end subroutine setwindow










      subroutine avgts(avgmodts,avgobsts,nwin)
      use fomod, only: modts,obsts,window,nloopcount
      use inout, only: nts
      implicit none
      integer :: nwin,low,hi,i,j
      real :: avgmodts(nwin,nts), avgobsts(nwin,nts)
      real :: owintot,owinn,mwintot,mwinn
      logical :: maskobs(nloopcount,nts),maskmod(nloopcount,nts)
      
      maskobs=obsts.ne.(-9999.99)
      maskmod=modts.ne.(-9999.99)

      do i=1,nwin
        low = (i-1)*window+1
        hi  = i*window

        do j=1,nts

          owintot=sum(obsts(low:hi,j),MASK=maskobs(low:hi,j))
          owinn =count(maskobs(low:hi,j))

          if (owinn .ne. 0.0) then
          avgobsts(i,j)=owintot/owinn
          else
          avgobsts(i,j)=-9999.99
          end if

          mwintot=sum(modts(low:hi,j),MASK=maskmod(low:hi,j))
          mwinn=count(maskmod(low:hi,j))

          if (mwinn .ne. 0.0) then
          avgmodts(i,j)=mwintot/mwinn
          else
          avgmodts(i,j)=-9999.99
          end if

        end do
      end do

      end subroutine avgts








      subroutine setmask(lmask)
      use fomod, only: l_mask
      implicit none
      real :: lmask
      if (lmask==1) then
      l_mask=.TRUE.
      else if (lmask==2) then
      l_mask=.FALSE.
      end if
      end subroutine setmask  
















                                                 
