!$TAF SUBROUTINE keep_init::init_time_fastopt  ACTIVE =
!$TAF SUBROUTINE keep_init::init_drive_fastopt ACTIVE =
!$TAF SUBROUTINE keep_init::init_veg_luke ACTIVE =
!$TAF SUBROUTINE keep_init::init_veg_vary_luke ACTIVE =
! module added to store all "initial data" read in from file for subsequent runs
!
!
     module keep_init

     use initial_mod
     use time_mod
     use readwrite_mod
     use snow_param
     use rad_param
     use csmin, only :      cs_min
     use seed, only :       frac_min,frac_seed
     use sigm, only :       pow
     use switches, only :   l_360
     use ancil_info                           
     use nstypes  
     use prognostics
     use p_s_parms 
     use nvegparm
     use pftparm
     use inout
     use coastal
     use time_loc
     use screen
     use c_z0h_z0m
     use trif
     use forcing
     use fluxes
     use aero
     use orog
     use u_v_grid
     use trifctl
     use top_pdm
     use drive_io_vars
     use route_mod	
     use veg_io_vars
     use offline_diag
     use spin_mod
     use fomod, only: iloopcount
     
! edithere begins
     use surf_param, only:  q10_leaf, q10_soil
! edithere ends
         
     implicit none        
!
! variables
!  
     real, allocatable :: frac_keep(:,:) 
     real, allocatable :: gs_keep(:)
     real, allocatable :: sthu_keep(:,:)
     real, allocatable :: sthf_keep(:,:)
     
     
! start of prognostic variables 
     real, allocatable :: nsnow_keep(:,:)
     real, allocatable :: sice_keep(:,:,:)
     real, allocatable :: sliq_keep(:,:,:)
     real, allocatable :: snowdepth_keep(:,:)
     real, allocatable :: tsnow_keep(:,:,:)
     real, allocatable :: rgrainl_keep(:,:,:)
     real, allocatable :: rho_snow_grnd_keep(:,:)
     real, allocatable :: rho_snow_keep(:,:,:)
     real, allocatable :: snow_soil_htf_keep(:,:)   
     real, allocatable :: canht_ft_keep(:,:)
     real, allocatable :: canopy_keep(:,:)
     real, allocatable :: canopy_gb_keep(:)
     real, allocatable :: cs_keep(:,:)
     real, allocatable :: di_keep(:,:)
     real, allocatable :: di_ncat_keep(:,:,:) 
     real, allocatable :: gc_keep(:,:)
     real, allocatable :: lai_keep(:,:)
     real, allocatable :: rgrain_keep(:,:)
     real, allocatable :: smc_keep(:)
     real, allocatable :: smcl_keep(:,:)
     real, allocatable :: snow_tile_keep(:,:)
     real, allocatable :: snow_grnd_keep(:,:)
     real, allocatable :: snow_mass_keep(:,:)
     real, allocatable :: snow_mass_sea_keep(:,:)     
     real, allocatable :: soot_keep(:)
     real, allocatable :: t_soil_keep(:,:)
     real, allocatable :: ti_keep(:,:)     
     real, allocatable :: tstar_tile_keep(:,:)
     real, allocatable :: z0msea_keep(:,:) 
!     real, allocatable :: routestore_keep(:,:) 
! end of prognostic variables

! start of trifctl variables      
     real, allocatable :: g_leaf_acc_keep(:,:)      
     real, allocatable :: npp_ft_acc_keep(:,:)
     real, allocatable :: g_leaf_phen_acc_keep(:,:)        
     real, allocatable :: resp_w_ft_acc_keep(:,:)        
     real, allocatable :: resp_s_acc_keep(:,:) 
     real, allocatable :: gpp_keep(:)
     real, allocatable :: npp_keep(:)
     real, allocatable :: resp_p_keep(:)
     real, allocatable :: g_leaf_keep(:,:)
     real, allocatable :: g_leaf_phen_keep(:,:)
     real, allocatable :: gpp_ft_keep(:,:)
     real, allocatable :: npp_ft_keep(:,:)
     real, allocatable :: resp_p_ft_keep(:,:)
     real, allocatable :: resp_s_keep(:,:)            
     real, allocatable :: resp_w_ft_keep(:,:)
     real, allocatable :: lai_phen_keep(:,:)
     real, allocatable :: c_veg_keep(:,:)
     real, allocatable :: cv_keep(:)
     real, allocatable :: g_leaf_day_keep(:,:)
     real, allocatable :: g_leaf_dr_out_keep(:,:)
     real, allocatable :: lit_c_keep(:,:)
     real, allocatable :: lit_c_mn_keep(:)
     real, allocatable :: npp_dr_out_keep(:,:)
     real, allocatable :: resp_w_dr_out_keep(:,:)
     real, allocatable :: frac_agr_keep(:)
! end of trifctl variables

     real, allocatable :: sea_frac_keep(:)            
     real, allocatable :: sice_frac_keep(:)            
     real, allocatable :: sice_frac_ncat_keep(:,:)     
     
! edithere begins      
     real :: q10_leaf_keep 
     real :: ts_bias_keep(nts)
     real :: ts_scal_keep(nts)                         
     real, allocatable :: nl0_keep(:)               
     real, allocatable :: alpha_keep(:)             
     real, allocatable :: f0_keep(:)
     real, allocatable :: tlow_keep(:)
     real, allocatable :: tupp_keep(:)
     real, allocatable :: dcatch_dlai_keep(:)	!Luke
     real, allocatable :: dqcrit_keep(:)	!Luke
     real, allocatable :: satcon_keep(:,:)
     real, allocatable :: b_keep(:,:)
     real, allocatable :: rootd_ft_keep(:)
     real, allocatable :: albsoil_keep(:)
     real, allocatable :: smvccl_keep(:,:)
     real, allocatable :: smvcwt_keep(:,:)

     real, allocatable :: sathh_keep(:,:)
     real, allocatable :: smvcst_keep(:,:)
     real, allocatable :: hcap_keep(:,:)
     real, allocatable :: hcon_keep(:,:)

     real :: q10_soil_keep			!Luke 

! edithere ends      

 
     logical :: spinup_keep                         
     logical :: spinend_keep                        
     real, allocatable :: spinValOld_keep(:,:,:)    

     contains



!------------------------------------------------------------------------
       subroutine keepinit(n)         ! to be called from initfunc
       implicit none
       integer n
       
       integer ierr
       

!
! allocate space for the reinitialisation variables
! (these will be used to store initial values)
! 

         allocate(frac_keep(LAND_PTS,NTYPE))
         allocate(gs_keep(land_pts))
	 allocate(sthu_keep(land_pts,sm_levels))
	 allocate(sthf_keep(land_pts,sm_levels))

! start of prognostics variables

         allocate(nsnow_keep(land_pts,ntiles))
         allocate(sice_keep(land_pts,ntiles,nsmax))
         allocate(sliq_keep(land_pts,ntiles,nsmax))
         allocate(snowdepth_keep(land_pts,ntiles))
         allocate(tsnow_keep(land_pts,ntiles,nsmax))
         allocate(rgrainl_keep(land_pts,ntiles,nsmax))
         allocate(rho_snow_grnd_keep(land_pts,ntiles))
         allocate(rho_snow_keep(land_pts,ntiles,nsmax),stat=ierr) ! here
         allocate(snow_soil_htf_keep(land_pts,ntiles))
         allocate(canht_ft_keep(land_pts,npft))
         allocate(canopy_keep(land_pts,ntiles))
         allocate(canopy_gb_keep(land_pts))
         allocate(cs_keep(land_pts,dim_cs1))
         allocate(di_keep(row_length,rows))
         allocate(di_ncat_keep(row_length,rows,nice))
         allocate(gc_keep(land_pts,ntiles))
         allocate(lai_keep(land_pts,npft))
         allocate(rgrain_keep(land_pts,ntiles))
         allocate(smc_keep(land_pts))
         allocate(smcl_keep(land_pts,sm_levels))
         allocate(snow_tile_keep(land_pts,ntiles))
         allocate(snow_grnd_keep(land_pts,ntiles))
         allocate(snow_mass_keep(row_length,rows))
         allocate(snow_mass_sea_keep(row_length,rows))
         allocate(soot_keep(row_length*rows))
         allocate(t_soil_keep(land_pts,sm_levels))
         allocate(ti_keep(row_length,rows))	 
         allocate(tstar_tile_keep(land_pts,ntiles))
         allocate(z0msea_keep(row_length,rows))
!         allocate(routestore_keep(nxroute,nyroute))
! end of prognostic variables

! start of trifctl variables
         allocate(g_leaf_acc_keep(land_pts,npft))                
         allocate(npp_ft_acc_keep(land_pts_trif,npft_trif))      
         allocate(g_leaf_phen_acc_keep(land_pts,npft))           
         allocate(resp_w_ft_acc_keep(land_pts_trif,npft_trif)) 
         allocate(resp_s_acc_keep(land_pts_trif,dim_cs1))        
         allocate(gpp_keep(land_pts))
         allocate(npp_keep(land_pts))
         allocate(resp_p_keep(land_pts))
         allocate(g_leaf_keep(land_pts,npft))
         allocate(g_leaf_phen_keep(land_pts,npft))
         allocate(gpp_ft_keep(land_pts,npft))
         allocate(npp_ft_keep(land_pts,npft))
         allocate(resp_p_ft_keep(land_pts,npft))
         allocate(resp_s_keep(land_pts_trif,dim_cs1))            
         allocate(resp_w_ft_keep(land_pts,npft))            
         allocate(lai_phen_keep(land_pts,npft))            
         allocate(c_veg_keep(land_pts,npft))            
         allocate(cv_keep(land_pts))            
         allocate(g_leaf_day_keep(land_pts,npft))            
         allocate(g_leaf_dr_out_keep(land_pts,npft))            
         allocate(lit_c_keep(land_pts,npft))            
         allocate(lit_c_mn_keep(land_pts))            
         allocate(npp_dr_out_keep(land_pts,npft))            
         allocate(resp_w_dr_out_keep(land_pts,npft))            
         allocate(frac_agr_keep(land_pts))            
! end of trifctl variables	 

         allocate(sea_frac_keep(ssi_pts))
         allocate(sice_frac_keep(ssi_pts))
         allocate(sice_frac_ncat_keep(ssi_pts,nice))  

! edithere begins      
         allocate(nl0_keep(npft))                                
         allocate(alpha_keep(npft))                              
         allocate(f0_keep(npft)) 
	 
                                
	                                
         allocate(tlow_keep(npft))                               
         allocate(tupp_keep(npft))
         allocate(dcatch_dlai_keep(npft))	!Luke                             
         allocate(dqcrit_keep(npft))		!Luke                                                         
	 allocate(satcon_keep(land_pts,0:sm_levels))
	 allocate(b_keep(land_pts,sm_levels))
         allocate(albsoil_keep(land_pts))
	 allocate(smvccl_keep(land_pts,sm_levels))
	 allocate(smvcwt_keep(land_pts,sm_levels))

	 allocate(sathh_keep(land_pts,sm_levels)) 
	 allocate(smvcst_keep(land_pts,sm_levels))                                
	 allocate(hcap_keep(land_pts,sm_levels))                                
	 allocate(hcon_keep(land_pts,0:sm_levels))

	 
	 allocate(rootd_ft_keep(1:npft))

! edithere ends      

         allocate( spinvalold_keep(nspinvar,npspinmax,nzspinmax) )
         if ( .not. allocated(spinvalold)) then
           allocate( spinvalold(nspinvar,npspinmax,nzspinmax) )
           spinvalold(:,:,:) = 0.
         endif                                             

         if ( .not. allocated(vegdatain)) then
             allocate( vegDataIn(nvegvar,land_pts,npft,vegTimeIndex(1):vegTimeIndex(2)) )
             vegdatain(:,:,:,:) = 0.
         endif
         if ( .not. allocated(vegunit)) then
           allocate( vegunit(nvegfilevar) )
           vegunit(:) = 0
         endif
         if ( .not. allocated(vegunitfile)) then
           allocate( vegunitfile(nvegfilevar) )
           vegunitfile(:) = 'blank'           
         endif                                                           

!
! now store values of variables in _keep variables
!
         frac_keep            = frac
         gs_keep              = gs
	 sthf_keep            = sthf
	 sthu_keep            = sthu	 	 
!
!   start of variables from prognostics
!	
         nsnow_keep           = nsnow
         sice_keep            = sice
         sliq_keep            = sliq
         snowdepth_keep       = snowdepth
         tsnow_keep           = tsnow
         rgrainl_keep         = rgrainl
         rho_snow_grnd_keep   = rho_snow_grnd
         rho_snow_keep        = rho_snow
         snow_soil_htf_keep   = snow_soil_htf
         canht_ft_keep        = canht_ft     
         canopy_keep          = canopy
         canopy_gb_keep       = canopy_gb
         cs_keep              = cs     
         di_keep              = di
         di_ncat_keep         = di_ncat
         gc_keep              = gc
         lai_keep             = lai
         rgrain_keep          = rgrain
         smc_keep             = smc
         smcl_keep            = smcl
         snow_tile_keep       = snow_tile
         snow_grnd_keep       = snow_grnd
         snow_mass_keep       = snow_mass
         snow_mass_sea_keep   = snow_mass_sea
         soot_keep            = soot
         t_soil_keep          = t_soil
         ti_keep              = ti
         tstar_tile_keep      = tstar_tile
         z0msea_keep          = z0msea
!         routestore_keep      = routestore
!
! end of variables from prognostics
!
!
! start of variables from trifctl
!
         g_leaf_acc_keep      = g_leaf_acc           
         npp_ft_acc_keep      = npp_ft_acc
         g_leaf_phen_acc_keep = g_leaf_phen_acc         
         resp_w_ft_acc_keep   = resp_w_ft_acc     
         resp_s_acc_keep      = resp_s_acc          
         gpp_keep             = gpp
         npp_keep             = npp	
         resp_p_keep          = resp_p
         g_leaf_keep          = g_leaf
         g_leaf_phen_keep     = g_leaf_phen
         gpp_ft_keep          = gpp_ft
         npp_ft_keep          = npp_ft
         resp_p_ft_keep       = resp_p_ft
         resp_s_keep          = resp_s 
         resp_w_ft_keep       = resp_w_ft
         lai_phen_keep        = lai_phen
         c_veg_keep           = c_veg
         cv_keep              = cv
         g_leaf_day_keep      = g_leaf_day
         g_leaf_dr_out_keep   = g_leaf_dr_out
         lit_c_keep           = lit_c 
         lit_c_mn_keep        = lit_c_mn
         npp_dr_out_keep      = npp_dr_out
         resp_w_dr_out_keep   = resp_w_dr_out
         frac_agr_keep        = frac_agr
!
! end of variables from trifctl
!             
         sea_frac_keep        = sea_frac
         sice_frac_keep       = sice_frac
         sice_frac_ncat_keep  = sice_frac_ncat	 
		 
! edithere begins      
         q10_leaf_keep        = q10_leaf      
         nl0_keep             = nl0                    
         alpha_keep           = alpha                 
         f0_keep              = f0  
               
         tlow_keep            = tlow                 
         tupp_keep            = tupp
         dcatch_dlai_keep     = dcatch_dlai	!Luke
         dqcrit_keep          = dqcrit		!Luke                 
	 satcon_keep          = satcon
	 b_keep               = b
	 rootd_ft_keep        = rootd_ft
         ts_bias_keep         = ts_bias
         ts_scal_keep         = ts_scal
         albsoil_keep         = albsoil
         smvccl_keep          = smvccl
         smvcwt_keep          = smvcwt
	 sathh_keep           = sathh
	 smvcst_keep          = smvcst
	 hcap_keep            = hcap
	 hcon_keep            = hcon

         q10_soil_keep        = q10_soil      !Luke

! edithere ends 

!         if (SpinUp) then
            spinup_keep            = spinup             
!            spinend_keep           = spinend            
            spinvalold_keep        = spinvalold  
!         end if


!
! now perturb parameter vector
      xoff(:) = 0.
! edithere begins
      xoff(1)      = q10_leaf    
      xoff(2)      = nl0(1)       
      xoff(3)      = nl0(2)       
      xoff(4)      = nl0(3)       
      xoff(5)      = nl0(4)                   
      xoff(6)      = nl0(5)       
      xoff(7)      = alpha(1)    
      xoff(8)      = alpha(2)     
      xoff(9)      = alpha(3)     
      xoff(10)     = alpha(4)    
      xoff(11)     = alpha(5)    
      xoff(12)     = f0(1)        
      xoff(13)     = f0(2)       
      xoff(14)     = f0(3)       
      xoff(15)     = f0(4)       
      xoff(16)     = f0(5)     
      xoff(17)     = tlow(1)       
      xoff(18)     = tlow(2)          
      xoff(19)     = tlow(3)        
      xoff(20)     = tlow(4)         
      xoff(21)     = tlow(5)        
      xoff(22)     = tupp(1)         
      xoff(23)     = tupp(2)           
      xoff(24)     = tupp(3)           
      xoff(25)     = tupp(4)            
      xoff(26)     = tupp(5)       
      xoff(27)     = lai(1,1)         
      xoff(28)     = lai(1,2)           
      xoff(29)     = lai(1,3)           
      xoff(30)     = lai(1,4)            
      xoff(31)     = lai(1,5) 
      xoff(32)     = canht_ft(1,1)         
      xoff(33)     = canht_ft(1,2)           
      xoff(34)     = canht_ft(1,3)           
      xoff(35)     = canht_ft(1,4)            
      xoff(36)     = canht_ft(1,5) 
      xoff(37)     = satcon(1,1)         
      xoff(38)     = satcon(1,2)           
      xoff(39)     = satcon(1,3)           
      xoff(40)     = satcon(1,4) 
           
      xoff(41)     = b(1,1)	 
      xoff(42)     = b(1,2)         
      xoff(43)     = b(1,3)           
      xoff(44)     = b(1,4) 

      xoff(45)     = cs(1,1) 
      xoff(46)     = rootd_ft(1) 
      xoff(47)     = rootd_ft(2) 
      xoff(48)     = rootd_ft(3) 
      xoff(49)     = rootd_ft(4) 
      xoff(50)     = rootd_ft(5)
      xoff(51)     = ts_bias(1)
      xoff(52)     = ts_bias(2)
      xoff(53)     = ts_bias(3)
      xoff(54)     = ts_bias(4)
      xoff(55)     = ts_scal(1)
      xoff(56)     = ts_scal(2)
      xoff(57)     = ts_scal(3)
      xoff(58)     = ts_scal(4)   
      xoff(59)     = albsoil(1)

      xoff(60)     = smvccl(1,1) / smvcst(1,1)		!Luke/Jupp
      xoff(61)     = smvccl(1,2) / smvcst(1,2)		!Defined soil moisture variables
      xoff(62)     = smvccl(1,3) / smvcst(1,3)		!so crit is a fraction of sat
      xoff(63)     = smvccl(1,4) / smvcst(1,4) 		!and wilt is a fraction of crit
      
      xoff(64)     = smvcwt(1,1) / smvccl(1,1)
      xoff(65)     = smvcwt(1,2) / smvccl(1,2)
      xoff(66)     = smvcwt(1,3) / smvccl(1,3)
      xoff(67)     = smvcwt(1,4) / smvccl(1,4)
      
      xoff(68)     = dcatch_dlai(1)	!Luke       
      xoff(69)     = dcatch_dlai(2)     !Luke     
      xoff(70)     = dcatch_dlai(3)     !Luke   
      xoff(71)     = dcatch_dlai(4)     !Luke    
      xoff(72)     = dcatch_dlai(5)     !Luke   
      xoff(73)     = dqcrit(1)          !Luke
      xoff(74)     = dqcrit(2)          !Luke
      xoff(75)     = dqcrit(3)          !Luke 
      xoff(76)     = dqcrit(4)          !Luke  
      xoff(77)     = dqcrit(5)          !Luke
      
      xoff(78)     = sathh(1,1)
      xoff(79)     = sathh(1,2)
      xoff(80)     = sathh(1,3)
      xoff(81)     = sathh(1,4)
      
      xoff(82)     = smvcst(1,1)
      xoff(83)     = smvcst(1,2)
      xoff(84)     = smvcst(1,3)
      xoff(85)     = smvcst(1,4)
      
      xoff(86)     = hcap(1,1)
      xoff(87)     = hcap(1,2)
      xoff(88)     = hcap(1,3)
      xoff(89)     = hcap(1,4)
      
      xoff(90)     = hcon(1,1)
      xoff(91)     = hcon(1,2)
      xoff(92)     = hcon(1,3)
      xoff(93)     = hcon(1,4)

      xoff(94)      = q10_soil    

         
! edithere ends      


       end subroutine keepinit
!------------------------------------------------------------------------


















!------------------------------------------------------------------------
       subroutine reinit           ! to be called from func
       implicit none
       if (echo) write(*,*) 'entering subroutine reinit right now'
!
! reset initial conditions as appropriate
!
! (some of these may be superfluous but we can keep them in for now)
! 	 
       if (echo) write(*,*) 'line 542'	 
         frac            = frac_keep
         gs              = gs_keep	 
	 sthf            = sthf_keep
	 sthu            = sthu_keep	

       if (echo) write(*,*) 'line 547'
 	 
!
!   start of variables from prognostics
!
         nsnow           = nsnow_keep
         sice            = sice_keep
         sliq            = sliq_keep
         snowdepth       = snowdepth_keep
         tsnow           = tsnow_keep
         rgrainl         = rgrainl_keep
         rho_snow_grnd   = rho_snow_grnd_keep
         rho_snow        = rho_snow_keep
         snow_soil_htf   = snow_soil_htf_keep	 
         canht_ft        = canht_ft_keep     
         canopy          = canopy_keep
         canopy_gb       = canopy_gb_keep
         cs              = cs_keep     
         di              = di_keep
         di_ncat         = di_ncat_keep
         gc              = gc_keep
         lai             = lai_keep
         rgrain          = rgrain_keep
         smc             = smc_keep
         smcl            = smcl_keep
         snow_tile       = snow_tile_keep
         snow_grnd       = snow_grnd_keep
         snow_mass       = snow_mass_keep
         snow_mass_sea   = snow_mass_sea_keep
         soot            = soot_keep
         t_soil          = t_soil_keep
         ti              = ti_keep
         tstar_tile      = tstar_tile_keep
         z0msea          = z0msea_keep
!         routestore      = routestore_keep
!
! end of variables from prognostics
!

!
! start of variables from trifctl
!
         g_leaf_acc      = g_leaf_acc_keep           
         npp_ft_acc      = npp_ft_acc_keep
         g_leaf_phen_acc = g_leaf_phen_acc_keep         
         resp_w_ft_acc   = resp_w_ft_acc_keep     
         resp_s_acc      = resp_s_acc_keep          
         gpp             = gpp_keep
         npp             = npp_keep	
         resp_p          = resp_p_keep
         g_leaf          = g_leaf_keep
         g_leaf_phen     = g_leaf_phen_keep
         gpp_ft          = gpp_ft_keep
         npp_ft          = npp_ft_keep
         resp_p_ft       = resp_p_ft_keep
         resp_s          = resp_s_keep 
         resp_w_ft       = resp_w_ft_keep
         lai_phen        = lai_phen_keep
         c_veg           = c_veg_keep
         cv              = cv_keep
         g_leaf_day      = g_leaf_day_keep
         g_leaf_dr_out   = g_leaf_dr_out_keep
         lit_c           = lit_c_keep 
         lit_c_mn        = lit_c_mn_keep
         npp_dr_out      = npp_dr_out_keep
         resp_w_dr_out   = resp_w_dr_out_keep
         frac_agr        = frac_agr_keep
!
! end of variables from trifctl
!             
         sea_frac        = sea_frac_keep
         sice_frac       = sice_frac_keep
         sice_frac_ncat  = sice_frac_ncat_keep	 		 
! edithere begins      
         q10_leaf        = q10_leaf_keep      
         nl0             = nl0_keep                    
         alpha           = alpha_keep                 
         f0              = f0_keep
	 
	 sathh           = sathh_keep
	 smvcst          = smvcst_keep
	 hcap            = hcap_keep
	 hcon            = hcon_keep
	                  
         tlow            = tlow_keep                 
         tupp            = tupp_keep
         dcatch_dlai     = dcatch_dlai_keep	!Luke
         dqcrit          = dqcrit_keep		!Luke                 
	 satcon          = satcon_keep
	 b               = b_keep
	 rootd_ft        = rootd_ft_keep
         ts_bias         = ts_bias_keep
         ts_scal         = ts_scal_keep
         albsoil         = albsoil_keep
         smvccl          = smvccl_keep
         smvcwt          = smvcwt_keep

         q10_soil        = q10_soil_keep	!Luke      

! edithere ends 
!         if (SpinUp) then
            spinup            = spinup_keep             
!            spinend           = spinend_keep            
            spinvalold        = spinvalold_keep  
!         end if	 
!
! following lines are culled from INIT.f
! and seem to be a minimal set of reinitialisation instructions
!	 
      iloopcount = 0 
       if (echo) write(*,*) 'about to call init_time_fastopt'
      call init_time_fastopt    ! reinitialise time data with file reads and allocations
       if (echo) write(*,*) 'about to call init_drive_fastopt'
      call init_drive_fastopt   ! reinitialise driving data    
      if (echo) write(*,*) 'about to rewind unit', driveunit
      rewind(drivefile)         ! rewind the driving data file
      if (echo) write(*,*) 'have rewound'

      IF ( vegVaryT ) THEN

      if (echo) write(*,*) 'about to rewind unit', vegunit
      rewind(vegfile)         ! rewind the vegfile
      if (echo) write(*,*) 'have rewound vegfile'

        call init_veg_luke

      if (echo) write(*,*) 'about to rewind unit', vegunit
      rewind(vegfile)         ! rewind the vegfile
      if (echo) write(*,*) 'have rewound vegfile'
      ENDIF

      call reinit_parms           ! reinitialise 
      IF ( .NOT.vegVaryT .AND. dumpFreq>1 ) CALL dump_io( .TRUE., dumpTypeArg='init' )	 			 
      if (echo) write(*,*) 'leaving subroutine reinit'
    end subroutine reinit
!------------------------------------------------------------------------




























!------------------------------------------------------------------------
       subroutine postkeep          ! to be called from postfunc
       implicit none

         deallocate(frac_keep)
         deallocate(gs_keep)
	 deallocate(sthf_keep)
	 deallocate(sthu_keep)

! start of prognostics variables

         deallocate(nsnow_keep)
         deallocate(sice_keep)
         deallocate(sliq_keep)
         deallocate(snowdepth_keep)
         deallocate(tsnow_keep)
         deallocate(rgrainl_keep)
         deallocate(rho_snow_grnd_keep)
         deallocate(rho_snow_keep)
         deallocate(snow_soil_htf_keep)
         deallocate(canht_ft_keep)
         deallocate(canopy_keep)
         deallocate(canopy_gb_keep)
         deallocate(cs_keep)
         deallocate(di_keep)
         deallocate(di_ncat_keep)
         deallocate(gc_keep)
         deallocate(lai_keep)
         deallocate(rgrain_keep)
         deallocate(smc_keep)
         deallocate(smcl_keep)
         deallocate(snow_tile_keep)
         deallocate(snow_grnd_keep)
         deallocate(snow_mass_keep)
         deallocate(snow_mass_sea_keep)
         deallocate(soot_keep)
         deallocate(t_soil_keep)
         deallocate(ti_keep)	 
         deallocate(tstar_tile_keep)
         deallocate(z0msea_keep)
!         deallocate(routestore_keep)
! end of prognostic variables




! start of trifctl variables
         deallocate(g_leaf_acc_keep)                
         deallocate(npp_ft_acc_keep)      
         deallocate(g_leaf_phen_acc_keep)           
         deallocate(resp_w_ft_acc_keep) 
         deallocate(resp_s_acc_keep)        
         deallocate(gpp_keep)
         deallocate(npp_keep)
         deallocate(resp_p_keep)
         deallocate(g_leaf_keep)
         deallocate(g_leaf_phen_keep)
         deallocate(gpp_ft_keep)
         deallocate(npp_ft_keep)
         deallocate(resp_p_ft_keep)
         deallocate(resp_s_keep)            
         deallocate(resp_w_ft_keep)            
         deallocate(lai_phen_keep)            
         deallocate(c_veg_keep)            
         deallocate(cv_keep)            
         deallocate(g_leaf_day_keep)            
         deallocate(g_leaf_dr_out_keep)            
         deallocate(lit_c_keep)            
         deallocate(lit_c_mn_keep)            
         deallocate(npp_dr_out_keep)            
         deallocate(resp_w_dr_out_keep)            
         deallocate(frac_agr_keep)            
! end of trifctl variables
	 

         deallocate(sea_frac_keep)
         deallocate(sice_frac_keep)
         deallocate(sice_frac_ncat_keep)
   

! edithere begins      
         deallocate(nl0_keep)                                
         deallocate(alpha_keep)                              
         deallocate(f0_keep)    
	 
	 deallocate(sathh_keep)
	 deallocate(smvcst_keep)
	 deallocate(hcap_keep)
	 deallocate(hcon_keep)
	                              
         deallocate(tlow_keep)                               
         deallocate(tupp_keep)
         deallocate(dcatch_dlai_keep)	!Luke
         deallocate(dqcrit_keep)	!Luke                               
	 deallocate(satcon_keep)
	 deallocate(b_keep)
	 deallocate(rootd_ft_keep)
         deallocate(albsoil_keep)
         deallocate(smvccl_keep)
         deallocate(smvcwt_keep)
         deallocate( spinvalold_keep )

! edithere ends      

                                                           
           
       end subroutine postkeep
!------------------------------------------------------------------------
       
     end module keep_init




































!
! This is a version of init_time, hacked to remove any read and deallocate statements
! TEJ, Jan 2009
!
!       init_time should be called *once* before multiple calls to JULES
! while init_time_fastopt should be called *between* multiple calls to JULES
! to reset the date and time information
!
!
! Read details of model timestep and run length.
 
  SUBROUTINE init_time_fastopt

  use file_utils, only: findTag
  use inout, only : echo,stdIn,npars

  use spin_mod, only : ispin,nspin,nspinFinal,nspinVarMax  &
             ,spinEnd,spinFail,spinTol,spinTolPercent,spinUp,spinVar,spinVarName
  use switches, only : l_360
  use time_loc, only : date,dateMainRun,dateNext,datePrev,dateRun,dateSpin,endMonth,endSec,endYear  &
                      ,stepFlag,time,timeNext,timePrev,timeRun,timeStep
  use time_mod, only : chhmmss_to_s,s_to_chhmmss,timeDate,timeDate_cmp
  use trifctl, only : phenol_period,triffid_period
  use timeConst, only : iSecInDay

!-------------------------------------------------------------------------------
  IMPLICIT NONE

!  real, external:: epsilon
!  integer, external:: isecinday
!  logical, external :: timedate_cmp
!  character(len=10), external :: s_to_chhmmss

  INTEGER ::  &!  local SCALARS
    i          !  work

  LOGICAL ::  &!  local SCALARS
    endrun     !  work (newTime needs an intent(out) logical variable, will be FALSE here)

  CHARACTER(len=8) ::  &!  local SCALARS
    cTimeRun(2)  !   times (hh:mm:ss)

  CHARACTER(len=10) ::  time_hms  ! Time in the form "hh:mm:ss H", used as a local variable
!                                 ! that is required because some compilers can not cope
!                                 ! with recursive write statements

!------------------------------------------------------------------------------------------
  if (echo) WRITE(*,"(50('-'),/,a)") 'init_time_fastopt'

!------------------------------------------------------------------------
! Locate the start of this section in input file.
!------------------------------------------------------------------------
!FastOpt  CALL findTag( stdIn,'init_time','>INIT_TIME' )

!------------------------------------------------------------------------
! Initialise.
!------------------------------------------------------------------------
  spinUp = .FALSE.
  spinEnd = .FALSE.
  stepFlag = 3  !  first timestep in "main" run

!------------------------------------------------------------------------
! Read data from JULES-in file.
!------------------------------------------------------------------------
!FastOpt  read(stdIn,*) timeStep
!FastOpt  read(stdIn,*) dateMainRun(1),ctimeRun(1),dateMainRun(2),cTimeRun(2)
!FastOpt  read(stdIn,*) dateSpin(1:2),nspin

!FastOpt  read(stdIn,*) spinFail

! For v2.0, only read the first two values of these spin up variables
! (these are for soil moisture and temperature), and ignore the third
! (routing store), since routing is not active in v2.0.
!  DO i=1,nspinVarMax    !   removed for v2.0
  DO i=1,2
!FastOpt    read(stdIn,*) spinVar(i),spinTolPercent(i),spinTol(i)
  ENDDO

! Ensure routing is not used for spin up (for v2.0).
  spinVar(3) = .FALSE.

!FastOpt  read(stdIn,*) PHENOL_PERIOD,TRIFFID_PERIOD
!FastOpt  read(stdIn,*) l_360

!------------------------------------------------------------------------
! Process time variables. Check that timestep is an integer number of 
! seconds.
!------------------------------------------------------------------------
  IF ( timeStep - REAL(INT(timeStep)) >= EPSILON(timeStep) )  &
      STOP'ERROR: INIT_TIME: timestep must be an integer number of seconds'

!------------------------------------------------------------------------
! Make sure one day is a multiple of timestep - makes life easier.
!------------------------------------------------------------------------
!write(*,*) 'keep_init isecinday',isecinday
!write(*,*) 'keep_init timestep',timestep

  IF ( MOD(iSecInDay,NINT(timeStep)) /= 0 ) STOP'ERROR: init_time: 24 hours must be a multiple of timestep'

!------------------------------------------------------------------------
! Convert input times (hh:mm:ss) to seconds.
!------------------------------------------------------------------------
  DO i=1,2
!FastOpt    timeRun(i) = chhmmss_to_s( cTimeRun(i),'init_time' )
  ENDDO

!------------------------------------------------------------------------
! At present, any spin up must be over times that either:
!    1. start at the same the "main" run, and end at or before the time
!       when the main run starts.
! or 2. end at the start time of the "main" run.
!
! For #1, the spin up is thus over times that are either the same as 
! those for the main run, or are a subset of those starting at same time.
! #1 means that spin up period is < or = to length of main run. For #2 
! the spin up can be over a period that is <,=, or > than main run.
!------------------------------------------------------------------------

! Check that run times are in chronological order.
  IF ( timeDate_cmp(timeRun(1),dateMainRun(1),'>=',timeRun(2),dateMainRun(2),'init_time') ) &
      STOP'ERROR: init_time: Start time/date for main run is >= end.'

! Set dates for integration, assuming no spin up.
  dateRun(:) = dateMainRun(:)

! Initially assume that all spin-up cycles will be required.
  nspinFinal = nspin

  ispin = 0 !FastOpt 
! Check any spinup period is correctly prescribed.
  IF ( nspin /= 0 ) THEN
    spinUp=.TRUE.
    ispin = 0   !  the call to newTime increments this

!------------------------------------------------------------------------
!   Check dates for spin up are reasonable.  Spinup must immediately 
!   precede the main run, or have same start date as run (eg if want to
!   spinup over same year). Spinup will be for a complete number of days.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!   First checking that spinup dates are successive.
!------------------------------------------------------------------------
    IF ( timeDate_cmp( 0,dateSpin(1),'>=',0,dateSpin(2),'init_time' ) ) STOP'ERROR: init_time: End date for spin up is <= start.'
!   Now check that spinup times agree with run times...
    IF ( dateSpin(2)/=dateMainRun(1) .AND. dateSpin(1)/=dateMainRun(1) ) &
        STOP'ERROR: init_time: Spinup must immediately precede main run OR start at same date/time as main run.'

!------------------------------------------------------------------------
!   Reset start date for integration.
!------------------------------------------------------------------------
    dateRun(1) = dateSpin(1)

!------------------------------------------------------------------------
!XX At present, output does not work particularly well for the case of a 
!XX single cycle of spin up followed by a main run with 
!XX dateMainRun(1)=dateSpin(2). Part of the problem is that endSec is set
!XX to T at end of spin up, whereas other bots of code want endSec=F for
!XX this special case.   This needs to be sorted out, possibly by always
!XX setting endSec=T at end of spin up, and always naming spin-up output 
!XX with .spin. in filenames. Any solution should probably avoid creating
!XX this "special case", since it has to be correctly accounted for at 
!XX all later stages in the code!  At present, if spin up ends mid-month,
!XX monthly output creates a file for the last month of spin up, then 
!XX endSec resets outFirstWrite, and then at start of main run a file 
!XX with same name is required...and goes wrong.  So, avoid all these 
!XX complications for now (although current code is likely OK for cases
!XX such as no output through spin up, or every time to separate file).
!------------------------------------------------------------------------
    IF ( nspin==1 .AND. dateMainRun(1)==dateSpin(2) ) THEN
      WRITE(*,*)'ERROR: init_time: precautionary error!'
      WRITE(*,*)'Output code is not robust/working correctly for nspin==1 .AND. dateMainRun(1)==dateSpin(2).'
      WRITE(*,*)'You can anyway do this integration by setting'
      WRITE(*,*)'dateMainRun(1) to current value of =dateSpin(1), and setting nspin=0'
      WRITE(*,*)'This means the "spin up" part of the run is now considered part of the "main" run.'
      WRITE(*,*)'Sorry about this - hopefully the code will be improved, one day...'
      STOP
    ENDIF

  ENDIF  !  nspin

!------------------------------------------------------------------------
! Set time to one timestep before start, so that first increment (call to
! newTime, below) takes to start time.
!------------------------------------------------------------------------
  IF ( .NOT. spinUp ) THEN
    date = dateRun(1)
  ELSE
    date = dateSpin(1)
  ENDIF
  time = timeRun(1)
  CALL timeDate( time,date,-1*NINT(timeStep),'sec',l_360,timePrev,datePrev,'init_time')
  date = datePrev
  time = timePrev

! Initialise past and next times.
  CALL timeDate( time,date,-1*NINT(timeStep),'sec',l_360,timePrev,datePrev,'init_time' )
  CALL timeDate( time,date,NINT(timeStep),'sec',l_360,timeNext,dateNext,'init_time' )

!------------------------------------------------------------------------
! Establish if we are at the start/end of a month/year. newTime will also
! advance time to the start time.  Set "end flags", so that newTime sets
! "new flags".
!------------------------------------------------------------------------
  endMonth = .TRUE.
  endYear = .TRUE.
  endSec = .TRUE.

!------------------------------------------------------------------------
! 1st argument to newTime (timestep number) is zero, to indicate that 
! none of the spin up code is to be activated.
!------------------------------------------------------------------------
  CALL newTime ( 0,endRun )

!------------------------------------------------------------------------
! Some output to screen.
!------------------------------------------------------------------------
  IF ( echo ) THEN

    IF ( nspin == 0 ) THEN
      WRITE(*,*) 'There is NO spin-up period.'
    ELSE
      time_hms = s_to_chhmmss( timeRun(1) )
      WRITE(*,"(a,i8,a,tr1,a)") 'Spin-up period starts ',dateSpin(1),' time=',time_hms
      time_hms = s_to_chhmmss( timeRun(2) )
      WRITE(*,"(a,i8,a,tr1,a)") 'Spin-up period ends   ',dateSpin(2),' time=',time_hms
      WRITE(*,*) 'Spin up period repeated up to ',ABS(nspin),' times.'
      IF ( spinFail ) THEN
        WRITE(*,*)'If model is not then spun-up, run will stop'
      ELSE
        WRITE(*,*)'If model is not then spun-up, run will CONTINUE'
      ENDIF
      WRITE(*,*)'Variables and tolerances used in spin up: see later.'
    ENDIF  !  nspin
    
    time_hms = s_to_chhmmss( timeRun(1) )
    WRITE(*,"(a,i8,a,tr1,a)") 'Main run (after spin up) starts at ',dateMainRun(1),' time=',time_hms
    time_hms = s_to_chhmmss( timeRun(2) )
    WRITE(*,"(a,i8,a,tr1,a)") 'and ends at                         ',dateMainRun(2),' time=',time_hms

  ENDIF  !  echo

  END SUBROUTINE init_time_fastopt



!############################################################################### 
SUBROUTINE Init_Drive_fastopt


  use ancil_info, only :  &!
!  imported arrays with intent(out)
     z1_uv,z1_tq

  use drive_io_vars, only :  &
!   imported scalars with intent(out)
      driveDataPer,driveDataStep,driveDataStepInit,driveDate,driveDateInit  &
     ,driveEndTime,driveFile,driveFilePer,driveFileStep,driveFileTemplate  &
     ,driveFormat,driveResetStep,driveResetStepPrev  &
     ,driveTemplateT,driveTemplateDate,driveTemplateUnits,driveTemplateTime  &
     ,driveTemplateV,driveTime,driveTimeInit  &
     ,ioPrecipType,io_rad_type &
     ,ndriveDataTime,ndriveFileTime,ndriveUnit  &
     ,ndriveHeaderField,ndriveHeaderFile,ndriveHeaderTime  &
     ,ndriveVar,ndriveVarIn,ndriveVarMax  &
     ,nfieldDriveFile,noNewLineDrive,notNextDrive,tForCRain,tForSnow    &
!   imported arrays with intent(out)  &
     ,driveFileDate,driveFileName,driveFileTime  &
     ,driveVarInterp,driveVarInSort  &
     ,driveTimeIndex,driveUnit,driveFileOnUnit,driveUnitUse  &
     ,driveVarFlag,driveVarNameSDF,driveVarStash,driveVarUse

  use file_utils, only : &
!   imported procedures
     closeFile,fileUnit,findTag,openFile


  use inout, only :  &
!  imported scalar parameters
     formatAsc,formatBin,formatNc,formatLen,formatPP,periodAnn,periodMon,periodOneFile  &
    ,stdIn,tagAscBin,tagNc  &
!  imported scalars with intent(in)
    ,echo,nxIn,nyIn
  
  use misc_utils, only :  &
!  imported procedures
     check_template,read_list,repeatVal

  use spin_mod, only :  &
!   imported scalars with intent(in)
     spinUp

  use switches, only :  &
!   imported scalars with intent(in)
     l_360,l_point_data,routeOnly

  use time_loc, only :  &
!   imported scalars with intent(in)
     timeStep

  use time_mod, only :  &
!   imported procedures
     chhmmss_to_s,timeDate,timeDate_cmp

  use timeConst, only : &
     iSecInDay

  use update_mod, only :  &
!   imported procedures
     data_init,calc_reset_step

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    i,ierr,ipos,iunit,ivar   &!  work
   ,ndriveVarNeed            &!  number of driving variables for chosen configuration
!                                  Only used for diagnostic message to screen.
   ,t                        &!  work
   ,tmpDataStepMax            !  work

  INTEGER ::  &!  local ARRAYS
    tmpFlag(ndriveVarMax)   !  values of driveVarFlag as read in

  REAL ::  &!  local scalars
    z1_uv_val,z1_tq_val    ! values of z1_uv and z1_tq that are used at all points on grid

  LOGICAL ::  &!  local SCALARS
    haveAve     &!  T means that one or more variable in input data is a time-average
!                     that places extra restriction on timestep
   ,readList    &!  T means read a list of file names, F means read a single file name
   ,templateT,templateV  &!  work
   ,ioWindSpeed           !  T means that the windspeed is input
!                            F means 2 components of wind are input

  LOGICAL ::  &! local ARRAYS
    done(ndriveVarMax)  &!  work
   ,needTime(-1:2)       !  flag indicating which times of data are required
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)

  CHARACTER(len=8) ::  &!  local SCALARS
    cdriveTime   !  time string

  CHARACTER(len=LEN(driveVarNameSDF)) ::  &!  local ARRAYS
    driveVarName(ndriveVarMax)   &!  names of all possible forcing variables
   ,driveVarFileName(ndriveVarMax)  &!  the part of the file name that changes between
!                  files that hold different variables, if driveTemplateV=TRUE.
   ,tmpName(ndriveVarMax)        &!  work: names of forcing variables, as read in.
   ,tmpNameSDF(ndriveVarMax)     &!  work: names of forcing variables for SDF files
   ,tmpNameFile(ndriveVarMax)     !  work: names used in file names, as read in

  CHARACTER(len=LEN(driveVarInterp)) ::  &!  local ARRAYS
    tmpInterp(ndriveVarMax)         !  work: values of driveVarInterp, as read in

  CHARACTER(len=150) ::  &!  local SCALARS
    listFile  !  work

  CHARACTER(len=70) ::  &!  local arrays
    varDesc(ndriveVarMax)  !  a description of each possible driving variable
!------------------------------------------------------------------------------------------


  if (echo) WRITE(*,"(50('-'),/,a)") 'init_drive_fastopt'   

!-------------------------------------------------------------------------------
! Initialise driving data. Read in past values that are needed on first timestep.
! The 2nd dummy argument is dataStepMax, which at present is only really needed
! for vegetation data (i.e. not here), but making it an optional argument has
! its disadvantages - so tmpDataStepMax is passed here (but not expected to change!).
! Admittedly this is rather untidy.....
!-------------------------------------------------------------------------------
   tmpDataStepMax = driveDataPer
   if (echo) write(*,*) 'calling data_init...',driveFilePer   
   templateT = driveTemplateT  !FastOpt templateT was uninitialised and causing errors - please verify!
   CALL data_init( driveDataPer,driveFilePer,driveTemplateDate,driveTemplateTime  &
    ,tmpDataStepMax,driveDataStep,driveDataStepInit  &
    ,driveDate,driveDateInit  &
    ,driveFile,driveFileStep,driveResetStep,driveResetStepPrev  &
    ,driveTime,driveTimeInit  &
    ,driveFileDate,driveFileTime,driveFileTemplate,driveTemplateUnits,driveFileName  &
    ,driveTimeIndex,.FALSE.,driveEndTime,templateT,notNextDrive,'drive' )
   if (echo) write(*,*) '...done'   

!-------------------------------------------------------------------------------
! Deal with the possibility that driving data may need to be "recycled" to cope 
! with spin up.
! Arguments: zero=timestep.
!-------------------------------------------------------------------------------

  if (echo) write(*,*) 'keep_init 855'
  IF ( spinUp ) CALL calc_reset_step( 0,driveResetStep,driveResetStepPrev,'drive' )
  if (echo) write(*,*) 'keep_init 857'

END SUBROUTINE init_drive_fastopt


!Luke
! subroutine init_veg_luke
! Driver to read parameters for vegetation surface types (PFTs).
 
  SUBROUTINE init_veg_luke()

  USE switches, ONLY :&
!  imported scalars with intent(in)
     routeOnly

  IMPLICIT NONE

!------------------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

! Read values for prescribed variables that may be functions of time and/or space.
  CALL init_veg_vary_luke()

  END SUBROUTINE init_veg_luke

! subroutine init_veg_vary
! Read details of prescribed vegetation variables that are a function of space and/or time.

  SUBROUTINE init_veg_vary_luke

  USE file_utils, ONLY:  &
!  imported procedures
    closeFile,fileUnit,findTag,openFile

  USE inout, ONLY :   &
!  imported scalar parameters
    formatAsc,formatBin,formatLen,formatNc,formatPP,npftInFile  &
   ,periodAnn,periodMon,periodOneFile,stdIn  &
   ,tagAscBin,tagNc  &
!  imported scalars with intent(in)
   ,echo,nxIn,nyIn

  USE misc_utils, ONLY:  &
!  imported procedures
     check_template,repeatVal

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     npft

  USE pftparm, ONLY :   &
!  imported arrays with intent(in)  &
     pftname

  USE prognostics, ONLY :  &
!   imported arrays with intent(out)
     canht_ft,lai

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
     spinUp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
    l_360,l_triffid
 
  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
    timeStep

  USE time_mod, ONLY :  &
!  imported procedures
    chhmmss_to_s,dateToBits,timeDate,timeDate_cmp,timeDate_diff

  USE timeConst, ONLY : &
    iSecInDay

  USE update_mod, ONLY :  &
!  imported procedures
    data_init,calc_reset_step,veg_update

  USE veg_io_vars, ONLY :  &
!  imported scalar parameters
    nvegVarMax,vegVarFlagPT,vegVarFlagPTX,vegVarFlagPX  &
!  imported scalars with intent(out)
   ,nfieldVegFile,notNextVeg,nvegDataTime,nvegFileTime,nvegFileVar  &
   ,nvegHeaderField,nvegHeaderFile,nvegHeaderTime,nvegVar  &
   ,noNewLineVeg,varNumCanht,varNumLAI,varNumRootd,vegClim  &
   ,vegDataPer,vegDataStep,vegDataStepInit,vegDataStepMax,vegDate,vegDateInit  &
   ,vegEndTime,vegFile,vegFilePer,vegFileStep,vegFileTemplate,vegFormat  &
   ,vegResetStep,vegResetStepPrev,vegTemplateDate,vegTemplateT,vegTemplateTime  &
   ,vegtemplateUnits,vegTemplateV,vegTime  &
   ,vegTimeInit,vegUpdatePer,vegUpdateStep,vegUpdateStepMax,vegUpdateStepInit,vegVaryT  &
!  imported arrays with intent(out)
   ,vegDataIn,vegFileDate,vegFileName  &
   ,vegFileTime,vegVarFlag,vegVarInterp  &
   ,vegVarName,vegVarNameFile,vegTimeIndex,vegUnit,vegUnitFile  &
   ,vegVarPos,vegVarStash
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &! local SCALARS
    i,ierr,ierrSum,ipos,iunit,ivar,j,jvar  &!  work
   ,nsec,nmin,nhr,ndy,nmon,nyr    &!  work
   ,t,tmpDay,tmpMonth,tmpYear     &!  work
   ,vegDataDay         !  day of month to which veg data refer to (monthly data only)

  INTEGER ::  &!  local ARRAYS
    ival(1)                   &!  work
   ,tmpPos(nvegVarMax)         !  work: temporary version of vegVarPos

  LOGICAL ::  &!  local SCALARS
    haveAve     &!  T means that one or more variable in input data is a time-average
!                     that places an extra restriction on timestep
   ,readList    &!  T means read a list of file names, F means read a single file name
   ,templateT,templateV   !  work

  LOGICAL ::  &! local ARRAYS
    done(nvegVarMax)  &!  work
   ,needTime(-1:2)     !  flag indicating which times of data are required
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)

  CHARACTER(len=8)  ::  &!  local SCALARS
    cvegTime   !  time string

  CHARACTER(len=LEN(vegVarName)) ::  &!  local ARRAYS
    vegName(nvegVarMax)     &!  work: names of veg variables
!                                     Only used to identify variables in this routine
   ,tmpName(nvegVarMax)     &!  work: temporary version of vegName
   ,tmpNameFile(nvegVarMax) &!  work: temporary version of vegVarNameFile
   ,tmpVegName(nvegVarMax)   !  work: temporary version of vegName

  CHARACTER(len=LEN(vegVarFlag)) ::  &!  local arrays
    tmpFlag(nvegVarMax)     !  work: temporary version of vegVarFlag

  CHARACTER(len=LEN(vegVarInterp)) ::  &!  local ARRAYS
    tmpInterp(nvegVarMax)         !  work: temporary version of vegVarInterp

  CHARACTER(len=150) ::  &!  local SCALARS
    line     &!  work
   ,listFile  !  work
!-------------------------------------------------------------------------------

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_veg_vary_luke'

!---------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Initialise veg data.
  IF ( vegVaryT ) THEN

!   Read in past values that are needed on first timestep.

  if ( echo ) write(*,*) 'calling data_init...',vegFilePer !Luke

    CALL data_init( vegDataPer,vegFilePer,vegTemplateDate,vegTemplateTime  &
        ,vegDataStepMax,vegDataStep,vegDataStepInit  &
        ,vegDate,vegDateInit  &
        ,vegFile,vegFileStep,vegResetStep,vegResetStepPrev  &
        ,vegTime,vegTimeInit  &
        ,vegFileDate,vegFileTime,vegFileTemplate,vegTemplateUnits,vegFileName  &
        ,vegTimeIndex,vegClim,vegEndTime,templateT,notNextVeg,'veg'  &
        ,vegUpdatePer,vegUpdateStepMax,vegUpdateStep,vegUpdateStepInit )
  if ( echo ) write(*,*) '..done',vegFilePer !Luke
!   Deal with the possibility that veg data may need to be "recycled" to cope with spin up.
!   Arguments: zero=timestep.
    IF ( spinUp ) CALL calc_reset_step( 0,vegResetStep,vegResetStepPrev,'veg' )

  ELSE

!   Read in a spatial field.
!   1st argument (a_step)=0, 2nd (next)=.TRUE. to read next data, 3rd should read
!   data into the only time level of vegDataIn.
    CALL veg_update( 0,.TRUE.,vegTimeIndex(1) )

  ENDIF

!-------------------------------------------------------------------------------
! If veg is not time-varying, don't need units or variables again.
! Also do optional print out.
  IF ( .NOT. vegVaryT ) THEN

!   Close files and deallocate.
    DO i=1,nvegFileVar
      CALL closeFile( vegUnit(i),vegFormat )
    ENDDO
    ierrSum = 0
!    DEALLOCATE( vegVarPos,vegVarFlag,vegVarInterp,vegVarName,vegVarNameFile, stat=ierr ); ierrSum=ierrSum+ierr
!    DEALLOCATE( vegFileName,vegFileDate,vegFileTime, stat=ierr ); ierrSum=ierrSum+ierr
!    DEALLOCATE( vegDataIn,vegUnit,vegUnitFile, stat=ierr ); ierrSum=ierrSum+ierr
!    IF ( ierrSum/=0 ) WRITE(*,"(50('#'),/,a,50('#'))")'WARNING: init_veg_vary: could not deallocate'

!   Optional print to screen for fields that are not functions of time.
    IF ( echo ) THEN

!     First, print entire fields.
      DO ivar=1,nvegVar
        IF ( ivar == varNumCanht ) THEN
          DO i=1,npft
            WRITE(*,*)'Canopy height for ',TRIM(pftName(i))
            WRITE(*,*) canht_ft(:,i)
          ENDDO
        ELSEIF ( ivar == varNumLAI ) THEN
          DO i=1,npft
            WRITE(*,*)'LAI for ',TRIM(pftName(i))
            WRITE(*,*) lai(:,i)
          ENDDO
        ENDIF
      ENDDO  !  ivar

!     Now, print summary info. Note that this will include locations where frac=0.
      DO ivar=1,nvegVar
        IF ( ivar == varNumCanht ) THEN
          DO i=1,npft
            WRITE(*,"(a,2(a,f6.2))") 'Range of canopy height for '  &
               ,TRIM(pftName(i)),MINVAL(canht_ft(:,i)),' to ',MAXVAL(canht_ft(:,i))
          ENDDO
        ELSEIF ( ivar == varNumLAI ) THEN
          DO i=1,npft
            WRITE(*,"(a,2(a,f6.2))") 'Range of LAI for '  &
              ,TRIM(pftName(i)),MINVAL(lai(:,i)),' to ',MAXVAL(lai(:,i))
          ENDDO
        ENDIF
      ENDDO  !  ivar

    ENDIF   !  echo

  ENDIF   !  vegVaryT

  END SUBROUTINE init_veg_vary_luke

