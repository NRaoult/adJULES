!**************************************************************         
!   Copyright (C) 2000, 2001                                            
!   FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR                  
!          , Hamburg, Germany                                           
                                                                        
                                                                        
                                                                        
                                                                        
      subroutine doit( n, m ) 
      use inout, only : weight,xoff  
                                                                        
      implicit none 
      integer m,n 
      real  x(n)
      real  y(m), eps 
      integer i 
      
      x =0.
      xoff=0.
      y=0.                                                                  
                                                                        
                                                                        
      eps = 1.e-2 
                                                                        
!==============================================================         
! initialize the model                                                
! and set the first guess of the control variables                      
!==============================================================         
      
      
      call initfunc( n, x ) 

      write(*,*) 'weight matrix is:' 
      write(*,*) weight              
                                                                        
!==============================================================         
! compute function twice                                                
!==============================================================         
      call func( n, x, m, y ) 
      print 40,'after 1st function call :' 
      print 30,'i','y(i)' 
      do i=1,m 
         print 20, i,y(i) 
      enddo 
      print 10,'sum = ',sum(y) 
                                                                        
      y=0. 
      call func( n, x, m, y ) 
      print 40,'after 2nd function call :' 
      print 30,'i','y(i)' 
      do i=1,m 
         print 20, i,y(i) 
      enddo 
      print 10,'sum = ',sum(y) 
                                                                        
      y=0. 
      call func( n, x, m, y ) 
      print 40,'after 3rd function call :' 
      print 30,'i','y(i)' 
      do i=1,m 
         print 20, i,y(i) 
      enddo 
      print 10,'sum = ',sum(y)                                                                         
                                                                        
      x = x + eps 
      y=0. 
      call func( n, x, m, y ) 
      print 50,'with x perturbed by :',eps 
      print 30,'i','y(i)' 
      do i=1,m 
         print 20, i,y(i) 
      enddo 
      print 10,'sum = ',sum(y) 
                                                                        
   10 format (a10,2x,e24.17) 
   20 format (i10,2x,e24.17) 
   30 format (a10,2x,a24) 
   40 format (a) 
   50 format (a,2x,e12.5) 
                                                                        
!==============================================================         
! postprocessing set                                                    
!==============================================================         
      call postfunc( n, x, m, y) 
      END                                           
