!**************************************************************         
!     file: 	        prgfunc.f90                                        
!     purpose: 	        driver to evaluate function                     
!     Copyright (C) 2000 - 2006                                         
!     FastOpt, Drs. Ralf Giering und Thomas Kaminski GbR                
!     , Hamburg, Germany                                                
!                                                                       
!     Email: info@FastOpt.de                                            
!                                                                       
!     All rights reserved.                                              
!**************************************************************         
      program main 
      use mo_n 
      implicit none                            

                                                                        
      call setfunc( n, m ) 



      call doit( n, m )                                                                   
      end program main 
                                                                        
                                                                      
      subroutine doit( n, m ) 
      use inout, only : weight,nts                                                                         

      implicit none 
      integer :: n, m 
      real :: x(n), y(m)
      call initfunc( n, x) 

      write(*,*) 'weight matrix is:' 
      write(*,*) weight 

      call func( n, x, m, y ) 
      call postfunc( n, x, m, y) 
      end subroutine doit                                          
