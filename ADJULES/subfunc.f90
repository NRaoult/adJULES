!**************************************************************         
!     file: 	        subfunc.f90                                        
!     purpose: 	        driver subroutine to evaluate function                                                                
!**************************************************************         
      subroutine subfunc 
      use mo_n 
      implicit none
                                                                        
      call setfunc( n, m ) 

      call rjulesdoit( n, m )                                                                   
      end subroutine subfunc 
                                                                        
                                                                      
      subroutine rjulesdoit( n, m ) 
      use inout, only : weight,nts ,xoff                                                                        

      implicit none 
      integer :: n, m 
      real :: x(n), y(m)
      call initfunc( n, x)
      write(*,*) 'called initfunc' 
      call setsmooth(2.0)
      write(*,*) 'called setsmooth' 
      call func( n, x, m, y )
      write(*,*) 'called func' 
      call postfunc( n, x, m, y)
      write(*,*) 'called postfunc'  
      end subroutine rjulesdoit   
      
                                             
