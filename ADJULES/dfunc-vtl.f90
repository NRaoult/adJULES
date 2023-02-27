!===============================================
! this subroutine evaluates the gradient using the vector tangent
!===============================================

      subroutine dfuncvtl( x, grady )
      use mo_n
      use mo_h
      use inout, only : echo, nts, weight
      implicit none
      real    x(n), x_tl(n,n), y(m), y_tl(n,m) 
      real    grady(n)
      character*23 outfile,itername
      integer i, icount
      save icount
      call setfunc(n,m) ! added by Jupp to set m,n    

! call tangent 
      icount=icount+1
      x_tl  = 0.0

      do i=1, n
         x_tl(i,i) = 1.
      enddo

      y_tl = 0.0
      y    = 0.0

      call func_tl( n, n, x, x_tl, n, m, y, y_tl, n) 
      call flush(6)
      grady(1:n) = y_tl(1:n,1)  
      
      
!      write(*,*) 'grady is', m,n,y_tl(1,1),y_tl(2,1),y_tl(3,1)
!      write(*,*) 'grady is', m,n,grady(1),grady(2),grady(3)
!      print*, 'done'

! avoid stop due to overflow
      do i = 1,n
         if(.not.grady(i).lt.huge(grady(i))) then
            grady(i) = huge(grady(i))/1.e50
            print*, 'PRGDFPMIN: avoiding derivative overflow for component ',i
         endif
         if(.not.-grady(i).lt.huge(grady(i))) then
            grady(i) = -huge(grady(i))/1.e50
            print*, 'PRGDFPMIN: avoiding derivative underflow for component ',i
         endif
      enddo
! store current point 
!      ih = ih + 1
!      do i = 1,n
!         hx(i,ih)=x(i)
!         hg(i,ih)=grady(i)
!      enddo
!      hf(ih)=y(1)
!      if (echo) write(*,*) 'Opti: ',ih,hf(ih),sqrt(sum(hg(:,ih)**2))
! save current point 
  if (echo)       write(itername,'(i18)') ih 
  if (echo)       outfile='output_FO/opti-'//trim(adjustl(itername))//'.b'
  if (echo)       OPEN (39,file=outfile,status='unknown',form='unformatted')
  if (echo)       WRITE(39) x
  if (echo)       CLOSE(39)
  if (echo)       outfile='output_FO/opti-'//trim(adjustl(itername))//'.dat'
  if (echo)       OPEN (39,file=outfile,status='unknown',form='formatted')
  if (echo)       WRITE(39,'(10(1x,e12.6))') x
  if (echo)       CLOSE(39)
    end	
