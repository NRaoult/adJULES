!-*-f90-*-
   module observation
    use inout, ONLY : nts, weight, stdIn
    use fomod, only : obsts,modts,l_mask,nloopcount

!Luke
    use forcing, only : sw_down
!Luke

        implicit none
        integer, parameter :: umodts = 2001 ! unit for modelled timeseries
	integer, parameter :: uobsts = 2002 ! unit for observed timeseries
        integer, parameter :: uavmodts = 2003 ! unit for averaged modelled timeseries
	integer, parameter :: uavobsts = 2004 ! unit for averaged observed timeseries
        real :: cost
	integer icall ! number of calls
      contains 
      
      
        subroutine initmodts
          integer :: recl
          real :: rdef
!          inquire(iolength=recl) rdef
          recl=8 !Jupp hack
	  close(umodts) ! testing
          open (umodts,file='modts.da',ACCESS='DIRECT',RECL=recl)
          cost = 0.
        end subroutine initmodts
	

        subroutine initobsts
          integer :: recl
          real :: rdef
!          inquire(iolength=recl) rdef
          recl=8 !Jupp hack
	  close(uobsts) ! testing
          open (uobsts,file='obsts.da',ACCESS='DIRECT',RECL=recl)
          cost = 0.
        end subroutine initobsts


        subroutine initavmodts
          integer :: recl
          real :: rdef
!          inquire(iolength=recl) rdef
          recl=8 !Jupp hack
	  close(uavmodts) ! testing
          open (uavmodts,file='avmodts.da',ACCESS='DIRECT',RECL=recl)
          cost = 0.
        end subroutine initavmodts
	

        subroutine initavobsts
          integer :: recl
          real :: rdef
!          inquire(iolength=recl) rdef
          recl=8 !Jupp hack
	  close(uavobsts) ! testing
          open (uavobsts,file='avobsts.da',ACCESS='DIRECT',RECL=recl)
          cost = 0.
        end subroutine initavobsts






        subroutine writemodts(step)
          integer :: step,i
	  do i=1,nts
             if (l_mask) then
               if (sw_down(1,1) < 20.0) modts(step, i)=-9999.99 !Luke mask nighttime fluxes
             end if
             write(unit=umodts,REC=i+nts*(step-1)) real(modts(step,i),8)
	  end do
        end subroutine writemodts



        subroutine writeavts(nwin,avmodts,avobsts)
          integer :: nwin,i,j
          real, intent(in) :: avmodts(nwin,nts),avobsts(nwin,nts)
	  do i=1,nwin
            do j=1,nts
              write(unit=uavmodts,REC=j+nts*(i-1)) real(avmodts(i,j),8)
              write(unit=uavobsts,REC=j+nts*(i-1)) real(avobsts(i,j),8)
	    end do
          end do
        end subroutine writeavts
 


	
	
        subroutine addcost(step,modts) ! read obsts from file and calculate cost
          integer :: step,i,j
          real :: modts(nts), obsts(nts)
          real :: temp 
          do i=1,nts
!             write(*,*) 'step is ',step
             read(unit=uobsts,REC=i+nts*(step-1)) temp
             obsts(i) = temp
	  end do
     

          do i=1,nts
             temp = 0.0
	     do j=1,nts
	        if (obsts(j) .ne. -9999.99) then                   ! this tests for NaNs...
!                   write(*,*) j,obsts(j),modts(j)
		   temp = temp + weight(i,j) * (obsts(j)-modts(j)) ! more robust than isnan()
	        end if
	     end do
             cost = cost + (temp)**2 
	  end do
        end subroutine addcost


	
        subroutine getobs ! read obsts from file
          integer :: step,i,j
          real :: obsts_read(nts)
          do i=1,nts
             do step=1,nloopcount
             read(unit=uobsts,REC=i+nts*(step-1)) obsts_read(i)
             obsts(step,i)=obsts_read(i)
             if (l_mask) then
               if (sw_down(1,1) < 20.0) obsts(step,i)=-9999.99 !Luke mask nighttime fluxes
             end if
            end do
	  end do

        end subroutine getobs



	
        subroutine postmodts
          close (umodts) 
        end subroutine postmodts
	
        subroutine postobsts
          close (uobsts) 
        end subroutine postobsts

        subroutine postavmodts
          close (uavmodts) 
        end subroutine postavmodts
	
        subroutine postavobsts
          close (uavobsts) 
        end subroutine postavobsts
	
		
      end module observation
