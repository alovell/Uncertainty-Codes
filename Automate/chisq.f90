   subroutine chisquare(nsets,ndata,parm,nfit,chi)
   implicit none
   real*8, intent(inout) :: chi
   integer, intent(in) :: nsets,ndata(3),nfit
   real*8, intent(in) :: parm(nfit)
   integer i,j
   real*8 err
   
   ! what if we just calculate the residuals instead?
   open(unit=1,file="parmfile.txt",status='replace')
   write(1,21) parm
   close(1)
   
   !run fdjacobian for error.dat
   call system("../bin/residuals -x_in parmfile.txt -fd_h 0.001 -search_file search-90Zrd.in > fdout.out")
   
   open(unit=2,file="errors.dat")
   do i=1,nsets
      do j=1,nfit
         read(2,*) ! do nothing, headers with parameters names
      enddo
      do j=1,ndata(i+1)
         read(2,*) err
	 chi = chi + err**2
      enddo 
   enddo
   
   close(2) 
   
21 format(8f12.6)
   end subroutine