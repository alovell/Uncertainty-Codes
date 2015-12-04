   program unweightjac
   
   implicit none
   real*8 sig(33,33),jac(33,6),uwjac(33,6),temp(2)
   integer i,j
   
   ! read in Jacobian
   open(unit=2,file="jacobian.txt")
   read(2,*) ! "Mat Object" line
   read(2,*) ! "type" line
   read(2,*) jac
   close(2)
   !print *, jac
   
   ! read in covariance matrix (for sigma)
   sig = 0.d0
   open(unit=1,file="experr.txt")
   do i=1,33
      read(1,*) temp(1),temp(2),sig(i,i)
      sig(i,i) = sig(i,i)*temp(2)
   enddo 
   close(1)
   
   ! multiply sigma and jacobian
   uwjac = matmul(sig,jac)
   print *, uwjac
   
   end program