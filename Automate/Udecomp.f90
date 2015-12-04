   ! calculates the upper diagonal matrix U for a symmetric matrix, cov
   ! such that transpose(U)U = cov
   
   subroutine udecomp(cov,nf)
   implicit none
   integer, intent(in) :: nf
   real*8, intent(inout) :: cov(nf,nf)
   integer i,j,k
   real*8 temp,orig(nf,nf)
   
   !print *, cov
   
   do i=1,nf
      ! calculate diagonal elements
      temp = 0
      do k=1,i-1
         temp = temp + cov(k,i)**2
	 !print *, i,i,"covki",cov(k,i),temp
      enddo 
      cov(i,i) = sqrt(cov(i,i)-temp)
      !print *, cov(i,i)
      do j=i+1,nf
         ! calculate off-diagonal elements
         temp = 0
	 do k=1,i-1
	    temp = temp + cov(k,i)*cov(k,j)
	    !print *, i,j,"covki",cov(k,i),"covkj",cov(k,j),temp
	 enddo
	 cov(i,j) = (cov(i,j) - temp)/cov(i,i) 
	 cov(j,i) = 0.d0
	 !print *, i,j,"covij",cov(i,j)
      enddo 
   enddo 
   
   !do i=1,nf
   !   print *, (cov(i,j),j=1,nf)
   !enddo 
   
   !orig = matmul(transpose(cov),cov)
   !print *, cov
   !print *, orig
   
   end subroutine udecomp