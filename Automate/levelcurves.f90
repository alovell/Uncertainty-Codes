   program levelcurves
   implicit none
   integer i,j
   integer ndraws,nsets,nfit,ndata(3),npulls
   real*8 s2,pi,y1,y2
   real*8 newchi
   real*8, allocatable :: bestfit(:)
   real*8, allocatable :: residual(:,:)
   real*8, allocatable :: covariance(:,:)
   real*8, allocatable :: pulls(:,:)
   
   npulls = 200   
   nsets = 1
   nfit = 6
   ndata(1) = 0
   ndata(2) = 29
   ndata(3) = 0
   pi = 4.d0*atan(1.d0)
   
   ! set up the multivariate gaussian
   
   ! calculate s^2
   ! read in residuals - errors.dat
   allocate(residual(3,ndata(2)+ndata(3)))
   open(unit=28,file="errors.dat")
   do j=1,nsets
   do i=1,nfit
      read(28,*)
   enddo 
   do i=1,ndata(j+1)
      read(28,*) residual(3,i+ndata(j))
   enddo 
   enddo 
   close(28)
   !print *, "read errors.dat"
   
   s2 = 0.d0
   do i=1,ndata(2)+ndata(3)
      s2 = s2 + residual(3,i)**2
   enddo 
   s2 = s2/(ndata(2)+ndata(3)-nfit)
   print *, "s2 value is ", s2   
   
   ! calculate s^2 * covariance matrix
   
   allocate (covariance(nfit,nfit))
   open(unit=29,file="correlation.txt")
   read(29,*) ! "Mat Object" line
   read(29,*) ! "type" line
   read(29,*) covariance 
   close(29)
      
   ! construct/call multivariate Gaussian
   ! need best fit and covariance
   ! have to make the covariance matrix symmetric
   covariance = 0.5d0 * (covariance + transpose(covariance))
   ! s2*covariance (to draw from)
   covariance = s2*covariance
   
   ! make covariance into upper triangular
   call udecomp(covariance,nfit)

   allocate (pulls(nfit,2))
   
   call init_random_seed()
   
   allocate (bestfit(nfit))
   
   ! read best fit information      
   open(unit=31,file="x.in")
   read(31,*) bestfit
   close(31)
         
   ! don't append to files from the last run
   open(unit=32,file="levelpulls.txt",status="replace")
   close(32)
   
   ! draw until 200 draws are within the level curves  
   open(unit=32,file="levelpulls.txt")
   ndraws = 0
   !print *, npulls/2
   do while (ndraws < npulls)
      ! make two draws
      call random_number(pulls)
      do i=1,nfit
         !do j=1,2/2
	    y1 = sqrt(-2.d0*log(pulls(i,1)))*cos(2.d0*pi*pulls(i,2))
   	    y2 = sqrt(-2.d0*log(pulls(i,1)))*sin(2.d0*pi*pulls(i,2))
   	    pulls(i,1) = y1
   	    pulls(i,2) = y2
	    !print *, "i=",i,y1,y2
	 !enddo 
      enddo

      ! create the parameters to remake the inputfile
      do i=1,2
         do j=1,nfit
            pulls(j,i) = bestfit(j) + pulls(j,i)
	    !print *, pulls(j,i),i,j
	 enddo 
      enddo 
      ! remake the input file
      open(unit=31,file="newfit.in",status="replace")
      close(31)
      call makeinputfile(pulls(:,1),nfit)
      ! run fresco
      call system("fresco < newfit.in > newfit.out")
      
      ! calculate new chi square value
      newchi = 0.d0
      call chisquare(nsets,ndata,pulls(:,1),nfit,newchi)
      newchi = newchi/(ndata(2) + ndata(3) - nfit)
      print *, newchi
      
      ! check if the new chi square is within \pm 1 of the old chi square
      if (abs(s2 - newchi) <= 1.d0) then
         ndraws = ndraws + 1
	 print *, newchi
	 print *, pulls(:,1)
	 write(32,"(8f15.8)") (pulls(j,1),j=1,nfit)
      end if 
      
      ! do above for the second set of parameters
      open(unit=31,file="newfit.in",status="replace")
      close(31)
      call makeinputfile(pulls(:,2),nfit)
      call system("fresco < newfit.in > newfit.out")
      newchi = 0.d0
      call chisquare(nsets,ndata,pulls(:,2),nfit,newchi)
      newchi = newchi/(ndata(2) + ndata(3) - nfit)
      print *, newchi
      if (abs(s2 - newchi) <= 1.d0) then
         ndraws = ndraws + 1
	 print *, newchi
	 print *, pulls(:,2)
	 write(32,"(8f15.8)") (pulls(j,2),j=1,nfit)
      end if 

   enddo
   close(32)
   
   ! calculate the confidence bands
   print *, "calculating confidence bands"
   call confidence(nfit,npulls,"levelpulls.txt")
   
   end program