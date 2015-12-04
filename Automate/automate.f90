   ! this program is designed to automate the fitting/analysis process
   ! needs to be updated to deal with multiple data sets
   
   program greasedlightnin
   implicit none
   integer ndata(3),io,nfit,i,j,k,nsets
   real*8 s2,rand,temp,y1,y2,pi,w
   character (len=20) fdata,fsearch,fsfres,fbf,fpulls
   integer npulls, run(3)
   real*8, allocatable :: bestfit(:)
   real*8, allocatable :: residual(:,:)
   real*8, allocatable :: covariance(:,:)
   real*8, allocatable :: pulls(:,:)
   
   ! read input file
   open(5)
   read(5,*) fdata
   read(5,*) fsearch
   read(5,*) fsfres
   read(5,*) fbf
   read(5,*) npulls
   read(5,*) fpulls
   read(5,*) run(1)
   read(5,*) run(2)
   read(5,*) run(3)
   !print*, fdata,fsearch,fsfres,fbf,npulls,fpulls,run(1),run(2),run(3)
   close(5)
   
   ! count the number of data points
   ndata = 0
   !open(unit=22,file="xexp.txt")
   !open(unit=22,file=fdata)
   !do 
   !   read(22,*,iostat=io)
   !   if (io /= 0) exit
   !   ndata = ndata + 1
   !enddo
   !close(22)
   
   ! get the number of fitted parameters
   !open(unit=23,file="search.in")
   open(unit=23,file=fsearch)
   read(23,*) ! file names
   read(23,*) nfit
   close(23) 
   
   nsets = 1
   ndata(1) = 0
   ndata(2) = 25
   ndata(3) = 0
   
   !print *, ndata,nfit,nsets
   
   ! allocate best fit array
   allocate (bestfit(nfit)) 
   allocate (residual(3,ndata(2)+ndata(3)))
   
   ! run sfresco?
   if (run(1)==0) then 
   
   ! find the inital fit
   ! make sure that PETSC_DIR and PETSC_ARCH are set
   ! needs sfresco.in, search.in, frescoinput.in
   call system("../bin/sfresco -tao_gatol 1e-12 -tao_frtol 0 -tao_fatol 0 -tao_grtol 0 -tao_pounders_delta 0.01 < " // fsfres // " > sfresco.out")
   ! outputs F and X files
   ! X gives parameter values
   ! F gives the residuals
   
   ! put the best fit data into x.in
   ! runs no matter what
   open(unit=24,file="X****.txt")
   read(24,*) ! header line
   read(24,*) (bestfit(i),i=1,nfit)
   close(24)
   !open(unit=27,file="x.in")
   open(unit=27,file=fbf)
   write(27,"(6f15.11)") (bestfit(i),i=1,nfit)
   close(27)       
   
   print *, "sfresco ran!  Output in sfresco.out  Fitting parameters in search.plot"
   
   end if   
   
   ! calculate the jacobian?
   if (run(2) == 0) then
   
   ! run fdjacobian to calculate the jacobian
   ! needs search.in and file with the best fit
   !call system("fdjacobian -x_in " // fbf // " -fd_h fdh -search_file " // fsearch // "> fdout.out")
    call system("fdjacobian -x_in " // fbf // " -fd_h 0.0001 -search_file " // fsearch // "> fdout.out")
   ! outputs correlation.txt and jacobian.txt
   
   print *, "Jacobian calculated.  Output in jacobian.txt and correlation.txt"
   
   end if 
   
   ! make the draws?
   if (run(3)==0) then
   
   ! calculate s^2
   ! read in residuals - errors.dat
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
   
   s2 = 0
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
   
   !npulls = 200 ! number of draws
   allocate (pulls(nfit,npulls))
   
   call init_random_seed()
   
   ! draw from multivariate normal distribution by calculating
   ! x = mu + R^T Y
   ! where R^T R = Covariance
   ! mu is the best fit vector
   ! Y is a nfit x npulls matrix of random numbers
   ! with guidance from 
   ! http://people.sc.fsu.edu/~jburkardt/f_src/normal_dataset/normal_dataset.html
   
   ! create the Y matrix
   !call random_seed()
   !call random_number(pulls)
   !do i=1,nfit
   !   do j=1,npulls
   !      pulls(i,j) = sqrt(-2.d0*log(pulls(i,j)))
   !	 call random_number(rand)
   !	 if (rand > 0.5) then
   !	    pulls(i,j) = -pulls(i,j)
   !	 end if 
   !   enddo 
   !enddo 
   
   ! better way to create Y matrix
   ! http://www.design.caltech.edu/erik/Misc/Gaussian.html
   ! define pi
   pi = 4.d0*atan(1.d0)
   call random_number(pulls)
   !do i=1,npulls
   !print *, pulls(2,i)
   !enddo
   !print *, pulls(4,:)
   do i=1,nfit
      do j=1,npulls/2
         y1 = sqrt(-2.d0*log(pulls(i,2*j-1)))*cos(2.d0*pi*pulls(i,2*j))
   	 y2 = sqrt(-2.d0*log(pulls(i,2*j-1)))*sin(2.d0*pi*pulls(i,2*j))
   	 pulls(i,2*j-1) = y1
   	 pulls(i,2*j) = y2
	 !print *, y1,y2
      enddo 
   enddo 
   
   !print *, " "
   !do i=1,npulls
   !print *, pulls
   !enddo 
   
   !do i=1,nfit
   !do j=1,npulls/2
   !do while (w > 1.d0)
   !   call random_number(pulls(i,2*j-1))
   !   call random_number(pulls(i,2*j))
      !print *, pulls(i,2*j-1),pulls(i,2*j),i,j
      !pulls(i,2*j-1) = 2*pulls(i,2*j-1) - 1
      !pulls(i,2*j) = 2*pulls(i,2*j) - 1
      !print *, pulls(i,2*j-1),pulls(i,2*j),i,j
   !   y1 = 2*pulls(i,2*j-1) - 1
   !   y2 = 2*pulls(i,2*j) - 1
      !print *, y1,y2, i,j
      !w = pulls(i,2*j-1)**2 + pulls(i,2*j)**2
   !   w = y1**2 + y2**2
      !print *, w
   !enddo 
   !w = sqrt((-2.d0*log(w))/w)
   !print *, y1,y2,i,j
   !pulls(i,2*j-1) = pulls(i,2*j-1)*w
   !pulls(i,2*j) = pulls(i,2*j)*w
   !pulls(i,2*j-1) = y1*w
   !pulls(i,2*j) = y2*w
   !enddo 
   !enddo
   
   pulls = matmul(transpose(covariance),pulls)
   
   if (run(1)==1) then
      open(unit=31,file=fbf)
      read(31,*) bestfit
      close(31)
   end if

      
   open(unit=30,file=fpulls)
   do i=1,npulls
      write(30,"(8f15.8)") (bestfit(j) + pulls(j,i),j=1,nfit)
   enddo 
   close(30)
 
   ! define confidence bands
   ! needs the file 200pulls.txt
   !call system("findconfidence")  
   !call confidence(nfit,npulls,"200pulls.txt")
   call confidence(nfit,npulls,fpulls)
   ! creates the files angle_el*, angle_tr*, elastic.201, transfer.202
   ! band_el.txt and band_tr.txt
   
   ! calculate the best fit curves
   open(unit=31,file="newfit.in",status="replace")
   close(31)
   call makeinputfile(bestfit,nfit)
   call system("fresco < newfit.in > newfit.out")
   ! now the best fit calculations will be in fort.201 and fort.202
   ! for elastic and transfer respectively
   
   print *, npulls,"pulls drawn.  Elastic band in band_el.txt, best fit in fort.201.  Transfer band in band_tr.txt, best fit in fort.202."
   
   end if
   
   
   end program greasedlightnin