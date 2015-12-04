   ! calculates the correlation matrix between all of the angles
   ! values in matrix values(pull#,angle0-180)
   program anglecorr
   implicit none
   integer i,j,jcount,k
   integer nparm,nsets,totdata
   real*8 angle
   character*50 sfile
   integer, allocatable :: ndata(:)
   real*8, allocatable :: angcorr(:,:)
   real*8, allocatable :: fitparm(:)
   
   ! need to run res to get error.dat for each of the pulls in pulls.txt
   
   ! get necessary information to set everything up
   ! number of fit parameters, number of data sets, points in each data set
   write(*,*) "How many parameters were fit?"
   read(*,*) nparm
   write(*,*) "How many data sets were fit?"
   read(*,*) nsets
   allocate (ndata(nsets))
   do i=1,nsets
      write(*,*) "How many data points in set",i,"?"
      read(*,*) ndata(i)
   enddo 
   write(*,*) "What is the name of the search file?"
   read(*,*) sfile
   
   ! set up matrix of values
   totdata = 0
   do i=1,nsets
      totdata = totdata + ndata(i) ! total number of data points
   enddo 
   
   allocate (angcorr(202,totdata))
   
   ! read in the experimental values
   open(unit=2,file="xexp.txt")
   do i=1,totdata
      read(2,*) angle,angcorr(1,i),angcorr(2,i)
   enddo 
   close(2)
   
   ! read in errors.data for each pull
   allocate (fitparm(nparm))
   open(unit=4,file="pulls.txt")
   do k=1,200
      read(4,*) fitparm
      !print *, fitparm
      
      call sfresco_printj(sfile)
      
      open(unit=3,file="errors.dat")
      jcount=1
      do i=1,nsets
         do j=1,nparm
            read(3,*) ! parameter values
         enddo 
         do j=jcount,jcount-1+ndata(i)
            read(3,*) angcorr(3,j)
	    !print *, j, angcorr(3,j)
         enddo 
         jcount = jcount + ndata(i)
      enddo 
      close(3)
   enddo
   close(4)
   
   ! formating
7  format(202f12.6)
8  format(8f12.6)
   
   end program