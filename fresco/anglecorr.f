   ! calculates the correlation matrix between all of the angles
   ! really prints a matrix that can be used to calculate the correlations
      program anglecorr
      implicit none
      integer i,j,jcount,k,io
      integer nparm,nsets,totdata
      real*8 angle
      character*50 sfile
      integer, allocatable :: ndata(:)
      real*8, allocatable :: angcorr(:,:)
      real*8, allocatable :: fitparm(:,:)
   
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
      allocate (fitparm(nparm,200))
      open(unit=4,file="pulls.txt")
      do k=1,200
         read(4,*) fitparm(:,k)
         !print *, fitparm(:,k)
      enddo
      close(4)
      
      do k=1,200	 
	 open(unit=10,file="ffile.txt",status='replace')
	 write(10,8) fitparm(:,k)
	 close(10)
      
         call sfresco_printj(sfile)
	 
	 close(16)
      
         open(unit=3,file="fort.16")
	 rewind(3)
	 j=1
	 do 
	    read(3,*,iostat=io) angle,angcorr(k+2,j)
	    if (io > 0) then
	       ! wrong character, continue
	       continue
	    else if (io < 0) then
	       ! end of file
	       exit
	    else
	       ! number, increment j
	       !print *, j, angcorr(k+2,j)
	       j = j + 1
	    end if 
	 enddo
	 close(3)
	 
	 !do i=1,9
	 !   read(3,*) ! xmgrace overall headers
	 !enddo 
         !jcount=1
         !do i=1,nsets
         !   do j=1,4
	 !      print *, j
         !      read(3,*) ! xmgrace headers
         !   enddo 
         !   do j=jcount,jcount-1+ndata(i)
	 !      print *, j
         !      read(3,*) angle,angcorr(k+2,j)
	 !      !print *, j, angcorr(k+2,j)
         !   enddo 
         !   jcount = jcount + ndata(i)
	 !   read(3,*) !END
         !enddo 
         !close(3)
      enddo
      
      ! print the values to a txt file
      open(unit=11,file="angcorr.txt",status='replace')
      do i=1,totdata
         write(11,7) angcorr(:,i)
      enddo 
      close(11)
   
   ! formating
7     format(202f12.6)
8     format(8f12.6)
   
      end program