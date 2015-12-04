   subroutine confidence(nparm,npulls,pullfile)
   implicit none
   integer, intent(in) :: nparm,npulls
   character (len=20), intent(in) :: pullfile
   integer i,j,nshift,per
   real*8 parm(nparm),data(2),xsval(npulls),angle,temp
   character (len=50) line
   character (len=15) file
   
   !program to run fresco and get histograms for each angle
   !nparm = 6
   !npulls = 200
   
   !print *, pullfile
   
   ! new file for all of the elastic
   open(unit=8,file="elastic.201",status="replace")
   close(8)
   ! new file for all of the transfer
   open(unit=10,file="transfer.202",status="replace")
   close(10)
   ! new file for histograms per angle
   do i=0,180
      write(file,'("angle_el",I0)') i
      open(unit=11,file=file,status="replace")
      write(11,*) i
      close(11)
      write(file,'("angle_tr",I0)') i
      open(unit=12,file=file,status="replace")
      write(12,*) i
      close(12)
   enddo 
   
   do i=1,npulls
      ! make sure the temp input file is new
      open(unit=4,file="newfit.in",status="replace")
      close(4)
      open(unit=2,file=pullfile)
      read(2,*) parm
      !print *, parm(1),parm(2),parm(3),parm(4),parm(5),parm(6)
      
      ! remake input file
      call makeinputfile(parm)
      
      ! run fresco - not sure if this is the best way to do this
      call system("fresco < newfit.in > newfit.out")
      
      ! save the elastic and transfer
      open(unit=7,file="fort.201")
      open(unit=8,file="elastic.201",access="append")
      open(unit=9,file="fort.202")
      open(unit=10,file="transfer.202",access="append")
      write(8,*) "#",parm(1),parm(2),parm(3)
      write(8,*) "#",parm(4),parm(5),parm(6)
      write(10,*) "#",parm(1),parm(2),parm(3)
      write(10,*) "#",parm(4),parm(5),parm(6)
      do j=1,192
         read(7,'(A)') line
	 write(8,'(A)') line
	 read(9,'(A)') line
	 write(10,'(A)') line
      enddo 
      close(7)
      close(8)
      close(9)
      close(10)
   enddo 
   close(2)
   
   ! histogram the "data" for each angle (elastic)
   open(unit=8,file="elastic.201")
   rewind(8)
   do i=1,npulls
      !print *, i
      do j=1,12
         !print *, j
         read(8,*) 
      enddo 
      do j=0,180
         write(file,'("angle_el",I0)') j
	 !print *, file
	 read(8,*) data(1),data(2)
	 open(unit=11,file=file,access="append")
	 write(11,*) data(2)
	 close(11)
      enddo 
      read(8,*) 
   enddo
   close(8)
   
   ! sort the array and "throw out" the 2.5% highest and 2.5% lowest
   !define the 2.5% highest/lowest
   per = nint(npulls*0.025)
   !print *, per
   
   open(unit=12,file="band_el.txt")
   do i=0,180
      !print *, i
      write(file,'("angle_el",I0)') i
      open(unit=11,file=file)
      read(11,*) angle
      do j=1,npulls
         read(11,*) xsval(j)
      enddo
      ! sort xsval
      nshift = 0
      do while (nshift/=npulls-1)
         nshift = 0
	 !print *, nshift
         do j=1,npulls-1
	    if (xsval(j) > xsval(j+1)) then
	        temp = xsval(j)
		xsval(j) = xsval(j+1)
		xsval(j+1) = temp
		nshift = 0
		!print *, nshift
	    else 	
		nshift=nshift + 1
	    end if 
	 enddo 
      enddo 
      
      !do j=1,npulls
      !   print *, xsval(j)
      !enddo 
      close(11)
      write(12,*) i,xsval(per),xsval(npulls-per)
   enddo 
   close(12)
   
   ! histogram the "data" for each angle (transfer)
   open(unit=8,file="transfer.202")
   rewind(8)
   do i=1,npulls
      !print *, i
      do j=1,12
         !print *, j
         read(8,*) 
      enddo 
      do j=0,180
         write(file,'("angle_tr",I0)') j
	 !print *, file
	 read(8,*) data(1),data(2)
	 open(unit=11,file=file,access="append")
	 write(11,*) data(2)
	 close(11)
      enddo 
      read(8,*) 
   enddo
   close(8)
   
   ! sort the array and "throw out" the five highest and five lowest
   open(unit=12,file="band_tr.txt")
   do i=0,180
      !print *, i
      write(file,'("angle_tr",I0)') i
      open(unit=11,file=file)
      read(11,*) angle
      do j=1,npulls
         read(11,*) xsval(j)
      enddo
      ! sort xsval
      nshift = 0
      do while (nshift/=npulls-1)
         nshift = 0
	 !print *, nshift
         do j=1,npulls-1
	    if (xsval(j) > xsval(j+1)) then
	        temp = xsval(j)
		xsval(j) = xsval(j+1)
		xsval(j+1) = temp
		nshift = 0
		!print *, nshift
	    else 	
		nshift=nshift + 1
	    end if 
	 enddo 
      enddo 
      
      !do j=1,npulls
      !   print *, xsval(j)
      !enddo 
      close(11)
      write(12,*) i,xsval(per+1),xsval(npulls-per)
   enddo 
   close(12)
   
   end subroutine confidence