      program grid
      implicit none
      integer ip,nparm,stat,it,ib,iq,is,nset,id
      character*50 grid_file,sfile
      character*50,allocatable::cfile(:)
      real*8 chi
      real*8, dimension(:,:), allocatable::corr_matrix
      real*8, dimension(:), allocatable::sig
      real*8, dimension(:), allocatable::start
      real*8, dimension(:), allocatable::fin
      real*8, dimension(:), allocatable::step
      real*8, dimension(:), allocatable::bestfitparm
      real*8, dimension(:), allocatable::parm
      integer, allocatable :: ndat(:)
      real*8, allocatable :: wij(:,:,:)
      logical file_exists
      
      write(*,*) "How many parameters were fit?"
      read(*,*) nparm
      write(*,*) "How many data sets?"
      read(*,*) nset
      allocate(ndat(nset))
      allocate(cfile(nset))
      do ip=1,nset
         write(*,*) "Number of data points in set", ip, "?"
	 read(*,*) ndat(ip)
         write(*,*) "Name of model correlation file (or not)?"
         read(*,*) cfile(ip)
	 ! check that file exists
	 inquire(file=cfile(ip),exist=file_exists)
	 if (file_exists .neqv. .TRUE.) then
	    print *, cfile(ip), "does not exist, no model correlations"
	    cfile(ip) = "not"
	 end if
	 ! if no mc file, type "not"
      enddo 
      write(*,*) "Name of the search file?"
      read(*,*) sfile
      ! check that sfile exists
      inquire(file=sfile,exist=file_exists)
      if (file_exists .neqv. .TRUE.) then
         print *, sfile, "does not exist; try again"
	 stop
      end if 
      !print *, nparm
      !nparm=8
      
      allocate(corr_matrix(nparm,nparm))
      allocate(sig(nparm))
      allocate(start(nparm))
      allocate(fin(nparm))
      allocate(step(nparm))
      allocate(bestfitparm(nparm))
      allocate(parm(nparm))
      allocate(wij(40,40,nset))
      grid_file='parm_grid.txt'
      
      ! read in correlation matrix
      open(unit=502,file="correlation.txt")
      read(502,*) !don't want first line
      read(502,*) !don't want second line
      do ip=1,nparm
        read(502,*,iostat=stat) corr_matrix
	if(stat /= 0) exit
      enddo 
      close(502)
      !print *, corr_matrix

      ! read in best fit parameters
      open(unit=503,file="x.in")
      !read(503,*) !header
      read(503,*) bestfitparm
      close(503)
      
      ! read in the wij matrix if there is one
      wij = 0.d0
      do id=1,nset
         if (cfile(id)=="not") then
	    do ip=1,ndat(id)
	       wij(ip,ip,id) = 1.d0
	    enddo
	 else
	    open(unit=503,file=cfile(id))
	    do ip=1,ndat(id)
	       !do iq=1,ndat(id)
	          !print *, ip,iq
	          read(503,*) (wij(ip,iq,id),iq=1,ndat(id))
		  !print *, wij(ip,:,id)
	       !enddo
	    enddo
	    close(503)
	 end if 
      enddo 
      
      ! calculate sigma values
      do ip=1,nparm
         sig(ip)=sqrt(corr_matrix(ip,ip))
	 ! grid based on standard deviation
	 start(ip)=bestfitparm(ip)-sig(ip)
	 fin(ip)=bestfitparm(ip)+sig(ip)
         step(ip)=2*sig(ip)/20.0
         ! grid based on 20% of value
         !start(ip)=bestfitparm(ip)-0.2*bestfitparm(ip)
         !fin(ip)=bestfitparm(ip)+0.2*bestfitparm(ip)
         !step(ip)=(fin(ip)-start(ip))/20.d0
         !step(ip) = 4*sig(ip)
         !print *, sig(ip),start(ip),fin(ip),step(ip)
         parm(ip)=bestfitparm(ip)
      enddo 
      !print *, "sig:",sig
      !print *, "parm:",parm
      !print *, "start:",start
      !print *, "fin:",fin
      !print *, "step:",step

      ! create grid of points and print them to a file
      ! will probably be read by bash
      open(504,file=grid_file,access="append")
      do it=1,nparm
	do ib=it+1,nparm
	   write(504,*) it,ib
	   print *, it,ib
	   !write(504,*) 
	   parm(it)=start(it)
	   !do is=1,21
	   do while (parm(it) <= fin(it))
	     parm(ib)=start(ib)
             !do ip=1,11
	     do while (parm(ib) <= fin(ib))
	       ! calculate the chi square value here
	       call chisq_calc(parm,nparm,sfile,ndat,nset,chi,wij)
	       !print *, nparm
	       write(504,20) parm,chi
               parm(ib)=parm(ib)+step(ib)
             enddo
	     parm(it)=parm(it)+step(it)
	     parm(ib)=bestfitparm(ib)
	     !print *, "ib,nparm", ib,nparm
	   enddo
	enddo 
	parm(it)=bestfitparm(it)
	!print *, "it,nparm",it,nparm
      enddo 
      close(504)
            
!20    format(4x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8)
! need to make this line more fluid (not hard coded)
20     format(10f12.6)
      end program grid