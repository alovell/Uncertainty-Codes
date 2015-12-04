      subroutine chisq_calc(parm,nparm,sfile,ndat,nset,fval,wij)
      implicit none
      integer, intent(in) :: nparm,nset
      integer, intent(in) :: ndat(nset)
      real*8, intent(in) :: parm(nparm)
      real*8, intent(in) :: wij(40,40,nset)
      character*50, intent(in) :: sfile
      real*8, intent(inout) :: fval
      integer i,j,k
      real*8 err(40)
      logical ok
      
      ! get the parameter values from parm_grid.txt
      ! put the current vec(x) values in xfile.txt
      !print *, parm
      open(unit=12,file="ffile.txt",status='replace')
      write(12,21) parm
      close(12)
	 
      ! run fdjacobian to get the errors
      ! give errors.dat

      !call system("../../bin/residuals -x_in ffile.txt -fd_h 0.001 
      !&        -search_file " // sfile // " > fdout.out")
      
      call sfresco_printj(sfile)
      
      fval=0.d0
      !reads in errors.dat and calculates a chi square value
      ! check to see if errors.dat is already open
      ! close it if it already is
      inquire(file="errors.dat",opened=ok)
      if (ok) close(400)
      
      open(unit=10,file="errors.dat")
      
      !print *, nset,nparm

      fval = 0.d0
      do j=1,nset
      !print *, "j=",j
         do k=1,nparm
	    read(10,*) !do nothing for header lines
	    !print *, "k=",k
	 enddo 
         i=0
         !print *, i,ndat(j)
         do while (i < ndat(j))
            !print *, "i=",i
            read(10,*) err(i)
            !fval = fval + err**2
            !print *, fval, err
            i=i+1
         enddo
	 do i=1,ndat(j)
	 do k=1,ndat(j)
	    fval=fval+wij(i,k,j)*err(i)*err(j)
	    print *, fval
	 enddo
	 enddo
      enddo 
      fval = sqrt(fval)
      !print *, fval
      
      close(10)

21    format(8f12.6)
      
      end subroutine chisq_calc