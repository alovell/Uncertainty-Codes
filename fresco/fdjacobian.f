      program sfresco_printj
      implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"
#include "finclude/petsctao.h"
      PetscBool      isset
      PetscErrorCode ierr
      integer   MAXNX,MAXNF,ip
      parameter (MAXNF=300,MAXNX=12)
      character*50 search_file
      character*50 xfile
      character*50 check
      character*50 chi_file
      PetscInt  nx,nf,i
      PetscReal startx(maxnx),h,one,Cscale
      PetscReal x(maxnx),f(maxnf)
      PetscReal x_v(0:1)
      PetscOffset x_i
      Vec       taoX,taoF
      Mat       taoJ,taoC,taoCinv,ident
      MatFactorInfo factorinfo
      Tao       tao 
      PetscViewer cviewer,jviewer
      external  FormFunction
      real*8 V,r,a,sig1,sig2,sig3,d1,d2,sta1,sta2,sta3,end1,end2,end3
      real*8 chitot,temp(12),step1,step2,step3,parm(12)
      real*8 expdat(40,3),thdat(40,2)
      real*8, dimension(:,:), allocatable::corr_matrix 
      real*8, dimension(:), allocatable::sig
      real*8, dimension(:), allocatable::start
      real*8, dimension(:), allocatable::finish
      real*8, dimension(:), allocatable::step
      real*8, dimension(:), allocatable::bestfitparm
      integer l,j,k,nlines,nemp,nfull,stat,test
      integer pm,pn
      logical lopen
      
      allocate (sig(6))
      allocate (start(6))
      allocate(finish(6))
      allocate(step(6))
      allocate(bestfitparm(6))

!     counts the number of data points
!     l=0 is a condition so that this only happens the first time (or when check=false)
      l=0
      nlines=0
      !open(502,file="xexp.txt",status="old")
      !do
      !   read(502,*,end=200)
!	 nlines=nlines+1
!      end do
!      close(502)
      
 200  continue
      search_file='search-102.in'
      xfile = 'x.in'
      chi_file='chi.txt'
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-search_file',       &
     &     search_file,isset,ierr)
      CHKERRQ(ierr)

      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-x_in',               &
     &     xfile,isset,ierr)
      CHKERRQ(ierr)
      
      !new condition if you want to compute the chi square for multiple points (then -condit true)
      !or just for one value
      check = 'false'
      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-condit',               &
     &     check,isset,ierr)
      CHKERRQ(ierr)
      
      !file to write chisq values to
      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-chi_file',       &
     &     chi_file,isset,ierr)
      CHKERRQ(ierr)

      h=1.0e-4
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-fd_h',h,isset       &
     &   ,ierr)
      CHKERRQ(ierr)
      
      !factor by which to scale the covariant matrix elements
      Cscale=1.0
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-scale',Cscale,       &
     &   isset,ierr)
      CHKERRQ(ierr)
          
      call sfresco_init(search_file,MAXNX,nx,nf,startx)
      
      call VecCreateSeq(PETSC_COMM_SELF,nx,taoX,ierr)
      CHKERRQ(ierr)
      call VecCreateSeq(PETSC_COMM_SELF,nf,taoF,ierr)
      CHKERRQ(ierr)
      
      call VecGetArray(taoX,x_v,x_i,ierr)
      open(111,file=xfile)
      read(111,*) (x_v(x_i+i),i=0,nx-1)
!      read(111,*) x_v(x_i+1)
!      read(111,*) x_v(x_i+2)
      rewind(111)
      close(111)
      print *,'computing fd jacobian at X=',(x_v(x_i+i),i=0,nx-1)

      call VecRestoreArray(taoX,x_v,x_i,ierr)


      call TaoCreate(PETSC_COMM_SELF,tao,ierr)
      CHKERRQ(ierr)
      call TaoSetType(tao,TAOPOUNDERS,ierr)
      CHKERRQ(ierr)
      call TaoSetInitialVector(tao,taoX,ierr)
      CHKERRQ(ierr)
      call TaoSetSeparableObjectiveRoutine(tao,taoF,FormFunction,
     & PETSC_NULL_OBJECT,ierr)
      CHKERRQ(ierr)
      call TaoSetFromOptions(tao,ierr)
      CHKERRQ(ierr)
      
      call MatCreateSeqDense(PETSC_COMM_SELF,nf,nx,PETSC_NULL_SCALAR,     &
     &   taoJ,ierr)
      CHKERRQ(ierr)
      call ComputeFDJacobian(tao,taoX,taoF,h,taoJ,ierr)
      CHKERRQ(ierr)
      call MatTransposeMatMult(taoJ,taoJ,MAT_INITIAL_MATRIX,              &
     &     PETSC_DEFAULT_REAL,taoCinv,ierr)
      CHKERRQ(ierr)
      call MatDuplicate(taoCinv,MAT_DO_NOT_COPY_VALUES,taoC,ierr)
      CHKERRQ(ierr)
      call MatDuplicate(taoCinv,MAT_DO_NOT_COPY_VALUES,ident,ierr)
      CHKERRQ(ierr)
      call MatZeroEntries(ident,ierr)
      CHKERRQ(ierr)
      one=1.0
      call MatShift(ident,one,ierr)
      CHKERRQ(ierr)

      if (check .ne. 'false') then
!     open the experimental cross section file, read points into expdat
      open(502,file="xexp.txt")
      rewind(502)
      do j=1,nlines
         read(502,*) expdat(j,1), expdat(j,2), expdat(j,3)
      end do
      close(502)
!      do k=1,nlines
!      print *, expdat(k,2)
!      end do
!     check if fort.16 is open and close it
!     read in theory points
      inquire(file="fort.16",opened=lopen)
      if (lopen) close(16)
      open(503,file="fort.16",status="old")
      do nemp=1,4
         read(503,*)
      end do
      do nfull=1,nlines
         read(503,*) thdat(nfull,1), thdat(nfull,2)
      end do
      close(503)
!      do k=1,nlines
!         print *, thdat(k,2)
!      end do
!     calculate chi square value (weighted w/ exp. error)
      chitot=0
      do k=1,nlines
         chitot=chitot+((expdat(k,2)-thdat(k,2))**2)/(expdat(k,3))**2
!         print *, expdat(k,2),thdat(k,2),expdat(k,3),chitot
      end do
      
      open(501,file=xfile)
      read(501,*) temp(1), temp(2), temp(3), temp(4), temp(5), temp(6)
      close(501)
      open(505,file=chi_file,access="append")
      write(505,*) (temp(ip),ip=1,6), chitot/6.00000
      close(505)
      end if
      
      !if computing chi sq for multiple points then don't calculate
      !the Jacobian - some give singular values
      if (l .eq. 1) then 
      go to 201
      end if
      
      call MatLUFactor(taoCinv,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,        &
     &     factorinfo,ierr)
      CHKERRQ(ierr)
      call MatMatSolve(taoCinv,ident,taoC,ierr)
      CHKERRQ(ierr)
      call PetscViewerASCIIOpen(PETSC_COMM_SELF,'jacobian.txt',            &
     &   jviewer,ierr)
      CHKERRQ(ierr)
      call PetscViewerASCIIOpen(PETSC_COMM_SELF,'correlation.txt',          &
     &   cviewer,ierr)
      print *
      print *,'Jacobian Matrix (saved in jacobian.txt)'
      call MatView(taoJ,jviewer,ierr)
      CHKERRQ(ierr)
      call MatView(taoJ,PETSC_VIEWER_STDOUT_SELF,ierr)
      CHKERRQ(ierr)
      print *
      print *,'Correlation Matrix inv(JT*J) (saved in correlation.txt)'
      call MatView(taoC,cviewer,ierr)
      CHKERRQ(ierr)
      call MatView(taoC,PETSC_VIEWER_STDOUT_SELF,ierr)
      CHKERRQ(ierr)
      

      
!      if (check .eq. 'false') then
!      call PetscFinalize()
!      end if


!     my additions to this program - calculating the chi squared function
      if (check .ne. 'false') then
      if(l .eq. 0)then
         allocate(corr_matrix(6,6))
         open(unit=502,file="bestfit.txt")
         open(unit=500,file="correlation.txt")
	 read(502,*)   !header line
         read(502,*) bestfitparm
         read(500,*)
         read(500,*)
	 do ip=1,6
           read(500,*,iostat=stat)corr_matrix
	   if (stat /= 0) exit
	 enddo 
	 close(500)
!         print *, V, r, a, sig1, sig2, sig3
         do ip=1,6
	   sig(ip) = sqrt(Cscale*corr_matrix(ip,ip))
	   start(ip) = bestfitparm(ip) - sig(ip)
	   finish(ip) = bestfitparm(ip) + sig(ip)
	   step(ip) = sig(ip)/7.0
	   parm(ip)=bestfitparm(ip)
	 enddo
	 
	 open(unit=510,file="bounds.txt")
	 write(510,*) (start(ip),',',finish(ip),ip=1,6)
	 close(510)

	 close(502)
	 pm=1
	 pn=2
	 open(507,file=chi_file,access="append")
	 write(507,*) pm,pn
	 close(507)
	 parm(pm)=start(pm)
	 parm(pn)=start(pn)      
         l=1
      end if 
      
  201 open(501,file=xfile,status="REPLACE")
      
      do while (pm<=6)
      do while (pn<=6)
        do while (parm(pm)<=finish(pm))
          do while (parm(pn)<=finish(pn))
            write(501,*) (parm(ip),ip=1,6)
            close(501)
            parm(pn)=parm(pn)+step(pn)
            go to 200
          enddo 
          write(501,*) (parm(ip),ip=1,6)
          close(501)
          parm(pn)=start(pn)
          parm(pm)=parm(pm)+step(pm)
          go to 200
        enddo 
	parm(pm)=start(pm)
        parm(pn)=bestfitparm(pn)
        pn=pn+1
        parm(pn)=start(pn)
        open(507,file=chi_file,access="append")
        write(507,*) pm,pn
        close(507)
      enddo 
        parm(pm)=bestfitparm(pm)
        pm=pm+1
	pn=pm+1
        parm(pm)=start(pm)
        open(507,file=chi_file,access="append")
        write(507,*) pm,pn
        close(507)
      enddo 
      
      !do ip=1,6
      !  do while (finish(ip) >= start(ip))
         !  write(501,*) (finish(j),j=1,6)
	 ! close(501)
	 ! finish(ip)=finish(ip)-step(ip)
	 ! go to 200
	!enddo 
	!finish(ip)=bestfitparm(ip)+sig(ip)
      !enddo
!      do while (end1 >= sta1)
!         do while (end2 >= sta2)
!	    do while (end3 >= sta3)
!	       write(501,*) end1, end2, end3
!	       close(501)
!	       end3=end3-step3
!	       go to 200
!	    end do
!	    end3=a+1*sig3
!	    end2=end2-step2
!            write(501,*) end1, end2, end3
!	    close(501)
!	    go to 200
!	 end do 
!	 end2=r+1*sig2
!	 end1=end1-step1
!	 write(501,*) end1, end2, end3
!	 close(501)
!	 go to 200
!      end do
!      print *, V+2*sig1,r+2*sig2,a+2*sig3
      end if
      
      call PetscViewerDestroy(jviewer,ierr)
      call PetscViewerDestroy(cviewer,ierr)
      call MatDestroy(taoJ,ierr)
      CHKERRQ(ierr)
      call MatDestroy(taoC,ierr)
      CHKERRQ(ierr)
      call MatDestroy(taoCinv,ierr)
      CHKERRQ(ierr)
      call MatDestroy(ident,ierr)
      CHKERRQ(ierr)
      call VecDestroy(taoX,ierr)
      CHKERRQ(ierr)
      call VecDestroy(taoF,ierr)
      CHKERRQ(ierr)
      call TaoDestroy(tao,ierr)
      CHKERRQ(ierr)
      call PetscFinalize ()
      
      stop
      end program sfresco_printj



      subroutine FormFunction(tao,x,f,dummy,ierr)
      implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petsctao.h"

      PetscErrorCode ierr
      PetscInt       dummy
      Vec            x,f
      Tao            tao
      PetscInt       nx,nf
      PetscScalar    x_v(0:1),f_v(0:1)
      PetscOffset    x_i,f_i

      call VecGetSize(x,nx,ierr)
      CHKERRQ(ierr)
      call VecGetSize(f,nf,ierr)
      CHKERRQ(ierr)
      call VecGetArray(x,x_v,x_i,ierr)
      CHKERRQ(ierr)
      call VecGetArray(f,f_v,f_i,ierr)
      CHKERRQ(ierr)
      
      call sfresco_wrapper(x_v(x_i),nx,f_v(f_i),nf)

      call VecRestoreArray(x,x_v,x_i,ierr)
      CHKERRQ(ierr)
      call VecRestoreArray(f,f_v,f_i,ierr)
      CHKERRQ(ierr)

      ierr = 0
      end subroutine FormFunction
      
      subroutine ComputeFDJacobian(tao,taoX,taoF,h,taoJ,ierr)
      implicit none
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petsctao.h"

      Tao tao
      Vec taoX,taoF
      PetscReal h
      Mat taoJ
      PetscErrorCode ierr

      integer   MAXNX,MAXNF
      parameter (MAXNF=300,MAXNX=12)
      Vec xwork,fplus,fminus
      PetscOffset  x_i,f_i
      PetscScalar x_v(0:1),f_v(0:1)
      PetscScalar oneoverh,minus1
      PetscInt m,n,i,j,dummy
      PetscInt jstart,jend
      PetscInt xstart,xend
      PetscInt fstart,fend
      PetscInt maxf
      PetscInt findices(0:MAXNF)


      call MatGetOwnershipRange(taoJ,jstart,jend,ierr)
      CHKERRQ(ierr)
      call VecGetSize(taoX,n,ierr)
      CHKERRQ(ierr)
      call VecGetSize(taoF,m,ierr)
      CHKERRQ(ierr)
      call VecDuplicate(taoX,xwork,ierr)
      CHKERRQ(ierr)
      call VecDuplicate(taoF,fplus,ierr)
      CHKERRQ(ierr)
      call VecDuplicate(taoF,fminus,ierr)
      CHKERRQ(ierr)
      call VecGetOwnershipRange(fplus,fstart,fend,ierr)
      CHKERRQ(ierr)
      call VecGetOwnershipRange(xwork,xstart,xend,ierr)
      CHKERRQ(ierr)
      call VecCopy(taoX,xwork,ierr)
      CHKERRQ(ierr)

      if ((fstart .ne. jstart) .or. (fend .ne. jend)) then
         print *,'Matrix J and Vector F incompatible'
         stop
      end if
      do j=fstart,fend-1
         findices(j) = j
      end do
      

!!      /* the jth column of the jacobian is approximated by central difference
!!     J(x)_j ~= (F(x_j+h/2) - F(x_j-h/2))/ h
!!      */

      call FormFunction(tao,xwork,taof,dummy,ierr)
      do i=0,n-1
         call VecCopy(taoX,xwork,ierr)
         CHKERRQ(ierr)
         if (i .ge. xstart .and. i.lt.xend) then
            call VecGetArray(xwork,x_v,x_i,ierr)
            CHKERRQ(ierr)
            x_v(x_i + i - xstart) = x_v(x_i + i - xstart) + h/2.0
            call VecRestoreArray(xwork,x_v,x_i,ierr)
            CHKERRQ(ierr)
         end if
         call FormFunction(tao,xwork,fplus,dummy,ierr)
         CHKERRQ(ierr)
         if (i .ge. xstart .and. i .lt. xend) then
            call VecGetArray(xwork,x_v,x_i,ierr)
            CHKERRQ(ierr)
            x_v(x_i + i - xstart) = x_v(x_i + i - xstart) - h
            call VecRestoreArray(xwork,x_v,x_i,ierr)
            CHKERRQ(ierr)
         end if
         call FormFunction(tao,xwork,fminus,dummy,ierr)
         CHKERRQ(ierr)
         minus1=-1.0
         oneoverh = 1.0/h
         call VecAXPY(fplus,minus1,fminus,ierr)
         CHKERRQ(ierr)
         call VecScale(fplus,oneoverh,ierr)
         CHKERRQ(ierr)
         
         call VecGetArray(fplus,f_v,f_i,ierr)
         CHKERRQ(ierr)
         call MatSetValues(taoJ,fend-fstart,findices,1,i,f_v(f_i),         &
     &        INSERT_VALUES,ierr)
         CHKERRQ(ierr)
         call VecRestoreArray(fplus,f_v,f_i,ierr)
         CHKERRQ(ierr)
      end do

      call MatAssemblyBegin(taoJ,MAT_FINAL_ASSEMBLY,ierr)
      CHKERRQ(ierr)
      call MatAssemblyEnd(taoJ,MAT_FINAL_ASSEMBLY,ierr)
      CHKERRQ(ierr)
      call VecDestroy(fplus,ierr)
      CHKERRQ(ierr)
      call VecDestroy(fminus,ierr)
      CHKERRQ(ierr)
      call VecDestroy(xwork,ierr)
      CHKERRQ(ierr)
      ierr = 0
      end subroutine computeFDJacobian
