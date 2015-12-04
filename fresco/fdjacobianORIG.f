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
      integer   MAXNX,MAXNF
      parameter (MAXNF=300,MAXNX=12)
      character*50 search_file
      character*50 xfile
      PetscInt  nx,nf,i
      PetscReal startx(maxnx),h,one
      PetscReal x(maxnx),f(maxnf)
      PetscReal x_v(0:1)
      PetscOffset x_i
      Vec       taoX,taoF
      Mat       taoJ,taoC,taoCinv,ident
      MatFactorInfo factorinfo
      Tao       tao 
      PetscViewer cviewer,jviewer
      external  FormFunction

      
      search_file='search-rv1.in'
      xfile = 'x.in'
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-search_file',       &
     &     search_file,isset,ierr)
      CHKERRQ(ierr)

      call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-x_in',               &
     &     xfile,isset,ierr)
      CHKERRQ(ierr)

      h=1.0e-4
      call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-fd_h',h,isset       &
     &   ,ierr)
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
      
      call PetscFinalize()

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
