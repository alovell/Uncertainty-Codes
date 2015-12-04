   ! creates a pull matrix from draws that are uniform across parameter space
   program pspulls
   implicit none
   integer nparm,i,j
   real*8, allocatable :: parm_lim(:,:)
   real*8, allocatable :: temp(:)
   real*8, allocatable :: parm(:)
   
   ! read in the parameter limits
   ! put a check on the number of parameters listed vs. # limits
   read(5,*) nparm
   allocate (parm_lim(2,nparm))
   do i=1,nparm
      read(5,*) parm_lim(1,i),parm_lim(2,i)
      !print *, parm_lim(1,i),parm_lim(2,i), i
   enddo 
   
   !print *, "out of loop"
   
   ! randomly generate numbers and map them to the parameter limits
   allocate (temp(nparm))
   allocate (parm(nparm))
   call random_seed()
   open(unit=2,file="pulls.txt",status='replace')
   do i=1,200
      call random_number(temp)
      parm = 0.d0
      do j=1,nparm
         parm(j) = (parm_lim(2,j)-parm_lim(1,j))*temp(j) + parm_lim(1,j)
      enddo
      write(2,7) parm
   enddo
   close(2)
   
7  format(12f12.6) 
   
   end program