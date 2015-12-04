   ! calculates the correlation matrix for a given covariance matrix
   ! for the covariance matrix located in correlation.txt
   ! (yes, that's weird; Jason misnamed it)
   ! executable is corremat
   
   program calc_corrmat
   implicit none
   integer size,i
   real*8, allocatable :: cov_matrix(:,:)
   real*8, allocatable :: sig_matrix(:,:)
   real*8, allocatable :: corr_matrix(:,:)
   real*8, allocatable :: temp_matrix(:,:)
   
   ! find out how big the covariance matrix is
   write(*,*) "How many parameters were fit for the covariance matrix?"
   read(*,*) size
   
   ! allocate matrices
   allocate(cov_matrix(size,size))
   allocate(sig_matrix(size,size))
   allocate(corr_matrix(size,size))
   allocate(temp_matrix(size,size))
   
   ! read in covariance matrix
   open(unit=1,file="correlation.txt")
   read(1,*) ! Mat Object
   read(1,*) ! type
   read(1,*) cov_matrix
   close(1)
   !print *, cov_matrix
   
   ! set up sigma matrix
   sig_matrix = 0.d0
   do i=1,size
      sig_matrix(i,i) = 1.d0/(sqrt(cov_matrix(i,i)))
   enddo 
   !print *, sig_matrix
   
   ! multiply sig*cov*sig = corr
   temp_matrix = matmul(sig_matrix,cov_matrix)
   corr_matrix = matmul(temp_matrix,sig_matrix)
   
   print *, corr_matrix
   
   end program