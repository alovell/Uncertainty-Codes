   subroutine makeinputfile(parm,nf)
   implicit none
   integer, intent(in) :: nf
   real*8, intent(in) :: parm(nf)
   character (len=100) line
   integer i
   
   ! make input file - replace parameters with fitted parameters
   ! will have to be read line by line, more or less
   open(unit=3,file="48Cain.in")
   rewind(3)
   open(unit=4,file="newfit.in",access="append")
   do i=1,16
      read(3,'(A)') line
      !print *, line
      write(4,'(A)') line
   enddo 
   read(3,'(A)') line
   write(line,'(" &POT kp=1 type=1 p1=" F7.3 " p2=" F8.5 " p3=" F9.6 "  p4=0.0 p5=1.270 p6=0.62  /")') parm(1),parm(3),parm(5)
   write(4,'(A)') line
   read(3,'(A)') line
   !write(line, '(" &pot kp=1 type=11 p2=" F7.3 "  /")') parm(5)
   write(4,'(A)') line
   read(3,'(A)') line
   write(line,'(" &POT kp=1 type=2 p1=0 p2=0 p3=0  p4=" F7.3 " p5=" F8.5 " p6=0.195695  /")') parm(2),parm(4)
   write(4,'(A)') line
   read(3,'(A)') line
   !write(line,'(" &pot kp=1 type=11 p2=" F7.3 "  /")') parm(6)
   write(4,'(A)') line
   do i=21,24
      read(3,'(A)') line
      write(4,'(A)') line
   enddo 
   close(3)
   close(4)
   
   end subroutine makeinputfile