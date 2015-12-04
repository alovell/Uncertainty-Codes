   subroutine makeinputfile(parm,nf)
   implicit none
   integer, intent(in) :: nf
   real*8, intent(in) :: parm(nf)
   character (len=100) line
   real*8 parm1,parm2,parm3,parm4,parm5,parm6
   integer i
   
   ! define parameters that don't change
   ! this was for a specific run - 12C-n @ 17.29 MeV
   !parm1 = 54.198708
   !parm2 = 5.091754
   !parm3 = 1.332515
   !parm4 = 1.274099
   !parm5 = 0.6d0
   !parm6 = 0.65d0
   
   ! make input file - replace parameters with fitted parameters
   ! will have to be read line by line, more or less
   open(unit=3,file="90Zrd.in")
   open(unit=4,file="newfit.in",access="append")
   do i=1,15
      read(3,'(A)') line
      !print *, line
      write(4,'(A)') line
   enddo 
   read(3,'(A)') line
   write(line,'(" &POT kp=1 type=1 p1=" F7.3 " p2=" F8.5 " p3=" F8.5 "  p4=0.0 p5=1.29 p6=0.58  /")') parm(1),parm(3),parm(5)
   write(4,'(A)') line
   !read(3,'(A)') line
   !write(line, '(" &pot kp=1 type=11 p2=" F7.3 "  /")') parm(5)
   !write(4,'(A)') line
   read(3,'(A)') line
   write(line,'(" &POT kp=1 type=2 p1=0 p2=0 p3=0  p4=" F7.3 " p5=" F8.5 " p6=" F8.5 "  /")') parm(2),parm(4),parm(6)
   write(4,'(A)') line
   !read(3,'(A)') line
   !write(line,'(" &pot kp=1 type=11 p2=" F7.3 "  /")') parm(6)
   !write(4,'(A)') line
   do i=18,21
      read(3,'(A)') line
      write(4,'(A)') line
   enddo 
   close(3)
   close(4)
   
   end subroutine makeinputfile