	subroutine sfresco_wrapper(x,nx,f,nf)
	use parameters
	use factorials
	use drier
	use io
	use searchpar
	use searchdata
	implicit none
	integer nx,nf
	real*8 x(nx),f(nf)
	integer counter,id,il,ip

	do ip=1,nx
	   srch_value(ip) = x(ip)
	enddo

	call fr

	counter=1
	do id=1,datasets
	   do il=1,datalen(id)
	      f(counter) = data_error(il,id)
!	      print *, f(counter),data_error(il,id)
	      counter = counter+1
	   enddo
	enddo



	end subroutine sfresco_wrapper
