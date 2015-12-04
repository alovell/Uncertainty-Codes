! sfresco_init - initializes sfresco with given search file(and max dimension),
!                returns problem size and default start point

	subroutine sfresco_init(search_file,maxnx,nx,nf,startx)
	use parameters
	use factorials
	use drier
        use io
	use searchpar
	use searchdata
	use fresco1, only:lmax,jbord,jump,thmin,thmax,thinc,elab,rterms,
     x    peli=>pel,exli=>exl,labi=>lab,lini=>lin,lexi=>lex
	implicit none
	integer kp,pline,col,nafrac,channel,q,par,nparameters,ndatasets,
     x          dataset,term,type,maxrank,points,pel,exl
	integer nx,nf,maxnx
        real*8 startx(maxnx)
	real*8 jtot,energy,afrac,potential,datanorm,width
	real*8 value,step,valmin,valmax,srch_error(mvars),esave(maxen),
     x         stepi,error,plist(1)
        real*8, allocatable :: twij(:,:)
	character*50 input_file,output_file,data_file,
     x       dat_file(mds),plot_file,tag
	character*(*) search_file
	character*3 cmlab
	character*8 stype(0:4),errortype
	character*10 name
	character*100 cmd
	character*20 cfile,corfile(40)
        character psign(3)
	logical xvals,abserr,lab,undef,noerror,ranks(0:4),nopot,nodat,op
        logical loge
        data psign / '-','?','+' /
	real a,angle,delta,factor,err,val,x,xmin
	integer i0,ic,id,ia,ie,idir,inf,ipiscale,k,koe
	integer labe,lex,lin,nangl,neng,ien,ip,iscale,kind
	
	namelist /variable/ name,kind,step,valmin,valmax, nul,
     x                      kp,pline,col,potential, dataset,datanorm,
     x                      nafrac,afrac, energy,jtot,par,channel,width,
     x                      term,nopot
	namelist /data/ type,data_file,points,xmin,delta,lab,energy,
     X          idir,iscale,abserr,ic,ia,k,q,angle,jtot,par,channel,
     x          pel,exl,labe,lin,lex,value,error,cfile
	
        data stype/ ' ','in Fm**2','in barns', 'in mb', 'in micbn'/
! DATA:
!    type      = 0       angular distribution for fixed energy
!              = 1       excitation & angular cross sections 
!              = 2       excitation cross section for fixed angle
!              = 3       excitation total cross section
!			  (ic=0: ia=0 is fusion; ia=1 is reaction xs)
!              = 4       excitation phase shift for fixed partial wave
!              = 5 	 search factor for bound state
!              = 6 	 value and error for a search parameter 

!    idir      = 0       cross-section data are given in absolute units.
!              = 1       cross-section data are ratio to rutherford.
!              = 2       cross sections are given in absolute units but will
!                             converted to ratio to rutherford.
!              =-1       cross sections are given as spectroscopic factors but will
!                             converted to absolute

!    iscale    =-1       dimensionless (eg ratio to rutherford if idir>0)
!              = 0       absolute cross-section units are fermi-squared/sr.
!              = 1       absolute scale is barn/sr.    (MeV-barn for idir=-1)
!              = 2       absolute scale is mb/sr. (MeV-mb for idir=-1)
!              = 3       absolute scale is micro-b/sr. (MeV-microbarn for idir=-1)
	
    	character*70 TMP,TMPN
     	common/cfs/ TMP,TMPN  
!				Change stdout recl on some machines
!	call stdout(6,140)
	inf = 0
	nul = -124578   ! 'undefined'
	fine = 10000.    ! this will be the penalty for FAILs, eg in EIGCC
	interactive = .false.
	number_calls = 0

	written = .false.; written(3) = .true.
	!write(6,1001) 
! 1001   format('SFRESCO - FRES 2.9: Search Coupled Reaction Channels'
!     x          /)
!
!####### Read in search input
!	write(6,*) 'Please give search file name'
!	read(5,*,end=999,err=1) search_file
	!write(6,1002) search_file
! 1002   format(' SFresco search file = ',a50/)
        open(303,file=search_file,status='old')

	read(303,*) input_file,output_file,nparameters,ndatasets
        if(nparameters>mvars) then
          write(6,*) 'ONLY ROOM FOR ',mvars,' SEARCH VARIABLES'
          stop
          endif
        if(ndatasets>mds) then
          write(6,*) 'ONLY ROOM FOR ',mds,' DATASETS'
          stop
          endif
	nvars=nparameters ; datasets=ndatasets
	!write(6,1003) input_file,output_file,nparameters,ndatasets
! 1003   format(' Fresco input file = ',a50/,
!     x         '   and output file = ',a50/,
!     X         ' Search on',i3,' variables for ',i3,' datasets,')
	
	
!
!####### Read in specification of search parameters
 !      	if(nparameters>0) write(6,'(/''  Define SEARCH VARIABLES:''/)')
        i0 = ichar('0')
	undef = .false.
	rterms=.false.
	stepi = 0.01
	srch_datanorm(:) = 1.
        do ip=1,nparameters
        name='Var'//char(mod(ip/10,10)+i0)//char(mod(ip,10)+i0)
        kind=0;valmin=0;valmax=0; kp=0;pline=0;col=0;
        potential=nul; step=stepi; width=0
        nafrac=0;afrac=nul; energy=nul;jtot=0;par=0;channel=1;term=1
        dataset=1;datanorm=1.0;    nopot=.false.
        read(303,nml=variable)

        srch_kind(ip) = kind; srch_name(ip) = name
        srch_minvalue(ip)=valmin; srch_maxvalue(ip)=valmax
	if(abs(step-stepi)<1e-10.and.step>abs(width) ! make default steps small for small widths
     x             .and.abs(width)>1e-20) step=step*abs(width)
        srch_step(ip)=step
	srch_error(ip)=0  ! initially

         if(kind==1) then !   potential parameter
 	  !write(6,1010) ip,name,kp,pline,col
!1010	  format('   Variable',i3,'=',a10,' is potential KP=',i2,
!     X       ', line=',i2,' col=',i2)
          srch_value(ip) = potential; 
          srch_kp(ip) = kp; srch_pline(ip)=pline; srch_col(ip)=col

         else if(kind==2) then ! spectroscopic amplitude
          !write(6,1012) ip,name,nafrac
!1012	  format('   Variable',i3,'=',a10,' is Afrac #',i3)
          srch_value(ip) = afrac; srch_nafrac(ip) = nafrac; 

         else if(kind==3) then ! R-matrix energy
 	  !write(6,1013) ip,name,term,jtot,psign(par+2),nopot
!1013	  format('   Variable',i3,'=',a10,' is energy of R-matrix term'
!     X       ,i3,' at J/pi =',f5.1,a1,' [NoPot=',L1,']')
          srch_value(ip) = energy; srch_rterm(ip) = term
	  srch_jtot(ip) = jtot; srch_par(ip) = par
	  srch_nopot(ip) = nopot
	  rterms=.true.

         else if(kind==4) then ! R-matrix width
 	  !write(6,1014) ip,name,term,channel
!1014	  format('   Variable',i3,'=',a10,' is width of R-matrix term',
!     X       i3,' in channel',i3)
          srch_value(ip) = width; srch_rterm(ip) = term
	  srch_r_ch(ip) = channel

         else if(kind==5) then ! dataset normalisation
          !write(6,1018) ip,name,dataset
!1018	  format('   Variable',i3,'=',a10,
!     X       ' is normalisation for dataset ',i3)
          srch_value(ip) = datanorm; srch_datanorm(ip) = dataset 
        endif

         if(abs(srch_value(ip)-nul)>.001) then
          !write(6,1019) srch_value(ip),step,valmin,valmax
!1019	  format('     value ',f8.4,10x,' step ',f7.4,
!     X           ', min,max ',2f8.4/)
	 else
          !write(6,1020)  step,valmin,valmax
!1020	  format('     value from Fresco input',
!     X           ', step ',f7.4,', min,max ',2f8.4/)
	  undef = .true.
	 endif
       enddo
	  peli=0; exli=0; labi=0; lini=0; lexi=0  ! we have not read in fresco input yet!
   
!
!####### Read in experimental data sets
       	ndof = 0
	ranks(:) = .false.; maxrank=-1
	num_energies=0
        do id=1,ndatasets
	  type=0; angle=0; jtot=-1;par=0;channel=1
          data_file="="; xmin=0;delta=-1;idir=0;iscale=-1;lab=.false.
          ic=1;ia=1;k=0;q=0; abserr=.false.;  points=-1; energy=-1
	  pel=peli; exl=exli; labe=labi; lin=lini; lex=lexi; cfile="not"
	  !write(6,*) ' Read definition of data set ',id
          read(303,nml=data)
	  corfile(id)=cfile
          if(idir.eq.1) iscale=-1
          errortype='relative'; if(abserr) errortype='absolute'
          cmlab='CM '; if(lab) cmlab='LAB'
	  dat_file(id) = data_file
          !write(6,1025) id,data_file,type,stype(iscale+1),
      !x      errortype,cmlab,ic,ia,idir
!1025      format('  Read DATA SET ',i2,' from file ',a50/
!     x     '    of type ',i2,'  ',a8,' with ',a8,' errors,',
!     x     1x,a3,' for state/excit',2i3,',  idir=',i2/)
	  ranks(k) = .true.
	  maxrank = max(maxrank,k)
	  !if(type==0.and.energy<0) write(6,1029) k,q
	  !if(type==0.and.energy>0) write(6,1030) k,q,energy
	  !if(type==1) write(6,1031) k,q
	  !if(type==2) write(6,1032) k,q,angle,cmlab
	  !if(type==3) write(6,1033) 
	  !if(type==4) write(6,1034) jtot,psign(par+2),channel
	  !if(type==5) write(6,1035) 
	  !if(type==6) write(6,1036) par,value,error,abserr
!1029	  format('  Angular cross section T',2i1)
!1030	  format('  Angular cross section T',2i1,
!     x ' for energy',f8.3,' MeV')
!1031	  format('  Excitation and Angular cross sections T',2i1)
!1032	  format('  Excitation cross section T',2i1,
!     x ' for angle',f8.3,' deg ',a3)
!1033	  format('  Excitation total cross section')
!1034	  format('  Excitation phase shift in channel',f5.1,a1,' #',i2)
!1035	  format('  Target search parameters for bound states')
!1036	  format('  Value and error for search parameter',i3,' are'/
!     x        2f10.5,' (abserr=',l1,')')
	  !if(type<0.or.type>6) write(6,*) 'Unrecognised data type ',type
	  data_type(id) = type
	if(type==6) then
	      datavals(1,id)=value
	      if(abserr) error = error*value
 	      dataerr(1,id)=error
 		datalen(id) = 1
 	      ip = 2
	else
              factor=1.0
	      if(iscale<0) factor=1. 
	      if(iscale==0) factor=10.
	      if(iscale>0) factor=1000.0/10.0**(3*(iscale-1)) 
	  inf=306
	  if(data_file=="=") inf=303
	  if(data_file=="<") inf=5
          if(inf==306) then
             !write(6,*) ' To open file for dataset ',id,': ',data_file
             open(inf,file=data_file,status='old')
             endif
          xvals = delta<=0.; x = xmin
          if(points<0) points=99999
	  !write(6,*) ' Read data set ',id
          do 10 ip=1,points
	    if(type==1) then
              read(inf,end=111,fmt=*,err=11) x,a,val,err
            else if(xvals.or.type==5) then
              read(inf,end=111,fmt=*,err=11) x,val,err
            else
              read(inf,end=111,fmt=*,err=11) val,err
              x = x + delta            
            endif
              if(x<0) go to 11
	      val =factor*val
	      if(abserr) then
	         err =factor*err
	        else
	         err = val * err
	        endif  ! Now all errors are absolute (and mb, except for r/ruth)
	      datavals(ip,id)=val
 	      dataerr(ip,id)=abs(err)

            if(type==0) then        ! angular distribution for fixed energy
              datangles(ip,id) = x
	      data_energies(ip,id) = energy
             else if(type==1) then  ! excitation and angular cross sections
 	      data_energies(ip,id) = x
              datangles(ip,id) = a
             else if(type==2) then  ! excitation cross section for fixed angle
 	      data_energies(ip,id) = x
	      datangles(ip,id) = angle
             else if(type==3) then  ! excitation total cross section 
 	      data_energies(ip,id) = x
             else if(type==4) then  ! excitation phase shift
 	      data_energies(ip,id) = x
             else if(type==5) then  ! bound state search
 	      bs_no(ip,id) = nint(x)
	     endif
  10      continue
          ip = points+1
	  go to 111
  11      backspace inf
	endif
 111      datalen(id) = ip-1
          if(inf==306) then
!             write(6,*) ' Close file for dataset ',id,': ',data_file
             close(inf)
             endif
	  !close(303)
  	  ndof = ndof + datalen(id)
          data_idir(id)=idir; data_idir1(id)=idir; data_lab(id)=lab; 
          data_ic(id) = ic; data_ia(id) = ia; 
          data_rank_k(id) = k; data_rank_q(id) = q; 
	  data_ch(id) = channel; data_jtot(id)=jtot; data_par(id)=par
	  if(data_type(id)/=6) then
          !if(energy<0) write(6,*) ' ',datalen(id),' data points:'
          !if(energy>0) write(6,*) ' ',datalen(id),' data points',
      !x              ' for lab energy ',real(energy)
     	  endif
	  if(datalen(id)==0) then
	    write(6,*) '   NO DATA POINTS !! Stop'
	    stop
	    endif
            if(pel.le.0) pel = 1
            if(exl.le.0) exl = 1
            if(labe.eq.0) labe = pel
            if(lin.eq.0) lin = 1
            if(lex.le.0) lex = 1
	    data_pel(id) = pel
	    data_exl(id) = exl
	    data_labe(id) = labe
	    data_lin(id) = lin
	    data_lex(id) = lex
            if(pel+exl+labe+lin+lex>5) then
	      !write(6,1250) pel,exl,labe,lin,lex
! 1250       format('     Incoming partition',I3,' in state #',I2,
!     X             ',   Lab Energy for part.',I3,' Nucleus',I2,
!     X             ' in Excitation pair',I2/)
              endif
	 if(type==0) then
          !write(6,*) ' Angle '//cmlab//' Datum     Absolute error'
	  do ip=1,datalen(id)
	  !write(6,12) datangles(ip,id),datavals(ip,id),dataerr(ip,id)
 ! 12 	   format(1x,f8.3,2g12.4)
	  enddo
	 else if(type==1) then
          !write(6,*) '   Energy Angle '//cmlab//
      !x               ' Datum     Absolute error'
	  do ip=1,datalen(id)
	  !write(6,13) data_energies(ip,id),datangles(ip,id),
      !x                datavals(ip,id),dataerr(ip,id)
 ! 13 	   format(1x,2f8.3,2g12.4)
	  enddo
	 else if(type==5) then
          !write(6,*) '   Bound state   Target     Absolute error'
	  do ip=1,datalen(id)
	  !write(6,131) bs_no(ip,id),datavals(ip,id),dataerr(ip,id)
 ! 131 	   format(1x,i8,3x,2f12.4)
	  enddo
	 else if(type==6) then
          !write(6,*) '   Search param  Target     Absolute error'
	  do ip=1,datalen(id)
	  !write(6,132) data_par(id),srch_name(par),
      !x                 datavals(ip,id),dataerr(ip,id)
!  132 	   format(1x,i2,':',a8,2f12.4)
	  enddo
	 else
          !write(6,*) '  Energy   Datum     Absolute error'
	  do ip=1,datalen(id)
	  !write(6,12) data_energies(ip,id),datavals(ip,id),
      !x                dataerr(ip,id)
	  enddo
	 endif
          close(1)
	  neng=1
	if(type<5) then
	  if(energy<0) neng=0
	  if(type>0) neng=datalen(id)
	  do ip=1,neng
	    if(type>0) energy=data_energies(ip,id)
	  ien=0
	  do ie=1,num_energies
	   if(abs(energy-energy_list(ie))<1e-5) then
 	     ien=ie; goto 14 ! found existing energy
	   endif
	  enddo
	  ien = num_energies+1 ! list new energy
	  num_energies = ien
	  if(ien>maxen) then
	    write(6,*) 'Should increase PARAMATER maxen!'
	    stop
	    endif
	  energy_list(ien) = energy
   14     continue
   	  enddo  
	endif

        enddo     
        close(2) 
	close(303) 
        !write(6,*) 
	
	! calculate the error matrix for each data set
	! W = (C + E)^(-1)
	! C is read in (nml=data), E has exp. errors squared
	! read in C for each data set
	!d = 1.d0
	do id=1,datasets
	   if (corfile(id) .ne. "not") then
	      ! read in C from corfile
	      open(unit=700,file=corfile(id))
	      allocate (twij(datalen(id),datalen(id)))
	      read(700,*) twij
	      close(700)
	   else
	      ! if no corfile, C=0
	      allocate (twij(datalen(id),datalen(id)))
	      twij = 0.d0
	   end if 
	   ! add the experimental errors in quadrature
	   do ip=1,datalen(id)
	      twij(ip,ip) = twij(ip,ip) + dataerr(ip,id)**2
	   enddo
	   ! invert the matrix W
	   !print *, twij
	   call matinv(twij,datalen(id),d)
	   !print *, twij
	   ! put twij into wij for easy transfer
	   do ip=1,datalen(id)
	      do iq=1,datalen(id)
	         wij(ip,iq,id) = twij(ip,iq)
	      enddo 
	   enddo
	   !wij(:,:,id) = twij
	   deallocate (twij)
	   !print *, wij(:,:,id)
	enddo 

!			Pre-read input to find SOME array-parameter limits
!
!####### Read in main fresco input
	ki = 306
	ko = 3
	koe = 307
!	call machine(mach)
	mach = 1
        open(ki,file=input_file,status='old')
	open(ko,form='formatted',delim='apostrophe')
	open(koe,file=trim(output_file)//'-init')
	!write(6,*) ' FRESCO output to ',trim(output_file)//'*'
        !write(6,*) 

	call freadf(ki,ko,koe,TMP,lmax,jbord,jump)
	close(ki); close(koe)
        NANGL = INT((abs(THMAX) - THMIN)/abs(THINC) + 1 + 0.5)
        allocate (theoryplot(max(mdl,NANGL),datasets))
	gettheoryplot = .false.
	  if(num_energies==0.and.elab(1)>0.) then
	    num_energies = 1
	    energy_list(1) = elab(1)
	    endif
	  if(num_energies>0) then
	     !write(6,15) (energy_list(ien),ien=1,num_energies)
!   15	  format(' Calculate scattering at energies'/(1x,10f8.3))
             else
	     write(6,*) ' No scattering energies'
	     endif
        !write(6,*) 
	
	ki = 3
	ko = 308
	open(ko,file=output_file)
	koi= ko
	written(ko) = .true.
!	rewind ki
	noerror = .true.
	final = .false.
	final = .true.
	MAXCH = 0; MAXICH = 0 !	no arrays allocated
       	
	nf = 0
	do id=1,ndatasets
	   nf = nf + datalen(id)
	enddo
	nx = nparameters
	do ip = 1,nx
	   startx(ip) = srch_value(ip)
	enddo
	if (allocated (theoryplot)) deallocate (theoryplot)
	end subroutine sfresco_init

