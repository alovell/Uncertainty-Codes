	module parameters
	integer mxx,mxp,mxpex,maxit,maxcpl,maxnnu,
     X		maxn,maxnln,maxnlo,msp,lmax1,mpair,nkmax,maxm
        integer   mloc,maxch,maxqrn,maxqrn2,mfnl,mclist,maxpw,maxb
	integer inh,maxich,maxcch,nfdec
!					derivative parameters:
	integer maxmul
	real*8 unitmass,finec,fmscal,coulcn,amu,hbc,etacns
	real*8 pmass,mun
	integer melfil,nfus,nfus1,kpmax,buttle
	integer maxnlr,maxnlc
	end module parameters

	module io
      	integer KI,KO,KOI
	logical written(399)
	end module io

	module factorials
	integer mfact
	real*8 pi,r4pi
	real*8,save,allocatable:: fact(:) 
!  	real*8 fact(1:2000) 
	end module factorials

	module drier
	integer mach
	logical dry,psiren
	real*8 fpmax,acc8
	end module drier
	
	module kcom
    	real*8 rintp,epc
    	integer nl,maxl,mlt,nln,nlo,nlc,
     X         n,mr,nbb(3),minl,jtest,lrec66,nlt,nlm,nnu
	end module kcom
	module kcom2
      	real*8,save,allocatable::     FNL(:,:),QRLN(:,:,:)
      	complex*16,save,allocatable:: FNC(:,:),QERN(:,:,:)
      	integer,save,allocatable::    WHOL(:),WHOI(:),WHERE(:,:)
      	real*8 z,eps,hfht,rinto
      	integer maxl1,minl1,minlrq,maxlrq,ib
	end module kcom2
	
	module trace
	integer chans,listcc,treneg,cdetr,smats,smatl,xstabl,nlpl,waves
     X         ,lampl,veff,kfus,wdisk,bpm,cdcc
	end module trace

	module fresco1
        character*100 headng
      	character*1 iso,rela
	real*8 hcm,rmatch,hnl,rnl,centre,hnn,rnn,rmin,rasym,accrcy,rsp,
     X   switch,ajswtch,sinjmax,jtmin,jtmax,absend,erange,dk,hktarg,
     X   thmin,thmax,thinc,cutl,cutr,cutc,ips,jbord(8),elab(5),jleast
 	integer nearfa,jump(7,3),koords,kqmax,pp,isocen,nnn,ngail,
     X		it0,iter,iblock,pade,mtmin,nlab(4),m,md,lmax,meigs,
     X		pel,exl,lin,lab,lex,nrbases,nrbmin,pcon,mrm,plane
	logical fatal,nosol,fcwfn,pralpha,symm,locfil,mcalls,rterms,
     x		allpart,cxwf,ccbins,nn2wf,sumccbins,sumkql
	real*8 rin,ebeta(2),rmatr,smallchan,smallcoup,gailacc,weak
	integer mlm,nlcn,mint,mintm2,nj,jset,pset,icutc,itmin,nsol,
     x          nforml1,nforml2,sumform,ompform,initwf
        character*10 ppkind(0:4)
        data ppkind / 'projectile','target','ejectile','residual',
     x		      'proj+eject' /
    	character*4 coords(4)
      	data coords / 'Mads','Tran','Recl','H-J ' /
     	character*11 machin(0:10)
      	data machin / 5*'Fortran 90','No PVM','f90+PVM3s','f90+PVM3',
     x 		'iPSC 860',2*'unknown '/
 	end module fresco1


	module gails
	implicit real*8(a-h,o-z)
	real*8 ryev
     	real*8, save,allocatable:: cf(:,:,:)
	integer numax
	end module gails

	module searchpar
 	integer nvars
	parameter (mvars=50)
	integer srch_kind(mvars)
!  Kind of search parameter: 0=ignore, 1=potential, 2=afrac, 3=R-mat energy,
!                            4=R-mat width, 5=dataset normalisation
!
	character*10 srch_name(mvars)
	real*8 srch_value(mvars),srch_step(mvars),nul,
     x         srch_minvalue(mvars),srch_maxvalue(mvars)
	integer srch_kp(mvars),srch_pline(mvars),srch_col(mvars)
	integer srch_nafrac(mvars),numafrac
	character*50 srch_afrac_overlap(mvars)
	real*8 srch_jtot(mvars),penalty,fine
	integer srch_par(mvars),srch_r_ch(mvars)
	integer srch_datanorm(mvars),srch_rterm(mvars)
	logical interactive,srch_nopot(mvars),final
	integer number_calls
	end module searchpar
	
	module searchdata
	parameter(mds=40, mdl=200, maxen=mdl*2)
	integer datasets,datalen(mds),chi_corr
	real*8 datangles(mdl,mds),datavals(mdl,mds),dataerr(mdl,mds)
	real*8 theoryvals(mdl,mds),data_energies(mdl,mds),data_jtot(mds)
        real*8,save,allocatable::  theoryplot(:,:)
	logical data_lab(mds),data_matched(mds),gettheoryplot
	integer data_type(mds),bs_no(mdl,mds)
	integer data_pel(mds),data_exl(mds),data_labe(mds),
     x 	        data_lin(mds),data_lex(mds)
	integer data_ic(mds),data_ia(mds),data_thfile(mds),ndof
	integer data_idir(mds),data_idir1(mds),data_rank_k(mds)
	integer data_ch(mds),data_par(mds),data_rank_q(mds)
        real*8 data_chisq(mds),ecmrat,etarat
	real*8 energy_list(maxen)
	real*8 data_error(mdl,mds)
	real*8 wij(mdl,mdl,mds)
	integer num_energies
	end module searchdata


