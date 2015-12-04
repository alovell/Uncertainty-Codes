#	Fresco: for Fortran
#
#	Ian J. Thompson,     I.Thompson@surrey.ac.uk
#	Aug 2006

FFLAGS = -g  -cpp -I$(PETSC_DIR)/include -I${PETSC_DIR$}/${PETSC_ARCH}/include  -O 		# generic
#FFLAGS = -O2 -Wall -ftrace=full	# generic
#FFLAGS = -O2 -Wtabs	# gfortran 
#FFLAGS = -O2  -Wall	# generic
#FFLAGS = -O2 -m32	# 32-bit machines generic
#FFLAGS = -O2 -m64	# 64-bit machines generic

#FFLAGS = -O2 -dcfuns -w				# for nag f95
#FFLAGS = -C=all -C=undefined -dcfuns -w -O0 -gline     # for nag f95 debug
#FFLAGS = -O2 -qmaxmem=16384                            #   optimisation for xlf
#FFLAGS = -Wall -g -fbounds-check -ftrace=full          # debug for g95
#FFLAGS = -g  -C  	 				# generic debug for f90
FC=ifort
F90C=$(FC)

include ${PETSC_DIR}/conf/variables


#               Select subroutine for CPU time:
#TIME = etime.o		# Comment out for g95 versions >= 0.91
#TIME = system_clock.o
TIME = cpu_time.o
#               Select subroutine for flush file buffers to disk
#FLUSH = flush.o
#FLUSH = 
VPATH=../bin
COMMON= globx7.o frxx0.o frxx1.o frxx2.o frxx3.o frxx4.o frxx5.o frxx6.o frxx7a.o \
        frxx7b.o frxx8.o frxx9.o frxx10.o frxx11.o frxx12.o \
        frxx13.o frxx16.o frxx17.o frxx18.o cdc.o   $(TIME) $(FLUSH)
	
CCOMMON= globx7.o frxx0.o cfrxx1.o frxx2.o frxx3.o cfrxx4.o cfrxx5.o frxx6.o frxx7a.o \
        frxx7b.o frxx8.o frxx9.o frxx10.o frxx11.o frxx12.o \
        frxx13.o frxx16.o frxx17.o frxx18.o cdc.o flib.o  $(TIME) $(FLUSH)

OBJ = $(COMMON) fresco.o 

target: fresco sfresco csfresco fdjacobian makegrid angcorr
install: fresco sfresco csfresco fdjacobian makegrid angcorr
#	cp -fp fresco sfresco $(HOME)/fresco/$(ARCH)

fresco: $(OBJ)
	$(FC) $(FFLAGS) -o ../bin/fresco $(OBJ)

OBJS = $(COMMON) sfresco.o minuit-cern.o
COBJS = $(CCOMMON) csfresco.o minuit-cern.o
#sfresco.o: sfresco.f
#	$(FC) $(FFLAGS)  -o sfresco sfresco.f
sfresco: $(OBJS)
	-@mkdir -p ../bin
	$(FLINKER)  -o ../bin/sfresco $(OBJS) $(PETSC_TAO_LIB)
	
csfresco: $(COBJS)
	-@mkdir -p ../bin
	$(FLINKER)  -o ../bin/csfresco $(COBJS) $(PETSC_TAO_LIB)

sfresco_tao: $(COMMON) sfresco_tao.o sfresco_wrapper.o sfresco_init.o
	-@mkdir -p ../bin
	${FLINKER} -o ../bin/sfresco_tao $(COMMON) sfresco_wrapper.o sfresco_init.o sfresco_tao.o $(PETSC_TAO_LIB)

fdjacobian: $(CCOMMON) fdjacobian.o sfresco_wrapper.o sfresco_init.o
	-@mkdir -p ../bin
	${FLINKER} -o ../bin/fdjacobian $(CCOMMON) sfresco_wrapper.o sfresco_init.o fdjacobian.o $(PETSC_TAO_LIB)

makegrid: $(COMMON) grid.o chisq_calc.o res.o sfresco_wrapper.o sfresco_init.o
	 -@mkdir -p ../bin
	 ${FLINKER} -o ../bin/makegrid $(COMMON) sfresco_init.o sfresco_wrapper.o res.o chisq_calc.o grid.o $(PETSC_TAO_LIB)
	 
angcorr: $(COMMON) anglecorr.o res.o sfresco_wrapper.o sfresco_init.o
	 -@mkdir -p ../bin
	 ${FLINKER} -o ../bin/angcorr $(COMMON) sfresco_init.o sfresco_wrapper.o res.o anglecorr.o $(PETSC_TAO_LIB)

sfresco_planes: $(COMMON) sfresco_planes.o sfresco_wrapper.o sfresco_init.o
	-@mkdir -p ../bin
	${FLINKER} -o ../bin/sfresco_planes $(COMMON) sfresco_wrapper.o sfresco_init.o sfresco_planes.o $(PETSC_TAO_LIB)

sfresco_eval: $(COMMON) sfresco_eval.o sfresco_wrapper.o sfresco_init.o
	-@mkdir -p ../bin
	${FLINKER} -o ../bin/sfresco_eval $(COMMON) sfresco_wrapper.o sfresco_init.o sfresco_eval.o $(PETSC_TAO_LIB)

rv1_R: $(COMMON) rv1_R.o sfresco_wrapper.o sfresco_init.o
	-@mkdir -p ../bin
	${FLINKER} -o ../bin/rv1_R $(COMMON) sfresco_wrapper.o sfresco_init.o rv1_R.o $(PETSC_TAO_LIB)

frxx1.o: globx7.o frxx1.f usescatwf.f veffpot.f bpmfus.f
frxx17.o: globx7.o frxx17.f usescatwf.f

.f.o:
	$(FC) $(FORM) -c $(FFLAGS)  $<


#########
clean:
	\rm -f *.o core  *.mod
