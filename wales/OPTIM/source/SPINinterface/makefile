#Makefile for SPINinterface

SPINOBJS = dipole_userpot.o init_module.o main_module.o math_module.o force_module.o rdm_module.o dipolsums.o 
LIBNAME = ../libSPIN.a

.SUFFIXES:
.SUFFIXES: .o .f90 .f
	
.f90.o:
	$(FC) $(FFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
../libSPIN.a: ${SPINOBJS}
	ar sr ${LIBNAME} ${SPINOBJS}
clean:
	rm -f ${LIBNAME} ${OBJS} *.mod ifc?????? *.il *.stb *.o

dipole_userpot.o: init_module.o main_module.o force_module.o math_module.o
init_module.o: main_module.o rdm_module.o
force_module.o: main_module.o math_module.o
math_module.o: main_module.o
