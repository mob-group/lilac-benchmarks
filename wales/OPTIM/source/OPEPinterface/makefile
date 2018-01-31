#Makefile for OPEPinterface

f77files = $(wildcard *.f)
f90files = $(wildcard *.f90 *.F90)
OPEPOBJS= $(f77files:.f=.o) $(f90files:.f90=.o) $(F90files:.F90=.o)

LIBNAME = ../libOPEP.a

.SUFFIXES:
.SUFFIXES: .o .f90 .f
	
.f90.o:
	$(FC) $(FFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
../libOPEP.a: ${OPEPOBJS}
	ar sr ${LIBNAME} ${OPEPOBJS}
clean:
	rm -f ${LIBNAME} ${OBJS} *.mod ifc?????? *.il *.stb *.o

# Dependencies
#
art.o : art.f90 art_initialise.o geometric_corrections.o lanczos.o art_main_module.o random.o defs.o  restart.o
art_initialise.o : art_initialise.f90 defs.o 
art_main_module.o : art_main_module.f90 defs.o 
art_read_parameters.o : art_read_parameters.f90 art_main_module.o defs.o 
calcforce.o : RNA-2006.o other-2010-RNA.o
find_saddle.o : find_saddle.f90 geometric_corrections.o art_main_module.o defs.o random.o 
geometric_corrections.o : geometric_corrections.f90 
lanczos.o : random.o geometric_corrections.o art_main_module.o defs.o 
matrix_mod.o : defs.o 
md_main_module.o : md_main_module.f90 geometric_corrections.o defs.o 
md_utils.o : defs.o geometric_corrections.o md_main_module.o 
md_initialise.o : restart.o md_utils.o md_main_module.o readtop.o
md_statistics.o : geometric_corrections.o md_main_module.o 
write_to_file.o : geometric_corrections.o defs.o 
md.o : md_utils.o md_statistics.o md_utils.o restart.o md_initialise.o write_to_file.o geometric_corrections.o md_main_module.o random.o defs.o 
md_read_parameters.o : md_read_parameters.f90 md_main_module.o 
OPEP_userpot.o : md_main_module.o md.o calcforce.o
other-2010-RNA.o : RNA-2006.o geometric_corrections.o
restart.o : restart.f90 geometric_corrections.o md_main_module.o fileio.o
mdmin.o : mdmin.f90 geometric_corrections.o defs.o 
read_parameters.o : read_parameters.f90 md_main_module.o random.o defs.o 
RNA-2006.o : fileio.o
#readtop1.o : readtop.f
saddle_converge.o : saddle_converge.f90 geometric_corrections.o lanczos.o art_main_module.o defs.o random.o 
simulateur.o : simulateur.f90 defs.o  md_main_module.o restart.o


