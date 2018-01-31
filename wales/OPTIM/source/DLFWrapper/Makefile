# This Makefile was adapted for use in OPTIM by 
# Judith Rommel (jbr36) and Joanne Carr in 2013

# COPYRIGHT
#
#  Copyright 2007 Johannes Kaestner (kaestner@theochem.uni-stuttgart.de),
#  Tom Keal (thomas.keal@stfc.ac.uk)
#
#  This file is part of DL-FIND.
#
#  DL-FIND is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as 
#  published by the Free Software Foundation, either version 3 of the 
#  License, or (at your option) any later version.
#
#  DL-FIND is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public 
#  License along with DL-FIND.  If not, see 
#  <http://www.gnu.org/licenses/>.

#for debugging the makefile
#FC=ifort

# pgf90, g95, ifort, gfortran
arch =$(FC)

PPFLAGS =  
F90EXT = f90
ALLOC = new
OBJDIR = $(arch)/objects
#PARA = 0
#$(MAKECMDGOALS)
VPATH = ../:ifort/objects/
LIBNAME = ../libdlfwrapper.a


########################################################################
##  Intel compiler
ifeq ($(arch),ifort)
  F90        = $(FC)
#  F90        = ifort 
  F90FLAGS   =  $(FFLAGS) -module $(OBJDIR)
#-C -ftrapuv -traceback -module $(OBJDIR)
  LDFLAGS    = -lmkl -lguide -lmkl_lapack -lpthread
  F90EXT     = F90
endif
########################################################################



########################################################################
##  Portland compiler (pgf)
ifeq ($(arch),pgf90)
  F90        = $(FC)
  MPI_F90    = mpif90
  F90FLAGS   =  $(FFLAGS) -module $(OBJDIR)
  LDFLAGS    = -llapack -lblas -lg2c
  MPI_LDFLAGS= -lacml -lmpi
  #PPFLAGS    = -DOLDALLOC
  F90EXT     = F90
  #ALLOC      = old
endif

########################################################################
##  g95 compiler
ifeq ($(arch),g95)
  F90        = $(FC)
  F90FLAGS   =  $(FFLAGS) -fmod=$(OBJDIR)
  #F90FLAGS   =  -C -ftrace=full -freal=nan -finteger=nan -flogical=none \
                -fpointer=none -fimplicit-none  -fbounds-check \
                -fintrinsic-extensions -std=f95 -ftr15581 -Wuninitialized \
                -Wprecision-loss -fmod=$(OBJDIR) 
  LDFLAGS    = 	-lm -llapack -lblas
  F90EXT     = F90
# the next line may be used to find memory leaks
#  ALLOC      = old
endif
########################################################################

########################################################################
##  gfortran (gcc) compiler
# at the moment, gfortran leads to numerically different results than 
# the others ...
ifeq ($(arch),gfortran)
  F90        = $(FC)
  F90FLAGS   =  $(FFLAGS) -M $(OBJDIR)
  LDFLAGS    = -lm -llapack -lblas
 #-lcblas -lf77blas -latlas -lg2c
  F90EXT     = F90
endif
########################################################################



.SUFFIXES:
.SUFFIXES: .f90 .o .mod

# modules used by more than one file
OBJ_MOD = $(OBJDIR)/dlf_stat_module.o $(OBJDIR)/dlf_global_module.o $(OBJDIR)/dlf_allocate.o $(OBJDIR)/dlf_util.o \
	$(OBJDIR)/dlf_sort.o

# Objects that work without numerical libraries (this distinction does not hold any more) 
OBJ_NOLAPACK = $(OBJDIR)/interfacetoOPTIM.o $(OBJDIR)/dl-find.o $(OBJDIR)/dlf_convergence.o $(OBJDIR)/dlf_checkpoint.o \
	$(OBJDIR)/dlf_svnversion.o  $(OBJDIR)/dlf_formstep.o $(OBJDIR)/dlf_scalestep.o  \
	$(OBJDIR)/dlf_lbfgs.o $(OBJDIR)/dlf_time.o \
	$(OBJDIR)/dlf_task.o $(OBJDIR)/dlf_neb.o $(OBJDIR)/dlf_qts.o $(OBJDIR)/dlf_dimer.o \
	$(OBJDIR)/dlf_conint.o $(OBJDIR)/dlf_coords.o

# Objects that need LAPACK linked in
OBJ_LAPACK = $(OBJDIR)/dlf_hdlc_matrixlib.o $(OBJDIR)/dlf_hdlc_primitive.o $(OBJDIR)/dlf_hdlc_constraint.o \
         $(OBJDIR)/dlf_hdlc_hdlclib.o $(OBJDIR)/dlf_hdlc_interface.o $(OBJDIR)/dlf_linalg.o

OBJ =  $(OBJ_NOLAPACK) $(OBJ_LAPACK) 

#find:	$(OBJ_MOD) $(OBJ) $(OBJDIR)/main.o $(OBJDIR)/dlf_serial.o
#	@echo "Linking ..."
#	$(F90) $(F90FLAGS) -o find.x $(OBJDIR)/main.o $(OBJ) $(OBJDIR)/dlf_serial.o $(OBJ_MOD)  $(LDFLAGS)
#	@cp find.x $(arch)/find_$(arch).x

obj:	$(OBJ_MOD) $(OBJ) $(OBJDIR)/dlf_serial.o

#all:	$(LIBNAME)

#$(LIBNAME):	${OBJ}
#	ar sr ${LIBNAME} ${OBJ}


#	echo 'subroutine dlf_svnversion(iver)' > dlf_svnversion.f90
#	echo '  character(*) :: iver' >> dlf_svnversion.f90
#	echo -n '  iver="' >> dlf_svnversion.f90
dlf_svnversion.f90: FORCE 
	cp dlf_svnversion.head dlf_svnversion.f90
	-svnversion -n . >> dlf_svnversion.f90
	echo '"' >> dlf_svnversion.f90
	echo "end subroutine dlf_svnversion" >> dlf_svnversion.f90

FORCE:

#
# The main rule for compiling f90 files:
#
$(OBJDIR)/%.o: %.f90
	@mkdir -p $(OBJDIR)
ifeq ($(ALLOC),old)
	sed -e"s/call allocate *(/ aLLocate(/g" \
	-e"/ aLLocate(/s/,/(/" \
	-e"/ aLLocate(/s/\(.*\))/\1))/" \
	-e"s/call deallocate *(/deallocate(/" $*.f90 > $(OBJDIR)/$*_d.$(F90EXT)
else
	cp $*.f90 $(OBJDIR)/$*_d.$(F90EXT)
endif

	$(F90) $(F90FLAGS) $(PPFLAGS) -c $(OBJDIR)/$*_d.$(F90EXT) -o $(OBJDIR)/$*.o
	@#rm $*_d.$(F90EXT)
	@echo

#
# remove $(OBJDIR)/dlf_mpi_d.$(F90EXT) if compiling the serial version, 
# and $(OBJDIR)/dlf_serial_d.$(F90EXT) if the parallel one.
#
# ifeq ($(findstring parallel,$(PARA)),parallel)
#	@rm -f $(OBJDIR)/dlf_serial_d.$(F90EXT)
#	$(MPI_F90) $(F90FLAGS) $(PPFLAGS) -c $(OBJDIR)/$*_d.$(F90EXT) -o $(OBJDIR)/$*.o
#	@#rm $*_d.$(F90EXT)
#	@echo 
#else
#	@rm -f $(OBJDIR)/dlf_mpi_d.$(F90EXT)
#endif


#
# Utilities
#
clean:
	rm -f $(OBJ_MOD) $(OBJ) $(OBJDIR)/*.mod $(OBJDIR)/interfacetoOPTIM.o $(OBJDIR)/*_d.?9? $(OBJDIR)/dlf_serial.o 

# remove all module directories and executables
clean_all:
	rm -r g95 pgf90 ifort gfortran

man:
	latex documentation.tex
	dvips documentation.dvi


#dependencies from g95 -M:
$(OBJDIR)/dlf_stat_module.o $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_parameter_module.mod: dlf_stat_module.f90
$(OBJDIR)/dlf_global_module.o $(OBJDIR)/dlf_global.mod: $(OBJDIR)/dlf_parameter_module.mod dlf_global_module.f90
$(OBJDIR)/dlf_allocate.o $(OBJDIR)/dlf_allocate.mod: $(OBJDIR)/dlf_parameter_module.mod dlf_allocate.f90
$(OBJDIR)/dlf_util.o $(OBJDIR)/dlf_store.mod $(OBJDIR)/bspline.mod: $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_parameter_module.mod dlf_util.f90
$(OBJDIR)/dl-find.o: $(OBJDIR)/dlf_store.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dl-find.f90
$(OBJDIR)/dlf_convergence.o $(OBJDIR)/dlf_convergence.mod: $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_convergence.f90
$(OBJDIR)/dlf_checkpoint.o $(OBJDIR)/dlf_checkpoint.mod: $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_parameter_module.mod $(OBJDIR)/dlf_global.mod dlf_checkpoint.f90
$(OBJDIR)/dlf_svnversion.o: dlf_svnversion.f90
$(OBJDIR)/dlf_formstep.o $(OBJDIR)/dlf_hessian.mod $(OBJDIR)/dlf_formstep_module.mod: $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_formstep.f90
$(OBJDIR)/dlf_scalestep.o $(OBJDIR)/dlf_linesearch.mod $(OBJDIR)/dlf_scalestep_module.mod: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_scalestep.f90
$(OBJDIR)/dlf_lbfgs.o $(OBJDIR)/newlbfgs_module.mod: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_lbfgs.f90
$(OBJDIR)/dlf_time.o $(OBJDIR)/dlf_time.mod: $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_time.f90
$(OBJDIR)/dlf_task.o $(OBJDIR)/dlf_task_module.mod: $(OBJDIR)/dlf_store.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_task.f90
$(OBJDIR)/dlf_neb.o $(OBJDIR)/dlf_neb.mod: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_stat.mod $(OBJDIR)/bspline.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_neb.f90
$(OBJDIR)/dlf_qts.o $(OBJDIR)/dlf_qts.mod: $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_neb.mod $(OBJDIR)/dlf_parameter_module.mod dlf_qts.f90
$(OBJDIR)/dlf_dimer.o $(OBJDIR)/dlf_dimer.mod: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_dimer.f90
$(OBJDIR)/dlf_conint.o $(OBJDIR)/dlf_conint.mod: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_hessian.mod $(OBJDIR)/dlf_stat.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_conint.f90
$(OBJDIR)/dlf_coords.o: $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_coords.f90
$(OBJDIR)/dlf_hdlc_matrixlib.o $(OBJDIR)/dlfhdlc_matrixlib.mod: $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_hdlc_matrixlib.f90
$(OBJDIR)/dlf_hdlc_primitive.o $(OBJDIR)/dlfhdlc_primitive.mod: $(OBJDIR)/dlfhdlc_matrixlib.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_hdlc_primitive.f90
$(OBJDIR)/dlf_hdlc_constraint.o $(OBJDIR)/dlfhdlc_constraint.mod: $(OBJDIR)/dlfhdlc_matrixlib.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_hdlc_constraint.f90
$(OBJDIR)/dlf_hdlc_hdlclib.o $(OBJDIR)/dlfhdlc_hdlclib.mod: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlfhdlc_constraint.mod $(OBJDIR)/dlfhdlc_primitive.mod $(OBJDIR)/dlfhdlc_matrixlib.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_hdlc_hdlclib.f90
$(OBJDIR)/dlf_hdlc_interface.o: $(OBJDIR)/dlf_checkpoint.mod $(OBJDIR)/dlf_allocate.mod $(OBJDIR)/dlfhdlc_hdlclib.mod $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_hdlc_interface.f90
$(OBJDIR)/dlf_linalg.o: $(OBJDIR)/dlfhdlc_matrixlib.mod $(OBJDIR)/dlf_parameter_module.mod dlf_linalg.f90
$(OBJDIR)/dlf_sort.o $(OBJDIR)/dlf_sort_module.mod: $(OBJDIR)/dlf_parameter_module.mod dlf_sort.f90
$(OBJDIR)/dlf_serial.o: $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod dlf_serial.f90
$(OBJDIR)/interfacetoOPTIM.o $(OBJDIR)/driver_module.mod: $(OBJDIR)/dlf_global.mod $(OBJDIR)/dlf_parameter_module.mod \
$(OBJDIR)/dlf_allocate.mod key.mod key.o commons.mod commons.o interfacetoOPTIM.f90

# They might look like:
#dlf_linalg.o: dlfhdlc_matrixlib.mod dlf_parameter_module.mod dlf_linalg.f90
#main.o driver_module.mod: dlf_parameter_module.mod main.f90
#$(OBJDIR)/main.o $(OBJDIR)/driver_module.mod: $(OBJDIR)/dlf_parameter_module.mod main.f90

