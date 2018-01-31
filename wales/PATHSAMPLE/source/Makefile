PROG =	../bin/PATHSAMPLE

OBJS =	commons.o nodes.o porfuncs.o utils.o dock.o main.o keywords.o setup.o input.o \
	KMCcommit.o mysystem.o sort.o KMC.o inertia.o eig.o amhdump.o charmmdump.o \
	cycle.o mywait.o submitoptimjob.o connectodata.o mindist.o unresdump.o \
	centre.o getnewpath.o addperm.o tssearch.o \
	mindouble.o tsdouble.o Pfold.o NGT.o GT.o getallpaths.o cycle2.o dprand.o sdprnd.o \
        Dijkstra.o getdpair.o Dijinit.o calcorder.o getrpair.o connectd.o getspair.o \
        getpair.o pairdouble.o donedouble.o newmindist.o mergedb.o getupair.o \
	GT2.o GT2data.o GT2FibonacciHeap.o GT2DLL.o GT2input.o GT2FreeMemory.o \
        minperm.o minpermdist.o rigidbodymod.o mathsconstants.o quaternionmatch.o \
        regroupfree.o regroup.o dsort.o savestate.o orderodata.o regroupfree2.o \
	probacc.o newconn.o getfreepair.o getfreebarrier.o kshortestpaths.o Dijinitfly.o \
	getallmin.o myorient.o getusepair.o NGTmem.o NGTrealloc.o NGTrenorm.o \
	NGTremoveid.o NGTremovei.o swapnode.o mymerge.o rbinertia.o rigidb.o diagonalise2.o reweight.o \
	rbperm.o minpermdistrbcom.o Cv.o DOS.o bulkmindist.o frictionfac.o make_conpot.o setup_sis.o \
        rateconst_setup.o checkspodata.o addminxyz.o printsummary.o \
	checkTS.o remove_unconnected.o getupairmetric.o \
	NGT_crstorage.o NGTremovei_crstorage.o NGTremoveid_crstorage.o NGTrenorm_crstorage.o \
        NGTrealloc_crstorage.o lopermdist.o getmetric.o rotations.o vec3.o newconnodata.o Gthomson.o \
	microtherm.o getbarrier2.o persistence.o regrouppersist.o shannon.o get_exact_persistence.o


# note that nag64/5.1-216 fails for large memory, but version 365 works
# large memory for pgi and ifort requires the -mcmodel flag!

# WARNING - points.min and point.ts created with ifort executables cannot be
# read by NAG or PGI executables, and vice versa, unless -assume byterecl is used

# Preprocessing
   DEFS =
   CPP = /lib/cpp
   CPFLAGS = -traditional -P
  
# NAG compiler 
# FC = nagfor
#  this line is for garden variety debugging 
# FFLAGS = -132 -maxcontin=3000 -C -g -kind=byte -mismatch_all -ieee=stop
#  this line is for thorough but slow debugging 
# FFLAGS = -132 -maxcontin=3000 -C=all -mtrace=all -gline -kind=byte -ieee=stop
# FFLAGS = -132 -maxcontin=3000 -mismatch_all -kind=byte -O0
# FFLAGS = -132 -maxcontin=3000 -mismatch_all -kind=byte -O4 -ieee=full 
# NOOPT= -132 -maxcontin=3000 -kind=byte -O0 -ieee=full
# NOOPT= -132 -maxcontin=3000 -kind=byte -O0 -ieee=stop
# NAGSWITCH = nag
# LIBS = libmyblas.a libmylapack.a libmyblas.a 
# LDFLAGS = -LBLAS -LLAPACK
# 
# pathscale compiler 
#  FC = pathf95
#  FFLAGS = -extend-source  -g -C
#  FFLAGS =  -extend-source -O3 
#  NOOPT= -extend-source -O0
#  NAGSWITCH = pathscale
#  LIBS = libmyblas.a libmylapack.a libmyblas.a 
#  LDFLAGS = -LBLAS -LLAPACK
# 
### Intel  
### Unformatted points files generated with ifort are incompatible with NAG and PGI
### without the -assume byterecl flag for ifort.
#
#  Without -heap-arrays ifort executables now generate a SEGV in Dijkstra.f90
#
FC = ifort
#  FFLAGS= -132 -C -g -heap-arrays -assume byterecl
# FFLAGS= -132 -g -debug all -check all -implicitnone -warn unused -fp-stack-check -heap-arrays -ftrapuv -check pointers -check bounds  \
# -assume byterecl
#
# FFLAGS= -132  -heap-arrays  -O0 -g -C -traceback -debug full -debug all -check uninit -stand f03 -assume realloc_lhs \
# -check all,noarg_temp_created -warn all -fstack-protector -assume protect_parens -implicitnone -fp-stack-check \
# -ftrapuv -check pointers -check bounds  -assume byterecl
FFLAGS= -132 -heap-arrays -O4 -assume byterecl
# FFLAGS= -132 -O0 -heap-arrays -assume byterecl
NOOPT= -O0  -assume byterecl # -mcmodel=large
NAGSWITCH=ifort
SEARCH_PATH =  -I..
LIBS = libmyblas.a libmylapack.a libmyblas.a 
LDFLAGS = -LBLAS -LLAPACK
# 
##### The Portland Group Compiler Technology Fortran 90 compiler 
#
# FC = pgf90
#
# FFLAGS= -Mextend -O3 -Munroll -Mnoframe -fast -Mipa=fast,inline
# FFLAGS= -Mextend -O3 -Munroll -Mnoframe -Mcache_align -Mflushz -mcmodel=medium # -C -g
# works for addperm3:
# FFLAGS= -Mextend -O1 -mcmodel=medium -Mnoframe -Mcache_align -Mflushz
# NOOPT = -O0 -Mextend 
# NAGSWITCH=pgi
# LIBS = -lblas -llapack -lblas
# LDFLAGS = -LBLAS -LLAPACK
#
##### end of The Portland Group Compiler Technology Fortran 90 compiler 
### Gfortran  
#
# FC = gfortran
# FFLAGS= -ffixed-line-length-132 -O0 
# FFLAGS= -ffixed-line-length-132 -O3
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv 
# FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
# NOOPT= -O0 -ffixed-line-length-132
# NAGSWITCH=gfortran
# SEARCH_PATH =  -I..
# LIBS = libmyblas.a libmylapack.a libmyblas.a 
# LDFLAGS = -LBLAS -LLAPACK
#
# 
###############################################################################################
# Pass the subversion revision number into the code to print it in the output
   DEFS+=-DSVNVERSION="'`./version.sh`'"
####################################### RULES AND TARGETS ###################################### 

# END OF COMPILER SPECIFIC STUFF

.SUFFIXES:
.SUFFIXES: .o .f .F .f90 .c

.f90.o:
	$(FC) $(FFLAGS) -c $<
.f.o:
	$(FC) $(FFLAGS) -c $<
.c.o:
	$(CC) -c $<
.F.f:
	$(CPP) $(CPFLAGS) $(DEFS) $< > $@
	
# first target encountered is the default build

default: $(PROG) 

$(PROG): $(OBJS) blas_lapack
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

blas_lapack: libmyblas.a libmylapack.a
libmyblas.a:
	cd ../../BLAS; make double FC="${FC}" FFLAGS="${FFLAGS}" BLAS_EXCLUDE_LIST="${BLAS_EXCLUDE_LIST}";\
        cp libmyblas.a ../PATHSAMPLE/source

libmylapack.a:
	cd ../../LAPACK; make selection FC="${FC}" FFLAGS="${FFLAGS}" NOOPT="${NOOPT}";\
        cp libmylapack.a ../PATHSAMPLE/source

porfuncs.f90: porfuncs.csh
	./porfuncs.csh ${NAGSWITCH} > porfuncs.f90

clean:
	rm -f $(PROG) $(OBJS) *.mod porfuncs.f90 main.f *.lst *.a
	cd ../../BLAS; make clean
	cd ../../LAPACK; make clean
$(OBJS): libmyblas.a libmylapack.a
KMC.o:     porfuncs.o commons.o
KMCcommit.o:      commons.o porfuncs.o
addperm.o:      porfuncs.o commons.o utils.o
centre.o:      commons.o
charmmdump.o:      commons.o utils.o
connectodata.o:    commons.o
checkspodata.o:    commons.o utils.o
newconnodata.o:    commons.o utils.o
addminxyzodata.o:    commons.o utils.o
dock.o:         porfuncs.o commons.o nodes.o
getallpaths.o:     porfuncs.o commons.o utils.o
getallmin.o:     porfuncs.o commons.o 
getnewpath.o:      porfuncs.o commons.o 
inertia.o:         commons.o rigidbodymod.o
keywords.o:      porfuncs.o nodes.o commons.o rigidbodymod.o
main.f: main.F
main.o:      commons.o porfuncs.o rigidbodymod.o dock.o
nodes.o:     porfuncs.o commons.o
setup.o:      porfuncs.o utils.o commons.o
setup_sis.o:	porfuncs.o utils.o commons.o
tssearch.o:   porfuncs.o commons.o
unresdump.o:  commons.o
Pfold.o:    commons.o porfuncs.o
cycle.o:    commons.o porfuncs.o
cycle2.o:    commons.o porfuncs.o utils.o savestate.o
mind.o: porfuncs.o
mindouble.o: commons.o
mysystem.o:   porfuncs.o
mywait.o:     porfuncs.o
nodes.o:           porfuncs.o
submitoptimjob.o:  porfuncs.o nodes.o commons.o
tsdouble.o:   commons.o
getdpair.o:   commons.o
diagonalise2.o: porfuncs.o commons.o
Dijkstra.o: commons.o porfuncs.o
Dijinit.o: commons.o porfuncs.o
Dijinitfly.o: commons.o porfuncs.o
connectd.o: commons.o
frictionfac.o: commons.o
getrpair.o: commons.o
getspair.o: commons.o porfuncs.o
getupair.o: commons.o porfuncs.o
getfreepair.o: commons.o porfuncs.o
getpair.o: commons.o porfuncs.o
pairdouble.o: commons.o 
donedouble.o: commons.o 
mergedb.o:     porfuncs.o commons.o utils.o
GT.o: porfuncs.o commons.o GT2.o GT2input.o GT2FreeMemory.o savestate.o
NGT.o: porfuncs.o commons.o savestate.o NGTmem.o
NGTrealloc.o: NGTmem.o
NGTrenorm.o: NGTmem.o commons.o
NGTremovei.o: commons.o NGTmem.o porfuncs.o
NGTremoveid.o: commons.o NGTmem.o
NGT_crstorage.o: porfuncs.o commons.o savestate.o NGTmem.o
NGTrealloc_crstorage.o: commons.o NGTmem.o
NGTrenorm_crstorage.o: NGTmem.o
NGTremovei_crstorage.o: commons.o NGTmem.o porfuncs.o
NGTremoveid_crstorage.o: commons.o NGTmem.o
GT2.o: commons.o GT2data.o GT2FibonacciHeap.o GT2DLL.o
GT2FibonacciHeap.o: GT2data.o
GT2DLL.o: GT2data.o
GT2input.o: GT2data.o GT2DLL.o
GT2FreeMemory.o: GT2data.o GT2DLL.o
Gthomson.o : commons.o
minpermdist.o: commons.o porfuncs.o rotations.o  Gthomson.o 
minpermdistrbcom.o: commons.o
newmindist.o:   commons.o rigidbodymod.o
rbinertia.o: commons.o
rbperm.o: commons.o rotations.o vec3.o
rotations.o: vec3.o
rigidb.o: commons.o
rigidbodymod.o: mathsconstants.o
regroupfree.o: commons.o savestate.o
regroupfree2.o: commons.o savestate.o
regrouppersist.o: commons.o savestate.o
regroup.o: commons.o
rotations.o: vec3.o
orderodata.o: commons.o
kshortestpaths: commons.o porfuncs.o
reweight.o: porfuncs.o commons.o
Cv.o: commons.o
DOS.o: commons.o
bulkmindist.o: commons.o
make_conpot.o: commons.o
rateconst_setup.o: commons.o porfuncs.o
checkTS.o: commons.o porfuncs.o
lopermdist.o: commons.o 
getmetric.o: commons.o porfuncs.o
myorient.o: commons.o 
microtherm.o: commons.o utils.o
getbarrier2.o: commons.o
persistence.o: commons.o
shannon.o: commons.o porfuncs.o savestate.o utils.o
get_exact_persistence.o: commons.o
printsummary.o:	commons.o porfuncs.o
