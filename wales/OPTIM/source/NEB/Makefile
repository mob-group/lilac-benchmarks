#   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
#   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
#   This file is part of NEB module. NEB module is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# Multiple directories make build for OPTIM
#
#    - source tree is organised into `root' directory and `module' subdirectories;
#    - Makefiles in module directories produce library files which are put into root directory;
#    - naming convention for library files: lib[module directory name].a
#
#    - resolving dependencies: 
#
#         1. top down
#            To specify that file `keywords.o' from root directory depends on file(s) in `newconnect' module directory
#            insert `in.o: libnc.a' into root Makefile;
#
#         2. bottom up
#            a. if `ncutils.o' file in `newconnect' directory depends on `mindist.o' file from root directory,
#               insert "libnc.a: mindist.o" into root Makefile
#            b. include the dependence in the corresponding module Makefile,
#               e.g., for the above case, "ncutils.o: mindist.o"
#                   
#              (root Makefile will recompile mind.o before descending into module directory,
#               module Makefile will do necessary recompilations in module directory)
#
#         3. module on module - similar to #2
#              a. insert "libnc.a: libnn.a" into root Makefile if newconnect module depends on newneb module
#              b. same as 2.b.
#              c. adjust VPATH variable so the modules can `see' each other
#
#    Semen Trygubenko  Tue Nov  9 23:24:36 GMT 2004
#

KEYOBJS = key_output.o key_sqvv.o key_grad.o key_tau.o key_minimizer.o key_lbfgs.o 
OBJS = nndata.o grad.o oldnebgradient.o lbfgs.o newneb.o output.o sqvv.o nnutils.o ntc.o key_neb.o tau.o unresnebbfgs.o \
charutils.o linkedlist.o ${KEYOBJS}
VPATH = ../:../CHARMM/
LIBNAME = ../libnn.a

.SUFFIXES:
.SUFFIXES: .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

../libnn.a: ${OBJS}
	ar sr ${LIBNAME} ${OBJS}

clean:
	rm -f ${LIBNAME} ${OBJS} *.mod ifc?????? *.il *.stb *.o

grad.o: key_grad.o key_tau.o nndata.o key_neb.o tau.o specialfuncts.o nnutils.o modunres.o key.o commons.o nndata.o key.o porfuncs.o nndata.o key_neb.o modunres.o key.o intcommons.o modcharmm.o
key_grad.o: charutils.o
key_lbfgs.o: charutils.o key.o
key_minimizer.o:          charutils.o key_lbfgs.o key_sqvv.o
key_neb.o: key.o key_grad.o key_output.o key_minimizer.o key_sqvv.o key_lbfgs.o key_tau.o charutils.o key.o
key_sqvv.o: charutils.o
lbfgs.o: key_minimizer.o key_lbfgs.o nndata.o key_neb.o nnutils.o grad.o output.o porfuncs.o modcharmm.o \
         key.o intcommons.o intcoords.o specialfuncts.o commons.o genrigid.o
newneb.o: porfuncs.o nndata.o key_neb.o lbfgs.o output.o nnutils.o key.o growstring.o gsdata.o modguess.o modmec.o intcommons.o intcoords.o specialfuncts.o ntc.o AMH/amhglobals.o sqvv.o unresnebbfgs.o
nnutils.o: nndata.o key_neb.o key.o commons.o \
           specialfuncts.o intcommons.o modcharmm.o modamber9.o porfuncs.o charutils.o AMH/amhglobals.o 
oldnebgradient.o: key_grad.o nndata.o key_neb.o tau.o nnutils.o key.o commons.o 
output.o: key.o gsdata.o key_output.o modcharmm.o nndata.o key_neb.o ntc.o charutils.o modguess.o modmec.o \
          linkedlist.o efol.o intcommons.o commons.o nndata.o linkedlist.o
sqvv.o: key_minimizer.o nndata.o key_sqvv.o key_neb.o grad.o nnutils.o charutils.o
tau.o: nndata.o key_neb.o modunres.o key_tau.o key.o commons.o
unresnebbfgs.o: key_minimizer.o key_lbfgs.o key_neb.o nndata.o nnutils.o grad.o output.o modunres.o porfuncs.o key.o
