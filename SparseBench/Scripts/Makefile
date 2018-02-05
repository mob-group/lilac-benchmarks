info :
	@ echo "make install, make clean"

PLAT = default_platform

include ../Make.$(PLAT)

install : lsq

lsq : lsq.o
	@ $(FC) -o lsq lsq.o

clean ::
	@ /bin/rm -f lsq
include ../Makefile.clean
