info:
	@ echo
	@ echo "Use:"
	@ echo "make install "
	@ echo "  possible: make OPT=<...> with OPT the name of"
	@ echo "      any of the subdirectories in SRC; default=reference"
	@ echo "  possible: make MACH=<...> with MACH any identifier;"
	@ echo "      object files will be saved in SRC/<OPT>/<MACH>"
	@ echo "  possible: make DUMP_FLAGS=<...> with"
	@ echo "      flags: -DDUMP, -DTRACE_VECTORS"
	@ echo "  possible: make F_OPT_FLAGS=<...> with -O default"
	@ echo " make clean : remove all executables and object files"
	@ echo " make total_clean : also remove crsmat and regmat files"
	@ echo

OPT = reference
PLAT = `Scripts/arch`
MACH = default_machine
install :
	@ echo ; \
	  if [ "$(OPT)" = "reference" ] ; \
	  then echo "Installing benchmark code" ; \
	  else echo "Installing benchmark code variant <$(OPT)>" ; fi ;\
	  if [ ! -d "SRC/$(OPT)/$(MACH)" ] ; \
	    then mkdir SRC/$(OPT)/$(MACH) ; fi ; \
	  if [ -f "Make.$(PLAT)" ] ; \
	    then USE_PLAT=$(PLAT) ; else USE_PLAT=default_platform ; fi ; \
	  cd SRC/$(OPT)/$(MACH) ; XDIR="../../.."; \
	  echo ; echo ".. in directory `pwd`" ; echo ; \
	  cp ../../Make.copy_up . ; \
	  cp ../../../Make.* . ; \
	   cp ../Makefile.archs Makefile ; \
	    cp ../../Makefile.clean . ; cp ../../Make.files . ; \
	     $(MAKE) PLAT=$$USE_PLAT install ; \
	       for file in bench_gen bench_sym crs_gen crs_sym reg_gen reg_sym ; do \
	         if [ -x $$XDIR/$$file ] ; then /bin/rm -f $$XDIR/$$file; fi; \
	         ln $$file $$XDIR ; done ; \
	  echo ; echo ".. installation complete " ; echo
total_test :
	@ for var in SRC/* ; do \
	  if [ -d $$var ] ; then \
	    make OPT=`echo $$var | sed 's/OPT\///'` install; \
	    Validate 10 20 ; \
	  fi ; done
plot_data : Out/$(MACH)/$(OPT)/plot_stamp
Out/$(MACH)/$(OPT)/plot_stamp : Out/$(MACH)/$(OPT)/test_stamp
	@ Scripts/plot_data.pl $(MACH) $(OPT)
	@ touch $@
Out/$(MACH)/$(OPT)/test_stamp :
	@ touch $@
clean_plots :
	@/bin/rm -f Out/*/*/*{grf,dat,plot,stamp}
total_clean :: clean
	@ /bin/rm -rf regmat* crsmat* \
	  Out/*
	@ echo "Cleaning up machine sources" ; \
	  cd SRC ; for dir in * ; do \
	  if [ -d $$dir ] ; \
	  then cd $$dir ; \
	       echo "in variant directory $$dir" ; \
	       for file in * ; \
	       do if [ -d $$file ] ; \
                  then echo ".. removing machine $$file" ; \
	               /bin/rm -rf $$file ; \
	          fi ; \
	       done ; cd .. ; \
	  fi ; done	
local_clean ::
	@ /bin/rm -f matdump* vecdump.m out *.log *.plot \
	  validation_results validation.log 1
clean :: local_clean
	@ /bin/rm -f bench_* crs_* reg_* plot.ps \
	  SRC/*/$(MACH)/bench* Out/$(MACH)/bench* \
	  Out/$(MACH)/validation_results Out/$(MACH)/validation.log \
	  Out/$(MACH)/Results.*

include Makefile.clean

