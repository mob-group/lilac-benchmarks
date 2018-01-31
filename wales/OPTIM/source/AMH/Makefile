ARCH = ar
ARCHFLAGS = cr
ARCHDELFLAGS = d
RANLIB = ranlib
VPATH = ..
AMHLIB = ../libamh.a

FMODULES_SOURCES = amhglobals.f altpot_interfaces.f90  amh_interfaces.f90  globals_alt.f90 key.f90 
FMODULES_OBJECTS = altpot_interfaces.o  amh_interfaces.o  globals_alt.o amhglobals.o key.o
FMODULES = altpot_interfaces.mod altpot_interfaces.mod amhglobals.mod globals_alt.mod key.mod 
F_RAW_SOURCES = $(wildcard *.f)
F_SOURCES = $(filter-out $(FMODULES_SOURCES),$(F_RAW_SOURCES))
F90_RAW_SOURCES = $(wildcard *.f90)
F90_SOURCES = $(filter-out $(FMODULES_SOURCES),$(F90_RAW_SOURCES))
F_OBJECTS = $(F_SOURCES:%.f=%.o)
F90_OBJECTS = $(F90_SOURCES:%.f90=%.o)

all: single
single: $(AMHLIB)

$(F_OBJECTS),$(F90_OBJECTS): $(FMODULES)

$(FMODULES_OBJECTS): $(FMODULES) 

${AMHLIB}: $(FMODULES_OBJECTS) $(F_OBJECTS) $(F90_OBJECTS)
	$(ARCH) $(ARCHFLAGS) $(AMHLIB) $(F_OBJECTS) $(F90_OBJECTS) $^

amhglobals.mod: amhglobals.f
	$(FC) $(FFLAGS)  -c $<  

amh_interfaces.mod: amh_interfaces.f90 amhglobals.mod
	$(FC) $(FFLAGS)  -c $< 

globals_alt.mod: globals_alt.f90 amhglobals.mod
	$(FC) $(FFLAGS)  -c $<  

altpot_interfaces.mod: altpot_interfaces.f90 amhglobals.mod globals_alt.mod
	$(FC) $(FFLAGS)  -c $< 

%.o: %.f
	$(FC) $(FFLAGS)  -c $< 

%.o: %.f90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm -f *.o $(F_OBJECTS) $(F90_OBJECTS) ${AMHLIB}  *.lst *mod

################### DEPENDENCIES ######################################

Q_bias_seg_a.o:      amhglobals.o
Q_bias_seg_b.o:      amhglobals.o
additive_ev.o:      amhglobals.o
dssp_hdrgn_eastwood.o:      amhglobals.o
dssp_hforce.o:      amhglobals.o
force.o:         amhglobals.o
force.o:         amh_interfaces.o
force.o:         globals_alt.o
gaspot.o:      amhglobals.o
gengrd.o:      amhglobals.o
gentab.o:      amhglobals.o
gentab.o:      amh_interfaces.o
getmem.o:      amhglobals.o
gomb.o:      amhglobals.o
harm_spring.o:      amhglobals.o
hforce.o:      amhglobals.o
initil.o:      amhglobals.o
intstr.o:      amhglobals.o
lookup.o:      amhglobals.o
lookup.o:      amh_interfaces.o
opnfil.o:      amhglobals.o
oxy.o:      amhglobals.o
oxy.o:      amh_interfaces.o
potind.o:      amhglobals.o
qchrgmk.o:      amhglobals.o
qchrgmk.o:      globals_alt.o
rama.o:      amhglobals.o
restrt.o:        amhglobals.o
restrt.o:        key.o
rndcol.o:      amhglobals.o
savstr.o:      amhglobals.o
scltab.o:      amhglobals.o
walesamh_initial.o:      amhglobals.o
zero_amh.o:      amhglobals.o
E_write.o:      amhglobals.o
Rg_bias.o:      amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_interfaces.o:             amhglobals.o
altpot_interfaces.o:             globals_alt.o
altpot_io.o:        globals_alt.o
altpot_io.o:        altpot_interfaces.o
altpot_io.o:        globals_alt.o
altpot_io.o:        altpot_interfaces.o
altpot_io.o:        amhglobals.o
altpot_io.o:        amhglobals.o
altpot_io.o:        globals_alt.o
altpot_io.o:        altpot_interfaces.o
altpot_io.o:        globals_alt.o
altpot_io.o:        amhglobals.o
altpot_io.o:        globals_alt.o
altpot_io.o:        amhglobals.o
altpot_io.o:        globals_alt.o
altpot_onebody.o:        globals_alt.o
altpot_onebody.o:        amhglobals.o
altpot_onebody.o:        globals_alt.o
altpot_onebody.o:        amhglobals.o
altpot_onebody.o:        globals_alt.o
altpot_onebody.o:        altpot_interfaces.o
altpot_onebody.o:        amhglobals.o
altpot_onebody.o:        globals_alt.o
altpot_onebody.o:        globals_alt.o
altpot_pairpotential.o:        amhglobals.o
altpot_pairpotential.o:        globals_alt.o
altpot_pairpotential.o:        amhglobals.o
altpot_pairpotential.o:        globals_alt.o
altpot_pairpotential.o:        altpot_interfaces.o
altpot_pairpotential.o:        amhglobals.o
altpot_pairpotential.o:        globals_alt.o
altpot_pairpotential.o:        amhglobals.o
altpot_pairpotential.o:        globals_alt.o
altpot_pairpotential.o:        amhglobals.o
altpot_pairpotential.o:        globals_alt.o
altpot_pairpotential.o:        altpot_interfaces.o
altpot_subroutines.o:        amhglobals.o
altpot_subroutines.o:        globals_alt.o
altpot_subroutines.o:        altpot_interfaces.o
altpot_subroutines.o:        amhglobals.o
altpot_subroutines.o:        globals_alt.o
altpot_subroutines.o:        altpot_interfaces.o
altpot_utilities.o:        amhglobals.o
altpot_utilities.o:        amhglobals.o
altpot_utilities.o:        globals_alt.o
altpot_utilities.o:        globals_alt.o
altpot_utilities.o:        amhglobals.o
altpot_utilities.o:        globals_alt.o
altpot_utilities.o:        amhglobals.o
altpot_utilities.o:        amhglobals.o
altpot_utilities.o:        amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
amh_interfaces.o:             amhglobals.o
contact_P_AP.o:      amhglobals.o
ev_gamma.o:      amhglobals.o
ev_set_up.o:      amhglobals.o
hdrgn.o:      amhglobals.o
non_add_contact.o:      amhglobals.o
num_to_char.o:      amhglobals.o
rep_bias.o:      amhglobals.o
rep_contact.o:      amhglobals.o
walesamh_interface.o:     amhglobals.o
