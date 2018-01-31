# This file sets up the package CHARMMXX.
# Compile the charmm libraries and specify the list of CHARMM libraries that
# should be linked to in the CHARMM GMIN/OPTIM target.
#

option(WITH_DFTB "Enable compilation with DFTB code" OFF)

set(CHARMM_SIZE medium CACHE STRING "Choose the charmm build size, options are: huge, xxlarge, xlarge, large, medium, small, xsmall, reduce")
set(CHARMM_COMPILE_ARGS gnu ${CHARMM_SIZE})

if(COMPILER_SWITCH MATCHES "pgi")
  set(CHARMM_COMPILE_ARGS ${CHARMM_COMPILE_ARGS} keepo keepf PGF90 OPT)
elseif(COMPILER_SWITCH MATCHES "ifort")
  set(CHARMM_COMPILE_ARGS ${CHARMM_COMPILE_ARGS} keepo keepf IFORT OPT)
else(COMPILER_SWITCH MATCHES "pgi")
  message( FATAL_ERROR "CHARMM35 can only be compiled with ifort or pgi" )
endif(COMPILER_SWITCH MATCHES "pgi")

if (WITH_CHARMM35)
  find_package(CHARMM35)
endif ()
if(WITH_CHARMM36)
  find_package(CHARMM36)
endif ()
