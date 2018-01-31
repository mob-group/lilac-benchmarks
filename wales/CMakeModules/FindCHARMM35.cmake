# This file sets up the package CHARMM35.
# Compile the charmm libraries and specify the list of CHARMM libraries that
# should be linked to in the CHARMM GMIN/OPTIM target.
#
#

# set the root of the charmm source directory
set(CHARMM_SOURCE_ROOT ${SVN_ROOT}/CHARMM35)

# set the out-of-source directory where charmm will be built
unset(CHARMM_BUILD_ROOT CACHE)
set(CHARMM_BUILD_ROOT ${CMAKE_BINARY_DIR}/C35 CACHE TYPE STRING)
file(MAKE_DIRECTORY ${CHARMM_BUILD_ROOT})

# Note: CHARMM compilation will fail if the file name of the ${CHARMM_BUILD_ROOT} is longer
# than 49 characters (for chm_host = gnu).  (50 fails, 48 works)
# This is because the charmm pre-processor sets a hard limit on the length of 
# the path to include files (fcm files) of 60 in total. If the build path is 
# too long, checking the original charmm directory to see if that can be used 
# instead.
string(LENGTH ${CHARMM_SOURCE_ROOT} howlong)
string(LENGTH ${CHARMM_BUILD_ROOT} howlongtwo)
set(EXTRACHARMMARGS "")
if(${howlongtwo} GREATER 49)
  if(${howlong} GREATER 49)
     message(FATAL_ERROR "Paths ${CHARMM_BUILD_ROOT} and ${CHARMM_SOURCE_ROOT} are both too long for charmm35 to deal with. One of them should be 49 characters or less")
  else(${howlong} GREATER 49)
     set(EXTRACHARMMARGS fcmdir ${CHARMM_SOURCE_ROOT})
  endif(${howlong} GREATER 49)
endif(${howlongtwo} GREATER 49)

set(LIBDIR ${CHARMM_BUILD_ROOT}/lib/gnu)

# Set the CHARMM_LIBS variable. These need to be specified manually because
# only certain of them should be linked to in the CHARMM GMIN target.
unset(CHARMM_LIBS CACHE)
set(CHARMM_LIBS 
  ${LIBDIR}/help.o ${LIBDIR}/iniall.o ${LIBDIR}/miscom.o ${LIBDIR}/usersb.o
  ${LIBDIR}/squantm.a ${LIBDIR}/io.a ${LIBDIR}/cheq.a ${LIBDIR}/moldyn.a
  ${LIBDIR}/emap.a ${LIBDIR}/solvation.a ${LIBDIR}/dimb.a ${LIBDIR}/graphics.a
  ${LIBDIR}/pipf.a ${LIBDIR}/correl.a ${LIBDIR}/pert.a ${LIBDIR}/mmff.a
  ${LIBDIR}/machdep.a ${LIBDIR}/shapes.a ${LIBDIR}/energy.a ${LIBDIR}/dynamc.a
  ${LIBDIR}/util.a ${LIBDIR}/gener.a ${LIBDIR}/flucq.a ${LIBDIR}/molvib.a
  ${LIBDIR}/vibran.a ${LIBDIR}/cadint.a ${LIBDIR}/manip.a ${LIBDIR}/minmiz.a
  ${LIBDIR}/quantum.a ${LIBDIR}/gamint.a ${LIBDIR}/nbonds.a ${LIBDIR}/zerom.a
  ${LIBDIR}/mscale.a ${LIBDIR}/prate.a ${LIBDIR}/adumb.a ${LIBDIR}/rxncor.a
  ${LIBDIR}/image.a ${LIBDIR}/mndint.a ${LIBDIR}/mc.a ${LIBDIR}/sccdftbint.a
  ${LIBDIR}/misc.a ${LIBDIR}/gukint.a ${LIBDIR}/mbond.a ${LIBDIR}/cff.a
  ${LIBDIR}/squantm.a ${LIBDIR}/io.a ${LIBDIR}/cheq.a ${LIBDIR}/moldyn.a
  ${LIBDIR}/emap.a ${LIBDIR}/solvation.a ${LIBDIR}/dimb.a ${LIBDIR}/graphics.a
  ${LIBDIR}/pipf.a ${LIBDIR}/correl.a ${LIBDIR}/pert.a ${LIBDIR}/mmff.a
  ${LIBDIR}/machdep.a ${LIBDIR}/shapes.a ${LIBDIR}/energy.a ${LIBDIR}/dynamc.a
  ${LIBDIR}/util.a ${LIBDIR}/gener.a ${LIBDIR}/flucq.a ${LIBDIR}/molvib.a
  ${LIBDIR}/vibran.a ${LIBDIR}/cadint.a ${LIBDIR}/manip.a ${LIBDIR}/minmiz.a
  ${LIBDIR}/quantum.a ${LIBDIR}/gamint.a ${LIBDIR}/nbonds.a ${LIBDIR}/zerom.a
  ${LIBDIR}/mscale.a ${LIBDIR}/prate.a ${LIBDIR}/adumb.a ${LIBDIR}/rxncor.a
  ${LIBDIR}/image.a ${LIBDIR}/mndint.a ${LIBDIR}/mc.a ${LIBDIR}/sccdftbint.a
  ${LIBDIR}/misc.a ${LIBDIR}/gukint.a ${LIBDIR}/mbond.a ${LIBDIR}/cff.a
  CACHE TYPE STRING) 

if(WITH_DFTB)
  set(CHARMM_COMPILE_ARGS ${CHARMM_COMPILE_ARGS} T)
  unset(CHARMM_LIBS CACHE)
  set(CHARMM_LIBS 
    ${LIBDIR}/help.o ${LIBDIR}/iniall.o ${LIBDIR}/miscom.o ${LIBDIR}/usersb.o
    ${LIBDIR}/squantm.a ${LIBDIR}/io.a ${LIBDIR}/cheq.a ${LIBDIR}/moldyn.a
    ${LIBDIR}/emap.a ${LIBDIR}/solvation.a ${LIBDIR}/dimb.a
    ${LIBDIR}/graphics.a ${LIBDIR}/pipf.a ${LIBDIR}/correl.a ${LIBDIR}/pert.a
    ${LIBDIR}/mmff.a ${LIBDIR}/machdep.a ${LIBDIR}/shapes.a ${LIBDIR}/energy.a
    ${LIBDIR}/dynamc.a ${LIBDIR}/util.a ${LIBDIR}/gener.a ${LIBDIR}/flucq.a
    ${LIBDIR}/molvib.a ${LIBDIR}/vibran.a ${LIBDIR}/cadint.a ${LIBDIR}/manip.a
    ${LIBDIR}/minmiz.a ${LIBDIR}/quantum.a ${LIBDIR}/gamint.a
    ${LIBDIR}/nbonds.a ${LIBDIR}/zerom.a ${LIBDIR}/mscale.a ${LIBDIR}/prate.a
    ${LIBDIR}/adumb.a ${LIBDIR}/rxncor.a ${LIBDIR}/image.a ${LIBDIR}/mndint.a
    ${LIBDIR}/mc.a ${LIBDIR}/sccdftb.a ${LIBDIR}/sccdftbint.a ${LIBDIR}/misc.a
    ${LIBDIR}/gukint.a ${LIBDIR}/mbond.a ${LIBDIR}/cff.a ${LIBDIR}/squantm.a
    ${LIBDIR}/io.a ${LIBDIR}/cheq.a ${LIBDIR}/moldyn.a ${LIBDIR}/emap.a
    ${LIBDIR}/solvation.a ${LIBDIR}/dimb.a ${LIBDIR}/graphics.a
    ${LIBDIR}/pipf.a ${LIBDIR}/correl.a ${LIBDIR}/pert.a ${LIBDIR}/mmff.a
    ${LIBDIR}/machdep.a ${LIBDIR}/shapes.a ${LIBDIR}/energy.a
    ${LIBDIR}/dynamc.a ${LIBDIR}/util.a ${LIBDIR}/gener.a ${LIBDIR}/flucq.a
    ${LIBDIR}/molvib.a ${LIBDIR}/vibran.a ${LIBDIR}/cadint.a ${LIBDIR}/manip.a
    ${LIBDIR}/minmiz.a ${LIBDIR}/quantum.a ${LIBDIR}/gamint.a
    ${LIBDIR}/nbonds.a ${LIBDIR}/zerom.a ${LIBDIR}/mscale.a ${LIBDIR}/prate.a
    ${LIBDIR}/adumb.a ${LIBDIR}/rxncor.a ${LIBDIR}/image.a ${LIBDIR}/mndint.a
    ${LIBDIR}/mc.a ${LIBDIR}/sccdftb.a ${LIBDIR}/sccdftbint.a ${LIBDIR}/misc.a
    ${LIBDIR}/gukint.a ${LIBDIR}/mbond.a ${LIBDIR}/cff.a 
    CACHE TYPE STRING)
endif(WITH_DFTB)

mark_as_advanced(CHARMM_LIBS)
mark_as_advanced(CHARMM_BUILD_ROOT)

# jmc> Adding ${EXTRACHARMMARGS} to the end of the args of install.com to 
# potentially tell it to use the fcm (Fortran common) files in the svn charmm source 
# directory rather than the copies in the 
# the build directory, for cases where the build path is too long for the charmm pre-processor
# but the svn path is OK.
set(CHARMM_COMPILE_ARGS ${CHARMM_COMPILE_ARGS} ${EXTRACHARMMARGS})

message("${SVN_ROOT}/CMakeModules/FindCHARMM35.cmake: compiling CHARMM35 with arguments = ${CHARMM_COMPILE_ARGS}")                      
message("${SVN_ROOT}/CMakeModules/FindCHARMM35.cmake: charmm libs = ${CHARMM_LIBS}")

# Define the command that builds the CHARMM libraries.  If the build is
# successful then file CHARM_WAS_BUILT_FILE is created.  This file is
# persistant, so CHARMM need not be rebuild every time the makefile is called.
# Do not make this file a dependency, use CHARMM_WAS_BUILT as defined below.
# We are building out-of-source, so copy the relevant source files to the build
# directory first.
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CHARMM_WAS_BUILT_FILE
   COMMAND rm -rf build tool source install.com clean.csh
   COMMAND mkdir build
   COMMAND cp -r ${CHARMM_SOURCE_ROOT}/build/UNX build/UNX
   COMMAND cp -r ${CHARMM_SOURCE_ROOT}/tool tool
   COMMAND cp -r ${CHARMM_SOURCE_ROOT}/source source
   COMMAND cp ${CHARMM_SOURCE_ROOT}/install.com .
   COMMAND cp ${CHARMM_SOURCE_ROOT}/clean.csh .
   COMMAND ./clean.csh
   COMMAND ./install.com ${CHARMM_COMPILE_ARGS} > build.log && touch ${CMAKE_CURRENT_BINARY_DIR}/CHARMM_WAS_BUILT_FILE
   WORKING_DIRECTORY ${CHARMM_BUILD_ROOT}
)


# NOTE: this doesn't work
#set(ALL_CHARMM_LIBS ${CHARMM_LIBS})
#add_custom_target( charmmlib
#  COMMAND ar crsv libcharmm_merged.a ${ALL_CHARMM_LIBS}
#  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CHARMM_WAS_BUILT_FILE
#  )
#unset(CHARMM_LIBS CACHE)
#SET(CHARMM_LIBS_NEW ${CMAKE_CURRENT_BINARY_DIR}/libcharmm_merged.a CACHE TYPE STRING)

# Create target CHARMM_WAS_BUILT that depends on CHARMM_WAS_BUILT_FILE.  Any
# target that depends on CHARMM can list CHARMM_WAS_BUILT as a dependency.
# This is to avoid problems with parallel compilation.
add_custom_target(CHARMM_WAS_BUILT DEPENDS CHARMM_WAS_BUILT_FILE)
