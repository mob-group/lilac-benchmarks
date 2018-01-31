###########################################################
#
# MAIN CMAKE FILE FOR OPTIM 
# ========================
#
# Usage: 
#   
#   $ mkdir ~/svn/OPTIM/builds/pgf90
#   $ cd ~/svn/OPTIM/builds/pgf90
#   $ FC=pgf90 cmake ../../source 
#   $ make -j8 <target>    
#              -- if <target> not specified, all active ones are built  
#
# Good to know: 
#   1. VERBOSE=1 make   
#         -- gives the full compilation command for each file, useful for
#         debugging  
#
#   2. To compile with MPI 
#         FC=mpif90 cmake ../source -DCOMPILER_SWITCH=pgi 
#         make OPTIM 
#       
#   3. GUI for CMakeCache.txt  
#         ccmake .    (after running cmake command) 
#         e.g.: 
#              set WITH_AMBER12 to ON -> Press [c] -> Press [e] -> Press [g]  
#           equivalent command line version
#            $ cmake . -DWITH_AMBER12=yes 
# 
#   4. Out-of-source build only 
#         do not run cmake in ~svn/OPTIM/source    
#
#   5. Tested with 
#        $ cmake --version 
#
#   6. To do a clean build, first delete CMakeCache.txt, and if that doesn't
#   help, delete the whole directory!
# 
############################################################

cmake_minimum_required(VERSION 2.8)
# 2.8 required particularly for the charmm-interfaced targets.

if (OPTIM_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (OPTIM_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(OPTIM)

set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "dw34@cam.ac.uk")
get_filename_component(SVN_ROOT ${OPTIM_SOURCE_DIR}/../.. ABSOLUTE)

EXECUTE_PROCESS(COMMAND "date" "+%d%m" OUTPUT_VARIABLE CURRENT_DATE)
if ($ENV{USER} MATCHES "wales" AND CURRENT_DATE MATCHES "0104")
  message( FATAL_ERROR "OPTIM 9000: I'm sorry, David. I'm afraid I can't do that." )
endif ($ENV{USER} MATCHES "wales" AND CURRENT_DATE MATCHES "0104")

if (CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(CMAKE_C_FLAGS_RELEASE "-O1 -DNDEBUG")
endif (CMAKE_C_COMPILER_ID MATCHES "Intel")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${SVN_ROOT}/CMakeModules)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Release Debug Debug_slow RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)
find_package(MYBLAS)
find_package(MYLAPACK)

message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
set(CMAKE_Fortran_MODULE_DIRECTORY ${OPTIM_BINARY_DIR}/modules)
include_directories(${CMAKE_SOURCE_DIR})

file(GLOB DUMMY_AMH     amhdummy.f90)
file(GLOB DUMMY_CHARMM  chdummy.f)
file(GLOB DUMMY_AMBER12 ${SVN_ROOT}/AMBER12/amber12_interface.F90)
file(GLOB DUMMY_AMBER9  amber9dummy.f90)
file(GLOB DUMMY_USERPOT dummy_userpot.f90)
file(GLOB DUMMY_UNRES   unresdummy.f)
file(GLOB DUMMY_DLF     dlfdummy.f90)
file(GLOB DUMMY_CPPNEB  dummy_cppneb.f90)
file(GLOB DUMMY_QUIP  optim_quipdummy.f90)
file(GLOB DUMMY_CUDA modcudadummylbfgs.f90 modcudadummybfgsts.f90)
file(GLOB DUMMY_OPEP    OPEPinterface/opep_interface.F90)

set(ALL_DUMMIES ${DUMMY_AMH}
                ${DUMMY_CHARMM}
                ${DUMMY_AMBER12}
                ${DUMMY_AMBER9}
                ${DUMMY_CUDA}
                ${DUMMY_USERPOT}
                ${DUMMY_UNRES} 
                ${DUMMY_DLF} 
                ${DUMMY_CPPNEB}
                ${DUMMY_QUIP}
                ${DUMMY_OPEP} )

# Glob all the sources
file(GLOB OPTIM_LIB_SOURCES *.f
                            *.f90
                            *.F
                            *.F90
                            NEB/*.f90
                            CONNECT/*.f90 
                            ALIGN/*.f90
                            AMH/amhglobals.f )
#                           sparse/*.f90
#                           sparse/*.c ) 

file(GLOB NOT_OPTIM_SOURCES getparams.f
                            dsygvx.f
                            amhdummy.f
                            #			    OPTIM.f
                            porfuncs.f90
                            header.f90
                            optim_quip_wrapper.f90 
                            modcudalbfgs.f90
                            modcudabfgsts.f90)                           

option(WITH_ALIGN "Enable ALIGN subroutines" ON)
if(WITH_ALIGN)
  file(GLOB DUMMY_ALIGN dummy_aligndecide.f90)
  list(REMOVE_ITEM OPTIM_LIB_SOURCES ${DUMMY_ALIGN})
else()
  file(GLOB ALIGN_SOURCES ALIGN/*.f90)
  list(REMOVE_ITEM OPTIM_LIB_SOURCES ${ALIGN_SOURCES})
endif(WITH_ALIGN)


# Due to a compiler bug in ifort 13.1.3, we can't use -O3 for genrigid.f90
# Investigations continue...
# There also seems to be a bug in the compiler which sometimes causes TSLOCATOR to go into an infinite loop instead of returning.
# This is removed by using O0 optimisation for NEB/output.f90.
if( ${COMPILER_SWITCH} STREQUAL "ifort" )
  SET_SOURCE_FILES_PROPERTIES( genrigid.f90 PROPERTIES COMPILE_FLAGS -O2 )
  SET_SOURCE_FILES_PROPERTIES( NEB/output.f90 PROPERTIES COMPILE_FLAGS -O0 )
endif ( ${COMPILER_SWITCH} STREQUAL "ifort" )
# And nag has an internal error of some kind with optimisation for dbpg.f90
if( ${COMPILER_SWITCH} STREQUAL "nag" )
  SET_SOURCE_FILES_PROPERTIES( dbpg.f90 
                               dmblpy.f90 PROPERTIES COMPILE_FLAGS -O1 )
endif ( ${COMPILER_SWITCH} STREQUAL "nag" )

list(REMOVE_ITEM OPTIM_LIB_SOURCES ${NOT_OPTIM_SOURCES}
                                   ${ALL_DUMMIES} )

# Add generated source files 
set(OPTIM_LIB_SOURCES ${OPTIM_LIB_SOURCES}
                      header.f90
                      porfuncs.f90 )

# porfuncs.f90
add_custom_command (
   OUTPUT porfuncs.f90 
   COMMAND bash ${OPTIM_SOURCE_DIR}/porfuncs.sh ${COMPILER_SWITCH} > porfuncs.f90
)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES porfuncs.f90)

# header.f90 
add_custom_command (
  OUTPUT header.f90 
  COMMAND bash ${OPTIM_SOURCE_DIR}/header.sh > header.f90
)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES header.f90)

# Collect all the extra source files into one place. It's just dummies at the moment.
set(ALL_EXTRA_SOURCES ${ALL_DUMMIES})

########################
# LIBRARIES START HERE #
########################

# Make an extras library for dummies, this is necessary because some of the
# files in optimlib use modules made by extralib
add_library(extralib ${ALL_EXTRA_SOURCES})
set_module_dir(extralib)
set_target_properties(extralib PROPERTIES LINKER_LANGUAGE "Fortran") 
set_target_properties(extralib PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA;__SPARSE")

# Make an optim library
add_library(optimlib ${OPTIM_LIB_SOURCES})
set_module_dir(optimlib)
set_module_depends(optimlib extralib)
set_target_properties(optimlib PROPERTIES LINKER_LANGUAGE "Fortran") 

add_subdirectory(Bowman)

# See notes in GMIN CMakeLists to explain this block. No MPI compilation for OPTIM, so simpler structure here.
if(WITH_ALIGN)
  include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
  set(FFTW_COMPILER "gfortran" CACHE TYPE STRING)
  mark_as_advanced(FORCE FFTW_COMPILER)
  ExternalProject_Add(fftw
    SOURCE_DIR ${SVN_ROOT}/MYFFTW
    PREFIX ${SVN_ROOT}/MYFFTW/local_build
    CONFIGURE_COMMAND ${SVN_ROOT}/MYFFTW/configure F77=${FFTW_COMPILER} --prefix=${SVN_ROOT}/MYFFTW/local_build/install --disable-doc
    BUILD_COMMAND make
    INSTALL_COMMAND ${MAKE_INSTALL})
  set_target_properties(fftw PROPERTIES Fortran_MODULE_DIRECTORY "${PREFIX}")
  set_module_depends(optimlib fftw)
  target_link_libraries(optimlib ${SVN_ROOT}/MYFFTW/local_build/install/lib/libfftw3.a)
endif(WITH_ALIGN)

include_directories(libmbpol)
add_subdirectory(libmbpol)
set_module_depends(optimlib mbpollib)
# target_link_libraries(optimlib libmbpol)

# Link CHOLMOD
#add_subdirectory(${SVN_ROOT}/SuiteSparse SuiteSparse)
#target_link_libraries(optimlib cholmod)
#include_directories(${SVN_ROOT}/SuiteSparse/SuiteSparse/CHOLMOD/Include)
#include_directories(${SVN_ROOT}/SuiteSparse/SuiteSparse/SuiteSparse_config)
#add_definitions(-D__SPARSE)

add_definitions(-D_SVN_ROOT_='${SVN_ROOT}')


##########################
# EXECUTABLES START HERE #
##########################

# standard OPTIM
add_executable(OPTIM getparams.f)
set_target_properties(OPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
set_module_depends(OPTIM optimlib
                         extralib )
target_link_libraries(OPTIM optimlib
                            extralib
                            ${MYLAPACK_LIBS}
                            ${BOWMAN_LIBS} 
			    mbpollib)
# AMHOPTIM
option(WITH_AMH "Enable AMHOPTIM compilation" OFF)
if(WITH_AMH)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
# AMH has issues with the e_write subroutine. The amh_interfaces_mod has the interface for e_write, but
# does not define it and so causes problems in linking. 
# At the moment, I just link in the dummy module instead.
#  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_AMH})
  add_subdirectory(${OPTIM_SOURCE_DIR}/AMH)
  add_executable(AMHOPTIM getparams.f ${EXTRA_SOURCES})
  set_module_dir(AMHOPTIM)
  set_module_depends(AMHOPTIM optimlib mbpollib)
  set_module_dir(amh)
  set_module_depends(amh optimlib
                         extralib )
  set_target_properties(AMHOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(AMHOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(AMHOPTIM optimlib 
                                 amh
                                 ${MYLAPACK_LIBS}
                                 ${BOWMAN_LIBS} 
				 mbpollib)
endif(WITH_AMH)

# A9OPTIM  
option(WITH_AMBER9 "Enable Amber OPTIM with AMBER9 compilation" OFF)
if(WITH_AMBER9)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_AMBER9})
  find_package(NAB)
  add_subdirectory(AMBER) # temporary folder to copy files from sander;  
  add_executable(A9OPTIM getparams.f
                         ${EXTRA_SOURCES} )
  set_module_dir(A9OPTIM)
  set_module_depends(A9OPTIM optimlib)
  set_module_dir(AMBER_LIB)
  set_module_depends(AMBER_LIB optimlib)
  set_target_properties(A9OPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(A9OPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA;__SPARSE")
  target_link_libraries(A9OPTIM optimlib
                                AMBER_LIB
                                NAB_LIB
                                ${MYLAPACK_LIBS}
                                ${BOWMAN_LIBS} 
				mbpollib)
endif(WITH_AMBER9)

# For A12OPTIM or CUDAOPTIM, we need to build and add the AMBER 12 interface.
# We also need to install networkx, for the chirality and cis/trans scripts.
if(WITH_AMBER12 OR WITH_CUDA)
  add_subdirectory(${SVN_ROOT}/AMBER12 AMBER12_interface)
#  find_program(PYTHON "python")
#  if (PYTHON)
#    add_custom_command(TARGET optimlib
#                       POST_BUILD
#                       COMMAND easy_install --user networkx)
#  else(PYTHON)
#    message(FATAL_ERROR "CMake cannot find Python.")
#  endif(PYTHON)
endif()

# A12OPTIM  
option(WITH_AMBER12 "Enable Amber OPTIM with AMBER12 compilation" OFF)
# A12OPTIM and CUDAOPTIM both use the same Amber interface. Thus, we can only build the CPU version or the
# GPU version at any given time.
if(WITH_AMBER12 AND NOT WITH_CUDA)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  add_executable(A12OPTIM getparams.f
                          ${EXTRA_SOURCES})
  set_module_dir(A12OPTIM)
  set_module_depends(A12OPTIM optimlib
                              amber12_base )
  set_target_properties(A12OPTIM PROPERTIES LINKER_LANGUAGE "Fortran") 
  set_target_properties(A12OPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_CUDA") 
  target_link_libraries(A12OPTIM optimlib
                                 amber12_base
                                 ${MYLAPACK_LIBS}
                                 ${BOWMAN_LIBS} 
				 mbpollib)
endif(WITH_AMBER12 AND NOT WITH_CUDA)

# CUDAOPTIM
option(WITH_CUDA "Enable CUDAOPTIM compilation" OFF)
if(WITH_CUDA)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES} modcudabfgsts.f90 modcudalbfgs.f90)
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_CUDA})
  add_subdirectory(${SVN_ROOT}/CUDA CUDAinterface)
  add_executable(CUDAOPTIM getparams.f ${EXTRA_SOURCES})
  set_module_dir(CUDAOPTIM)
  set_module_depends(CUDAOPTIM optimlib amber12_base)
  set_target_properties(CUDAOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(CUDAOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};CUDA")
  target_link_libraries(CUDAOPTIM optimlib CUDAinterface amber12_base ${MYLAPACK_LIBS} ${BOWMAN_LIBS} mbpollib)
endif(WITH_CUDA)

#SPIN
option(WITH_SPIN "Enable SPINOPTIM compilation" OFF)
if(WITH_SPIN)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_USERPOT})
  add_subdirectory(SPINinterface)
  add_executable(SPINOPTIM getparams.f
                           ${EXTRA_SOURCES})
  set_module_dir(SPINOPTIM)
  set_module_depends(SPINOPTIM optimlib)
  set_target_properties(SPINOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(SPINOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(SPINOPTIM optimlib
                                  SPIN
                                  ${MYLAPACK_LIBS}
                                  ${BOWMAN_LIBS} 
				  mbpollib)
endif(WITH_SPIN)

#OPEP
option(WITH_OPEP "Enable OPEPOPTIM compilation" OFF)
if(WITH_OPEP)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  add_subdirectory(OPEPinterface)
  add_executable(OPEPOPTIM getparams.f
                           ${EXTRA_SOURCES} )
  set_module_dir(OPEPOPTIM)
  set_module_depends(OPEPOPTIM optimlib)
  set_module_dir(OPEP)
  set_module_depends(OPEP optimlib)
  set_target_properties(OPEPOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(OPEPOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA;__OPEP")
  target_link_libraries(OPEPOPTIM optimlib
                                  OPEP
                                  ${MYLAPACK_LIBS}
                                  ${BOWMAN_LIBS} 
				  mbpollib)
endif(WITH_OPEP)



# QUIP
option(WITH_QUIP "Enable QUIPOPTIM compilation" OFF)
if(WITH_QUIP)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES} optim_quip_wrapper.f90)
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_QUIP})
  add_executable(QUIPOPTIM getparams.f ${EXTRA_SOURCES})
  set_module_dir(QUIPOPTIM)
  set_module_depends(QUIPOPTIM optimlib)
  set_target_properties(QUIPOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
 set_target_properties(QUIPOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
# set_target_properties(QUIPOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA;BLAS_LIBRARIES=${BLAS};LAPACK_LIBRARIES=${LAPACK}")
# set_target_properties(QUIPOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA;BLAS_LIBRARIES=$BLAS;LAPACK_LIBRARIES=$LAPACK")
  target_link_libraries(QUIPOPTIM optimlib /home/wales/QUIP/build/linux_x86_64_gfortran/libquip.a ${MYLAPACK_LIBS} ${BOWMAN_LIBS} mbpollib)
endif(WITH_QUIP)


option(WITH_PLUGIN "Enable shared library plugin support" OFF)
if(WITH_PLUGIN)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_USERPOT})
  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
  add_subdirectory(pluginloader)
  add_executable(PLOPTIM getparams.f
                           ${EXTRA_SOURCES} )
  set_module_dir(PLOPTIM)
  set_module_depends(PLOPTIM optimlib)
  set_target_properties(PLOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(PLOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(PLOPTIM optimlib
                                pluginloader
                                ${MYLAPACK_LIBS}
                                ${BOWMAN_LIBS} 
				mbpollib)
endif(WITH_PLUGIN)

option(WITH_CPPNEB "Build with new C++ nudged elastic band support" OFF)
if(WITH_CPPNEB)
#  SET(CMAKE_CXX_FLAGS -std=c++0x)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_CPPNEB})
  add_subdirectory(cppneb)
  add_dependencies(cppneb optimlib)
  target_link_libraries(optimlib cppneb)
  set_module_dir(cppneb)
  set_module_depends(cppneb optimlib)
  add_executable(OPTIM_CPPNEB getparams.f
                              ${EXTRA_SOURCES} )
  set_target_properties(OPTIM_CPPNEB PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(OPTIM_CPPNEB PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(OPTIM_CPPNEB optimlib ${MYLAPACK_LIBS} ${BOWMAN_LIBS} mbpollib)
endif(WITH_CPPNEB)

option(WITH_CHARMM35 "Enable C35OPTIM compilation" OFF)
if(WITH_CHARMM35)
  if(WITH_CHARMM36)
    message( FATAL_ERROR "cannot build CHARMM35 and CHARMM36 at the same time" )
  endif(WITH_CHARMM36)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_CHARMM})
  add_subdirectory(${OPTIM_SOURCE_DIR}/CHARMMinterface)
  add_executable(C35OPTIM getparams.f
                          myblas.f
                          mylapack.f
                          ${EXTRA_SOURCES})
  set_module_dir(C35OPTIM)
  set_module_depends(C35OPTIM optimlib)
  set_module_dir(charmm35_interface_lib)
  set_module_depends(charmm35_interface_lib optimlib)
  set_target_properties(C35OPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(C35OPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  # add the CHARMM library as a dependency
  add_dependencies(C35OPTIM CHARMM_WAS_BUILT)
  # Link in the CHARMM libraries and the CHARMM35 OPTIM interface
  target_link_libraries(C35OPTIM optimlib
                                 mylapack
                                 myblascharmm
                                 charmm35_interface_lib
                                 ${CHARMM_LIBS} 
                                 ${BOWMAN_LIBS} 
				 mbpollib)
endif(WITH_CHARMM35)

# C36OPTIM
option(WITH_CHARMM36 "Enable C36OPTIM compilation" OFF)
if(WITH_CHARMM36)
  if(WITH_CHARMM35)
    message( FATAL_ERROR "cannot build CHARMM35 and CHARMM36 at the same time" )
  endif(WITH_CHARMM35)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  LIST(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_CHARMM})
  add_subdirectory(${OPTIM_SOURCE_DIR}/CHARMMinterface)
  add_executable(C36OPTIM getparams.f
                          myblas.f
                          mylapack.f
                          ${EXTRA_SOURCES})
  set_module_dir(C36OPTIM)
  set_module_depends(C36OPTIM optimlib)
  set_module_dir(charmm36_interface_lib)
  set_module_depends(charmm36_interface_lib optimlib)
  add_dependencies(C36OPTIM CHARMM_WAS_BUILT)
  target_link_libraries(C36OPTIM optimlib
                                 mylapack
                                 myblascharmm
                                 charmm36_interface_lib
                                 ${CHARMM_LIBS}
                                 ${BOWMAN_LIBS} 
				 mbpollib)
endif(WITH_CHARMM36)

option(WITH_DLF "Enable DLFOPTIM compilation" OFF)
if(WITH_DLF)
  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_DLF})
  add_subdirectory(DLFWrapper)
  # include_directories(libmbpol)
  add_executable(DLFOPTIM ${EXTRA_SOURCES} getparams.f)
  set_module_dir(DLFOPTIM)
  set_module_depends(DLFOPTIM optimlib)
  set_target_properties(DLFOPTIM PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(DLFOPTIM PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(DLFOPTIM optimlib
                                 ${MYLAPACK_LIBS}
                                 ${BOWMAN_LIBS}
                                 DLF_LIB
				 mbpollib)
endif(WITH_DLF)



# Unres interface here for possible future reference. NOT to be routinely used
#option(WITH_UNRES "Enable UNOPTIM compilation" OFF)
#if(WITH_UNRES)
#  SET(EXTRA_SOURCES ${ALL_EXTRA_SOURCES})
#  list(REMOVE_ITEM EXTRA_SOURCES ${DUMMY_UNRES})
#  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/UNRESinterface)
#  add_executable(UNOPTIM ${EXTRA_SOURCES})
#  target_link_libraries(UNOPTIM optimlib ${MYLAPACK_LIBS} unresinterface, ${BOWMAN_LIBS})
#endif(WITH_UNRES)


# Test for duplicates.
# In the standard build we first merge the object files together in a library
# (optimlib) to avoid duplicate compilation for multiple targets.  A
# drawback of this is that duplicate symbols (duplicate names of subroutines,
# functions, etc.) are simply ignored.  This is a separate target which builds
# OPTIM directly from the sources, specifically to check for duplicate
# symbols.
option(WITH_DUPLICATE_TESTING "Compile OPTIM in such a way that duplicate symbols are flagged" OFF)
if(WITH_DUPLICATE_TESTING)
  add_executable(OPTIM_TEST_DUPLICATES getparams.f
                                       ${OPTIM_LIB_SOURCES}
                                       ${ALL_EXTRA_SOURCES} )
  set_module_dir(OPTIM_TEST_DUPLICATES)
  set_module_depends(OPTIM_TEST_DUPLICATES BOWMAN_LIB)
  set_target_properties(OPTIM_TEST_DUPLICATES PROPERTIES LINKER_LANGUAGE "Fortran")
  set_target_properties(OPTIM_TEST_DUPLICATES PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};DUMMY_AMBER12;DUMMY_CUDA")
  target_link_libraries(OPTIM_TEST_DUPLICATES ${MYLAPACK_LIBS}
                                              BOWMAN_LIB
					      mbpollib
					      ${SVN_ROOT}/MYFFTW/local_build/install/lib/libfftw3.a)
#                                             cholmod ) << This isn't linked to OPTIM yet
  add_dependencies(OPTIM_TEST_DUPLICATES optimlib)
endif(WITH_DUPLICATE_TESTING)
