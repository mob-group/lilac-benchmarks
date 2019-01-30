###########################################################
#
# MAIN CMAKE FILE FOR PATHSAMPLE 
# ========================
#
# Usage: 
#   
#   $ mkdir ~/svn/PATHSAMPLE/builds/pgf90 
#   $ cd ~/svn/PATHSAMPLE/builds/pgf90 
#   $ FC=pgf90 cmake ../../source 
#   $ make -j8

cmake_minimum_required(VERSION 2.8)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(PATHSAMPLE)

set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "dw34@cam.ac.uk")
get_filename_component(SVN_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)

EXECUTE_PROCESS(COMMAND "date" "+%d%m" OUTPUT_VARIABLE CURRENT_DATE)
if ($ENV{USER} MATCHES "wales" AND CURRENT_DATE MATCHES "0104")
  message( FATAL_ERROR "PATHSAMPLE 9000: I'm sorry, David. I'm afraid I can't do that." )
endif ($ENV{USER} MATCHES "wales" AND CURRENT_DATE MATCHES "0104")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${SVN_ROOT}/CMakeModules)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Release Debug Debug_slow RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_language(CXX)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)

# get the svn version and add it to the -D definitions
execute_process(
  COMMAND bash version.sh SVNREV
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE SVNVERSION OUTPUT_STRIP_TRAILING_WHITESPACE
  )
set(SVNVERSION "'${SVNVERSION}'")
message("SVNVERSION ${SVNVERSION}")
add_definitions(-DSVNVERSION=${SVNVERSION})


message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/modules)

if(DEFINED ENV{SPMV_VERSION})
  if("$ENV{SPMV_VERSION}" STREQUAL "mkl")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lmkl-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "mkl-slow")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L$ENV{MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -lmkl-slow-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "gpu")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L$ENV{CUDA_ROOT}/lib64 -lcublas -lcusparse -lcudart -lgpu-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "unsafe")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L$ENV{CUDA_ROOT}/lib64 -lcublas -lcusparse -lcudart -lunsafe-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl00")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl00-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl01")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl01-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl10")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl10-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl11")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl11-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl00-slow")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl00-slow-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl01-slow")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl01-slow-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl10-slow")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl10-slow-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl11-slow")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lopencl11-slow-spmv")
  elseif("$ENV{SPMV_VERSION}" STREQUAL "sparsex")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lsparsex-spmv")
  else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lnative-spmv")
  endif()
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lnative-spmv")
endif()

set(CMAKE_Fortran_FLAGS "-L$ENV{SPMV_ROOT}/lib ${CMAKE_Fortran_FLAGS}")
message(STATUS "${CMAKE_Fortran_FLAGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -D_GNU_SOURCE")

include_directories(${CMAKE_SOURCE_DIR})

file(GLOB PATHSAMPLE_MAIN main.F)
file(GLOB PATHSAMPLE_SOURCES *.f *.f90 *.F *.F90)
file(GLOB NOT_PATHSAMPLE_SOURCES
    main.f
    main.F
    porfuncs.f90
)
if (NOT (NOT_PATHSAMPLE_SOURCES STREQUAL ""))
  list(REMOVE_ITEM PATHSAMPLE_SOURCES ${NOT_PATHSAMPLE_SOURCES})
endif (NOT (NOT_PATHSAMPLE_SOURCES STREQUAL ""))


# ---- add generated source files 
set(PATHSAMPLE_SOURCES ${PATHSAMPLE_SOURCES} porfuncs.f90)

# ---- porfuncs.f90 TODO: shuold it build every time?
#message("compiler switch ${COMPILER_SWITCH}")
add_custom_command (
   OUTPUT porfuncs.f90 
   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/porfuncs.sh ${COMPILER_SWITCH} > porfuncs.f90)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES porfuncs.f90)

# ----- Libraries 
# ----- generate a pathsample library to avoid recompiling the whole code for every executable
add_library(pathsamplelib ${PATHSAMPLE_SOURCES})
set_target_properties(pathsamplelib PROPERTIES
  LINK_FLAGS efwkpo
)

find_package(MYBLAS)
find_package(MYLAPACK)

add_definitions(-D_SVN_ROOT_='${SVN_ROOT}')

# --- Finally, the execs

# standard PATHSAMPLE
add_executable(PATHSAMPLE ${PATHSAMPLE_MAIN})
target_link_libraries(PATHSAMPLE pathsamplelib ${MYLAPACK_LIBS})
set_target_properties(PATHSAMPLE PROPERTIES LINKER_LANGUAGE "Fortran") 

if(DEFINED ENV{SPMV_VERSION})
  if("$ENV{SPMV_VERSION}" STREQUAL "mkl")
    target_link_libraries(PATHSAMPLE mkl-spmv)
  elseif("$ENV{SPMV_VERSION}" STREQUAL "gpu")
    target_link_libraries(PATHSAMPLE gpu-spmv)
  elseif("$ENV{SPMV_VERSION}" STREQUAL "unsafe")
    target_link_libraries(PATHSAMPLE unsafe-spmv)
  elseif("$ENV{SPMV_VERSION}" STREQUAL "opencl")
    target_link_libraries(PATHSAMPLE opencl-spmv)
  elseif("$ENV{SPMV_VERSION}" STREQUAL "clgpu")
    target_link_libraries(PATHSAMPLE clgpu-spmv)
  else()
    target_link_libraries(PATHSAMPLE native-spmv)
  endif()
else()
  target_link_libraries(PATHSAMPLE native-spmv)
endif()
