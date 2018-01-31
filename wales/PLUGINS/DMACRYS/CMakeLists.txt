# So far build is with cmake.
# run run_cython.sh to generate all interface files
cmake_minimum_required(VERSION 2.8)
find_package(PythonLibs)

include_directories(../include)
 
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message("sources ${sources}")

include_directories(include ${PYTHON_INCLUDE_PATH})

enable_language(Fortran)
enable_language(C)
enable_language(CXX)

file(GLOB PLDMACRYS_SOURCES *.f90 *.f *.c ./dmacrys/*.f90 ./dmacrys/*.c)
add_library(pldmacrys SHARED ${PLDMACRYS_SOURCES})
set_target_properties(pldmacrys PROPERTIES PREFIX "")

