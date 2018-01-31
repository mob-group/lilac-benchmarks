cmake_minimum_required(VERSION 2.6)

if (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")
  message( FATAL_ERROR "CMake should not be built in the source directory" )
endif (CMAKE_CURRENT_SOURCE_DIR MATCHES "^${CMAKE_BINARY_DIRECTORY}$")

project(disconnectionDPS)

get_filename_component(SVN_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../.. ABSOLUTE)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
# CMake modules implement some additional functionality to check
# for external dependencies
set(CMAKE_MODULE_PATH ${SVN_ROOT}/CMakeModules)

# if not specified by user, the standard build type is release
if(NOT CMAKE_BUILD_TYPE)
  #Release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# the FORTRANCOMPILER check is in svn/CMakeModules
# it adds a compiler switch and sets flags
find_package(FORTRANCOMPILER)

# add the executables
add_executable(disconnectionDPS disconnectionDPS.f90)
add_executable(manipulate manipulate.f90)
