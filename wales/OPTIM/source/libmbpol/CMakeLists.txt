cmake_minimum_required(VERSION 2.6)

project(libmbpol)
# Francesco Paesani's MBPOL flexible water potential
set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "cv320@cam.ac.uk")

# optimisation of polynomials takes forever
# would be better as a CMAKE variable

if(${COMPILER_SWITCH} STREQUAL "g++")
  set(CMAKE_CXX_FLAGS_RELEASE -O0 -fPIC)
else()
  set(CMAKE_CXX_FLAGS_RELEASE -O0)
endif()

file(GLOB LIBMBPOL_SOURCES *.cpp *.h)
add_library(mbpollib ${LIBMBPOL_SOURCES})
