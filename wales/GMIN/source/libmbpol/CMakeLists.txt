cmake_minimum_required(VERSION 2.8)

# Silence warnings about PROJECT_VERSION not being set
cmake_policy(SET CMP0048 OLD)

project(libmbpol)
# Francesco Paesani's MBPOL flexible water potential
set(PROJECT_VERSION "0.0")
set(PROJECT_CONTACT "jdf43@cam.ac.uk")

# optimisation of polynomials takes forever
# would be better as a CMAKE variable

if(${COMPILER_SWITCH} STREQUAL "g++")
  set(CMAKE_CXX_FLAGS_RELEASE -O0 -fPIC)
else()
  set(CMAKE_CXX_FLAGS_RELEASE -O0)
endif()

file(GLOB LIBMBPOL_SOURCES *.cpp *.h)
add_library(libmbpol ${LIBMBPOL_SOURCES})
