# Sketched-out cmake file for unres interface

# need to add the unres source directory below (not in the repository!)
include_directories(${CMAKE_BINARY_DIR} <path-to-unres-source>)
file(GLOB UNRES_SOURCES *.f *.f90 <path-to-unres-source>/*.f)

add_library(unresinterface ${UNRES_SOURCES})
ADD_DEPENDENCIES(unresinterface optimlib)
