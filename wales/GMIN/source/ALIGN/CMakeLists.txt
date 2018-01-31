file(GLOB ALIGN_SOURCES *.f90)
add_library(ALIGN ${ALIGN_SOURCES})
add_dependencies(ALIGN gminlib)
