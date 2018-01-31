# routines for myblas
# for now we add an extra library which is compiled within the project.
# The only reason for doing that is to ensure that compile flags are identical as
# in the rest of the project. It should in principle work without this trick
# for BLAS/MYBLAS (see BLAS/LAPACK  system libraries). My feeling is that 
# the compiler flag issue is a CHARMM artifact (vr274)

option(WITH_MYLAPACK "Compile own lapack (needed for charmm, can cause problems with gfortran 4.7)" ON)
if(WITH_MYLAPACK)
    file(GLOB MYLAPACK_SOURCES ${SVN_ROOT}/LAPACK/*.f)
    add_library(mylapack ${MYLAPACK_SOURCES})

    message("${SVN_ROOT}/CMakeModules/FindMYLAPACK.cmake: creating LAPACK library.")                      

    set(MYLAPACK_LIBS mylapack ${MYBLAS_LIBS})
    mark_as_advanced(MYLAPACK_LIBS)
else(WITH_MYLAPACK)
    find_package(LAPACK)
    message("using installed lapack: ${LAPACK_LIBRARIES}")
    set(MYLAPACK_LIBS ${LAPACK_LIBRARIES} ${MYBLAS_LIBS})
endif(WITH_MYLAPACK)
