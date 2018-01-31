# routines for myblas
# for now we add an extra library which is compiled within the project.
# The only reason for doing that is to ensure that compile flags are identical as
# in the rest of the project. It should in principle work without this trick
# for BLAS/MYBLAS (see BLAS/LAPACK  system libraries). My feeling is that 
# the compiler flag issue is a CHARMM artifact (vr274)

file(GLOB MYBLAS_SOURCES ${SVN_ROOT}/BLAS/*.f)
file(GLOB NOT_MYBLAS_SOURCES 
  ${SVN_ROOT}/BLAS/zhpmv.f
  ${SVN_ROOT}/BLAS/ztpsv.f
  ${SVN_ROOT}/BLAS/ztrsv.f
  ${SVN_ROOT}/BLAS/zgbmv.f
  ${SVN_ROOT}/BLAS/zgerc.f
  ${SVN_ROOT}/BLAS/zgemv.f
  ${SVN_ROOT}/BLAS/ztbsv.f
  ${SVN_ROOT}/BLAS/ztrsm.f
  ${SVN_ROOT}/BLAS/zher2k.f
  ${SVN_ROOT}/BLAS/zher2.f
  ${SVN_ROOT}/BLAS/zgemm.f
  ${SVN_ROOT}/BLAS/zhbmv.f
  ${SVN_ROOT}/BLAS/zhpr2.f
  ${SVN_ROOT}/BLAS/zhpr.f
  ${SVN_ROOT}/BLAS/ztbmv.f
  ${SVN_ROOT}/BLAS/zherk.f
  ${SVN_ROOT}/BLAS/ztrmm.f
  ${SVN_ROOT}/BLAS/ztpmv.f
  ${SVN_ROOT}/BLAS/zher.f
  ${SVN_ROOT}/BLAS/ztrmv.f
  ${SVN_ROOT}/BLAS/zhemm.f
  ${SVN_ROOT}/BLAS/zhemv.f
# charmm specific 
  ${SVN_ROOT}/BLAS/dnrm2.f 
  ${SVN_ROOT}/BLAS/daxpy.f 
  ${SVN_ROOT}/BLAS/dcopy.f 
  ${SVN_ROOT}/BLAS/ddot.f 
)

list(REMOVE_ITEM MYBLAS_SOURCES ${NOT_MYBLAS_SOURCES})

add_library(myblascharmm ${MYBLAS_SOURCES})

message("${SVN_ROOT}/CMakeModules/FindMYBLAScharmm.cmake: creating CHARMM specific BLAS library.") 

SET(MYBLAS_LIBS myblascharmm)
MARK_AS_ADVANCED(MYBLAS_LIBS)
