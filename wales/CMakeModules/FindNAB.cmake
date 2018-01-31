
# set the root of the nab source directory
set(NAB_SOURCE_ROOT ${SVN_ROOT}/NAB)

# set the out-of-source directory where nab will be built
set(NAB_BUILD_ROOT ${CMAKE_CURRENT_BINARY_DIR}/NAB)
file(MAKE_DIRECTORY ${NAB_BUILD_ROOT})


# Copy the NAB sources to the build directory and build.
# If successfull create file nab_binaries_built
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/nab_binaries_built
  #OUTPUT libnab.a DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/nab_binaries_built
  COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/nab_binaries_built
  COMMAND rm -rf src ucpp-1.3 byacc semantics Makefile config.h bin cifparse
  COMMAND cp -r ${NAB_SOURCE_ROOT}/Makefile .
  COMMAND cp -r ${NAB_SOURCE_ROOT}/config.h .
  COMMAND cp -r ${NAB_SOURCE_ROOT}/src src
  COMMAND cp -r ${NAB_SOURCE_ROOT}/ucpp-1.3 ucpp-1.3
  COMMAND cp -r ${NAB_SOURCE_ROOT}/byacc byacc
  COMMAND cp -r ${NAB_SOURCE_ROOT}/semantics semantics
  COMMAND cp -r ${NAB_SOURCE_ROOT}/cifparse cifparse
  COMMAND mkdir bin
  COMMAND make clean
  COMMAND cd ucpp-1.3/ && make -j1 install NABHOME="${NAB_BUILD_ROOT}"
  COMMAND cd byacc/ && make -j1 install NABHOME="${NAB_BUILD_ROOT}"
  COMMAND make -j1 install NABHOME="${NAB_BUILD_ROOT}" LIBDIR="${NAB_BUILD_ROOT}/lib" && echo "nab built" > ${CMAKE_CURRENT_BINARY_DIR}/nab_binaries_built
  #COMMAND cp lib/libnab.a ${CMAKE_CURRENT_BINARY_DIR}/libnab.a && echo "nab built" > ${CMAKE_CURRENT_BINARY_DIR}/nab_binaries_built
  WORKING_DIRECTORY ${NAB_BUILD_ROOT}
)

add_custom_target (nab_build DEPENDS nab_binaries_built)
add_library(NAB_LIB STATIC IMPORTED DEPENDS nab_build)
set_target_properties( NAB_LIB PROPERTIES IMPORTED_LOCATION ${NAB_BUILD_ROOT}/lib/libnab.a )
add_dependencies(NAB_LIB nab_build)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES nab_binaries_built)
