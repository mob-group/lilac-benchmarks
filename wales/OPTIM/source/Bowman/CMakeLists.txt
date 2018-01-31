#
# build commands for OPTIM with Bowman
#
# Bowman (as the OPTIM interface is written) can't be compiled as a pure
# library because the the optim sources depend on bowmanwater.mod.  Thus the
# choice of whether to use fake or real Bowman must be decided at compile time,
# not at link time.  If WITH_BOWMAN is on then target OPTIM will be compiled
# against the full Bowman library, a separete target for JBOPTIM would simply
# be a duplicate

# Create a library for the dummy bowman sources
option(WITH_BOWMAN "Compile OPTIM with Bowman sources" OFF)
if (NOT WITH_BOWMAN)
  # create the fake bowman library

  file(GLOB sources fake/*.f90)
  add_library(BOWMAN_LIB ${sources})

  set(BOWMAN_LIBS BOWMAN_LIB PARENT_SCOPE)
  add_dependencies(optimlib BOWMAN_LIB)
  set_module_dir(BOWMAN_LIB)
  set_module_depends(optimlib BOWMAN_LIB)


else (NOT WITH_BOWMAN)
  # create the real Bowman library

  if(NOT COMPILER_SWITCH MATCHES "ifort")
    message( FATAL_ERROR "Bowman can only be compiled with ifort" )
  endif(NOT COMPILER_SWITCH MATCHES "ifort")

  # locate the precompiled librarires
  file(GLOB precompiled_lib_files libs/*a)
  message("lib files ${precompiled_lib_files}")

  # build the Bowman library
  file(GLOB sources src/*f90)
  message("sources ${sources}")
  add_library(BOWMAN_LIB ${sources})
  set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS "-r8 -O -u -I${CMAKE_CURRENT_SOURCE_DIR}/mod")

  # both the built and precompiled libraries should be linked against
  set(BOWMAN_LIBS BOWMAN_LIB ${precompiled_lib_files} PARENT_SCOPE)
  add_dependencies(optimlib BOWMAN_LIB)
  set_module_dir(BOWMAN_LIB)
  set_module_depends(optimlib BOWMAN_LIB)

    
endif (NOT WITH_BOWMAN)

