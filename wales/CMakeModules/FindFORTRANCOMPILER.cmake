# setup fortran compiler

set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")

enable_language(Fortran)

GET_FILENAME_COMPONENT(FC_PROGNAME ${CMAKE_Fortran_COMPILER} NAME)
message("FC_PROGNAME = ${FC_PROGNAME}")

# if no compiler switch is given try to automatically determine it
if(NOT COMPILER_SWITCH)
   GET_FILENAME_COMPONENT(FC_PROGNAME ${CMAKE_Fortran_COMPILER} NAME)

   if(FC_PROGNAME MATCHES "pgf90" OR FC_PROGNAME MATCHES "pgfortran")
      set(COMPILER_SWITCH "pgi" CACHE TYPE STRING)
   elseif(FC_PROGNAME MATCHES "gfortran")
      set(COMPILER_SWITCH "gfortran" CACHE TYPE STRING)
   elseif(FC_PROGNAME MATCHES "nagfor")
      set(COMPILER_SWITCH "nag" CACHE TYPE STRING)
   elseif(FC_PROGNAME MATCHES "ifort")
      set(COMPILER_SWITCH "ifort" CACHE TYPE STRING)
   else()
      message(FATAL_ERROR "ERROR: could not determine compiler switch from fortran compiler name (${FC_PROGNAME})\nplease specify with -DCOMPILER_SWITCH=<pgi/gfortran/nag/ifort>")
   endif()
endif()

message("Compiler switch = ${COMPILER_SWITCH}")

# check if the compiler is mpif90 
# if("${CMAKE_Fortran_COMPILER}" MATCHES "mpif90$") 
#  message("MPI build") 
# endif("${CMAKE_Fortran_COMPILER}" MATCHES "mpif90$") 

# now setup some special variables
if(NOT COMPILER_FLAGS_WERE_SET)
   message("Setting initial values for compiler flags")
   if(COMPILER_SWITCH MATCHES "pgi")
      set (CMAKE_Fortran_FLAGS "-Mextend" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -Munroll -Mnoframe" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG "-Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Melf -Mpgicoff -traceback" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG_SLOW "-Mextend -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mpgicoff -traceback" CACHE TYPE STRING FORCE)
      set (FORTRAN_FREEFORM_FLAG "-Mfree" CACHE TYPE STRING)
   elseif(COMPILER_SWITCH MATCHES "gfortran")
      set (CMAKE_Fortran_FLAGS "-ffixed-line-length-200 -ffree-line-length-0" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE TYPE STRING FORCE)
#      set (CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG "-g -fbounds-check -Wuninitialized -O -ftrapv -fno-automatic" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG_SLOW "${CMAKE_Fortran_FLAGS_DEBUG} -fimplicit-none" CACHE TYPE STRING FORCE)
      set (FORTRAN_FREEFORM_FLAG "-ffree-form" CACHE TYPE STRING)
   elseif(COMPILER_SWITCH MATCHES "nag")
      set (CMAKE_Fortran_FLAGS "-132 -kind=byte -maxcontin=6000" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_RELEASE "-mismatch_all -O4" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG "-g -mismatch_all -ieee=stop" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG_SLOW "-C=all -mtrace=all -gline -g -mismatch_all -ieee=stop" CACHE TYPE STRING FORCE)
      set (FORTRAN_FREEFORM_FLAG "-free" CACHE TYPE STRING) # js850> is this ever used?
   elseif(COMPILER_SWITCH MATCHES "ifort")
      set (CMAKE_Fortran_FLAGS "-132 -heap-arrays -assume byterecl" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_RELEASE "-O3" CACHE TYPE STRING FORCE)
# Warnings about temporary argument creation and edit descriptor widths are disabled with the final flags.
      set (CMAKE_Fortran_FLAGS_DEBUG "-C -g -traceback -debug full -check all,noarg_temp_created -diag-disable 8290,8291" CACHE TYPE STRING FORCE)
      set (CMAKE_Fortran_FLAGS_DEBUG_SLOW "-debug all -check all,noarg_temp_created -implicitnone -warn unused -fp-stack-check -ftrapuv -check pointers -check bounds" CACHE TYPE STRING FORCE)
      set (FORTRAN_FREEFORM_FLAG "-free" CACHE TYPE STRING)
   else()
      message(FATAL_ERROR "unknown comiler switch: ${COMPILER_SWITCH}")
   endif()
    SET(COMPILER_FLAGS_WERE_SET yes CACHE TYPE INTERNAL)
endif(NOT COMPILER_FLAGS_WERE_SET)

mark_as_advanced(FORTRAN_FREEFORM_FLAG COMPILER_SWITCH COMPILER_FLAGS_WERE_SET CMAKE_Fortran_FLAGS_DEBUG_SLOW)

set(CMAKE_GEN_FILES ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY})

function(set_module_dir target_name)
# Sets the module output directory for ${target_name} to ${YOUR_BUILD_DIRECTORY}/CMakeFiles/${target_name}.dir/modules
# e.g. GMIN/build/CMakeFiles/gminlib.dir/modules for ${target_name} = gminlib
#
# Also adds its own modules directory to the include directory for that target.
  set_target_properties(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_GEN_FILES}/${target_name}.dir/modules)
  if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
      get_target_property(MOD_COMPILE_FLAGS ${target_name} COMPILE_FLAGS)
      if(NOT MOD_COMPILE_FLAGS)
          set(MOD_COMPILE_FLAGS "")
      endif(NOT MOD_COMPILE_FLAGS)
      set_target_properties(${target_name} PROPERTIES COMPILE_FLAGS "${MOD_COMPILE_FLAGS} -I${CMAKE_GEN_FILES}/${target_name}.dir/modules") 
#      message("Compile flags for ${target_name}: ${MOD_COMPILE_FLAGS}")
  else(${CMAKE_VERSION} VERSION_LESS 2.8.12)
      target_include_directories(${target_name} PUBLIC ${CMAKE_GEN_FILES}/${target_name}.dir/modules)
  endif(${CMAKE_VERSION} VERSION_LESS 2.8.12)
endfunction(set_module_dir)

function(set_module_depends target_using_mods)
# Adds the modules made by target ${target_making_mods} to the include directory for
# ${target_using_mods}.
  foreach(target_making_mods ${ARGN})
    get_target_property(TARGET_MOD_LOC ${target_making_mods} Fortran_MODULE_DIRECTORY)
#    message("${target_making_mods} mod location: ${TARGET_MOD_LOC}")
    if(${CMAKE_VERSION} VERSION_LESS 2.8.12)
        get_target_property(MOD_COMPILE_FLAGS ${target_using_mods} COMPILE_FLAGS)
        if(NOT MOD_COMPILE_FLAGS)
            set(MOD_COMPILE_FLAGS "")
        endif(NOT MOD_COMPILE_FLAGS)
#        message("${target_making_mods} mod location now: ${TARGET_MOD_LOC}")
        set_target_properties(${target_using_mods} PROPERTIES COMPILE_FLAGS "${MOD_COMPILE_FLAGS} -I${TARGET_MOD_LOC}") 
        get_target_property(MOD_COMPILE_FLAGS ${target_using_mods} COMPILE_FLAGS)
#        message("depends: Compile flags for ${target_using_mods}: ${MOD_COMPILE_FLAGS}")
    else(${CMAKE_VERSION} VERSION_LESS 2.8.12)
        target_include_directories(${target_using_mods} 
                                   PUBLIC ${TARGET_MOD_LOC})
    endif(${CMAKE_VERSION} VERSION_LESS 2.8.12)
    add_dependencies(${target_using_mods} ${target_making_mods})
  endforeach()
endfunction(set_module_depends)

