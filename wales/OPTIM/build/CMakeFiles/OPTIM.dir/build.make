# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/bruce/code/sparse/wales/OPTIM/source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bruce/code/sparse/wales/OPTIM/build

# Include any dependencies generated for this target.
include CMakeFiles/OPTIM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/OPTIM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/OPTIM.dir/flags.make

CMakeFiles/OPTIM.dir/getparams.f.o: CMakeFiles/OPTIM.dir/flags.make
CMakeFiles/OPTIM.dir/getparams.f.o: /home/bruce/code/sparse/wales/OPTIM/source/getparams.f
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bruce/code/sparse/wales/OPTIM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/OPTIM.dir/getparams.f.o"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/bruce/code/sparse/wales/OPTIM/source/getparams.f -o CMakeFiles/OPTIM.dir/getparams.f.o

CMakeFiles/OPTIM.dir/getparams.f.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/OPTIM.dir/getparams.f.i"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/bruce/code/sparse/wales/OPTIM/source/getparams.f > CMakeFiles/OPTIM.dir/getparams.f.i

CMakeFiles/OPTIM.dir/getparams.f.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/OPTIM.dir/getparams.f.s"
	/usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/bruce/code/sparse/wales/OPTIM/source/getparams.f -o CMakeFiles/OPTIM.dir/getparams.f.s

CMakeFiles/OPTIM.dir/getparams.f.o.requires:

.PHONY : CMakeFiles/OPTIM.dir/getparams.f.o.requires

CMakeFiles/OPTIM.dir/getparams.f.o.provides: CMakeFiles/OPTIM.dir/getparams.f.o.requires
	$(MAKE) -f CMakeFiles/OPTIM.dir/build.make CMakeFiles/OPTIM.dir/getparams.f.o.provides.build
.PHONY : CMakeFiles/OPTIM.dir/getparams.f.o.provides

CMakeFiles/OPTIM.dir/getparams.f.o.provides.build: CMakeFiles/OPTIM.dir/getparams.f.o


# Object files for target OPTIM
OPTIM_OBJECTS = \
"CMakeFiles/OPTIM.dir/getparams.f.o"

# External object files for target OPTIM
OPTIM_EXTERNAL_OBJECTS =

OPTIM: CMakeFiles/OPTIM.dir/getparams.f.o
OPTIM: CMakeFiles/OPTIM.dir/build.make
OPTIM: liboptimlib.a
OPTIM: libextralib.a
OPTIM: libmylapack.a
OPTIM: libmyblas.a
OPTIM: libBOWMAN_LIB.a
OPTIM: libmbpollib.a
OPTIM: /home/bruce/code/sparse/wales/MYFFTW/local_build/install/lib/libfftw3.a
OPTIM: CMakeFiles/OPTIM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bruce/code/sparse/wales/OPTIM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable OPTIM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OPTIM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/OPTIM.dir/build: OPTIM

.PHONY : CMakeFiles/OPTIM.dir/build

CMakeFiles/OPTIM.dir/requires: CMakeFiles/OPTIM.dir/getparams.f.o.requires

.PHONY : CMakeFiles/OPTIM.dir/requires

CMakeFiles/OPTIM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/OPTIM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/OPTIM.dir/clean

CMakeFiles/OPTIM.dir/depend:
	cd /home/bruce/code/sparse/wales/OPTIM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bruce/code/sparse/wales/OPTIM/source /home/bruce/code/sparse/wales/OPTIM/source /home/bruce/code/sparse/wales/OPTIM/build /home/bruce/code/sparse/wales/OPTIM/build /home/bruce/code/sparse/wales/OPTIM/build/CMakeFiles/OPTIM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/OPTIM.dir/depend

