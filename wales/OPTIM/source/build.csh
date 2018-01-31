#!/bin/csh -f
#
#   OPTIM installation script
#   Copyright (C) 2003-2008 Semen A. Trygubenko and David J. Wales
#   This file is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# set parameters
set defaulttarget = optim
set defaultcompiler = ifort 
set defaultlevel = opt
set defaultexetype = dynamic
set defaultarchitecture = intel
set defaultchmsize = medium
set defaultprofile = 0
set defaultchmqm = none
set programname = OPTIM
set EmptyString = --debug=v
set SingleQuote = "'"
set DoubleQuote = '"'

# set defaults for variables
set USER = `echo $USER`
set MACHINE = `uname -n`
set OS = `uname -s`
set QMTYPE = ""
set CTYPE = c35
# default CHARMM directories
set c36src = $PWD"/../../CHARMM36"
set c35src = $PWD"/../../CHARMM35"
set c31src = $PWD"/../../CHARMM31"
#set amhsrc = "/home/$USER/svn/OPTIM/source/AMH"
set amhsrc = "AMH"
set target = $defaulttarget
set compiler = $defaultcompiler
set level = $defaultlevel
set exetype = $defaultexetype
set architecture = $defaultarchitecture
set chmsize = $defaultchmsize
set chmqm = $defaultchmqm
set optsrc = `pwd`
set quiet = $EmptyString
set profile = defaultprofile
set exemessage = 1

# process keyword-based arguments
set readnext = 0
foreach opt ( $argv )
 if ( $readnext == 1 ) then
   set readnext = 0
   if ( $optsrc == FILL ) then
     set optsrc = $opt
   else if ( $c31src == FILL ) then
     set c31src = $opt
   else if ( $c35src == FILL ) then
     set c35src = $opt
   else if ( $c36src == FILL ) then
     set c36src = $opt
   endif
 endif
 if ( $opt == optsrc ) then
   set readnext = 1
   set optsrc = FILL
 else if ( $opt == c31src ) then
   set readnext = 1
   set c31src = FILL
 else if ( $opt == c35src ) then
   set readnext = 1
   set c35src = FILL
 else if ( $opt == c36src ) then
   set readnext = 1
   set c36src = FILL
 endif
end

# make sure we are in the source directory. After this point all the references can be relative to the source dir
cd $optsrc

# read in the last arguments list into lastargs; if unavailable, initialise lastargs to an empty string
if ( -e lastargs ) then
  set lastargs = `cat lastargs`
else
  set lastargs = $EmptyString
endif

# check if we can make use of the lastargs file
if (! -e 'lastargs' && "$argv" == $EmptyString ) then # we'd better display help message
  set arguments = help
else if ( "$argv" == last ) then # file 'lastargs' exists and we are instructed to read it
  echo "build.csh> lastargs = $lastargs"
  set arguments = "$lastargs"
else
  set arguments = "$argv"
endif

# process the command line input
set readcompiler = 0
set readlevel = 0
set readtarget = 0
set readexetype = 0
set readarchitecture = 0
set readchmsize = 0
set readchmqm = 0
set readchmtype = 0
set readnext = 1
foreach opt ( $arguments )
  # keyword arguments and their values need to be ignored
  if ( $readnext != 1 ) then
    set readnext = 1
  else if ( $opt == 'optsrc' || $opt == 'c31src' || $opt == 'c35src' || $opt == 'c36src' ) then
    set readnext = 0
  else if ( $opt == opt || $opt == debug || $opt == debugslow || $opt == noopt ) then
    if ( $readlevel == 1 ) then
      echo "build.csh> ERROR: Two or more mutually exclusive compiler options 1 were specified"
      echo "build.csh> ERROR: First two conflicting arguments are '"$level"' and '"$opt"'"
      exit
    endif
    set level = $opt
    if ( $level == $defaultlevel ) then
      echo "build.csh> '$defaultlevel' is the default and can be omitted"
    endif
    set readlevel = 1
  else if ( $opt == optim || $opt == coptim || $opt == unoptim || $opt == amb9optim || $opt == amhoptim || $opt == dlfoptim || $opt == jboptim || $opt == quipoptim) then
    if ( $readtarget == 1 ) then
      echo "build.csh> ERROR: Two or more mutually exclusive targets were specified"
      echo "build.csh> ERROR: First two conflicting arguments are '"$target"' and '"$opt"'"
      exit
    endif
    set target = $opt
    if ( $target == $defaulttarget ) then
      echo "build.csh> '$defaulttarget' is the default and can be omitted"
    endif
    set readtarget = 1
  else if ( $opt == pgi || $opt == ifort || $opt == nag || $opt == ibm || $opt == g95 || $opt == gfortran || $opt == pathscale ) then
    if ( $readcompiler == 1 ) then
      echo "build.csh> ERROR: Two or more mutually exclusive compiler types were specified"
      echo "build.csh> ERROR: First two conflicting arguments are '"$compiler"' and '"$opt"'"
      exit
    endif
    set compiler = $opt
    if ( $compiler == $defaultcompiler ) then
      echo "build.csh> '$defaultcompiler' is the default and can be omitted"
    endif
    set readcompiler = 1
  else if ( $opt == dynamic || $opt == static || $opt == nobflag ) then
    if ( $readexetype == 1 ) then
      echo "build.csh> ERROR: Two mutually exclusive compiler options 2 were specified"
      echo "build.csh> ERROR: Conflicting arguments are '"$exetype"' and '"$opt"'"
      exit
    endif
    set exetype = $opt
    set readexetype = 1
    if ( $exetype == $defaultexetype ) then
      echo "build.csh> '$defaultexetype' is the default and can be omitted"
    endif
  else if ( $opt == intel || $opt == amd || $opt == mips ) then
    if ( $readarchitecture == 1 ) then
      echo "build.csh> ERROR: Two mutually exclusive compiler options 3 were specified"
      echo "build.csh> ERROR: Conflicting arguments are '"$architecture"' and '"$opt"'"
      exit
    endif
    set architecture = $opt
    set readarchitecture = 1
    if ( $architecture == $defaultarchitecture ) then
      echo "build.csh> '$defaultarchitecture' is the default and can be omitted"
    endif
  else if ( $opt == xxlarge || $opt == large || $opt == medium || $opt == small || $opt == reduce ) then
    if ( $readchmsize == 1 ) then
      echo "build.csh> ERROR: Two mutually exclusive CHARMM options were specified"
      echo "build.csh> ERROR: Conflicting arguments are '"$chmsize"' and '"$opt"'"
      exit
    endif
    set readchmsize = 1
    set chmsize = $opt
    if ( $chmsize == $defaultchmsize ) then
      echo "build.csh> '$defaultchmsize' is the default and can be omitted"
    endif
  else if ( $opt == dftb ) then
    if ( $readchmqm == 1 ) then
      echo "build.csh> ERROR: Two mutually exclusive CHARMM QM methods specified"
      exit
    endif
    set readchmqm = 1
    set chmqm = $opt
  else if ( $opt == c31 || $opt == c35 || $opt == c36 ) then
    if ( $readchmtype == 1 ) then
      echo "build.csh> ERROR: Two mutually exclusive CHARMM versions specified"
      exit
    endif
    set readchmtype = 1
    set CTYPE = $opt
  else if ( $opt == last ) then
    echo "build.csh> ERROR: Argument 'last' cannot be accompanied by other arguments"
    exit
  else if ( $opt == quiet ) then
    set quiet = "-s"
  else if ( $opt == profile || $opt == noprofile ) then
    if ( $opt == profile ) then
      set profile = 1
    else
      set profile = 0
    endif
    if ( $profile == $defaultprofile ) then
      echo "build.csh> '$opt' is the default and can be omitted"
    endif
  else if ( $opt == help ) then
    echo
    echo "     OPTIM installation script"
    echo "     Copyright (C) 2003-2008 Semen A. Trygubenko and David J. Wales"
    echo
    echo "         Usage: ./build.csh [arg1] [arg2] ... [argN]"
    echo
    echo "         Arguments can be one or more keywords listed below or a valid target. The order in which the arguments"
    echo "         are given is unimportant. Keywords within the same section are mutually exclusive. Default value is"
    echo "         denoted by <>. All keywords must be given in lowercase. Valid keywords are:"
    echo
    echo "         Section            Keyword             Explanation"
    echo
    echo "         Targets            <optim>             Build OPTIM"
    echo "                            coptim              Build OPTIM and link against CHARMM c35 or c36 source directory"
    echo "                            amhoptim            Build OPTIM with Associative Memory Hamiltonian"
    echo "                            unoptim             Build OPTIM and link against UNRES source directory"
    echo "                            amb9optim           Build OPTIM with the interfaced AMBER 9 and NAB libraries"
    echo "                            jboptim             Build OPTIM and link against Bowman source directory"
    echo "                            dlfoptim            Build OPTIM and link against DL-find source directory"
    echo "                            quipoptim           Build OPTIM and link against QUIP source directory"
    echo "                            atarget             Any other valid target 'atarget' understood by OPTIM makefile"
    echo
    echo "         Compilers          <pgi>               The Portland Group Compiler Technology Fortran 90 compiler"
    echo "                            ifort               Intel(R) Fortran Compiler, Version 8.0 or later"
    echo "                            pathscale           pathscale compiler"
    echo "                            g95                 Open source Fortran 95 compiler"
    echo "                            nag                 NAGWare Fortran 95 compiler"
    echo "                            ibm                 IBM compiler"
    echo
    echo "         Compiler options 1 <opt>               Switch on compiler optimization"
    echo "                            noopt               Switch off compiler optimization"
    echo "                            debug               Switch off compiler optimization, include debug flags"
    echo "                            debugslow           Switch off compiler optimization, include debug flags with slow execution"
    echo
    echo "         Compiler options 2 <dynamic>           Produce dynamic executable"
    echo "                            static              Produce static executable"
    echo "                            nobflag             Do not pass binding flag to the compiler (= use compiler default)"
    echo
    echo "         Compiler options 3 <intel>             Optimize for Intel architecture"
    echo "                            amd                 Optimize for AMD architecture"
    echo "                            mips                Optimize for SiCortex"
    echo
    echo "         Compiler options 4 profile             Profile the executable"
    echo "                            <noprofile>"
    echo
    echo "         Miscellaneous      help                Prints this message"
    echo "                            examples            Prints several examples of usage of this script"
    echo "                            last                Use the same arguments as in the last invocation of this script"
    echo "                            quiet               Passes -s option to make"
    echo
    echo "         CHARMM size        xxlarge, large, <medium>, small, reduce"
    echo
    echo "         CHARMM QM method   dftb                Use the SCC-DFTB routines for QM or QM/MM input"
    echo
    echo "         CHARMM version     <c35>               Use the CHARMM35 source code"
    echo "                            c36                 Use the CHARMM36 source code no eef1 interface"
    echo "                            c31                 Use the CHARMM31 source code"
    echo
    echo "         Some keywords drawn from different sections are incompatible. Certain keyword combinations are not"
    echo "         supported by this script, compiler or OPTIM source. Consult this script's source code for details."
    echo "         Known keyword dependencies are given below:"
    echo
    echo "         amd, intel: opt, pgi"
    echo "         CHARMM size, CHARMM type: coptim"
    echo "         jboptim: ifort"
    echo
    echo "         Following keywords require an argument that must follow immediately after the keyword separated by one"
    echo "         or more spaces:"
    echo
    echo "         Keyword Default_argument          Explanation"
    echo "         optsrc $optsrc          OPTIM source directory"
    echo "         c31src $c31src      CHARMM31 source directory"
    echo "         c35src $c35src      CHARMM35 source directory"
    echo "         c36src $c36src      CHARMM36 source directory"
    echo
    exit
  else if ( $opt == examples ) then
    echo
    echo "     OPTIM installation script"
    echo "     Copyright (C) 2003-2008 Semen A. Trygubenko and David J. Wales"
    echo
    echo "         Examples of usage:"
    echo
    echo "         1. To build coptim with Portland and AMD hardware:"
    echo
    echo "               ./build.csh coptim amd pgi"
    echo
    echo "         2. To clean up the source directory without deleting any of the executable files in OPTIM/bin/ run:"
    echo
    echo "               ./build.csh clean"
    echo
    echo "         3. To install optim on a Windows machine with cygwin using nag execute:"
    echo
    echo "               ./build.csh nag nobflag"
    echo
    echo "         because at the moment of writing this (8 18:39:45 BST 2005) nag compiler for Windows does"
    echo "         not support -B option?"
    exit
  else if ( $opt == cleanc31 ) then
    echo Executing "cd $c31src && ./clean.csh"
    cd $c31src && ./clean.csh
    exit
  else if ( $opt == cleanc35 ) then
    echo Executing "cd $c35src && ./clean.csh"
    cd $c35src && ./clean.csh
    exit
  else if ( $opt == cleanc36 ) then
    echo Executing "cd $c36src && ./clean.csh"
    cd $c36src && ./clean.csh
    exit
  else
    echo "build.csh> '"$opt"' is neither a valid option nor a top-level target - will be passed to make"
    set target = $opt
  endif
end

echo "build.csh> optsrc = $optsrc"

# set compiler name, compiler flags and search path
echo "build.csh> compiler = $compiler"
if ( $compiler == pgi ) then
  if($CTYPE == c36 || $CTYPE == C36) then
    set FC = "mpif90"
  #  set FC = "pgf95"
  else
    set FC = "pgf90"
  endif
  set GENFLAGS = "-Mextend"
  set FREEFORMAT_FLAG="-Mfree"
  set EXTRA_FLAGS="-module"
  if ( $profile == 1 ) then
    set GENFLAGS = "-Mprof=func $GENFLAGS"
  endif
  set NOOPT = "$GENFLAGS -O0 "
  if ( $level == opt ) then
    set EXTRA_CFLAGS="-O3"
    if ( $architecture == intel ) then
      if ( $target == coptim ) then 
         set FFLAGS = "$GENFLAGS -O3 -Munroll -Mnoframe"
      else
         set FFLAGS = "$GENFLAGS -O3 -Munroll -Mscalarsse -Mnoframe -Mvect=sse -Mcache_align -Mflushz "
      endif
      set NOOPT = "$NOOPT "
    else if ( $architecture == amd ) then
      if ( $target == coptim ) then 
         set FFLAGS = "$GENFLAGS -Mextend -O3 -Munroll -Mnoframe"
      else 
# -Mvect=sse is currently incompatible between clust head node and other nodes !!!
         set FFLAGS = "$GENFLAGS -O3 -Munroll -Mscalarsse -Mnoframe -Mcache_align -Mflushz"
      endif
      set NOOPT = "$NOOPT "
    endif
  else if ( $level == noopt ) then
    set FFLAGS = "$NOOPT"
    set EXTRA_CFLAGS="-O0"
  else if ( $level == debug ) then
    set FFLAGS = "$NOOPT -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mpgicoff -traceback "
    set EXTRA_CFLAGS="-g"
  endif
  if ( $exetype == static ) then
     set FFLAGS = "-Bstatic $FFLAGS"
  else if ( $exetype == dynamic ) then
     set FFLAGS = "-Bdynamic $FFLAGS"
  endif
  set NOOPT2 = "$FFLAGS"
  if($CTYPE == c36 || $CTYPE == C36) then
   set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB -Impi"
  else
   set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB "
  endif
  set SWITCH = "pgi"
else if ( $compiler == ifort ) then
  set FC = "ifort"
  set GENFLAGS = "-132  -heap-arrays "
  set FREEFORMAT_FLAG="-free"
  set EXTRA_FLAGS="-I"
  if ( $profile == 1 ) then
    echo 'WARNING: Profile option was NOT tested with ifort compiler!'
    set GENFLAGS = "-prof_gen -prof_file=summary $GENFLAGS"
  endif
  set NOOPT = "$GENFLAGS -O0"
  if ( $level == opt ) then
    # set FFLAGS = "$GENFLAGS -O3 -axW -Vaxlib"
    set FFLAGS = "$GENFLAGS -O3 -Vaxlib"
    set EXTRA_CFLAGS="-O3"
  else if ( $level == noopt ) then
    set FFLAGS = "$NOOPT"
    set EXTRA_CFLAGS="-O0"
  else if ( $level == debug ) then
#   set FFLAGS = "$NOOPT -g -C -traceback -debug full -check uninit"
    set FFLAGS = "$NOOPT -g -C -traceback -debug full -check uninit -stand f03 -assume realloc_lhs -check all,noarg_temp_created -traceback -warn all -fstack-protector -assume protect_parens -implicitnone"
    set EXTRA_CFLAGS="-g"
  endif
  if ( $exetype == static ) then
     set FFLAGS = "-static $FFLAGS"
  endif
  set NOOPT2 = "$FFLAGS"
  set SEARCH_PATH = "-INEB -ICONNECT -IAMH -I.. -I../NEB"
  set SWITCH = "ifort"
else if ( $compiler == nag ) then
  set FC = "nagfor"
  set GENFLAGS = "-132 -maxcontin=3000 -kind=byte -ieee=full "
  set FREEFORMAT_FLAG="-free"
  set EXTRA_FLAGS="-I"
  if ( $profile == 1 ) then
    echo 'WARNING: Profile option was NOT tested with nag compiler!'
    set GENFLAGS = "-pg $GENFLAGS"
  endif
  set NOOPT = "$GENFLAGS -O0"
  if ( $level == opt ) then
#   set FFLAGS = "$GENFLAGS -O4 -mismatch_all "
    set FFLAGS = "$GENFLAGS -O4 "
    set EXTRA_CFLAGS="-O3"
  else if ( $level == noopt ) then
#   set FFLAGS = "$NOOPT -mismatch_all"
    set FFLAGS = "$NOOPT "
    set EXTRA_CFLAGS="-O0"
  else if ( $level == debug ) then
    set GENFLAGS = "-132 -maxcontin=3000 -kind=byte "
    set NOOPT = "$GENFLAGS -O0"
#   set FFLAGS = "$NOOPT -g -C -mismatch_all -ieee=stop "
    set FFLAGS = "$NOOPT -g -C -ieee=stop "
    set EXTRA_CFLAGS="-O0 -g"
  else if ( $level == debugslow ) then 
    set GENFLAGS = "-132 -maxcontin=3000 -kind=byte "
    set NOOPT = "$GENFLAGS -O0"
#   set FFLAGS = "$NOOPT -g -C=all -mtrace=all -gline -mismatch_all -ieee=stop "
    set FFLAGS = "$NOOPT -g -C=all -mtrace=all -gline -ieee=stop "
    set EXTRA_CFLAGS="-O0 -g"
  endif
  if ( $exetype == static ) then
     set FFLAGS = "-Bstatic $FFLAGS"
  else if ( $exetype == dynamic ) then
     set FFLAGS = "-Bdynamic $FFLAGS"
  endif
  set NOOPT2 = "$GENFLAGS -O1 -mismatch_all "
  set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB"
  set SWITCH = "nag"
else if ( $compiler == pathscale ) then
  set FC = "pathf95"
  set GENFLAGS = "-extend-source"
  set FREEFORMAT_FLAG="-freeform"
  set EXTRA_FLAGS="-I"
  set EXTRA_CFLAGS="-G0"
  if ( $profile == 1 ) then
    echo 'WARNING: Profile option was NOT tested with pathscale compiler!'
    set GENFLAGS = "-pg $GENFLAGS"
  endif
  set NOOPT = "$GENFLAGS -O0 -G0"
  if ( $level == opt ) then
    set FFLAGS = "$GENFLAGS -O3 -G0"
    set EXTRA_CFLAGS="-G0 -O3"
  else if ( $level == noopt ) then
    set FFLAGS = "$NOOPT"
    set EXTRA_CFLAGS="-G0 -O0"
  else if ( $level == debug ) then
    set FFLAGS = "$NOOPT -g -C -G0"
    set EXTRA_CFLAGS="-G0 -g"
  endif
  if ( $exetype == static ) then
     set FFLAGS = "-Bstatic $FFLAGS"
  else if ( $exetype == dynamic ) then
     set FFLAGS = "$FFLAGS"
  endif
  set NOOPT2 = "$FFLAGS"
  set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB "
  set SWITCH = "pathscale"
else if ( $compiler == gfortran ) then
  set FC = "gfortran"
#  set GENFLAGS = "-ffixed-line-length-132"
  set GENFLAGS = "-ffixed-line-length-132 -ffree-line-length-none"
  set NOOPT = "$GENFLAGS -O0"
  set FREEFORMAT_FLAG="-free "
  set EXTRA_FLAGS=" "
  if ( $level == opt ) then
    set FFLAGS = "$GENFLAGS -O2"
    set EXTRA_CFLAGS="-O2"
  else if ( $level == noopt ) then
    set FFLAGS = "$NOOPT"
    set EXTRA_CFLAGS="-O0"
  else if ( $level == debug ) then
    set FFLAGS = "$NOOPT -g -C -fbounds-check -Wuninitialized -O -ftrapv"
#    set FFLAGS = "$NOOPT -g -C -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic"
    set EXTRA_CFLAGS="-g"
  endif
  if ( $exetype == static ) then
     set FFLAGS = "-Bstatic $FFLAGS"
  else if ( $exetype == dynamic ) then
     set FFLAGS = "-Bdynamic $FFLAGS"
  endif
  set NOOPT2 = "$FFLAGS"
  set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB"
  set SWITCH = "gfortran"
else if ( $compiler == ibm ) then
  echo "build.csh> This is untested - please correct the compiler flags below if needed and let us know about the changes!"
  set FC = "xlf90"
  set GENFLAGS = "-qfixed=132"
  set FREEFORMAT_FLAG=" "
  set EXTRA_FLAGS=" "
  if ( $profile == 1 ) then
    echo 'WARNING: Profile option was NOT tested with ibm compiler!'
  endif
  set NOOPT = "$GENFLAGS -O0"
  if ( $level == opt ) then
    set FFLAGS = "$GENFLAGS -O3 -qstrict -qxlf77=leadzero -qmaxmem=8192 -qarch=pwr3  -qtune=pwr3 -qnolm"
  else if ( $level == noopt ) then
    set FFLAGS = "$NOOPT"
  else if ( $level == debug ) then
    set FFLAGS = "$NOOPT -g -C"
  endif 
  if ( $exetype == static ) then
     set FFLAGS = "-bstatic $FFLAGS"
  endif
  set NOOPT2 = "$FFLAGS"
  set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB"
  set SWITCH = "ibm"
else if ( $compiler == g95 ) then
  set EXTRA_CFLAGS="-g"
  set FC = "gfortran"
  set GENFLAGS = "-ffixed-line-length-132"
  set FREEFORMAT_FLAG="-ffree-form"
  # set EXTRA_FLAGS="-fmod="
  # ss2029 
  echo 'changing EXTRA_FLAGS from -fmod to -fintrinsic-modules-path'
  set EXTRA_FLAGS="-fintrinsic-modules-path "
  if ( $profile == 1 ) then
    echo 'WARNING: Profile option was NOT tested with g95 compiler!'
  endif
  set NOOPT = "$GENFLAGS -O0"
  if ( $level == opt ) then
    set FFLAGS = "$GENFLAGS -O3"
  else if ( $level == noopt ) then
    set FFLAGS = "$NOOPT"
  else if ( $level == debug ) then
    set FFLAGS = "$NOOPT -g -C"
  endif 
  if ( $exetype == static ) then
     echo 'Unsupported'
     exit
     set FFLAGS = "-bstatic $FFLAGS"
  endif
  set NOOPT2 = "$FFLAGS"
  set SEARCH_PATH = "-INEB -ICONNECT -I.. -I../NEB"
  set SWITCH = "gfortran"
endif

# print compiler flags info
if ( $level == opt ) then
  echo 'build.csh> compiler option 1 = opt'
else if ( $level == noopt ) then
  echo 'build.csh> compiler option 1 = noopt'
else if ( $level == debug ) then
  echo 'build.csh> compiler option 1 = debug'
else if ( $level == debugslow ) then
  echo 'build.csh> compiler option 1 = debugslow'
endif
if ( $exetype == static ) then
  echo 'build.csh> exetype = static'
else if ( $exetype == dynamic ) then
  echo 'build.csh> exetype = dynamic'
else
  echo 'build.csh> exetype = compiler default'
endif

# get OPTIM version number
#if (! -e VERSION ) then
#  echo "build.csh> ERROR: File 'VERSION' was not found - corrupted source directory or invalid path"
#  exit
#else
#  set version = `cat VERSION`
#  setenv VERSION $version
#  echo "build.csh> version = $version"
#endif

# generate the name of the executable and appropriate for it make options
if (! -e ../bin ) mkdir ../bin
if (! -e ../bin/$compiler ) mkdir ../bin/$compiler
if ( $target == optim ) then
  set EXE = "../bin/$compiler/$programname"
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == amb9optim ) then
  set EXE = "../bin/$compiler/A9$programname"
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set NABHOME = "../../NAB"
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == quipoptim ) then
  set EXE = "../bin/$compiler/QUIP$programname"
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == coptim ) then
  if ( $CTYPE == c31 ) then
    set EXE = "../bin/$compiler/C31$programname"
    set CTYPE = C31
    set chsrc = $c31src
    echo "build.csh> chsrc reset to " $c31src
    set FCMDIR = $c31src"/source/fcm"
    echo "build.csh> Assuming the platform is 'gnu' when compiling c31"
    set PREFLX = $c31src"/tool/prefx_gnu"
    set PREFDIR = $c31src"/build/gnu"
    set LIBDIRCH = $c31src"/lib/gnu"
  else if ( $CTYPE == c36 ) then 
    set EXE = "../bin/$compiler/C36$programname"
    set CTYPE = C36
    set chsrc = $c36src
    echo "build.csh> chsrc reset to " $c36src
    set FCMDIR = $c36src"/source/fcm"
    echo "build.csh> Assuming the platform is 'gnu' when compiling c36"
    set PREFLX = $c36src"/tool/prefx_gnu"
    set PREFDIR = $c36src"/build/gnu"
    set LIBDIRCH = $c36src"/lib/gnu"
  else
    set EXE = "../bin/$compiler/C35$programname"
    set CTYPE = C35
    set chsrc = $c35src
    echo "build.csh> chsrc reset to " $c35src
    set FCMDIR = $c35src"/source/fcm"
    echo "build.csh> Assuming the platform is 'gnu' when compiling c35"
    set PREFLX = $c35src"/tool/prefx_gnu"
    set PREFDIR = $c35src"/build/gnu"
    set LIBDIRCH = $c35src"/lib/gnu"
  endif
  set BLAS_EXCLUDE_LIST = "dnrm2.o daxpy.o dcopy.o ddot.o"
  set SRCCH = "charmm_main.src energy.src"
else if ( $target == unoptim ) then
  set EXE = "../bin/$compiler/UN$programname"
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == amhoptim ) then
  set EXE = "../bin/$compiler/AMH$programname"
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == dlfoptim ) then
  if ( $compiler != ifort && $compiler != pgi && $compiler != g95 && $compiler != gfortran) then
    echo "build.csh> ERROR: can only compile DLFWrapper with ifort, pgi, g95, and gfortran"
	 exit
  endif
  set EXE = "../bin/$compiler/DLF$programname."
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == jboptim ) then
  if ( $compiler != ifort ) then
    echo "build.csh> ERROR: can only compile Bowman potential with ifort"
	 exit
  endif
  set EXE = "../bin/$compiler/JB$programname."
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set chsrc = ""
else if ( $target == clean || $target == cleanall ) then
  set EXE = $target
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = "charmm_main.src energy.src"
  set LIBDIRCH = ""
  set chsrc = ""
else
  set EXE = $target
  set BLAS_EXCLUDE_LIST = ""
  set CTYPE = ""
  set PREFLX = ""
  set FCMDIR = ""
  set PREFDIR = ""
  set SRCCH = ""
  set LIBDIRCH = ""
  set exemessage = 0
  set chsrc = ""
endif

if ( $target == jboptim ) then
	set JBMODS = "-IBowman"
else
	set JBMODS = "-IBowman/fake"
endif

echo "build.csh> target = $target"

set make = "make "

# need to determine if we could use some of the previous build
if ( "$arguments" == "$lastargs" ) then # full compatability of current build and the previous one
else
  # determine last target and compiler
  set LastTarget = $defaulttarget
  set LastCompiler = $defaultcompiler
  set LastQM = $defaultchmqm
  set LastChtype = ""
  foreach opt ( $lastargs )
    if ( $opt == optim || $opt == coptim || $opt == unoptim || $opt == jboptim ) then
      set LastTarget = $opt
    endif
    if ( $opt == pgi || $opt == ifort || $opt == nag || $opt == ibm || $opt == g95 || $opt == pathscale ) then
      set LastCompiler = $opt
    endif
  # determine last QM method used 
    if ( $opt == dftb ) then
      set LastQM = $opt 
    endif
    if ( $opt == c31 || $opt == c35 || $opt == c36 ) then
      set LastChtype = $opt
    endif
  end
  # compare current and previous targets
  if ( "$LastTarget" == "$target" ) then
  else
     if ( "$LastTarget" == coptim && "$target" != coptim || "$LastTarget" != coptim && "$target" == coptim ) then
       echo "build.csh> Need to relink BLAS library because BLAS_EXCLUDE_LIST has changed"
       cd ../../BLAS && make cleanlib && cd ../OPTIM/source
     endif
  endif
  # compare current and previous compilers
  if ( "$LastCompiler" == "$compiler" ) then
  else # we need to rebuild from scratch
     echo "build.csh> Last target was built with a different compiler - rebuilding from scratch"
     $make clean
  endif
  if ( $target == coptim ) then
     # compare CHARMM QM methods
     if ( ("$LastQM" == dftb && "$chmqm" == dftb) || ("$LastQM" == none && "$chmqm" == none) ) then
     else
       echo "build.csh> Last target used a different CHARMM QM method - rebuilding CHARMM"
       echo "build.csh> cd ${chsrc} && ./clean.csh && cd $optsrc"
       cd ${chsrc} && ${chsrc}/clean.csh && cd $optsrc
       echo "build.csh> cd ${optsrc}/CHARMMinterface && make clean && cd $optsrc"
       cd ${optsrc}/CHARMMinterface && make clean && cd $optsrc
     endif
     # compare CHARMM versions
     if ( `expr substr "$LastChtype" 2 2` == `expr substr "$CTYPE" 2 2` ) then
     else
       echo "build.csh> Cleaning the CHARMMinterface directory because the CHARMM version may have changed"
       cd CHARMMinterface && make clean && cd ../
     endif
  endif
endif

# all set up and ready to run make
set success = 0
if ( $target == coptim ) then # need to compile CHARMM first
   if (! -e $chsrc ) then
     echo "build.csh> ERROR: CHARMM source directory $chsrc does not exist!"
     exit
   else
     echo "build.csh> chsrc = $chsrc"
   endif
   set command = "install.com gnu"
   set command = "$command $chmsize keepo keepf"
   if ( $chmqm == dftb ) then
     if ( $CTYPE == C31 || $CTYPE == C36 ) then
       echo "build.csh> ERROR: DFTB and CHARMM31 are incompatible!"
       echo "build.csh> ERROR: DFTB and CHARMM36 not interfaced yet!"
       exit
     endif
     set command = "$command T"
# QMTYPE variable passed to Makefile to include SCC-DFTB libraries 
     set QMTYPE = DFTB
   endif
   if ( $compiler == pgi ) then
     if($CTYPE == C36) then #compile.csh adjustments
     	 #export MPI_INCLUDE=/usr/local/shared/ubuntu-12.04/x86_64/openmpi-1.4.5-pgi12/include
		 #export MPI_LIB=/usr/local/shared/ubuntu-12.04/x86_64/openmpi-1.4.5-pgi12/lib	
	     #set command = "install.com gnu xxlarge PGF95 keepo keepf" 
	     #set command = "install.com gnu xxlarge PGF95 keepo keepf E MPICH"  #mpif90
         set command = "compile.csh"
     else 
	     set command = "$command PGF90"
	     if ( $level == opt ) then
	       set command = "$command OPT"
	     endif
	     if ( $architecture == amd ) then
	       set command = "$command AMD"
	     endif
	 endif
   else if ( $compiler == ifort ) then
     set command = "$command ifort"
     if($CTYPE == C36) then
     echo "build.csh> ERROR: CHARMM36 only with PGI95 so far"
     endif
   else
     echo "build.csh> ERROR: This script cannot compile CHARMM directory with the selected compiler!"
     exit
   endif
#   echo "build.csh> cd ${chsrc} && ${chsrc}/$command && cd $optsrc"
   echo "build.csh> cd ${chsrc} && ./$command && cd $optsrc"
   echo
#   cd ${chsrc} && ${chsrc}/$command && cd $optsrc
   cd ${chsrc} && ./$command && cd $optsrc
   echo "build.csh> tail of" $PREFDIR/gnu.log
   tail -2 $PREFDIR/gnu.log
   echo " "
else if ( $target == amhoptim ) then
   if (! -e $amhsrc ) then
     echo "build.csh> ERROR: AMH source directory $amhsrc does not exist!"
     exit
   else
     echo "build.csh> amhsrc = $amhsrc"
   endif
#   echo "build.csh> cd ${amhsrc} && make  && cd $optsrc"
#   cd ${amhsrc} && make  && cd $optsrc
   echo "  "
endif
$make $quiet $EXE FC=$FC NOOPT="$NOOPT" NOOPT2="$NOOPT2" SEARCH_PATH="$SEARCH_PATH $JBMODS" FFLAGS="$FFLAGS" FREEFORMAT_FLAG="$FREEFORMAT_FLAG" EXTRA_FLAGS="$EXTRA_FLAGS" EXTRA_CFLAGS="$EXTRA_CFLAGS" $target="$EXE" LIBDIRCH="$LIBDIRCH" \
BLAS_EXCLUDE_LIST="$BLAS_EXCLUDE_LIST" CTYPE="$CTYPE" QMTYPE="$QMTYPE" PREFLX="$PREFLX" FCMDIR="$FCMDIR" PREFDIR="$PREFDIR" SRCCH="$SRCCH" \
CHSRC="$chsrc" SWITCH="$SWITCH" >& make.out && set success = 1

# save arguments that were passed to this script into lastargs file
if ( "$argv" == 'last' ) then
else if ( $target == 'clean' || $target == 'cleanall' ) then
else
  echo "$argv" > lastargs
endif

if ( $success == 1 ) then
  if ( $exemessage == 1 ) then
    echo "build.csh> OPTIM executable = $optsrc/$EXE"
  endif
  echo 'build.csh> Success. See make.out for detailed compiler output'
else
  echo 'build.csh> Failure - refer to make.out'
  echo 'build.csh> grep for error messages in make.out gives:'
  grep -i error make.out
endif
