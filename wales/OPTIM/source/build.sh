#!/bin/bash
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
# For syntax checking
# set -u
# set -x

# set parameters
defaulttarget=optim
defaultcompiler=ifort 
defaultlevel=opt
defaultexetype=dynamic
defaultarchitecture=intel
defaultchmsize=medium
defaultprofile=0
defaultchmqm=none
programname=OPTIM
EmptyString=""
SingleQuote="'"
DoubleQuote='"'

# set defaults for variables
USER=`echo $USER`
MACHINE=`uname -n`
OS=`uname -s`
QMTYPE=""
CTYPE=c35
# default CHARMM directories
c35src=$PWD"/../../CHARMM35"
c31src=$PWD"/../../CHARMM31"
#set amhsrc = "/home/$USER/svn/OPTIM/source/AMH"
amhsrc="AMH"
target=$defaulttarget
compiler=$defaultcompiler
level=$defaultlevel
exetype=$defaultexetype
architecture=$defaultarchitecture
chmsize=$defaultchmsize
chmqm=$defaultchmqm
optsrc=`pwd`
quiet=$EmptyString
profile=defaultprofile
exemessage=1

# process keyword-based arguments
readnext=0
for opt in $@; do 
  if [ "$readnext" = "1" ]; then
    readnext=0
    if [ "$optsrc" = "FILL" ]; then
      optsrc=$opt
    elif [ "$c31src" = "FILL" ]; then
      c31src=$opt
    elif [ "$c35src" = "FILL" ]; then
      c35src=$opt
    fi
  fi
  if [ "$opt" = "optsrc" ]; then
    readnext=1
    optsrc="FILL"
  elif [ "$opt" = "c31src" ]; then
    readnext=1
    c31src="FILL"
  elif [ "$opt" = "c35src" ]; then
    readnext=1
    c35src="FILL"
  fi
done

# make sure we are in the source directory. After this point all the references can be relative to the source dir
cd $optsrc

# read in the last arguments list into lastargs; if unavailable, initialise lastargs to an empty string
if [ -e lastargs ]; then
  lastargs=`cat lastargs`
else
  lastargs=$EmptyString
fi

# check if we can make use of the lastargs file
if [[ ! -e 'lastargs' && "$@" = "$EmptyString" ]]; then # we'd better display help message
  arguments=help
elif [[ "$@" = "last" ]]; then # file 'lastargs' exists and we are instructed to read it
  echo "build.sh> lastargs = $lastargs"
  arguments="$lastargs"
else
  arguments="$@"
fi

# process the command line input
readcompiler=0
readlevel=0
readtarget=0
readexetype=0
readarchitecture=0
readchmsize=0
readchmqm=0
readchmtype=0
readnext=1
for opt in $arguments; do 
  # keyword arguments and their values need to be ignored
  if [ "$readnext" != 1 ]; then
    readnext=1
  elif [[ "$opt" = "optsrc" || "$opt" = "c31src" || "$opt" = "c35src" ]]; then
    readnext=0
  elif [[ "$opt" = "opt" || "$opt" = "debug" || "$opt" = "debugslow" || "$opt" = "noopt" ]]; then
    if [ "$readlevel" = 1 ]; then
      echo "build.sh> ERROR: Two or more mutually exclusive compiler options 1 were specified"
      echo "build.sh> ERROR: First two conflicting arguments are '$level' and '$opt'"
      exit
    fi
    level=$opt
    if [ "$level" = "$defaultlevel" ]; then
      echo "build.sh> '$defaultlevel' is the default and can be omitted"
    fi
    readlevel=1
  elif [[ "$opt" = "optim" || "$opt" = "coptim" || "$opt" = "unoptim" || "$opt" = "amb9optim" || "$opt" = "amhoptim" || "$opt" = "dlfoptim" || "$opt" = "jboptim" || "$opt" = "quipoptim" ]]; then
    if [ "$readtarget" = 1 ]; then
      echo "build.sh> ERROR: Two or more mutually exclusive targets were specified"
      echo "build.sh> ERROR: First two conflicting arguments are '$target' and '$opt'"
      exit
    fi
    target=$opt
    if [ "$target" = "$defaulttarget" ]; then
      echo "build.sh> '$defaulttarget' is the default and can be omitted"
    fi
    readtarget=1
  elif [[ "$opt" = "pgi" || "$opt" = "ifort" || "$opt" = "nag" || "$opt" = "ibm" || "$opt" = "g95" || "$opt" = "gfortran" || "$opt" = "pathscale" ]]; then
    if [ "$readcompiler" = 1 ]; then
      echo "build.sh> ERROR: Two or more mutually exclusive compiler types were specified"
      echo "build.sh> ERROR: First two conflicting arguments are '$compiler' and '$opt'"
      exit
    fi
    compiler=$opt
    if [ "$compiler" = "$defaultcompiler" ]; then
      echo "build.sh> '$defaultcompiler' is the default and can be omitted"
    fi
    readcompiler=1
  elif [[ "$opt" = "dynamic" || "$opt" = "static" || "$opt" = "nobflag" ]]; then
    if [ "$readexetype" = 1 ]; then
      echo "build.sh> ERROR: Two mutually exclusive compiler options 2 were specified"
      echo "build.sh> ERROR: Conflicting arguments are '$exetype' and '$opt'"
      exit
    fi
    exetype=$opt
    readexetype=1
    if [ "$exetype" = "$defaultexetype" ]; then
      echo "build.sh> '$defaultexetype' is the default and can be omitted"
    fi
  elif [[ "$opt" = "intel" || "$opt" = "amd" || "$opt" = "mips" ]]; then
    if [ "$readarchitecture" = 1 ]; then
      echo "build.sh> ERROR: Two mutually exclusive compiler options 3 were specified"
      echo "build.sh> ERROR: Conflicting arguments are '$architecture' and '$opt'"
      exit
    fi
    architecture=$opt
    readarchitecture=1
    if [ "$architecture" = "$defaultarchitecture" ]; then
      echo "build.sh> '$defaultarchitecture' is the default and can be omitted"
    fi
  elif [[ "$opt" = "xxlarge" || "$opt" = "large" || "$opt" = "medium" || "$opt" = "small" || "$opt" = "reduce" ]]; then
    if [ "$readchmsize" = "1" ]; then
      echo "build.csh> ERROR: Two mutually exclusive CHARMM options were specified"
      echo "build.csh> ERROR: Conflicting arguments are '$chmsize' and '$opt'"
      exit
    fi
    readchmsize=1
    chmsize=$opt
    if [ "$chmsize" = "$defaultchmsize" ]; then
      echo "build.sh> '$defaultchmsize' is the default and can be omitted"
    fi
  elif [ "$opt" = "dftb" ]; then
    if [ "$readchmqm" = 1 ]; then
      echo "build.sh> ERROR: Two mutually exclusive CHARMM QM methods specified"
      exit
    fi
    readchmqm=1
    chmqm=$opt
  elif [[ "$opt" = "c31" || "$opt" = "c35" ]]; then
    if [ "$readchmtype" = "1" ]; then
      echo "build.sh> ERROR: Two mutually exclusive CHARMM versions specified"
      exit
    fi
    readchmtype=1
    CTYPE=$opt
  elif [ "$opt" = "last" ]; then
    echo "build.sh> ERROR: Argument 'last' cannot be accompanied by other arguments"
    exit
  elif [ "$opt" = "quiet" ]; then
    quiet="-s"
  elif [[ "$opt" = "profile" || "$opt" = "noprofile" ]]; then
    if [ "$opt" = "profile" ]; then
      profile=1
    else
      profile=0
    fi
    if [ "$profile" = "$defaultprofile" ]; then
      echo "build.sh> '$opt' is the default and can be omitted"
    fi
  elif [ "$opt" = "help" ]; then
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
    echo "                            coptim              Build OPTIM and link against CHARMM c35 source directory"
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
    echo
    exit
  elif [ "$opt" = "examples" ]; then
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
  elif [ "$opt" = "cleanc31" ]; then
    echo 'Executing "cd $c31src && ./clean.csh"'
    cd $c31src && ./clean.csh
    exit
  elif [ "$opt" = "cleanc35" ]; then
    echo 'Executing "cd $c35src && ./clean.csh"'
    cd $c35src && ./clean.csh
    exit
  else
    echo "build.sh> '$opt' is neither a valid option nor a top-level target - will be passed to make"
    target=$opt
  fi
done

echo "build.sh> optsrc = $optsrc"

# set compiler name, compiler flags and search path
echo "build.sh> compiler = $compiler"
if [ "$compiler" = "pgi" ]; then
  FC="pgf90"
  GENFLAGS="-Mextend"
  FREEFORMAT_FLAG="-Mfree"
  EXTRA_FLAGS="-module"
  if [ "$profile" = 1 ]; then
    GENFLAGS="-Mprof=func $GENFLAGS"
  fi
  NOOPT="$GENFLAGS -O0 "
  if [ "$level" = "opt" ]; then
    EXTRA_CFLAGS="-O3"
    if [ "$architecture" = "intel" ]; then
      if [ "$target" = "coptim" ]; then 
        FFLAGS="$GENFLAGS -O3 -Munroll -Mnoframe"
      else
        FFLAGS="$GENFLAGS -O3 -Munroll -Mscalarsse -Mnoframe -Mvect=sse -Mcache_align -Mflushz "
      fi
      NOOPT="$NOOPT "
    elif [ "$architecture" = "amd" ]; then
      if [ "$target" = "coptim" ]; then 
        FFLAGS="$GENFLAGS -Mextend -O3 -Munroll -Mnoframe"
      else 
# -Mvect=sse is currently incompatible between clust head node and other nodes !!!
        FFLAGS="$GENFLAGS -O3 -Munroll -Mscalarsse -Mnoframe -Mcache_align -Mflushz"
      fi
      NOOPT="$NOOPT "
    fi
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT"
    EXTRA_CFLAGS="-O0"
  elif [ "$level" = "debug" ]; then
    FFLAGS="$NOOPT -C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mpgicoff -traceback "
    EXTRA_CFLAGS="-g"
  fi
  if [ "$exetype" = "static" ]; then
     FFLAGS="-Bstatic $FFLAGS"
  elif [ "$exetype" = "dynamic" ]; then
     FFLAGS="-Bdynamic $FFLAGS"
  fi
  NOOPT2="$FFLAGS"
  SEARCH_PATH="-INEB -ICONNECT -I.. -I../NEB "
  SWITCH="pgi"
elif [ "$compiler" = "ifort" ]; then
  FC="ifort"
  GENFLAGS="-132  -heap-arrays "
  FREEFORMAT_FLAG="-free"
  EXTRA_FLAGS="-I"
  if [ "$profile" = 1 ]; then
    echo 'WARNING: Profile option was NOT tested with ifort compiler!'
    GENFLAGS="-prof_gen -prof_file=summary $GENFLAGS"
  fi
  NOOPT="$GENFLAGS -O0"
  if [ "$level" = "opt" ]; then
    # set FFLAGS = "$GENFLAGS -O3 -axW -Vaxlib"
    FFLAGS="$GENFLAGS -O3 -Vaxlib"
    EXTRA_CFLAGS="-O3"
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT"
    EXTRA_CFLAGS="-O0"
  elif [ "$level" = "debug" ]; then
#   set FFLAGS = "$NOOPT -g -C -traceback -debug full -check uninit"
    FFLAGS="$NOOPT -g -C -traceback -debug full -check uninit -stand f03 -assume realloc_lhs -check all,noarg_temp_created -traceback -warn all -fstack-protector -assume protect_parens -implicitnone"
    EXTRA_CFLAGS="-g"
  fi 
  if [ "$exetype" = "static" ]; then
    FFLAGS="-static $FFLAGS"
  fi
  NOOPT2="$FFLAGS"
  SEARCH_PATH="-INEB -ICONNECT -IAMH -I.. -I../NEB"
  SWITCH="ifort"
elif [ "$compiler" = "nag" ]; then
  FC="nagfor"
  GENFLAGS="-132 -maxcontin=3000 -kind=byte -ieee=full "
  FREEFORMAT_FLAG="-free"
  EXTRA_FLAGS="-I"
  if [ "$profile" = 1 ]; then
    echo 'WARNING: Profile option was NOT tested with nag compiler!'
    GENFLAGS="-pg $GENFLAGS"
  fi
  NOOPT="$GENFLAGS -O0"
  if [ "$level" = "opt" ]; then
    FFLAGS="$GENFLAGS -O4 -mismatch_all "
    EXTRA_CFLAGS="-O3"
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT -mismatch_all"
    EXTRA_CFLAGS="-O0"
  elif [ "$level" = "debug" ]; then
    GENFLAGS="-132 -maxcontin=3000 -kind=byte "
    NOOPT="$GENFLAGS -O0"
    FFLAGS="$NOOPT -g -C -mismatch_all -ieee=stop "
    EXTRA_CFLAGS="-O0 -g"
  elif [ "$level" = "debugslow" ]; then 
    GENFLAGS="-132 -maxcontin=3000 -kind=byte "
    NOOPT="$GENFLAGS -O0"
    FFLAGS="$NOOPT -g -C=all -mtrace=all -gline -mismatch_all -ieee=stop "
    EXTRA_CFLAGS="-O0 -g"
  fi
  if [ "$exetype" = "static" ]; then
     FFLAGS="-Bstatic $FFLAGS"
  elif [ "$exetype" = "dynamic" ]; then
     FFLAGS="-Bdynamic $FFLAGS"
  fi
  NOOPT2="$GENFLAGS -O1 -mismatch_all "
  SEARCH_PATH="-INEB -ICONNECT -I.. -I../NEB"
  SWITCH="nag"
elif [ "$compiler" = "pathscale" ]; then
  FC="pathf95"
  GENFLAGS="-extend-source"
  FREEFORMAT_FLAG="-freeform"
  EXTRA_FLAGS="-I"
  EXTRA_CFLAGS="-G0"
  if [ "$profile" = 1 ]; then
    echo 'WARNING: Profile option was NOT tested with pathscale compiler!'
    GENFLAGS="-pg $GENFLAGS"
  fi
  NOOPT="$GENFLAGS -O0 -G0"
  if [ "$level" = "opt" ]; then
    FFLAGS="$GENFLAGS -O3 -G0"
    EXTRA_CFLAGS="-G0 -O3"
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT"
    EXTRA_CFLAGS="-G0 -O0"
  elif [ "$level" = "debug" ]; then
    FFLAGS="$NOOPT -g -C -G0"
    EXTRA_CFLAGS="-G0 -g"
  fi
  if [ "$exetype" = "static" ]; then
    FFLAGS="-Bstatic $FFLAGS"
  elif [ "$exetype" = "dynamic" ]; then
    FFLAGS="$FFLAGS"
  fi
  NOOPT2="$FFLAGS"
  SEARCH_PATH="-INEB -ICONNECT -I.. -I../NEB "
  SWITCH="pathscale"
elif [ "$compiler" = "gfortran" ]; then
  FC="gfortran"
#  set GENFLAGS = "-ffixed-line-length-132"
  GENFLAGS="-ffixed-line-length-132 -ffree-line-length-none"
  NOOPT="$GENFLAGS -O0"
  FREEFORMAT_FLAG="-free "
  EXTRA_FLAGS=" "
  if [ "$level" = "opt" ]; then
    FFLAGS="$GENFLAGS -O2"
    EXTRA_CFLAGS="-O2"
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT"
    EXTRA_CFLAGS="-O0"
  elif [ "$level" = "debug" ]; then
    FFLAGS="$NOOPT -g -C -fbounds-check -Wuninitialized -O -ftrapv"
#    set FFLAGS = "$NOOPT -g -C -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic"
    EXTRA_CFLAGS="-g"
  fi
  if [ "$exetype" = "static" ]; then
    FFLAGS="-Bstatic $FFLAGS"
  elif [ "$exetype" = "dynamic" ]; then
    FFLAGS="-Bdynamic $FFLAGS"
  fi
  NOOPT2="$FFLAGS"
  SEARCH_PATH="-INEB -ICONNECT -I.. -I../NEB"
  SWITCH="gfortran"
elif [ "$compiler" = "ibm" ]; then
  echo "build.sh> This is untested - please correct the compiler flags below if needed and let us know about the changes!"
  FC="xlf90"
  GENFLAGS="-qfixed=132"
  FREEFORMAT_FLAG=" "
  EXTRA_FLAGS=" "
  if [ "$profile" = 1 ]; then
    echo 'WARNING: Profile option was NOT tested with ibm compiler!'
  fi
  NOOPT="$GENFLAGS -O0"
  if [ "$level" = "opt" ]; then
    FFLAGS="$GENFLAGS -O3 -qstrict -qxlf77=leadzero -qmaxmem=8192 -qarch=pwr3  -qtune=pwr3 -qnolm"
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT"
  elif [ "$level" = "debug" ]; then
    FFLAGS="$NOOPT -g -C"
  fi 
  if [ "$exetype" = "static" ]; then
    FFLAGS="-bstatic $FFLAGS"
  fi
  NOOPT2="$FFLAGS"
  SEARCH_PATH="-INEB -ICONNECT -I.. -I../NEB"
  SWITCH="ibm"
elif [ "$compiler" = "g95" ]; then
  EXTRA_CFLAGS="-g"
  FC="gfortran"
  GENFLAGS="-ffixed-line-length-132"
  FREEFORMAT_FLAG="-ffree-form"
  # set EXTRA_FLAGS="-fmod="
  # ss2029 
  echo 'changing EXTRA_FLAGS from -fmod to -fintrinsic-modules-path'
  EXTRA_FLAGS="-fintrinsic-modules-path "
  if [ "$profile" = 1 ]; then
    echo 'WARNING: Profile option was NOT tested with g95 compiler!'
  fi
  NOOPT="$GENFLAGS -O0"
  if [ "$level" = "opt" ]; then
    FFLAGS="$GENFLAGS -O3"
  elif [ "$level" = "noopt" ]; then
    FFLAGS="$NOOPT"
  elif [ "$level" = "debug" ]; then
    FFLAGS="$NOOPT -g -C"
  fi 
  if [ "$exetype" = "static" ]; then
    echo 'Unsupported'
    exit
    FFLAGS="-bstatic $FFLAGS"
  fi
  NOOPT2="$FFLAGS"
  SEARCH_PATH="-INEB -ICONNECT -I.. -I../NEB"
  SWITCH="gfortran"
fi

# print compiler flags info
echo "build.sh> compiler option 1 = $level"
if [[ "$exetype" = "static" || "$exetype" = "dynamic" ]]; then
  echo "build.sh> exetype = $exetype"
else
  echo "build.sh> exetype = compiler default"
fi

# get OPTIM version number
#if [ ! -e VERSION ]; then
#  echo "build.sh> ERROR: File 'VERSION' was not found - corrupted source directory or invalid path"
#  exit
#else
#  version=`cat VERSION`
#  export VERSION=$version
#  echo "build.sh> version = $version"
#fi

# generate the name of the executable and appropriate for it make options
if [ ! -e ../bin ]; then
  mkdir ../bin
fi
if [ ! -e ../bin/$compiler ]; then
  mkdir ../bin/$compiler
fi
if [ "$target" = "optim" ]; then
  EXE="../bin/$compiler/$programname"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [ "$target" = "amb9optim" ]; then
  EXE="../bin/$compiler/A9$programname"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  NABHOME="../../NAB"
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [ "$target" = "quipoptim" ]; then
  EXE="../bin/$compiler/QUIP$programname"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [ "$target" = "coptim" ]; then
  if [ "$CTYPE" = c31 ]; then
    EXE="../bin/$compiler/C31$programname"
    CTYPE="C31"
    chsrc="$c31src"
    echo "build.sh> chsrc reset to $c31src"
    FCMDIR="$c31src/source/fcm"
    echo "build.sh> Assuming the platform is 'gnu' when compiling c31"
    PREFLX="$c31src/tool/prefx_gnu"
    PREFDIR="$c31src/build/gnu"
    LIBDIRCH="$c31src/lib/gnu"
  else 
    EXE="../bin/$compiler/C35$programname"
    CTYPE="C35"
    chsrc="$c35src"
    echo "build.sh> chsrc reset to $c35src"
    FCMDIR="$c35src/source/fcm"
    echo "build.sh> Assuming the platform is 'gnu' when compiling c35"
    PREFLX="$c35src/tool/prefx_gnu"
    PREFDIR="$c35src/build/gnu"
    LIBDIRCH="$c35src/lib/gnu"
  fi
  BLAS_EXCLUDE_LIST="dnrm2.o daxpy.o dcopy.o ddot.o"
  SRCCH="charmm_main.src energy.src"
elif [ "$target" = "unoptim" ]; then
  EXE="../bin/$compiler/UN$programname"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [ "$target" = "amhoptim" ]; then
  EXE="../bin/$compiler/AMH$programname"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [ "$target" = "dlfoptim" ]; then
  if [ "$compiler" != "ifort" && "$compiler" != "pgi" && "$compiler" != "g95" && "$compiler" != "gfortran" ]; then
    echo "build.sh> ERROR: can only compile DLFWrapper with ifort, pgi, g95, and gfortran"
    exit
  fi
  EXE="../bin/$compiler/DLF$programname."
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [ "$target" = "jboptim" ]; then
  if [ "$compiler" != "ifort" ]; then
    echo "build.sh> ERROR: can only compile Bowman potential with ifort"
    exit
  fi
  EXE="../bin/$compiler/JB$programname."
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  chsrc=""
elif [[ "$target" = "clean" || "$target" = "cleanall" ]]; then
  EXE="$target"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH="charmm_main.src energy.src"
  LIBDIRCH=""
  chsrc=""
else
  EXE="$target"
  BLAS_EXCLUDE_LIST=""
  CTYPE=""
  PREFLX=""
  FCMDIR=""
  PREFDIR=""
  SRCCH=""
  LIBDIRCH=""
  exemessage=0
  chsrc=""
fi

if [ "$target" = "jboptim" ]; then
  JBMODS="-IBowman"
else
  JBMODS="-IBowman/fake"
fi

echo "build.sh> target = $target"

make="make "

# need to determine if we could use some of the previous build
if [ ! "$arguments" = "$lastargs" ]; then
  # determine last target and compiler
  LastTarget="$defaulttarget"
  LastCompiler="$defaultcompiler"
  LastQM="$defaultchmqm"
  LastChtype=""
  for opt in $lastargs; do 
    if [[ "$opt" = "optim" || "$opt" = "coptim" || "$opt" = "unoptim" || "$opt" = "jboptim" ]]; then
      LastTarget="$opt"
    fi
    if [[ "$opt" = "pgi" || "$opt" = "ifort" || "$opt" = "nag" || "$opt" = "ibm" || "$opt" = "g95" || "$opt" = "pathscale" ]]; then
      LastCompiler="$opt"
    fi
  # determine last QM method used 
    if [ "$opt" = "dftb" ]; then
      LastQM="$opt" 
    fi
    if [[ "$opt" = "c31" || "$opt" = "c35" ]]; then
      LastChtype="$opt"
    fi
  done
  # compare current and previous targets
  if [ ! "$LastTarget" = "$target" ]; then
    if [[ "$LastTarget" = "coptim" && "$target" != "coptim" || "$LastTarget" != "coptim" && "$target" == "coptim" ]]; then
      echo "build.sh> Need to relink BLAS library because BLAS_EXCLUDE_LIST has changed"
      cd ../../BLAS && make cleanlib && cd ../OPTIM/source
    fi
  fi
  # compare current and previous compilers
  if [ ! "$LastCompiler" = "$compiler" ]; then # we need to rebuild from scratch
    echo "build.sh> Last target was built with a different compiler - rebuilding from scratch"
    $make clean
  fi
  if [ "$target" = "coptim" ]; then
    # compare CHARMM QM methods
    if [ ! ["$LastQM" = "dftb" -a "$chmqm" = "dftb"] -o ["$LastQM" = "none" -a "$chmqm" = "none"] ]; then
      echo "build.sh> Last target used a different CHARMM QM method - rebuilding CHARMM"
      echo "build.sh> cd ${chsrc} && ./clean.csh && cd $optsrc"
      cd ${chsrc} && ${chsrc}/clean.csh && cd $optsrc
      echo "build.sh> cd ${optsrc}/CHARMMinterface && make clean && cd $optsrc"
      cd ${optsrc}/CHARMMinterface && make clean && cd $optsrc
    fi
    # compare CHARMM versions
    if [ ! `expr substr "$LastChtype" 2 2` = `expr substr "$CTYPE" 2 2` ]; then
      echo "build.sh> Cleaning the CHARMMinterface directory because the CHARMM version may have changed"
      cd CHARMMinterface && make clean && cd ../
    fi
  fi
fi

# all set up and ready to run make
success=0
if [ "$target" = "coptim" ]; then # need to compile CHARMM first
   if [ ! -e $chsrc ]; then
     echo "build.sh> ERROR: CHARMM source directory $chsrc does not exist!"
     exit
   else
     echo "build.sh> chsrc = $chsrc"
   fi
   command="install.com gnu"
   command="$command $chmsize keepo keepf"
   if [ "$chmqm" = "dftb" ]; then
     if [ "$CTYPE" = "C31" ]; then
       echo "build.sh> ERROR: DFTB and CHARMM31 are incompatible!"
       exit
     fi
     command="$command T"
# QMTYPE variable passed to Makefile to include SCC-DFTB libraries 
     QMTYPE="DFTB"
   fi
   if [ "$compiler" = "pgi" ]; then
     command="$command PGF90"
     if [ "$level" = "opt" ]; then
       command="$command OPT"
     fi
     if [ "$architecture" = "amd" ]; then
       command="$command AMD"
     fi
   elif [ "$compiler" = "ifort" ]; then
     command="$command ifort"
   else
     echo "build.sh> ERROR: This script cannot compile CHARMM directory with the selected compiler!"
     exit
   fi
#   echo "build.sh> cd ${chsrc} && ${chsrc}/$command && cd $optsrc"
   echo "build.sh> cd ${chsrc} && ./$command && cd $optsrc"
   echo
#   cd ${chsrc} && ${chsrc}/$command && cd $optsrc
   cd ${chsrc} && ./$command && cd $optsrc
   echo "build.sh> tail of $PREFDIR/gnu.log"
   tail -2 $PREFDIR/gnu.log
   echo " "
elif [ "$target" = "amhoptim" ]; then
   if [ ! -e $amhsrc ]; then
     echo "build.sh> ERROR: AMH source directory $amhsrc does not exist!"
     exit
   else
     echo "build.sh> amhsrc = $amhsrc"
   fi
#   echo "build.csh> cd ${amhsrc} && make  && cd $optsrc"
#   cd ${amhsrc} && make  && cd $optsrc
   echo "  "
fi
$make $quiet $EXE FC="$FC" NOOPT="$NOOPT" NOOPT2="$NOOPT2" SEARCH_PATH="$SEARCH_PATH $JBMODS" FFLAGS="$FFLAGS" FREEFORMAT_FLAG="$FREEFORMAT_FLAG" EXTRA_FLAGS="$EXTRA_FLAGS" EXTRA_CFLAGS="$EXTRA_CFLAGS" $target="$EXE" LIBDIRCH="$LIBDIRCH" \
BLAS_EXCLUDE_LIST="$BLAS_EXCLUDE_LIST" CTYPE="$CTYPE" QMTYPE="$QMTYPE" PREFLX="$PREFLX" FCMDIR="$FCMDIR" PREFDIR="$PREFDIR" SRCCH="$SRCCH" \
CHSRC="$chsrc" SWITCH="$SWITCH" >& make.out && success=1

# save arguments that were passed to this script into lastargs file
if [[ ! ["$@" = "last" || "$target" = "clean" || "$target" = "cleanall"] ]]; then
  echo "$@" > lastargs
fi

if [ "$success" = 1 ]; then
  if [ "$exemessage" = 1 ]; then
    echo "build.sh> OPTIM executable = $optsrc/$EXE"
  fi
  echo "build.sh> Success. See make.out for detailed compiler output"
else
  echo "build.sh> Failure - refer to make.out"
  echo "build.sh> grep for error messages in make.out gives:"
  grep -i error make.out
fi
