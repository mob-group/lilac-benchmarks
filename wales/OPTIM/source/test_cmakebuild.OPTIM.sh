#!/bin/bash

build_dir=$(readlink -f ../build_test)
source_dir=${HOME}/svn/OPTIM/source
cmake=${HOME}/bin/cmake

mkdir ${build_dir} || exit 1
cd ${build_dir}

for compiler in gfortran pgf90 ifort; do
  # Check whether our compiler is loaded and try to load it, if not.
  ${compiler} --version > /dev/null 2>&1 || {
    module load ${compiler}
  }
  # If it still doesn't work, exit.
  ${compiler} --version > /dev/null 2>&1 || exit 1

  # Iterate through release and debug modes.
  for build_type in Release Debug
  do
    curdir=${build_dir}/${compiler}_${build_type}
    mkdir ${curdir}
    cd ${curdir}

    if [ ${compiler} == ifort ] || [ ${compiler} == pgf90 ]; then
      # Now build with AMBER9
      mkdir ./AMBER
      cd ./AMBER
      FC=${compiler} ${cmake} -DCMAKE_BUILD_TYPE=${build_type} -DWITH_DMAOPTIM=no -DWITH_CHARMM35=no -DWITH_AMBER9=yes -DWITH_AMH=no ${source_dir} > /dev/null 2>&1
      for target in OPTIM A9OPTIM; do
        echo -n "${compiler} ${build_type} ${target}: "
        make -j 4 ${target}  > /dev/null 2>&1
        if [ "$?" == "0" ]; then
          echo -e "\e[01;32mok\e[00m"
        else 
          echo -e "\e[00;31mfailed\e[00m"
        fi
      done
      # If building with ifort or pgf90, we can build CHARMM35.
      mkdir ../CHARMM35/
      cd ../CHARMM35/
      FC=${compiler} ${cmake} -DCMAKE_BUILD_TYPE=${build_type} -DWITH_DMAOPTIM=no -DWITH_CHARMM35=yes -DWITH_AMBER9=no -DWITH_AMH=no ${source_dir} > /dev/null 2>&1
      echo -n "${compiler} ${build_type} C35OPTIM: "
      make C35OPTIM > /dev/null 2>&1
      if [ "$?" == "0" ]; then
        echo -e "\e[01;32mok\e[00m"
      else 
        echo -e "\e[00;31mfailed\e[00m"
      fi
    else
      # Otherwise, don't build CHARMM35.
      FC=${compiler} ${cmake} -DCMAKE_BUILD_TYPE=${build_type} -DWITH_DMAOPTIM=no -DWITH_CHARMM35=no -DWITH_AMBER9=yes -DWITH_AMH=no ${source_dir} > /dev/null 2>&1
      for target in OPTIM A9OPTIM; do
        echo -n "${compiler} ${build_type} ${target}: "
        make -j 8 ${target}  > /dev/null 2>&1
        if [ "$?" == "0" ]; then
          echo -e "\e[01;32mok\e[00m"
        else 
          echo -e "\e[00;31mfailed\e[00m"
        fi
      done
    fi
  done
done
