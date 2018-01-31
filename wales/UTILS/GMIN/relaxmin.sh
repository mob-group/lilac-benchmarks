#!/bin/bash
#PBS -q h1
#PBS -j oe
#PBS -N relaxmin

cd $PBS_O_WORKDIR
TMP=/scratch/wales/$PBS_JOBID
mkdir -p $TMP

for file in * ; do
   cp $file $TMP
done

cd $TMP

OPTIM > output

cp min.data.info $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
rm -rf $TMP

echo
qstat -f ${PBS_JOBID}@volkhan
echo

