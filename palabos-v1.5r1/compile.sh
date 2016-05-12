#!/bin/bash
#PBS -N PALABOS
#PBS -m abe
#PBS -M jannick.seelen@gmail.com
#PBS -l nodes=2:ppn=48
#PBS -q rrr-pnr
#PBS -o "~/Palabos/run.log"
#PBS -j oe

set -e

cd ~/Palabos/palabos-v1.5r1/build/

echo "Loading mpi module"

module load "/etc/modulefiles/openmpi-x86_64"

echo "Running CMake"

cmake -DCMAKE_BUILD_TYPE=DEBUG ~/Palabos/palabos-v1.5r1

echo "Running make"

make -j10 -d
