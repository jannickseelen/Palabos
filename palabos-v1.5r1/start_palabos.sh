#!/bin/bash
#PBS -N PALABOS 
#PBS -m abe
#PBS -M jannick.seelen@gmail.com
#PBS -l nodes=2:ppn=48
#PBS -q rrr-pnr
#PBS -o "~/Palabos/run.log"
#PBS -j oe

set -e

cd ~/Palabos

git pull

cd ~/Palabos/palabos-v1.5r1/build/

echo "Removing previous build"

rm -rf cmake*
rm -rf CMake*
rm -rf *.cmake
rm -rf *.ninja
rm -rf ninja*
rm -rf *.dSYM
rm -rf src
rm -rf viscosityTest
rm -rf *.dylib
rm -rf dlib*

echo "Loading mpi module"

module load "/etc/modulefiles/openmpi-x86_64"

echo "Running CMake"

cmake -v ~/Palabos/palabos-v1.5r1

echo "Running make"

make -j48 -d

echo "Running MPI"

chmod +x viscosityTest

mpirun ./viscosityTest /input/parameters.xml
