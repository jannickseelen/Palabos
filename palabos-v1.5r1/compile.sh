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

cmake -DCMAKE_BUILD_TYPE=DEBUG -DDLIB_ENABLE_ASSERTS=ON -DDLIB_ENABLE_STACK_TRACE=ON -DDLIB_ISO_CPP_ONLY=ON ~/Palabos/palabos-v1.5r1

echo "Running make"

make -j10 --debug=b
