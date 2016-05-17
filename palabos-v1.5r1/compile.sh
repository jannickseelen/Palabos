#!/bin/bash
#PBS -N PALABOS
#PBS -m abe
#PBS -M jannick.seelen@gmail.com
#PBS -l nodes=2:ppn=48
#PBS -q rrr-pnr
#PBS -o "~/Palabos/run.log"
#PBS -j oe

set -e

cd ~

cd build

echo "Running CMake"

sudo cmake -DCMAKE_BUILD_TYPE=DEBUG -GNinja ~/Palabos/palabos-v1.5r1

echo "Running Ninja"

sudo ninja -v

cd ..
