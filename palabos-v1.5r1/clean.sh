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

sudo rm -rf cmake*
sudo rm -rf CMake*
sudo rm -rf *.cmake
sudo rm -rf *.ninja
sudo rm -rf ninja*
sudo rm -rf *.dSYM
sudo rm -rf src
sudo rm -rf viscosityTest
sudo rm -rf *.dylib
sudo rm -rf dlib*
sudo rm -rf Make*
sudo rm -rf *.so
sudo rm -rf *.txt
