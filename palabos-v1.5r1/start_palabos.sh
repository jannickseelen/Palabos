#!/bin/bash
#PBS -N PALABOS 
#PBS -m abe
#PBS -M jannick.seelen@gmail.com
#PBS -l nodes=1:ppn=48
#PBS -q rrr-pnr
#PBS -o "~/Palabos/run.log"
#PBS -j oe

cd ~/Palabos/palabos-v1.5r1/build

make clean

make

mpirun ./viscosityTest /input/parameters.xml
