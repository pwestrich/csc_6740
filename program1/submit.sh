#!/bin/bash
#PBS -N pmw_matrix
#PBS -W x=PARTITION:csc
#PBS -l nodes=1:ppn=36
#PBS -l mem=8192mb
#PBS -l walltime=3:00:00
#PBS -m bea
#PBS -M pmwestrich42@students.tntech.edu

#change to the program directory
cd ~/csc_6740/program1

#build the program
make clean
make -j32

#run the program
mpiexec -n 36 main 3600
