#!/bin/bash
#PBS -S /bin/sh
#PBS -N hw3
#PBS -A eecs587f16_flux
#PBS -l qos=flux
#PBS -l procs=36,walltime=0:05:00
#PBS -l pmem=4000mb
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
echo "I ran on:"
#cat $PBS_NODEFILE
# Let PBS handle your output
cd /home/binbinli/hw3
mpirun -np 16 ./test 1000 4000