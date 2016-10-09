#!/bin/bash
#PBS -S /bin/sh
#PBS -N hw3
#PBS -A eecs587f16_flux
#PBS -l qos=flux
#PBS -l procs=36,walltime=0:05:00
#PBS -l pmem=1000mb
#PBS -q flux
#PBS -m abe
#PBS -j oe
#PBS -V
echo "I ran on:"
#cat $PBS_NODEFILE
# Let PBS handle your output
cd /home/binbinli/hw3/EECS-587/
mpirun -np 1 ./test 1000 4000
mpirun -np 4 ./test 1000 4000
mpirun -np 16 ./test 1000 4000
mpirun -np 36 ./test 1000 4000
mpirun -np 1 ./test 1000 4000
mpirun -np 4 ./test 1000 4000
mpirun -np 16 ./test 1000 4000
mpirun -np 36 ./test 1000 4000

mpirun -np 1 ./test 2000 500
mpirun -np 4 ./test 2000 500
mpirun -np 16 ./test 2000 500
mpirun -np 36 ./test 2000 500
mpirun -np 1 ./test 2000 500
mpirun -np 4 ./test 2000 500
mpirun -np 16 ./test 2000 500
mpirun -np 36 ./test 2000 500
