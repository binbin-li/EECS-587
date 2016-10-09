#!/bin/bash

rm test
mpiCC -o test test.cpp -O3
#mpirun -np 1 ./test 1000 4000
