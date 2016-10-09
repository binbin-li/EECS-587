#!/bin/bash

rm test
mpic++ -o test test.cpp -O3
#mpirun -np 1 ./test 1000 4000
