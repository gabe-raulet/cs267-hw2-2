#!/bin/bash

CORRECTNESS_DIR=/Users/gabrielraulet/Dev/March2024/hw2-2/hw2-correctness
BUILD_DIR=/Users/gabrielraulet/Dev/March2024/hw2-2/build

mpirun -np 1 $BUILD_DIR/mpi -s 1 -n 1000 -o $BUILD_DIR/mpi.n1000.out
mpirun -np 1 $BUILD_DIR/mpi -s 1 -n 5000 -o $BUILD_DIR/mpi.n5000.out

python $CORRECTNESS_DIR/correctness-check.py $BUILD_DIR/mpi.n1000.out $CORRECTNESS_DIR/verf.n1000.out
python $CORRECTNESS_DIR/correctness-check.py $BUILD_DIR/mpi.n5000.out $CORRECTNESS_DIR/verf.n5000.out
