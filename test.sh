#!/bin/bash

#SRC_DIR=/pscratch/sd/g/gabeh98/March2024/cs267-hw2-2
SRC_DIR=/Users/gabrielraulet/Dev/March2024/cs267-hw2-2
CORRECTNESS_DIR=$SRC_DIR/hw2-correctness
BUILD_DIR=$SRC_DIR/build

mpirun -np 4 $BUILD_DIR/mpi -s 1 -n 1000 -o $BUILD_DIR/mpi.n1000.out

python $CORRECTNESS_DIR/correctness-check.py $BUILD_DIR/mpi.n1000.out $CORRECTNESS_DIR/verf.n1000.out
