#!/bin/bash

SRC_DIR=/Users/gabrielraulet/Dev/March2024/cs267-hw2-2
CORRECTNESS_DIR=$SRC_DIR/hw2-correctness
BUILD_DIR=$SRC_DIR/build

$BUILD_DIR/mpi -s 1 -n 1000 -o $BUILD_DIR/serial.n1000.out
python $CORRECTNESS_DIR/correctness-check.py $BUILD_DIR/serial.n1000.out $CORRECTNESS_DIR/verf.n1000.out

$BUILD_DIR/mpi -s 1 -n 5000 -o $BUILD_DIR/serial.n5000.out
python $CORRECTNESS_DIR/correctness-check.py $BUILD_DIR/serial.n5000.out $CORRECTNESS_DIR/verf.n5000.out

#mpirun -np 1 $BUILD_DIR/mpi -s 1 -n 1000 -o $BUILD_DIR/mpi.n1000.out
#rm $BUILD_DIR/mpi.n1000.out

#mpirun -np 4 $BUILD_DIR/mpi -s 1 -n 1000 -o $BUILD_DIR/mpi.n1000.out
#rm $BUILD_DIR/mpi.n1000.out

#mpirun -np 9 $BUILD_DIR/mpi -s 1 -n 1000 -o $BUILD_DIR/mpi.n1000.out
#rm $BUILD_DIR/mpi.n1000.out

