#!/bin/bash

SRC_DIR=/Users/gabrielraulet/Dev/March2024/cs267-hw2-2
BUILD_DIR=$SRC_DIR/build

rm -rf $BUILD_DIR
mkdir $BUILD_DIR
cd $BUILD_DIR
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++
#cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_CXX_FLAGS="-O0 -g"
make -j12
