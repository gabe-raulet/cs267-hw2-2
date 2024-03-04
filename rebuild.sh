#!/bin/bash

rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++
make -j12
