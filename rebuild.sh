#!/bin/bash

rm -rf build
mkdir build
cd build
#cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-O0 -g -fsanitize=address -fno-omit-frame-pointer"
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_CXX_FLAGS="-O0 -g"
make -j12
