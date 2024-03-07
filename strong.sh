#!/bin/bash

num_parts="1000 10000 100000 1000000 6000000"
threads="2 8 32"


for t in $threads
do
  for n in $num_parts
  do
    srun -N 1 --ntasks-per-node=$t ./build/mpi -s 1 -n $n > data_strong/t$t-$n.out
  done
done