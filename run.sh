#!/bin/bash

num_parts="6000000"
threads="64 16 4 1"


for t in $threads
do
  for n in $num_parts
  do
    srun -N 1 --ntasks-per-node=$t ./build/mpi -s 1 -n $n > data/t$t-$n.out
  done
done