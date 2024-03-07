#!/bin/bash

num_parts="1000 6000 10000 60000 100000 600000 1000000 6000000"
threads="1 2 4 8 16 32 64"
ratio="1000 4000 10000 40000"



for r in $ratio
do
  for t in $threads
  do
    let n=$r*$t
    srun -N 1 --ntasks-per-node=$t ./build/mpi -s 1 -n $n > data_weak/t$t-$n.out
  done
  let n=$r*128
  srun -N 2 --ntasks-per-node=64 ./build/mpi -s 1 -n $n > data_weak/t128-$n.out
done