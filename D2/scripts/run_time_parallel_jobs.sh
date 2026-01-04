#!/bin/bash

MAX_PROCS=128
INPUT="$1"    # Input from terminal
MODE="$2"    # Input from terminal

mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out # C++ code

# for NPROCS in 2 4 8 16 32 64 128
for NPROCS in 2 4
do
  for MODE in 0 1 2
  do
    echo "Submitting job for NPROCS=$NPROCS MODE=$MODE"
    qsub -v INPUT="$INPUT",NPROCS="$NPROCS",MODE="$MODE" scripts/run_time.pbs
  done
done
