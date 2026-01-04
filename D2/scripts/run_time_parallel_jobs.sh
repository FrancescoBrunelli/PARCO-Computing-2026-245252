#!/bin/bash

MAX_PROCS=128
INPUT="$1"    # Input from terminal
MODE="$1"    # Input from terminal

mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out # C++ code

# for NPROCS in 2 4 8 16 32 64 128
for NPROCS in 2 4
do
  echo "Submitting job for NPROCS=$NPROCS"
  qsub -v INPUT="$INPUT",NPROCS="$NPROCS",MODE="$MODE" scripts/run.pbs
done
