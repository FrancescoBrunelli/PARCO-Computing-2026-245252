#!/bin/bash

MAX_PROCS=128
INPUT="$1"    # Input from terminal
#MODE="$2"    # Input from terminal

# Load modules
module load gcc91
module load mpich-3.2.1--gcc-9.1.0

mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out # C++ code

for NPROCS in 2 4 8 16 32 64 128
do
  for MODE in 0 1 2
  do
    echo "Submitting job for NPROCS=$NPROCS MODE=$MODE"
    qsub -v INPUT="$INPUT",NPROCS="$NPROCS",MODE="$MODE" scripts/run_time.pbs
  done
done

NPROCS=1
MODE=3
echo "Submitting job for NPROCS=$NPROCS MODE=$MODE"
qsub -v INPUT="$INPUT",NPROCS="$NPROCS",MODE="$MODE" scripts/run_time.pbs
