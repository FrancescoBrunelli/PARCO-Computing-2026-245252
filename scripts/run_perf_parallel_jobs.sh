#!/bin/bash

MAX_THREADS=64
INPUT="$1"    # Input passato da terminale

g++ -std=c++11 -O3 -march=native -ftree-vectorize src/main.cpp src/matrix.cpp -o src/main.out

for NTHREADS in 1 2 4 8 16 32 64
do
    echo "Submitting job for NTHREADS=$NTHREADS"
    qsub -v INPUT="$INPUT",NTHREADS=$NTHREADS run_perf.pbs
done
