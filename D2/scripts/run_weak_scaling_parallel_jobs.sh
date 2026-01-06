#!/bin/bash
# Load modules
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load python-3.10.14_gcc91

mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out # C++ code

for NPROCS in 2 4 8 16 32 64 128
do
  python3 scripts/matrix_generator.py $NPROCS
done

for NPROCS in 2 4 8 16 32 64 128
do
  for MODE in 0 1 2
  do
    MATRIX="weak_scaling_NP${NPROCS}.mtx"
    echo "Submitting job for NPROCS=$NPROCS MODE=$MODE MATRIX=$MATRIX"
    qsub -v MATRIX="$MATRIX",NPROCS="$NPROCS",MODE="$MODE" scripts/run_weak_scaling.pbs
  done
done
