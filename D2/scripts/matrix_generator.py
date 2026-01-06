#!/usr/bin/env python3
import sys
import random

# Weak-scaling parameters
rows_per_rank = 20000
nnz_per_row   = 20
matrix_name   = "weak_scaling"

# Read number of MPI processes
if len(sys.argv) != 2:
    print("Usage: python matrix_generator.py <num_processes>")
    sys.exit(1)

P = int(sys.argv[1])

# Global matrix dimensions
N = rows_per_rank * P
total_nnz = rows_per_rank * nnz_per_row * P

print("Generating weak-scaling matrix")
print(f"  MPI processes      : {P}")
print(f"  Rows per rank      : {rows_per_rank}")
print(f"  Nonzeros per row   : {nnz_per_row}")
print(f"  Global size        : {N} x {N}")
print(f"  Total NNZ          : {total_nnz}")

random.seed(42)     # For reproducibility

filename = f"{matrix_name}_NP{P}.mtx"

with open(filename, "w") as f:
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write("% Weak-scaling sparse matrix for MPI SpMV\n")
    f.write(f"% MPI processes      : {P}\n")
    f.write(f"% Rows per rank      : {rows_per_rank}\n")
    f.write(f"% Nonzeros per row   : {nnz_per_row}\n")
    f.write(f"% Global size        : {N} x {N}\n")
    f.write(f"% Total NNZ          : {total_nnz}\n")
    f.write(f"{N} {N} {total_nnz}\n")

    # Generate entries row by row
    for row in range(N):
        for _ in range(nnz_per_row):
            col = random.randrange(N)
            val = random.random()
            f.write(f"{row+1} {col+1} {val:.6e}\n")

print(f"Matrix written to: {filename}")
