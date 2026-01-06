#!/usr/bin/env python3
import numpy as np
import sys

# Weak-scaling parameters
rows_per_rank = 20000     # rows handled by each MPI rank
nnz_per_row   = 20        # nonzeros per row (SpMV realistic)
matrix_name   = "weak_scaling"

# Read number of MPI processes
if len(sys.argv) != 2:
    print("Usage: python generate_weak_scaling_mtx.py <num_processes>")
    sys.exit(1)

P = int(sys.argv[1])

# Global matrix dimensions
N = rows_per_rank * P
total_nnz = rows_per_rank * nnz_per_row * P

print("Generating weak-scaling matrix for MPI SpMV")
print(f"  MPI processes        : {P}")
print(f"  Rows per rank        : {rows_per_rank}")
print(f"  Nonzeros per row     : {nnz_per_row}")
print(f"  Global size          : {N} x {N}")
print(f"  Total NNZ            : {total_nnz}")

# Generate COO entries
rng = np.random.default_rng()

# Each row has exactly nnz_per_row entries
row_indices = np.repeat(
    np.arange(N, dtype=np.int64),
    nnz_per_row
)

# Column indices uniformly distributed over the global matrix
col_indices = rng.integers(
    low=0,
    high=N,
    size=total_nnz,
    dtype=np.int64
)

# Random values
data = rng.random(total_nnz)

# Write Matrix Market file
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

    # (1-based indexing)
    for i, j, v in zip(row_indices, col_indices, data):
        f.write(f"{i+1} {j+1} {v:.6e}\n")

print(f"Matrix written to: {filename}")
