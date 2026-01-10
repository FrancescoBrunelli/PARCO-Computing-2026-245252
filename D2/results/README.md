# Results
This directory contains all the performance metrics collected from the experiments. Each `.csv` file contains:
- Time elapsed (expressed in milliseconds) during the multiplication between the CSR format of the matrix and a randomly generated vector
- Average number of nonzero (NNZ) elements among MPI-processes
- Maximum number of nonzero (NNZ) elements among MPI-processes
- Minimum number of nonzero (NNZ) elements among MPI-processes
- Number of floating point operations per second (FLOPS) (expressed in GFLOPS)
- Communication time (expressed in milliseconds)
- Peak memory footprint among MPI processes (expressed in MB)  
For each one of 100 executions.
The results are organized by **matrix** and decomposition strategy adopted.

---

## Directory Structure

``` text
results
├── README.md
├── matrix_name         # One for each matrix
│   ├── 1D_Block
│   ├── 1D_Cyclic
│   ├── 2D_Block
│   └── <matrix_name>_NP1.csv
└── weak_scaling
    ├── 1D_Block
    ├── 1D_Cyclic
    └── 2D_Block
```
