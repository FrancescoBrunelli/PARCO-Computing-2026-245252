# src

## Table of Contents
- [Introduction](#introduction)
- [Directory Structure](#directory-structure)
- [Main Components](#main-components)
- [Compilation](#compilation)
- [Change Parameters](#change-parameters)
- [Compatibility](#compatibility)
- [Notes](#notes)
- [Example Usage](#example-usage)

---

## Introduction
This directory contains the source files for the **Distributed Sparse Matrix-Vector Multiplication (SpMV)** project.
## Directory Structure
``` text
src
├── README.md
├── main.cpp
├── main.out
├── matrix.cpp
├── matrix.h
├── mpiSpMV.cpp
└── mpiSpMV.h
```

## Main Components
### `main.cpp`
- Calls the correct matrix decomposition strategy based on selected input mode

### `matrix.h` / `matrix.cpp`
Implement functions responsible for:
- matrix reading and loading
- COO sorting
- COO -> CSR conversion
- Provides both sequential (`CSRmul`), partial (`PartialCSRmul`) and parallel `P_CSRmul` SpMV routines

### `mpiSpMV.h` / `mpiSpMV.cpp`
Implement the three decomposition strategies:
- 1D Block
- 1D Cyclic
- 2D Cartesian Block

And a sequential MPI version of the SpMV computation
## Compilation
The code can be compiled using mpicxx with the following flags

``` sh
mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out
```

## Change Parameters
### Change Decomposition Strategy
To change the decomposition function adopted by the program is sufficient to call it using a different MODE parameter:
- MODE = 0 --> 1D Block Decomposition
- MODE = 1 --> 1D Cyclic Decomposition
- MODE = 2 --> 2D Cartesian Block Decomposition
- MODE = 3 --> Sequential execution

## Compatibility
The program needs a .mtx matrix as an input.

**NOTE:** Not all matrix formats are compatible.

**Supported storage formats:**
- general
- symmetric

**Supported datatypes:**
- integer
- real
- complex
- pattern

## Notes
- The input matrix must be in the **Matrix Market (.mtx)** format
- It is recommended to run the program from the main folder
- **Outputs:** 
	- Computation time (expressed in milliseconds)
	- Average number of nonzero (NNZ) elements among MPI-processes
	- Maximum number of nonzero (NNZ) elements among MPI-processes
	- Minimum number of nonzero (NNZ) elements among MPI-processes
	- Number of floating point operations per second (FLOPS) (expressed in GFLOPS)
	- Communication time (expressed in milliseconds)
	- Peak memory footprint among MPI processes (expressed in MB)

## Example Usage
``` sh
# Compile
mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out

# Run
mpirun -np 2 ./src/main.out matrix.mtx 1

```
