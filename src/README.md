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
This directory contains the source files for the **Sparse Matrix-Vector Multiplication (SpMV)** project.

The implementation supports both **sequential** and **parallel (OpenMP)** execution, allowing performance comparisons across different scheduling policies.

## Directory Structure
``` text
src
├── README.md
├── main.cpp
├── matrix.cpp
└── matrix.h
```

## Main Components
### `main.cpp`
- Parses input matrix via command line
- Stores the COO format and sorts it in row-major order
- Converts it into CSR
- Computes the product between CSR and a randomly generated vector
- Measures and prints the elapsed time of the multiplication in milliseconds

### `matrix.h` / `matrix.cpp`
Implement functions responsible for:
- matrix reading and loading
- COO sorting
- COO -> CSR conversion
- Provides both sequential (`CSRmul`) and parallel `P_CSRmul` SpMV routines

## Compilation
The code can be compiled using g++ with OpenMP and the following flags

``` sh
g++ -std=c++11 -O3 -march=native -ftree-vectorize -fopenmp main.cpp matrix.cpp -o main.out
```

## Change Parameters

### Switch Between Sequential and Parallel Execution
To switch between sequential and parallel execution is necessary to edit the `main.cpp` file:
- To make it run sequentially: comment `P_CSRmul(aRow, aCol, aVal, v, out);` and uncomment `CSRmul(aRow, aCol, aVal, v, out);`
- To make it run in parallel: comment `CSRmul(aRow, aCol, aVal, v, out);`  and uncomment `P_CSRmul(aRow, aCol, aVal, v, out);`

### Change Scheduling Policy
To change the OMP scheduling policy, edit the following line inside `P_CSRmul()` in `matrix.cpp`:
``` c++
#pragma omp parallel for schedule(guided)
```

**Available options:**
- `static, <chunk_size>`
- `guided, <chuck_size> (optional)`
- `dynamic, <chunk_size> (optional)`
- `runtime`

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
- The program (if run in parallel) will use the maximum amount of threads available unless otherwise specified with `OMP_NUM_THREADS`
- **Output:** elapsed time during the SpMV multiplication, expressed in milliseconds, printed to stdout

## Example Usage
``` sh
# Compile
g++ -std=c++11 -O3 -march=native -ftree-vectorize -fopenmp main.cpp matrix.cpp -o main.out

# Run
./main.out ../fd12.mtx

```
