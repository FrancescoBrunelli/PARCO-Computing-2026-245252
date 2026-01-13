# Scripts

## Table of Contents

- [Introduction](#introduction)
- [Setup Scripts](#setup-scripts)
- [Run Scripts](#run-scripts)
- [Python Scripts](#python-scripts)
- [Matrices used](#matrices-used)

---

## Introduction

This folder contains all the scripts needed to set up the environment, run different types of tests, process the data and generate plots based on different performance metrics

## Setup Scripts

Scripts in this section are responsible for setting up the environment and downloading required data

- `download_matrices.sh` - Downloads and extracts five sparse matrices
  
  Usage:
``` sh
# Make sure to make it executable first:
chmod +x ./scripts/download_matrices.sh

# Example usage of run_time_parallel_jobs.sh
./scripts/download_matrices.sh
```
  
  Alternatively you can download the matrices you prefer from [here](https://sparse.tamu.edu), make sure they are written in any of the different formats compatible with the program [check compatible formats](../src/README.md#compatibility).

## Run Scripts

Scripts used to run multiple tests on the compute nodes automatically

It is important to run them from the main folder of the repository

### `run_time.pbs`
Runs 100 executions of the program, saving: computation time, communication time, NNZ distribution statistics, number of floating point operations per second and memory footprint into a `.csv` file named: `<matrix_name>_NP<number_of_MPI_processes_used>`.
The `.csv` file is located at results/matrix_name/decomposition_used/

**Usage (after program compilation):**
``` sh
# Example usage of run_time.pbs
qsub -v INPUT=matrix.mtx,NPROCS=32,MODE=2 scripts/run_time.pbs
```
### `run_weak_scaling.pbs`
A similar version of the script `run_time.pbs` dedicated for weak_scaling experiments. It is meant to run by call from the bash script `run_weak_scaling_parallel_jobs.sh`

**Usage (after program compilation):**
``` sh
# Example usage of run_time.pbs
qsub -v INPUT=matrix.mtx,NPROCS=32,MODE=2 scripts/run_weak_scaling.pbs
```

### `run_time_parallel_jobs.sh`
Submits a total of 21 `run_time.pbs` jobs, one for each combination of number of MPI processes (from 2 to 128 in powers of two) and decomposition strategies (from 0 to 2)

**Usage:**
``` sh
# Make sure to make it executable first:
chmod +x ./scripts/run_time_parallel_jobs.sh

# Example usage of run_time_parallel_jobs.sh
./scripts/run_time_parallel_jobs.sh matrix.mtx
```
### `run_weak_scaling_parallel_jobs.sh`
Generates 7 pseudo-random sparse matrices (one for each amount of MPI processes from 2 to 128 in powers of two) by calling a dedicated python script: `matrix_generator.py`
Once each matrix is generated, this script will submit a total of 21 `run_weak_scaling.pbs` jobs (one for each combination of matrix and modes).

**Note:** the generation of the 7 matrices will be 20+ minutes long. Before running this script it is recommended to open an interactive session.

**Usage:**
``` sh
# Make sure to make it executable first:
chmod +x ./scripts/run_weak_scaling_parallel_jobs.sh

# Open an interactive session
qsub -I -q queueName -l select=1:ncpus=64: mem=64gb:walltime=01:00:00

# Example usage of run_perf_parallel_jobs.sh
./scripts/run_weak_scaling_parallel_jobs.sh
```
## Python Scripts

Organized as follows:
### `matrix_generator.py`
Generates a pseudo-random sparse matrix proportional to the number of processes provided as input. The matrix is in MatrixMarket format and is stored as: `weak_scaling_NP<nprocs>.mtx`
  
  Usage: `python3 scripts/matrix_generator.py <nprocs>`

### Plot Scripts
- `run_plot.py`: For each matrix in /results and for each decomposition strategy adopted, reads all the `.csv` files, extracts for each file the 90th percentile of: computation time, communication time, GFLOPS, memory footprint and for each one creates a linear plot for each decomposition strategy
  
  Usage: `python3 scripts/run_plot.py`
  
- `speedup_plot.py`: For each matrix in /results and for each decomposition strategy adopted, reads all the `.csv` files, extracts for each file the 90th percentile of the computation time, calculates the speedup factor using as baseline se sequential results stored in results/matrix_name/`<matrix_name>_NP1.csv` and creates an histogram plot for each decomposition strategy

  Usage: `python3 scripts/speedup_plot.py`
  
- `weak_scaling_plot.py`: Similarly to `run_plot.py` generates linear plots from weak_scaling results

  Usage: `python3 scripts/weak_scaling_plot.py`
  
- `efficiency_plot.py`: For each matrix in /results and for each decomposition strategy adopted generates the corresponding efficiency linear plot

  Usage: `python3 scripts/efficiency_plot.py`

## Matrices used
- **fd12:** 7.500x7.500 28.462nnz (129KB)
``` sh
wget https://suitesparse-collection-website.herokuapp.com/MM/Hohn/fd12.tar.gz
```
  
- **sinc12:** 7.500x7.500 283.992nnz (1,5MB)
``` sh
wget https://suitesparse-collection-website.herokuapp.com/MM/Hohn/sinc12.tar.gz
```

- **epb1:** 14.724x14.734 95.053nnz (935KB)
``` sh
wget https://suitesparse-collection-website.herokuapp.com/MM/Averous/epb1.tar.gz
```

- **c-46:** 14.913x14.913 130.397nnz (541KB)
``` sh
wget https://suitesparse-collection-website.herokuapp.com/MM/Schenk_IBMNA/c-46.tar.gz
```

- **human_gene2:** 14.340x14.340 18.068.388nnz (60MB)
``` sh
wget https://suitesparse-collection-website.herokuapp.com/MM/Belcastro/human_gene2.tar.gz
```
  