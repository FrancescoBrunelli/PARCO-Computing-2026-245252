# Scripts

## Table of Contents

- [Introduction](#introduction)
- [Setup Scripts](#setup-scripts)
- [Run Scripts](#run-scripts)
- [Plot Scripts](#plot-scripts)
- [Matrices used](#matrices-used)

---
## Introduction

This folder contains all the scripts needed to set up the environment, run different types of tests, process the data and generate plots based on different performance metrics

## Setup Scripts

Scripts in this section are responsible for setting up the environment and downloading required data

- `download_matrices.sh` - Downloads and extracts five sparse matrices
  Alternatively you can download the matrices you prefer from [here](https://sparse.tamu.edu), make sure they are written in any of the different formats compatible with the program [check compatible formats](../README.md).

## Run Scripts

Scripts used to run multiple tests on the compute nodes automatically
It is important to run them from the main folder of the repository

### `run_time.pbs`
Runs 100 executions of the program saving the runtime of each execution into a .txt file named: <matrix_name>_NT<number_of_threads_used>.
The .txt file is located at results/matrix_name/

**Usage (after program compilation):**
``` sh
# Example usage of run_time.pbs
qsub -v INPUT=matrix.mtx,NTHREADS=NTHREADS scripts/run_time.pbs
```

### `run_perf.pbs`
Runs 10 executions of perf on the program saving the output of each execution into a .txt file named: <matrix_name>_perf_NT<number_of_threads_used>.
The .txt file is located at results/<matrix_name>/

**Usage (after program compilation):**
``` sh
# Example usage of run_perf.pbs
qsub -v INPUT=matrix.mtx,NTHREADS=NTHREADS scripts/run_perf.pbs
```
### `run_cachegrind.pbs`
Runs 1 execution of cachegrind (valgrind) on the program saving the output into a .err file located at results/cachegrind_logs/<matrix_name>/

**Usage (after program compilation):**
``` sh
# Example usage of run_cachegrind.pbs
qsub -v INPUT=matrix.mtx scripts/run_cachegrind.pbs
```
### `run_time_parallel_jobs.sh`
Submits a total of 7 `run_time.pbs` jobs, one for each power of 2 (from 0 to 7) number of threads (from 1 to 64)

**Usage (after program compilation):**
``` sh
# Make sure to make it executable first:
chmod +x ./scripts/run_time_parallel_jobs.sh

# Example usage of run_time_parallel_jobs.sh
./scripts/run_time_parallel_jobs.sh matrix.mtx
```
### `run_perf_parallel_jobs.sh`
Submits a total of 7 `run_perf.pbs` jobs, one for each power of 2 (from 0 to 7) number of threads (from 1 to 64)

**Usage (after program compilation):**
``` sh
# Make sure to make it executable first:
chmod +x ./scripts/run_perf_parallel_jobs.sh

# Example usage of run_perf_parallel_jobs.sh
./scripts/run_perf_parallel_jobs.sh matrix.mtx
```
## Plot Scripts

Scripts used to process the output data and generate performance plots.


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
  