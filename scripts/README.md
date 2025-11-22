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
  
  Alternatively you can download the matrices you prefer from [here](https://sparse.tamu.edu), make sure they are written in any of the different formats compatible with the program [check compatible formats](../src/README.md#compatibility).

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
## Python Scripts

Organized as follows:
### CVS Scripts
Starting from a folder organization like the one of [results](../results/README.md), these scripts are responsible for the creation of the .cvs files that will be used for plotting

- `fetch_runtime_results.py`
  
  Usage: `python3 scripts/fetch_runtime_results.py`
- `fetch_perf_results.py`
  
  Usage: `python3 scripts/fetch_perf_results.py`
- `fetch_cachegrind_results.py`
  
  Usage: `python3 scripts/fetch_cachegrind_results.py`

### Plot Scripts
- `plot_runtime.py`: Fetches all runtime results for a given input matrix name, calculates the 90th percentile and creates a linear plot
  
  Usage: `python3 scripts/plot_runtime.py <matrix_name>`
- `plot_perf.py`: creates an histogram showing the percentages of different types of cache-miss for each scheduling policy given the amount of threads and memory (hardcoded)
  
  Usage: `python3 scripts/plot_perf.py`
- `plot_scaling.py`: creates an histogram showing the (strong) scaling for each scheduling type of a given input matrix name, comparing parallel executions to the sequential one
  
  Usage: `python3 scripts/plot_scaling.py <matrix_name>`

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
  