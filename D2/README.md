# PARCO-Computing-2026-245252
This deliverable evaluates distributed-memory SpMV using MPI across three decomposition strategies (1D block, 1D cyclic, 2D Cartesian block) on a 72-core Intel Xeon cluster. Strong scaling experiments on 5 real SuiteSparse matrices and weak scaling experiments on synthetic matrices show that 1D cyclic outperforms block strategies on irregular matrices due to better load balance, while block decompositions achieve better performance on regular matrices. 1D block decomposition achieves >50x speedup on human_gene2 (strong scaling, 32 processes); 2D block peaks at ~18 GFLOPS (weak scaling, 64 processes), including superlinear speedup attributed to improved cache utilization in parallel execution. Communication time exhibits latency-dominated behaviour at high process counts.

(See the [technical report](Brunelli_245252_D2.pdf) for full results)

## Table of Contents

- [Get Started](#get-started)
- [Compiler version and flags](#compiler-version-and-flags)
- [How to Compile](#how-to-compile)
- [How to Run](#how-to-run)
- [Input and Output](#input-and-output)
- [Change Parameters and Default values](#change-parameters-and-default-values)
- [Cluster-specific notes](#cluster-specific-notes)

---

## Get Started
First, clone this repository to your local machine:
``` sh
git clone https://github.com/FrancescoBrunelli/PARCO-Computing-2026-245252.git
```

After cloning, visit the [set up scripts](scripts/README.md#setup-scripts) section.

## Compiler version and flags
- Compiler: mpicxx (GCC 9.1.0)
- Flags: -std=c++11 -O3

## How to Compile
To compile the program type:
``` sh
# First, make sure to load all necessary modules
module load gcc91
module load mpich-3.2.1--gcc-9.1.0

mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out
```

## How to Run
### Run (locally):
``` sh
mpirun -np <nprocs> ./src/main.out <matrix.mtx> <mode (from 0 to 3)>
```

### Run (on the cluster)
To run the program on the cluster you can either:
#### Open an interactive session
``` sh
# Example:
qsub -I -q queueName -l select=1:ncpus=64: mem=64gb:walltime=01:00:00
```

And then run it in the same way you do it locally:
``` sh
# Example with 2 processes and MODE:1
mpirun -np 2 ./src/main.out matrix.mtx 1
```

OR
#### Run it in batch mode (using a .pbs script)
By creating a .pbs script. An example follows:
``` pbs
#!/bin/bash  
# Job name:  
#PBS -N SpMV 
# Output files:  
#PBS -o ./SpMV.out
# Error files:  
#PBS -e ./SpMV.err
# Queue name
#PBS -q short_cpuQ
# Set the maximum wall time
#PBS -l walltime=0:10:00
# Set the memory, cpus and resources
#PBS -l select=2:ncpus=64:mpiprocs=64:mem=64gb

# Load modules
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
# Select the working directory
cd "$PBS_O_WORKDIR"
# Compile and run
mpicxx -std=c++11 -O3 src/main.cpp src/matrix.cpp src/mpiSpMV.cpp -o src/main.out
mpirun -np 2 ./src/main.out matrix.mtx 1
```

And then type:
``` sh
qsub scriptName.pbs
```
### Input and Output
#### Input
The program will take as an input a Matrix Market sparse matrix (.mtx) and the decomposition mode to adopt. All inputs are provided via terminal as showed in the following example:
``` sh
mpirun -np 2 ./src/main.out matrix.mtx 1
```

#### Output
The output produced by the program is a `.csv` file containing:
- Computation time (expressed in milliseconds)
- Average number of nonzero (NNZ) elements among MPI-processes
- Maximum number of nonzero (NNZ) elements among MPI-processes
- Minimum number of nonzero (NNZ) elements among MPI-processes
- Number of floating point operations per second (FLOPS) (expressed in GFLOPS)
- Communication time (expressed in milliseconds)
- Peak memory footprint among MPI processes (expressed in MB)

## Change Parameters and Default values
### Default values
- By default, without modifying the source files, the program will submit 21 jobs: one for each combination of number of MPI processes (from 2 to 128 in powers of two) and modes (from 0 to 2).
- By default, the `.pbs` scripts will submit jobs to the `short_cpuQ` with: 128 CPUs, 128GB of RAM, 128 MPI processes and a walltime of 6 hours

#### Change decomposition function
To change the decomposition function adopted by the program is sufficient to call it using a different MODE parameter:
- MODE = 0 --> 1D Block Decomposition
- MODE = 1 --> 1D Cyclic Decomposition
- MODE = 2 --> 2D Cartesian Block Decomposition
- MODE = 3 --> Sequential execution

#### Change Memory Requested
To change the amount of RAM requested by the job when running in batch mode (via `qsub`), you need to edit the `.pbs` script and modify the line `#PBS -l select=2:ncpus=64:mpiprocs=64:mem=64gb` with the amount of memory you want, for example, to request 128 GB:
``` pbs
#PBS -l select=2:ncpus=64:mpiprocs=64:mem=128gb
```

#### Change Number of CPUs Requested
To change the amount of CPUs requested by the job when running in batch mode (via qsub), you need to edit the .pbs script and modify the line `#PBS -l select=2:ncpus=64:mpiprocs=64:mem=64gb` with the amount of CPUs you want, for example, to request 32 CPUs:
``` pbs
#PBS -l select=2:ncpus=32:mpiprocs=32:mem=64gb
# Remember to adjust the number of mpi processes accordingly
```
### Cluster-specific notes
**Modules to load:**
``` sh
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
```

**Recommended queue:** `short_cpuQ`

**Restrictions:**
- **Max concurrent jobs per user in queue:** 30
- **Max concurrent jobs per user in execution:** 30
- **Max wall time:** 6h
