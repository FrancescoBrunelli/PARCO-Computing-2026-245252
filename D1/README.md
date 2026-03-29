# PARCO-Computing-2026-245252
This deliverable evaluates shared-memory SpMV performance using OpenMP across three scheduling strategies (static, dynamic, guided) on a 72-core Intel Xeon cluster. Experiments on 5 real SuiteSparse matrices show that static scheduling (chunksize = 100) consistently outperforms dynamic and guided policies due to lower synchronization overhead, despite their better load balancing. Cache profiling via perf and cachegrind reveals that dynamic scheduling incurs the highest LLC-miss rate while guided achieves the lowest. NUMA-induced memory bandwidth limitations and scheduling overhead are identified as the main performance bottlenecks.

(See the [technical report](Brunelli-245252-D1.pdf) for full results)

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
- Compiler: g++ 9.1.0 (GCC)
- Flags: -std=c++11 -O3 -march=native -ftree-vectorize -fopenmp

## How to Compile
To compile the program type:
``` sh
g++ -std=c++11 -O3 -march=native -ftree-vectorize -fopenmp src/main.cpp src/matrix.cpp -o src/main.out
```

## How to Run
### Run (locally):
``` sh
./src/main.out matrix.mtx
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
./src/main.out matrix.mtx
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
#PBS -l select=1:ncpus=64:mem=64gb

# Load modules
module load gcc91
# Select the working directory
cd "$PBS_O_WORKDIR"
# Compile and run
g++ -std=c++11 -O3 -march=native -ftree-vectorize -fopenmp src/main.cpp src/matrix.cpp -o src/main.out
./src/main.out
```

And then type:
``` sh
qsub scriptName.pbs
```
### Input and Output
#### Input
The program will take as an input a Matrix Market sparse matrix (.mtx) provided via terminal as showed in the following example:
``` sh
./src/main.out matrix.mtx
```

#### Output
The output produced by the program is the time elapsed (expressed in milliseconds) during the multiplication between the CSR format of the matrix and a randomly generated vector

## Change Parameters and Default values
### Default values
- By default, without modifying the source files, the program will compute the SpMV multiplication using a guided scheduling option and using the maximum amount of threads available.
- By default, the .pbs scripts will submit jobs to the `short_cpuQ` with: 64 CPUs, 64GB of RAM and a walltime of 6 hours

#### Change scheduling
To change the OpenMP scheduling policy, modify the line: `#pragma omp parallel for schedule(guided)` inside the function `P_CSRmul` in `src/matrix.cpp`.

#### Change Memory Requested
To change the amount of RAM requested by the job when running in batch mode (via qsub), you need to edit the .pbs script and modify the line `#PBS -l select=1:ncpus=64:mem=64gb` with the amount of memory you want, for example, to request 128 GB:
``` pbs
#PBS -l select=1:ncpus=64:mem=128gb
```

#### Change Number of CPUs Requested
To change the amount of CPUs requested by the job when running in batch mode (via qsub), you need to edit the .pbs script and modify the line `#PBS -l select=1:ncpus=64:mem=64gb` with the amount of CPUs you want, for example, to request 32 CPUs:
``` pbs
#PBS -l select=1:ncpus=32:mem=64gb
```
### Cluster-specific notes
**Modules to load:**
``` sh
module load gcc91
module load perf
```

**Recommended queue:** `short_cpuQ`

**Restrictions:**
- **Max concurrent jobs per user in queue:** 30
- **Max concurrent jobs per user in execution:** 30
- **Max wall time:** 6h
