# PARCO Computing a.y. 2025/2026

This project evaluates the performance of Sparse Matrix-Vector Multiplication (SpMV) on a 72-core Intel Xeon cluster across two parallelization paradigms:

* **D1** — Shared-memory parallelism with OpenMP: comparison of static, dynamic and guided scheduling strategies across different thread counts and memory configurations
* **D2** — Distributed-memory parallelism with MPI: comparison of 1D block, 1D cyclic and 2D Cartesian block decomposition strategies across up to 128 processes

Both deliverables include technical reports, automation scripts and reproducible experiments on 5 real matrices from the SuiteSparse collection.
