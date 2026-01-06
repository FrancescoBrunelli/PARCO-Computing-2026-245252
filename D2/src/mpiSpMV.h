#ifndef D2_MPISPMV_H
#define D2_MPISPMV_H

long get_rss_kb(void);

int sequential_execution(int argc, char** argv);

int MPI_1D_Partitioning(int argc, char** argv);

int MPI_1D_CyclingPartitioning(int argc, char** argv);

int MPI_2D_Partitioning(int argc, char** argv);

#endif //D2_MPISPMV_H
