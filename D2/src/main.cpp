#include <iostream>
#include <mpi.h>
#include "mpiSpMV.h"

using namespace std;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	//sequential_execution(argc, argv);

	MPI_1D_Partitioning(argc, argv);

	MPI_1D_CyclingPartitioning(argc, argv);

	MPI_2D_Partitioning(argc, argv);

	MPI_Finalize();

	return 0;
}
