#include <iostream>
#include <mpi.h>
#include "mpiSpMV.h"

using namespace std;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	if (argc < 3) {
		fprintf(stderr, "Usage: %s [martix-market-filename] [mode: 0 for 1D-partitioning, 1 for 1D-cyclic partitioning, or 2 for 2D partitioning]\n", argv[0]);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	int mode = stoi(argv[2]);
	switch (mode) {
		case 0: {
			MPI_1D_Partitioning(argc, argv);
			break;
		}
		case 1: {
			MPI_1D_CyclingPartitioning(argc, argv);
			break;
		}
		case 2: {
			MPI_2D_Partitioning(argc, argv);
			break;
		}
		case 3: {
			sequential_execution(argc, argv);
			break;
		}
	}

	MPI_Finalize();
	return 0;
}
