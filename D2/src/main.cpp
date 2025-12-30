#include <iostream>
#include <ctime>
#include <immintrin.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <array>
#include <float.h>
#include <mpi.h>
#include "matrix.h"

using namespace std;

double toMilliseconds(const timespec& t) {
	return (double) t.tv_sec * 1000 + (double) t.tv_nsec / 1e6;
}

void sequential_execution(int argc, char** argv) {
	cout << "--- SEQUENTIAL EXECUTION ---" << endl;
	srand(time(NULL));
	timespec t0, t1;
	// Get input from command line
	if (argc < 2) {
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
	string filename = argv[1];

	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;

	// -------- Fetch Matrix using vector
	vector<tuple<int, int, float>> COOvect;
	fetch_COO(COOvect, filename, aRow, aCol, aVal);

	float* mat = printFullMatrix(COOvect);

	cout << "---- COO: ----" << endl;
	cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
	cout << "--- aRow: ---" << endl;
	printVector(aRow);
	cout << "--- aCol: ---" << endl;
	printVector(aCol);
	cout << "--- aVal: ---" << endl;
	printVector(aVal);

	COOtoCSR(aRow);

	cout << "---- CSR: ----" << endl;
	cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
	cout << "--- aRow: ---" << endl;
	printVector(aRow);
	cout << "--- aCol: ---" << endl;
	printVector(aCol);
	cout << "--- aVal: ---" << endl;
	printVector(aVal);

	/*
	double start = toMilliseconds(t0);
	double end = toMilliseconds(t1);
	cout << "CSR elapsed time: " << end - start << " ms" << endl;
	*/

	float* v = new float[col];
	float* out = new float[row];
	initV(v);
	cout << "--- V: ---" << endl;
	printV(v, col);

	clock_gettime(CLOCK_MONOTONIC, &t0);		// Get start time
	CSRmul(aRow, aCol, aVal, v, out);
	//P_CSRmul(aRow, aCol, aVal, v, out);
	clock_gettime(CLOCK_MONOTONIC, &t1);		// Get end time

	cout << "--- OUT: ---" << endl;
	printV(out, row);

	double start = toMilliseconds(t0);
	double end = toMilliseconds(t1);
	//cout << end - start << endl;
	cout << "Mul elapsed time: " << end - start << " ms" << endl;

	delete[] mat;
	delete[] v;
	delete[] out;
}

void MPI_1D_Partitioning(int argc, char** argv) {

}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	double t_start, t_end, time, total_time = 0;
	double min = DBL_MAX;
	double max = -DBL_MAX;
	int n;		// number of nnz in a certain chunk
	int chunk_size;
	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;
	float* v = nullptr;
	float* out = nullptr;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Request requests[5];
	MPI_Status status;

	if (rank == 0) {
		string filename = argv[1];
		cout << "Input matrix: " << filename << endl;
		sequential_execution(argc, argv);
	}

	/*
	start = MPI_Wtime();
	end = MPI_Wtime();
	total_time = end - start;
	*/

	// ---- MASTER CODE ----
	if (rank == 0) {
		cout << "--- MPI EXECUTION ---" << endl;
		// Get input from command line
		if (argc < 2) {
			fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
			exit(1);
		}
		string filename = argv[1];

		// -------- Fetch Matrix using vector
		vector<tuple<int, int, float>> COOvect;
		fetch_COO(COOvect, filename, aRow, aCol, aVal);

		v = new float[col];
		out = new float[row];
		initV(v);
	}

	// Provide each process the number of rows
	MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the number of cols
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the number of nnz elements
	MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank != 0) {
		// Initialize v for all workers
		v = new float[col];
	}

	// Provide each process the randomly generated vector
	MPI_Bcast(v, col, MPI_FLOAT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		/*
		// Provide each process the empty output vector
		MPI_Bcast(out, row, MPI_FLOAT, 0, MPI_COMM_WORLD);
		*/

		// Compute chunksize and send a chunk to each process
			// NOTE: 1D PARTITIONING
		int spare_rows = row % size;
		chunk_size = (row - spare_rows) / size;		// number of rows per chunk
		int first = 0;	// first row of the chunk
		int last;		// last row of the chunk
		int current_index = 0;		// index of the first element to parse
		int tmp;

		for (int i = 1; i < size; i++) {
			if (i <= spare_rows) {
				last = first + chunk_size + 1;
			} else {
				last = first + chunk_size;
			}
			tmp = last - first;		// Actual number of rows for the i-th process
			MPI_Send(&tmp, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

			n = count_if(aRow.begin() + current_index, aRow.end(),  [first, last] (int val) {
				return val >= first && val < last;
			});

			MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);	// Send message size
			MPI_Send(aRow.data() + current_index, n, MPI_INT, i, 2, MPI_COMM_WORLD);		// Send aRow
			MPI_Send(aCol.data() + current_index, n, MPI_INT, i, 3, MPI_COMM_WORLD);		// Send aCol
			MPI_Send(aVal.data() + current_index, n, MPI_FLOAT, i, 4, MPI_COMM_WORLD);	// Send aVal

			current_index += n + 1;
			first = last;
		}
	} else {
		// Receive number of rows the process has to work on
		MPI_Recv(&chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		out = new float[chunk_size];

		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	// Receive message size
		aRow.resize(n);
		aCol.resize(n);
		aVal.resize(n);

		// Receive COO-chunk
		MPI_Recv(aRow.data(), n, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);		// Receive aRow
		MPI_Recv(aCol.data(), n, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);		// Receive aCol
		MPI_Recv(aVal.data(), n, MPI_FLOAT, 0, 4, MPI_COMM_WORLD, &status);	// Receive aVal
	}
		// DEBUG PRINT CHECK
		if (rank == 3) {
			cout << "**** RANK: " << rank << " ****" << endl;
			cout << "rows: " << chunk_size << endl;
			cout << "nnz: " << n << endl;
			cout << "---- COO: ----" << endl;
			cout << "--- aRow: ---" << endl;
			printVector(aRow);
			cout << "--- aCol: ---" << endl;
			printVector(aCol);
			cout << "--- aVal: ---" << endl;
			printVector(aVal);
		}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		printf("Setup finished");
	}

	// --- CSR + SpMV Computation
	if (rank == 0) {
		time = 0;
	} else {
		COOtoCSR(aRow);
		t_start = MPI_Wtime();
		CSRmul(aRow, aCol, aVal, v, out);
		t_end = MPI_Wtime();
		printf("rank: %d has performed SpMV successfully\n", rank);
		time = t_end - t_start;

	}
	printf("-- %d\n", rank);
	MPI_Allreduce(&time, &total_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (rank == 0) {
		time = DBL_MAX;
	}
	MPI_Allreduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	if (rank == 0) {
		time = -DBL_MAX;
	}
	MPI_Allreduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	// Output printing
	bool token = true;
	if (rank == 0) {
		printf("Average time: %f\nMax execution time: %f\nMin execution time: %f\n", total_time / (size - 1), max, min);
		if (size > 1) {
			MPI_Send(&token, 1, MPI_C_BOOL, rank + 1, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(&token, 1, MPI_C_BOOL, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printV(out, chunk_size);
		fflush(stdout);
		if (rank < size - 1) {
			MPI_Send(&token, 1, MPI_C_BOOL, rank + 1, 0, MPI_COMM_WORLD);
		}
	}

	delete[] v;
	delete[] out;

	MPI_Finalize();

	return 0;
}
