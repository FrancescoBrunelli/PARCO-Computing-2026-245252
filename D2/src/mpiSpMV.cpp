#include <iostream>
#include <immintrin.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <array>
#include <float.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <mpi.h>
#include "matrix.h"
#include "mpiSpMV.h"

using namespace std;

int sequential_execution(int argc, char** argv) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
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

		/*
		float* mat = printFullMatrix(COOvect);

		cout << "---- COO: ----" << endl;
		cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
		cout << "--- aRow: ---" << endl;
		printVector(aRow);
		cout << "--- aCol: ---" << endl;
		printVector(aCol);
		cout << "--- aVal: ---" << endl;
		printVector(aVal);
		*/

		COOtoCSR(aRow);

		/*
		cout << "---- CSR: ----" << endl;
		cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
		cout << "--- aRow: ---" << endl;
		printVector(aRow);
		cout << "--- aCol: ---" << endl;
		printVector(aCol);
		cout << "--- aVal: ---" << endl;
		printVector(aVal);
		*/

		/*
		double start = toMilliseconds(t0);
		double end = toMilliseconds(t1);
		cout << "CSR elapsed time: " << end - start << " ms" << endl;
		*/

		float* v = new float[col];
		float* out = new float[row];
		initV(v);
		/*
		cout << "--- V: ---" << endl;
		printV(v, col);
		*/

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

		//delete[] mat;
		delete[] v;
		delete[] out;
	}
	return 0;
}

int MPI_1D_Partitioning(int argc, char** argv) {
	double t_start, t_end, time, comm_time = 0, total_time;
	int minNNZ, maxNNZ;
	double max;
	int n;		// number of nnz in a certain chunk
	int chunk_size;
	int local_row = 0;
	int local_col = 0;
	int local_nnz = 0;
	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;
	float* v = nullptr;
	float* out = nullptr;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	MPI_Comm workers_comm;
	int color;
	int COOaRow0 = MPI_UNDEFINED;

	if (size < 2) {
		if (rank == 0)
			fprintf(stderr, "Run with at least 2 processes.\n");
		MPI_Finalize();
		return 1;
	}

	if (rank == 0) {
		color = MPI_UNDEFINED;
	}
	else {
		color = 1;
	}

	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &workers_comm);

	// ---- MASTER CODE ----
	if (rank == 0) {
		//cout << "--- MPI EXECUTION [1D] ---" << endl;
		// Get input from command line
		string filename = argv[1];

		// -------- Fetch Matrix using vector
		vector<tuple<int, int, float>> COOvect;
		fetch_COO(COOvect, filename, aRow, aCol, aVal);

		v = new float[col];
		out = new float[row]{};
		initV(v);
	}

	t_start = MPI_Wtime();
	// Provide each process the number of rows
	MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the number of cols
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the number of nnz
	MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	if (rank != 0) {
		v = new float[col];			// Initialize v for all workers
		out = new float[row]{};
	}

	t_start = MPI_Wtime();
	// Provide each process the randomly generated vector
	MPI_Bcast(v, col, MPI_FLOAT, 0, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	if (rank == 0) {
		// Compute chunksize and send a chunk to each process
		int spare_rows = row % (size - 1);
		chunk_size = (row - spare_rows) / (size - 1);		// number of rows per chunk
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
			t_start = MPI_Wtime();
			MPI_Send(&tmp, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			comm_time += MPI_Wtime() - t_start;

			n = count_if(aRow.begin() + current_index, aRow.end(),  [first, last] (const int& val) {
				return val >= first && val < last;
			});

			t_start = MPI_Wtime();
			MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);	// Send message size
			MPI_Send(aRow.data() + current_index, n, MPI_INT, i, 2, MPI_COMM_WORLD);		// Send aRow
			MPI_Send(aCol.data() + current_index, n, MPI_INT, i, 3, MPI_COMM_WORLD);		// Send aCol
			MPI_Send(aVal.data() + current_index, n, MPI_FLOAT, i, 4, MPI_COMM_WORLD);	// Send aVal
			comm_time += MPI_Wtime() - t_start;

			current_index += n;
			first = last;
		}
	} else {
		t_start = MPI_Wtime();
		// Receive number of rows the process has to work on
		MPI_Recv(&local_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		comm_time += MPI_Wtime() - t_start;
		//out = new float[local_row];

		t_start = MPI_Wtime();
		MPI_Recv(&local_nnz, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	// Receive message size
		comm_time += MPI_Wtime() - t_start;
		aRow.resize(local_nnz);
		aCol.resize(local_nnz);
		aVal.resize(local_nnz);

		t_start = MPI_Wtime();
		// Receive COO-chunk
		MPI_Recv(aRow.data(), local_nnz, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);		// Receive aRow
		MPI_Recv(aCol.data(), local_nnz, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);		// Receive aCol
		MPI_Recv(aVal.data(), local_nnz, MPI_FLOAT, 0, 4, MPI_COMM_WORLD, &status);	// Receive aVal
		comm_time += MPI_Wtime() - t_start;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/*
	if (rank == 0) {
		printf("Setup finished\n");
	}
	*/

	// --- CSR + SpMV Computation
	if (rank != 0 && local_nnz > 0) {
		COOaRow0 = aRow.front();
		COOtoCSR(aRow);
	}

	if (rank != 0) {
		MPI_Barrier(workers_comm);		// Wait all workers are ready for SpMV
	}

	if (rank == 0 || local_nnz == 0) {
		time = -DBL_MAX;
	} else {
		t_start = MPI_Wtime();
		//CSRmul(aRow, aCol, aVal, v, out);
		PartialCSRmul(COOaRow0, aRow, aCol, aVal, v, out);
		t_end = MPI_Wtime();
		time = t_end - t_start;
	}

	t_start = MPI_Wtime();
	MPI_Allreduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	/*
	MPI_Allreduce(&time, &total_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (rank == 0 || local_nnz == 0) {
		time = DBL_MAX;
	}
	MPI_Allreduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	*/

	MPI_Allreduce(&local_nnz, &maxNNZ, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);		// Compute max nnz among processes

	if (local_nnz <= 0) {
		local_nnz = INT_MAX;
	}
	MPI_Allreduce(&local_nnz, &minNNZ, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);		// Compute min nnz among processes

	float* output = new float[row]{};
	t_start = MPI_Wtime();
	MPI_Allreduce(out, output, row, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	/*
	if (rank == 0) {
		// printf("Execution time: %f\n", max * 1000);
		// printf("AVG: %d\n", (nnz / (size - 1)));
		// printf("MAX: %d\n", maxNNZ);
		// printf("MIN: %d\n", minNNZ);
		// printf("GFLOPS: %f\n", (2.0 * nnz / max) / 1e9);
		cout << "**** OUT: ****" << endl;
		printV(output, row);
	}
	*/

	MPI_Allreduce(&comm_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("%f,%d,%d,%d,%f,%f\n", max * 1000, (nnz / (size - 1)), maxNNZ, minNNZ, (2.0 * nnz / max) / 1e9, total_time);
	}

	delete[] v;
	delete[] out;
	return 0;
}

int MPI_1D_CyclingPartitioning(int argc, char** argv) {
	double t_start, t_end, time, comm_time = 0, total_time;
	int minNNZ, maxNNZ;
	double max = -DBL_MAX;
	int n;		// number of nnz in a certain chunk
	int local_nnz = 0;
	vector<tuple<int, int, float>> COOvect;
	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;
	float* v = nullptr;
	float* out = nullptr;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int COOaRow0 = MPI_UNDEFINED;

	if (size < 2) {
		if (rank == 0)
			fprintf(stderr, "Run with at least 2 processes.\n");
		MPI_Finalize();
		return 1;
	}

	// ---- MASTER CODE ----
	if (rank == 0) {
		//cout << "--- MPI EXECUTION [1D-Cyclic] ---" << endl;
		// Get input from command line
		string filename = argv[1];

		// -------- Fetch Matrix using vector
		fetch_COO(COOvect, filename, aRow, aCol, aVal);

		v = new float[col];
		initV(v);
	}

	t_start = MPI_Wtime();
	// Provide each process the number of rows
	MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the number of cols
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the number of nnz
	MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	if (rank != 0) {
		v = new float[col];			// Initialize v for all workers
	}

	out = new float[row]{};

	t_start = MPI_Wtime();
	// Provide each process the randomly generated vector
	MPI_Bcast(v, col, MPI_FLOAT, 0, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	if (rank == 0) {
		// Compute chunksize and send a chunk to each process
		vector<int> aRowB[size];
		vector<int> aColB[size];
		vector<float> aValB[size];
		for_each(COOvect.begin(), COOvect.end(), [size, &aRowB, &aColB, &aValB] (const tuple<int, int, float>& val) {
			int r = get<0>(val);
			int c = get<1>(val);
			int worker_rank = r % size;
			aRowB[worker_rank].push_back(r);
			aColB[worker_rank].push_back(c);
			aValB[worker_rank].push_back(get<2>(val));
		});

		for (int i = 1; i < size; i++) {
			n = aRowB[i].size();

			t_start = MPI_Wtime();
			MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);	// Send worker number of nnz to receive
			MPI_Send(aRowB[i].data(), n, MPI_INT, i, 1, MPI_COMM_WORLD);		// Send aRow
			MPI_Send(aColB[i].data(), n, MPI_INT, i, 2, MPI_COMM_WORLD);		// Send aCol
			MPI_Send(aValB[i].data(), n, MPI_FLOAT, i, 3, MPI_COMM_WORLD);		// Send aVal
			comm_time += MPI_Wtime() - t_start;
		}

		local_nnz = aRowB[rank].size();
		aRow = aRowB[rank];
		aCol = aColB[rank];
		aVal = aValB[rank];

	} else {
		t_start = MPI_Wtime();
		MPI_Recv(&local_nnz, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);		// Receive number of nnz inside the chunk
		comm_time += MPI_Wtime() - t_start;
		aRow.resize(local_nnz);
		aCol.resize(local_nnz);
		aVal.resize(local_nnz);

		t_start = MPI_Wtime();
		MPI_Recv(aRow.data(), local_nnz, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	// Receive chunk aRow

		MPI_Recv(aCol.data(), local_nnz, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);	// Receive chunk aCol
		MPI_Recv(aVal.data(), local_nnz, MPI_FLOAT, 0, 3, MPI_COMM_WORLD, &status);	// Receive chunk aVal
		comm_time += MPI_Wtime() - t_start;
	}

	/*
	if (rank == 0) {
		printf("Setup finished\n");
	}
	*/

	// --- CSR + SpMV Computation
	if (local_nnz > 0) {
		COOaRow0 = aRow.front();
		COOtoCSR(aRow);
	}

	MPI_Barrier(MPI_COMM_WORLD);		// Wait all processes are ready for SpMV

	if (local_nnz == 0) {
		time = -DBL_MAX;
	} else {
		t_start = MPI_Wtime();
		PartialCSRmul(COOaRow0, aRow, aCol, aVal, v, out);
		t_end = MPI_Wtime();
		time = t_end - t_start;
	}

	t_start = MPI_Wtime();
	MPI_Allreduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	//MPI_Allreduce(&time, &total_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	/*
	if (local_nnz == 0) {
		time = DBL_MAX;
	}
	MPI_Allreduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	*/

	MPI_Allreduce(&local_nnz, &maxNNZ, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);		// Compute max nnz among processes

	if (local_nnz <= 0) {
		local_nnz = INT_MAX;
	}
	MPI_Allreduce(&local_nnz, &minNNZ, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);		// Compute min nnz among processes


	float* output = new float[row]{};
	t_start = MPI_Wtime();
	MPI_Allreduce(out, output, row, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	/*
	if (rank == 0) {
		// printf("Execution time: %f\n", max * 1000);
		// printf("AVG: %d\n", (nnz / size));
		// printf("MAX: %d\n", maxNNZ);
		// printf("MIN: %d\n", minNNZ);
		// printf("GFLOPS: %f\n", (2.0 * nnz / max) / 1e9);
		cout << "**** OUT: ****" << endl;
		printV(output, row);
	}
	*/

	MPI_Allreduce(&comm_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("%f,%d,%d,%d,%f,%f\n", max * 1000, (nnz / size), maxNNZ, minNNZ, (2.0 * nnz / max) / 1e9, total_time);
	}

	delete[] v;
	delete[] out;
	return 0;
}

int MPI_2D_Partitioning(int argc, char** argv) {
	double t_start, t_end, time, comm_time = 0, total_time;
	int local_row = 0;
	int local_col = 0;
	int local_nnz = 0;
	int minNNZ, maxNNZ;
	double max;
	int n;		// number of nnz in a certain chunk
	int row_size, col_size;
	vector<tuple<int, int, float>> COOvect;
	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;
	float* v = nullptr;
	float* out = nullptr;

	int rank, size;
	int recv_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	MPI_Comm workers_comm;
	int color;
	int COOaRow0 = MPI_UNDEFINED;

	if (size < 2) {
		if (rank == 0)
			fprintf(stderr, "Run with at least 2 processes.\n");
		MPI_Finalize();
		return 1;
	}

	if (rank == 0) {
		color = MPI_UNDEFINED;
	}
	else {
		color = 1;
	}

	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &workers_comm);

	int dims[2] = {0, 0};
	MPI_Dims_create(size - 1, 2, dims);	// Get number of row and col chunks
	MPI_Comm cart_comm = MPI_COMM_NULL;
	if (rank != 0) {
		int periods[2] = {0, 0};
		int reorder = 1;
		MPI_Cart_create(workers_comm, 2, dims, periods, reorder, &cart_comm);
	}
	int coords[2];
	int cart_rank;
	if (cart_comm != MPI_COMM_NULL) {
		MPI_Comm_rank(cart_comm, &cart_rank);
		MPI_Cart_coords(cart_comm, cart_rank, 2, coords);	// Each worker process gets its coordinates
	}

	// ---- MASTER CODE ----
	if (rank == 0) {
		//cout << "--- MPI EXECUTION [2D] ---" << endl;
		// Get input from command line
		string filename = argv[1];

		// -------- Fetch Matrix using vector
		fetch_COO(COOvect, filename, aRow, aCol, aVal);

		v = new float[col];
		out = new float[row]{};
		initV(v);
	}

	t_start = MPI_Wtime();
	// Provide each process the total number of rows
	MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the total number of cols
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Provide each process the total number of nnz
	MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	if (rank != 0) {
		v = new float[col];			// Initialize v for all workers
		out = new float[row]{};		// Create and initialize with 0s the output array
	}

	t_start = MPI_Wtime();
	// Provide each process the randomly generated vector
	MPI_Bcast(v, col, MPI_FLOAT, 0, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	MPI_Group world_group, cart_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	if (rank != 0) {
		MPI_Comm_group(cart_comm, &cart_group);
	}

	// Setup
	vector<int> cart_ranks(dims[0] * dims[1]);
	if (rank != 0) {
		if (cart_rank == 0) {
			int tmp_coords[2];
			int tmp_rank;
			for (int i = 0; i < dims[0]; i++) {
				for (int j = 0; j < dims[1]; j++) {
					tmp_coords[0] = i; tmp_coords[1] = j;
					MPI_Cart_rank(cart_comm, tmp_coords, &tmp_rank);
					int world_rank;
					MPI_Group_translate_ranks(cart_group, 1, &tmp_rank, world_group, &world_rank);
					cart_ranks[i * dims[1] + j] = world_rank;
				}
			}
			t_start = MPI_Wtime();
			MPI_Send(cart_ranks.data(), size - 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			comm_time += MPI_Wtime() - t_start;
		}
	}
	if (rank == 0) {
		//printf("Dims: (%d, %d)\n", dims[0], dims[1]);
		t_start = MPI_Wtime();
		MPI_Recv(cart_ranks.data(), size - 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);		// Receive cart_comm -> world_comm ranks mapping
		comm_time += MPI_Wtime() - t_start;
		int spare_rows = row % dims[0];
		int spare_cols = col % dims[1];
		row_size = (row - spare_rows) / dims[0];		// Number of rows per chunk
		col_size = (col - spare_cols) / dims[1];		// Number of cols per chunk
		int first_row = 0;	// first row of the chunk
		int last_row;		// last row of the chunk
		int first_col = 0;	// first col of the chunk
		int last_col;		// last col of the chunk
		int tmp;

		for (int i = 0; i < dims[0]; i++) {		// For: number of row chunks
			if (i < spare_rows) {
				last_row = first_row + row_size + 1;
			} else {
				last_row = first_row + row_size;
			}
			for (int j = 0; j < dims[1]; j++) {	// For: number of col chunks
				vector<int> aRowB;
				vector<int> aColB;
				vector<float> aValB;
				if (j < spare_cols) {
					last_col = first_col + col_size + 1;
				} else {
					last_col = first_col + col_size;
				}

				recv_rank = cart_ranks[i * dims[1] + j];

				for_each(COOvect.begin(), COOvect.end(), [first_row, first_col, last_row, last_col, &aRowB, &aColB, &aValB] (const tuple<int, int, float>& val) {
					int r = get<0>(val);
					int c = get<1>(val);
					if (r >= first_row && r < last_row && c >= first_col && c < last_col) {
						aRowB.push_back(r);
						aColB.push_back(c);
						aValB.push_back(get<2>(val));
					}
				});

				n = aValB.size();

				// Get and send number of rows and cols the i,j process has to work on
				tmp = last_row - first_row;
				t_start = MPI_Wtime();
				MPI_Send(&tmp, 1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD);
				comm_time += MPI_Wtime() - t_start;
				tmp = last_col - first_col;
				t_start = MPI_Wtime();
				MPI_Send(&tmp, 1, MPI_INT, recv_rank, 1, MPI_COMM_WORLD);

				MPI_Send(&n, 1, MPI_INT, recv_rank, 2, MPI_COMM_WORLD);		// Send number of nnz elements inside the chunk
				MPI_Send(aRowB.data(), n, MPI_INT, recv_rank, 3, MPI_COMM_WORLD);		// Send chunk aRow
				MPI_Send(aColB.data(), n, MPI_INT, recv_rank, 4, MPI_COMM_WORLD);		// Send chunk aCol
				MPI_Send(aValB.data(), n, MPI_FLOAT, recv_rank, 5, MPI_COMM_WORLD);		// Send chunk aVal
				comm_time += MPI_Wtime() - t_start;
				first_col = last_col;
			}
			first_col = 0;
			first_row = last_row;
		}
	} else {
		t_start = MPI_Wtime();
		MPI_Recv(&local_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);		// Receive number of rows of the chunk
		MPI_Recv(&local_col, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);		// Receive number of cols of the chunk
		MPI_Recv(&local_nnz, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);		// Receive number of nnz inside the chunk
		comm_time += MPI_Wtime() - t_start;
		aRow.resize(local_nnz);
		aCol.resize(local_nnz);
		aVal.resize(local_nnz);

		t_start = MPI_Wtime();
		MPI_Recv(aRow.data(), local_nnz, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);	// Receive chunk aRow

		MPI_Recv(aCol.data(), local_nnz, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);	// Receive chunk aCol
		MPI_Recv(aVal.data(), local_nnz, MPI_FLOAT, 0, 5, MPI_COMM_WORLD, &status);	// Receive chunk aVal
		comm_time += MPI_Wtime() - t_start;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/*
	if (rank == 0) {
		printf("Setup Finished\n");
	}
	*/
	/*
	// DEBUG PRINT
	if (rank != 0) {
		for (int i = 0; i < dims[0]; i++) {
			for (int j = 0; j < dims[1]; j++) {
				if (i == coords[0] && j == coords[1]) {
					printf("\t--- World Rank: %d with coordinates (%d, %d)\n", rank, coords[0], coords[1]);
					cout << "--- aRow: ---" << endl;
					printVector(aRow);
					cout << "--- aCol: ---" << endl;
					printVector(aCol);
					cout << "--- aVal: ---" << endl;
					printVector(aVal);
					fflush(stdout);
				}
				MPI_Barrier(cart_comm);
			}
		}
	}
	*/

	// --- CSR + SpMV Computation
	if (rank != 0 && local_nnz > 0) {
		COOaRow0 = aRow.front();
		COOtoCSR(aRow);
	}

	if (rank != 0) {
		MPI_Barrier(workers_comm);		// Wait all workers are ready for SpMV
	}

	if (rank == 0 || local_nnz == 0) {
		time = -DBL_MAX;
	} else if (local_nnz > 0) {
		t_start = MPI_Wtime();
		PartialCSRmul(COOaRow0, aRow, aCol, aVal, v, out);
		t_end = MPI_Wtime();
		time = t_end - t_start;
	}

	t_start = MPI_Wtime();
	MPI_Allreduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	/*
	MPI_Allreduce(&time, &total_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (rank == 0 || local_nnz == 0) {
		time = DBL_MAX;
	}
	MPI_Allreduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	*/


	/*
	// Print Check
	if (rank != 0) {
		for (int i = 0; i < dims[0]; i++) {
			for (int j = 0; j < dims[1]; j++) {
				if (i == coords[0] && j == coords[1] && nnz > 0) {
					printf("\t--- World Rank: %d with coordinates (%d, %d)\n", rank, coords[0], coords[1]);
					cout << "--- CSR aRow: ---" << endl;
					printVector(aRow);
					cout << "--- OUT: ---" << endl;
					printV(out, row);
					fflush(stdout);
				}
				MPI_Barrier(cart_comm);
			}
		}
	}
	*/

	MPI_Allreduce(&local_nnz, &maxNNZ, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);		// Compute max nnz among processes

	if (local_nnz <= 0) {
		local_nnz = INT_MAX;
	}
	MPI_Allreduce(&local_nnz, &minNNZ, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);		// Compute min nnz among processes

	float* output = new float[row]{};
	t_start = MPI_Wtime();
	MPI_Allreduce(out, output, row, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	comm_time += MPI_Wtime() - t_start;

	/*
	if (rank == 0) {
		// printf("Execution time: %f\n", max * 1000);
		// printf("AVG: %d\n", (nnz / (size - 1)));
		// printf("MAX: %d\n", maxNNZ);
		// printf("MIN: %d\n", minNNZ);
		// printf("GFLOPS: %f\n", (2.0 * nnz / max) / 1e9);
		cout << "**** OUT: ****" << endl;
		printV(output, row);
	}
	*/

	MPI_Allreduce(&comm_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if (rank == 0) {
		//printf("NPROCS = %d", size);
		printf("%f,%d,%d,%d,%f,%f\n", max * 1000, (nnz / (size - 1)), maxNNZ, minNNZ, (2.0 * nnz / max) / 1e9, total_time);
	}

	delete[] v;
	delete[] out;
	return 0;
}
