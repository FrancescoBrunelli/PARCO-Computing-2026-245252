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

int MPI_1D_Partitioning(int argc, char** argv) {
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

	/*
	start = MPI_Wtime();
	end = MPI_Wtime();
	total_time = end - start;
	*/

	if (size < 2) {
		if (rank == 0)
			fprintf(stderr, "Run with at least 2 processes.\n");
		MPI_Finalize();
		return 1;
	}

	if (rank == 0) {
		sequential_execution(argc, argv);
	}

	// ---- MASTER CODE ----
	if (rank == 0) {
		cout << "--- MPI EXECUTION ---" << endl;
		// Get input from command line
		if (argc < 2) {
			fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		string filename = argv[1];

		// -------- Fetch Matrix using vector
		vector<tuple<int, int, float>> COOvect;
		fetch_COO(COOvect, filename, aRow, aCol, aVal);

		v = new float[col];
		out = new float[row];
		initV(v);
	}

	// Provide each process the number of cols
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
			MPI_Send(&tmp, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

			n = count_if(aRow.begin() + current_index, aRow.end(),  [first, last] (const int& val) {
				return val >= first && val < last;
			});

			MPI_Send(&n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);	// Send message size
			MPI_Send(aRow.data() + current_index, n, MPI_INT, i, 2, MPI_COMM_WORLD);		// Send aRow
			MPI_Send(aCol.data() + current_index, n, MPI_INT, i, 3, MPI_COMM_WORLD);		// Send aCol
			MPI_Send(aVal.data() + current_index, n, MPI_FLOAT, i, 4, MPI_COMM_WORLD);	// Send aVal

			current_index += n;
			first = last;
		}
	} else {
		// Receive number of rows the process has to work on
		MPI_Recv(&chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		row = chunk_size;
		out = new float[chunk_size];

		MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	// Receive message size
		nnz = n;
		printf("rank: %d has nnz = %d\n", rank, nnz);
		aRow.resize(n);
		aCol.resize(n);
		aVal.resize(n);

		// Receive COO-chunk
		MPI_Recv(aRow.data(), n, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);		// Receive aRow
		MPI_Recv(aCol.data(), n, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);		// Receive aCol
		MPI_Recv(aVal.data(), n, MPI_FLOAT, 0, 4, MPI_COMM_WORLD, &status);	// Receive aVal
	}
	// DEBUG PRINT
	if (rank == 1) {
		cout << "COO: " << endl;
		printVector(aRow);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		printf("Setup finished");
	}

	// --- CSR + SpMV Computation
	if (rank == 0) {
		time = 0;
	} else if (nnz > 0) {
		COOtoCSR(aRow);
		printf("rank %d finished CSR conversion\n", rank);
		t_start = MPI_Wtime();
		CSRmul(aRow, aCol, aVal, v, out);
		t_end = MPI_Wtime();
		printf("rank: %d has performed SpMV successfully\n", rank);
		time = t_end - t_start;
	}

	// DEBUG PRINT CHECK
	if (rank == 1) {
		cout << "**** RANK: " << rank << " ****" << endl;
		cout << "rows: " << chunk_size << endl;
		cout << "nnz: " << n << endl;
		cout << "---- CSR: ----" << endl;
		cout << "--- aRow: ---" << endl;
		printVector(aRow);
		cout << "--- aCol: ---" << endl;
		printVector(aCol);
		cout << "--- aVal: ---" << endl;
		printVector(aVal);
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
	return 0;
}

int MPI_2D_Partitioning(int argc, char** argv) {
	double t_start, t_end, time, total_time = 0;
	double min = DBL_MAX;
	double max = -DBL_MAX;
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
	MPI_Request requests[5];
	MPI_Status status;
	MPI_Comm workers_comm;
	int color;

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

	if (rank == 0) {
		sequential_execution(argc, argv);
	}

	int dims[2] = {0, 0};
	MPI_Dims_create(size - 1, 2, dims);	// Get number of row and col chunks
	MPI_Comm cart_comm;
	if (rank != 0) {
		int periods[2] = {0, 0};
		int reorder = 1;
		MPI_Cart_create(workers_comm, 2, dims, periods, reorder, &cart_comm);
	}
	int coords[2];
	int cart_rank;
	if (rank != 0) {
		MPI_Comm_rank(cart_comm, &cart_rank);
		MPI_Cart_coords(cart_comm, cart_rank, 2, coords);	// Each worker process gets its coordinates
	}

	// ---- MASTER CODE ----
	if (rank == 0) {
		cout << "--- MPI EXECUTION ---" << endl;
		// Get input from command line
		if (argc < 2) {
			fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		string filename = argv[1];

		// -------- Fetch Matrix using vector
		fetch_COO(COOvect, filename, aRow, aCol, aVal);

		v = new float[col];
		out = new float[row];
		initV(v);
	}

	// Provide each process the number of cols
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank != 0) {
		// Initialize v for all workers
		v = new float[col];
	}

	// Provide each process the randomly generated vector
	MPI_Bcast(v, col, MPI_FLOAT, 0, MPI_COMM_WORLD);

	MPI_Group world_group, cart_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	MPI_Comm_group(cart_comm, &cart_group);

	// Setup
	if (rank == 0) {
		int spare_rows = row % dims[0];
		int spare_cols = col % dims[1];
		row_size = (row - spare_rows) / dims[0];		// Number of rows per chunk
		col_size = (col - spare_cols) / dims[1];		// Number of cols per chunk
		int first_row = 0;	// first row of the chunk
		int last_row;		// last row of the chunk
		int first_col = 0;	// first col of the chunk
		int last_col;		// last col of the chunk
		int tmp;
		int worker_world_rank;

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

				for_each(COOvect.begin(), COOvect.end(), [first_row, first_col, last_row, last_col, &aRowB, &aColB, &aValB] (const tuple<int, int, float>& val) {
					int r = get<0>(val);
					int c = get<1>(val);
					if (r >= first_row && r <= last_row && c >= first_col && c <= last_col) {
						aRowB.push_back(r);
						aColB.push_back(c);
						aValB.push_back(get<2>(val));
					}
				});

				n = aValB.size();

				/*
				n = count_if(COOvect.begin(), COOvect.end(),  [first_row, first_col, last_row, last_col, &aRowB, &aColB, &aValB] (const tuple<int, int, float>& val) {
					if ((std::get<0>(val) >= first_row && std::get<0>(val) < last_row) && (std::get<1>(val) >= first_col && std::get<1>(val) < last_col)) {
						aRowB.push_back(get<0>(val));
						aColB.push_back(get<1>(val));
						aValB.push_back(get<2>(val));
						return true;
					}
					return false;
				});		// Count number of nnz elements inside the chunk
				*/

				// Get rank of process at coordinate (i, j)
				coords[0] = i; coords[1] = j;
				MPI_Cart_rank(cart_comm, coords, &cart_rank);
				MPI_Group_translate_ranks(cart_group, 1, &cart_rank, world_group, &recv_rank);

				// Get and send number of rows and cols the i,j process has to work on
				tmp = last_row - first_row;
				MPI_Send(&tmp, 1, MPI_INT, recv_rank, 0, MPI_COMM_WORLD);
				tmp = last_col - first_col;
				MPI_Send(&tmp, 1, MPI_INT, recv_rank, 1, MPI_COMM_WORLD);

				MPI_Send(&n, 1, MPI_INT, recv_rank, 2, MPI_COMM_WORLD);		// Send number of nnz elements inside the chunk
				MPI_Send(aRowB.data(), n, MPI_INT, recv_rank, 3, MPI_COMM_WORLD);		// Send chunk aRow
				MPI_Send(aColB.data(), n, MPI_INT, recv_rank, 4, MPI_COMM_WORLD);		// Send chunk aCol
				MPI_Send(aValB.data(), n, MPI_FLOAT, recv_rank, 5, MPI_COMM_WORLD);		// Send chunk aVal


				first_col = last_col;
			}



			first_row = last_row;
		}


	} else {
		MPI_Recv(&row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);		// Receive number of rows of the chunk
		MPI_Recv(&col, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);		// Receive number of cols of the chunk
		MPI_Recv(&nnz, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);		// Receive number of nnz inside the chunk
		aRow.resize(nnz);
		aCol.resize(nnz);
		aVal.resize(nnz);

		MPI_Recv(aRow.data(), nnz, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);	// Receive chunk aRow
		MPI_Recv(aCol.data(), nnz, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);	// Receive chunk aCol
		MPI_Recv(aVal.data(), nnz, MPI_FLOAT, 0, 5, MPI_COMM_WORLD, &status);	// Receive chunk aVal
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		printf("Setup Finished\n");
	}

	// DEBUG PRINT
	if (rank != 0) {
		for (int i = 0; i < dims[0]; i++) {
			for (int j = 0; j < dims[1]; j++) {
				if (i == coords[0] && j == coords[1]) {
					printf("\t--- World Rank: %d with coordinates (%d, %d)", rank, coords[0], coords[1]);
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

	return 0;
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	//MPI_1D_Partitioning(argc, argv);

	MPI_2D_Partitioning(argc, argv);

	MPI_Finalize();

	return 0;
}
