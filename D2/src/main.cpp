#include <iostream>
#include <ctime>
#include <immintrin.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <array>
#include <mpi.h>
#include "matrix.h"

using namespace std;

double toMilliseconds(const timespec& t) {
	return (double) t.tv_sec * 1000 + (double) t.tv_nsec / 1e6;
}

void sequential_execution(int argc, char** argv) {
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

	/*
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

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, size;
	double start, end, total_time;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		cout << "--- SEQUENTIAL EXECUTION ---" << endl;
		sequential_execution(argc, argv);
		cout << endl;
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

		vector<int> aRow;
		vector<int> aCol;
		vector<float> aVal;

		// -------- Fetch Matrix using vector
		vector<tuple<int, int, float>> COOvect;
		fetch_COO(COOvect, filename, aRow, aCol, aVal);


		//float* mat = printFullMatrix(COOvect);

		/*
		cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
		cout << "--- aRow: ---" << endl;
		printVector(aRow);
		cout << "--- aCol: ---" << endl;
		printVector(aCol);
		cout << "--- aVal: ---" << endl;
		printVector(aVal);
		*/

		float* v = new float[col];
		float* out = new float[row];
		initV(v);

		cout << "--- V: ---" << endl;
		printV(v, col);

		//delete[] mat;
		delete[] v;
		delete[] out;
	}

	MPI_Finalize();

	return 0;
}
