#include <iostream>
#include <ctime>
#include <immintrin.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <array>
#include "matrix.h"

using namespace std;

double toMilliseconds(const timespec& t) {
	return (double) t.tv_sec * 1000 + (double) t.tv_nsec / 1e6;
}

int main(int argc, char** argv) {
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
	COOtoCSRvect(COOvect, filename, aRow, aCol, aVal);

	// -------- Fetch Matrix using multimap
	//multimap<array<int, 2>, float> COOmap;
	//fetch_matrix(filename, COOmap);
	//COOtoCSRmap(COOmap, aRow, aCol, aVal);

	//float* mat = printFullMatrix(COOmap);

	/*
	cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
	double start = toMilliseconds(t0);
	double end = toMilliseconds(t1);
	cout << "CSR elapsed time: " << end - start << " ms" << endl;
	*/

	float* v = new float[col];
	float* out = new float[row];
	initV(v);
	//cout << "--- V: ---" << endl;
	//printV(v, col);

	clock_gettime(CLOCK_MONOTONIC, &t0);		// Get start time
	//CSRmul(aRow, aCol, aVal, v, out);
	P_CSRmul(aRow, aCol, aVal, v, out);
	clock_gettime(CLOCK_MONOTONIC, &t1);		// Get end time

	 //cout << "--- OUT: ---" << endl;
	 //printV(out, row);

	double start = toMilliseconds(t0);
	double end = toMilliseconds(t1);
	cout << end - start << endl;
	//cout << "Mul elapsed time: " << end - start << " ms" << endl;

	//delete[] mat;
	delete[] v;
	delete[] out;
	return 0;
}
