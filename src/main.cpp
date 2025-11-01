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

	// -------- Fetch Matrix using multimap
	multimap<array<int, 2>, float> COOmap;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	fetch_matrix(filename, COOmap);

	//float* mat = printFullMatrix(COOmap);

	// PRINT CHECK
	for(const auto& item: COOmap) {
		cout << "(" << item.first[0] + 1 << ", " << item.first[1] + 1 << "): " << item.second << endl;
	}
	
	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;
	COOtoCSRmap(COOmap, aRow, aCol, aVal);
	clock_gettime(CLOCK_MONOTONIC, &t1);
	
	/*---
	// Fetch Matrix using vectors
	vector<int> aRow;
	vector<int> aCol;
	vector<float> aVal;
	fetch_matrix(filename, aRow, aCol, aVal);		// Read .mtx file and get the COO format of the sparse matrix
	
	COOtoCSR(aRow);
	---*/
	
	cout << "Rows: " << row << ", Cols: " << col << ", nnz: " << nnz << endl;
	double start = toMilliseconds(t0);
	double end = toMilliseconds(t1);
	cout << "CSR elapsed time: " << end - start << " ms" << endl;

	/*
	// PRINT CHECK FOR CSR CORRECTNESS
	cout << "--- CSR with multimap ---" << endl;
	cout << "aRow: ";
	printVector(aRow);
	cout << "aCol: "; 
	printVector(aCol);
	cout << "aVal: ";
	printVector(aVal);
	*/

	float* v = new float[col];
	float* out = new float[col];
	initV(v);

	clock_gettime(CLOCK_MONOTONIC, &t0);		// Get start time
	CSRmul(aRow, aCol, aVal, v, out);
	clock_gettime(CLOCK_MONOTONIC, &t1);		// Get end time
	
	start = toMilliseconds(t0);
	end = toMilliseconds(t1);
	cout << "Mul elapsed time: " << end - start << " ms" << endl;
	
	/*
	for(int i = 0; i < nnz; i++) {
		cout << "(" << aRow[i] << ", " << aCol[i] << ") = " << aVal[i] << endl;
	}
	*/
	
	
	
	//float* mat = new float[ROW * COL];
	//float* v = new float[COL];
	//initM(mat);
	//printM(mat);

	/*
	//	CHECK FOR COO
	vector<int> aRow1;
	vector<int> aCol1;
	vector<float> aVal1;
	cout << "--- COO: ---" << endl;
	COO(mat, aRow1, aCol1, aVal1);
	printVector(aRow1);
	printVector(aCol1);
	printVector(aVal1);
	
	// CHECK FOR CSR
	vector<int> aRow2;
	vector<int> aCol2;
	vector<float> aVal2;
	cout << "--- CSR: ---" << endl;
	CSR(mat, aRow2, aCol2, aVal2);
	printVector(aRow2);
	printVector(aCol2);
	printVector(aVal2);
	
	// CHECK FOR COO to CSR
	nnz = aVal1.size();
	row = ROW;
	col = COL;
	cout << "--- COO to CSR: ---" << endl;
	vector<int> aRow3 = aRow1;
	vector<int> aCol3 = aCol1;
	vector<float> aVal3 = aVal1;
	COOtoCSR(aRow3);
	printVector(aRow3);
	printVector(aCol3);
	printVector(aVal3);
	*/

	//delete[] mat;
	delete[] v;
	delete[] out;
	/*
	// CHECK FOR MV MUL
	float* out = new float[COL];
	initV(v);
	cout << "v:" << endl;
	printV(v);
	cout << "--- MV: ---" << endl;
	MVmul(mat, v, out);
	printV(out);
	
	cout << "--- CSR MV: ---" << endl;
	float* out2 = new float[COL];
	CSRmul(aRow2, aCol2, aVal2, v, out2);
	printV(out2);
	
	cout << "--- COO MV: ---" << endl;
	float* out1 = new float[COL];
	COOmul(aRow1, aCol1, aVal1, v, out1);
	printV(out1);
	
	cout << "--- post-COOtoCSR MV: ---" << endl;
	float* out3 = new float[COL];
	CSRmul(aRow3, aCol3, aVal3, v, out3);
	printV(out3);
	*/
	return 0;
}
