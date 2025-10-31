#include <iostream>
#include <ctime>
#include <immintrin.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <array>

using namespace std;

const int ROW = 5;
const int COL = 5;

int row;	// Number of rows
int col;	// Number of columns
int nnz;	// Number of non-zero elements

double toMilliseconds(const timespec& t) {
	return (double) t.tv_sec * 1000 + (double) t.tv_nsec / 1e6;
}

int randomNumber() {
	return rand() % 11;
}

void initM(float* mat) {
	for(int i = 0; i < ROW; i++) {
		for(int j = 0; j < COL; j++) {
			if (randomNumber() % 2 == 0) {
				mat[i * COL + j] = 0;
			} else {
				mat[i * COL + j] = randomNumber();
			}
		}
	}
}

void initV(float* v) {
	for(int i = 0; i < col; i++) {
		v[i] = randomNumber();
	}
}

void COO(const float* mat, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	for(int i = 0; i < ROW; i++) {
		for(int j = 0; j < COL; j++) {
			if(mat[i * COL + j] != 0) {
				aRow.push_back(i);
				aCol.push_back(j);
				aVal.push_back(mat[i * COL + j]);
			}
		}
	}
}

void CSR(const float* mat, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	int count = 0;
	aRow.push_back(0);
	for(int i = 0; i < ROW; i ++) {
		for(int j = 0; j < COL; j++) {
			if(mat[i * COL + j] != 0) {
				count++;
				aCol.push_back(j);
				aVal.push_back(mat[i * COL + j]);
			}
		}
		aRow.push_back(count);
	}
}

void SELL(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, int C) {		//Start from the output of the CSR format
	//C is the slice height
	int nnz = aVal.size();		// nnz = number of non zero elements in the matrix
	
}

void MVmul(const float* mat, const float* v, float* out) {
	int temp;
	for(int i = 0; i < ROW; i++) {
		temp = 0;
		for(int j = 0; j < COL; j++) {
			temp += mat[i * COL + j] * v[j];
		}
		out[i] = temp;
	}
}

void COOmul(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, const float* v, float* out) {
	int temp;
	for(int i = 0; i < aRow.size(); i++) {
		temp = 0;
		for(int j = 0; j < aRow.size(); j++) {
			if(aRow[j] == i) {
				temp += aVal[j] * v[aCol[j]];
			}
		}
		out[i] = temp;
	}
}

void CSRmul(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, const float* v, float* out) {
	int temp;
	int j = 0;
	for(int i = 0; i < aRow.size() - 1; i++) {
		temp = 0;
		for(; j < aRow[i+1]; j++) {
			temp += aVal[j] * v[aCol[j]];
		}
		out[i] = temp;
	};
}

void printM(const float* mat) {
	for(int i = 0; i < ROW; i++) {
		cout << "[\t";
		for(int j = 0; j < COL; j++) {
			cout << mat[i * COL + j] << "\t";
		}
		cout << "]\n";
	}
	cout << endl;
}

void printV(const float* v) {
	for(int i = 0; i < col; i++) {
		cout << "[" << v[i] << "] ";
	}
	cout << endl;
}

void printVector(const vector<float>& v) {
	for(const auto& item: v) {
		cout << "[" << item << "] ";
	}
	cout << endl;
}

void printVector(const vector<int>& v) {
	for(const auto& item: v) {
		cout << "[" << item << "] ";
	}
	cout << endl;
}

void fetch_matrix(const string filename, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	fstream file(filename);
	while(file.peek() == '%') {
		file.ignore(2048, '\n');
	}
	file >> row >> col >> nnz;
	aRow.reserve(nnz);
	aCol.reserve(nnz);
	aVal.reserve(nnz);
	int r, c;
	float v;
	for(int i = 0; i < nnz; i++) {
		file >> r >> c >> v;
		r--; c--;
		aRow.push_back(r);
		aCol.push_back(c);
		aVal.push_back(v);
	}
}

void fetch_matrix2(const string filename, multimap<array<int, 2>, float>& COOmap) {
	fstream file(filename);
	while(file.peek() == '%') {
		file.ignore(2048, '\n');
	}
	file >> row >> col >> nnz;
	int r, c;
	float v;
	for(int i = 0; i < nnz; i++) {
		file >> r >> c >> v;
		r--; c--;
		COOmap.insert({{r, c}, v});
	}
}

void COOtoCSRmap(const multimap<array<int, 2>, float>& COOmap, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	aRow.reserve(row);
	aCol.reserve(nnz);
	aVal.reserve(nnz);
	vector<int> tmp;
	int count = 0;
	aRow.push_back(0);
	
	for(const auto& item: COOmap) {
		aCol.push_back(item.first[1]);
		aVal.push_back(item.second);
	}
	
	for(int i = 0; i < row; i++) {
		for(const auto& item: COOmap) {
			if(item.first[0] == i) {
				count++;
			}
		}
		aRow.push_back(count);
	}
}

void COOtoCSR(vector<int>& aRow) {
	vector<int> tmp;
	int count = 0;
	tmp.push_back(0);
	for(int i = 0; i < row; i++) {
		for(int j = count; j < nnz; j++) {
			if(aRow[j] == i) {
				count++;
			}
		}
		tmp.push_back(count);
	}
	aRow = tmp;
}

int main(int argc, char** argv) {
	srand(time(NULL));
	// Get input from command line
	if (argc < 2) {
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
	string filename = argv[1];
	
	timespec t0, t1;
	
	// -------- Fetch Matrix using multimap
	multimap<array<int, 2>, float> COOmap;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	fetch_matrix2(filename, COOmap);
	
	/*
	// PRINT CHECK
	for(const auto& item: COOmap) {
		cout << "(" << item.first[0] << ", " << item.first[1] << "): " << item.second << endl;
	}
	*/
	
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
	cout << "aRow: ";
	printVector(aRow);
	cout << "\naCol: "; 
	printVector(aCol);
	cout << "\naVal: ";
	printVector(aVal);
	*/
	
	float* v = new float[col];
	float* out = new float[col];
	
	clock_gettime(CLOCK_MONOTONIC, &t0);		// Get start time
	CSRmul(aRow, aCol, aVal, v, out);
	clock_gettime(CLOCK_MONOTONIC, &t1);		// Get end time
	
	start = toMilliseconds(t0);
	end = toMilliseconds(t1);
	cout << "Mul elapsed time: " << end - start << " ms" << endl;
	
	return 0;
}
