#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <time.h>
#include <immintrin.h>
#include <vector>

using namespace std;

const int ROW = 5;
const int COL = 5;

double toSeconds(const timespec& t) {
	return (double) t.tv_sec + (double) t.tv_nsec / 1e9;
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
	for(int i = 0; i < COL; i++) {
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
	for(int i = 0; i < COL; i++) {
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

int main(int argc, char** argv) {
	srand(time(NULL));
	timespec t0, t1;
	float* mat = new float[ROW * COL];
	float* v = new float[COL];
	float* out = new float[COL];
	initM(mat);
	initV(v);
	printM(mat);
	printV(v);
	clock_gettime(CLOCK_MONOTONIC, &t0);
	clock_gettime(CLOCK_MONOTONIC, &t1);
	double start = toSeconds(t0);
	double end = toSeconds(t1);
	cout << "elapsed time: " << end - start << " seconds" << endl;
	
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
	
	/*
	// CHECK FOR MV MUL
	cout << "--- MV: ---" << endl;
	MVmul(mat, v, out);
	printV(out);
	
	cout << "--- CSR MV: ---" << endl;
	CSRmul(aRow2, aCol2, aVal2, v, out);
	printV(out);
	
	cout << "--- COO MV: ---" << endl;
	COOmul(aRow1, aCol1, aVal1, v, out);
	printV(out);
	*/
	
	return 0;
}
