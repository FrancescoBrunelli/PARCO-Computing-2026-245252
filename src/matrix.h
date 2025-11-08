#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <iostream>
#include <ctime>
#include <immintrin.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <array>
#include <sstream>
#include <algorithm>

using namespace std;

extern int row;	// Number of rows
extern int col;	// Number of columns
extern int nnz;	// Number of non-zero elements
extern string object, format, datatype, storage;		// Banner elements in the matrix file

int randomNumber(const int& min, const int& max);

void initM(float* mat);

void initV(float* v);

void COO(const float* mat, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal);

void CSR(const float* mat, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal);

void SELL(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, int C);

void MVmul(const float* mat, const float* v, float* out);

void COOmul(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, const float* v, float* out);

void CSRmul(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, const float* v, float* out);

void P_CSRmul(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, const float* v, float* out);

void printM(const float* mat);

void printV(const float* v, int dim);

void printVector(const vector<float>& v);

void printVector(const vector<int>& v);

void COOtoCSRmap(const multimap<array<int, 2>, float>& COOmap, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal);

void COOtoCSRvect(vector<tuple<int, int, float>>& COOvect, const string& filename, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal);

// COO to CSR using vector
void COOtoCSR(vector<int>& aRow);

void fetch_general_matrix_vector(fstream& file, vector<tuple<int, int, float>>& COOvect);

void fetch_symmetric_matrix_vector(fstream& file, vector<tuple<int, int, float>>& COOvect);

void read_banner(fstream& file);

float* printFullMatrix(const multimap<array<int, 2>, float>& COOmap);

void fetch_matrix(const string &filename, multimap<array<int, 2>, float>& COOmap);

void fetch_general_matrix(fstream& file, multimap<array<int, 2>, float>& COOmap);

void fetch_symmetric_matrix(fstream& file, multimap<array<int, 2>, float>& COOmap);

#endif //__MATRIX_H__
