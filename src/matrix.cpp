#include "matrix.h"

int row;	// Number of rows
int col;	// Number of columns
int nnz;	// Number of non-zero elements
string object, format, datatype, storage;		// Banner elements in the matrix file


int randomNumber(const int& min, const int& max) {
	return rand() % (max - min + 1) + min;
}

void initM(float* mat) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++) {
			if (randomNumber(1, 10) % 2 == 0) {
				mat[i * col + j] = 0;
			} else {
				mat[i * col + j] = randomNumber(1, 10);
			}
		}
	}
}

void initV(float* v) {
	for(int i = 0; i < col; i++) {
		v[i] = randomNumber(2, 2);
	}
}

void COO(const float* mat, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++) {
			if(mat[i * col + j] != 0) {
				aRow.push_back(i);
				aCol.push_back(j);
				aVal.push_back(mat[i * col + j]);
			}
		}
	}
}

void CSR(const float* mat, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	int count = 0;
	aRow.push_back(0);
	for(int i = 0; i < row; i ++) {
		for(int j = 0; j < col; j++) {
			if(mat[i * col + j] != 0) {
				count++;
				aCol.push_back(j);
				aVal.push_back(mat[i * col + j]);
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
	for(int i = 0; i < row; i++) {
		temp = 0;
		for(int j = 0; j < col; j++) {
			temp += mat[i * col + j] * v[j];
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
	for(int i = 0; i < aRow.size() - 1; i++) {
		temp = 0;
		for(int j = aRow[i]; j < aRow[i+1]; j++) {
			temp += aVal[j] * v[aCol[j]];
		}
		out[i] = temp;
	};
}

void P_CSRmul(const vector<int>& aRow, const vector<int>& aCol, const vector<float>& aVal, const float* v, float* out) {
	#pragma omp parallel for schedule(dynamic)
	//pragma omp parallel for schedule(guided)
	for(int i = 0; i < aRow.size() - 1; i++) {
		int temp = 0;
		for(int j = aRow[i]; j < aRow[i+1]; j++) {
			temp += aVal[j] * v[aCol[j]];
		}
		out[i] = temp;
	};
}

void printM(const float* mat) {
	for(int i = 0; i < row; i++) {
		cout << "[\t";
		for(int j = 0; j < col; j++) {
			cout << mat[i * col + j] << "\t";
		}
		cout << "]\n";
	}
	cout << endl;
}

void printV(const float* v, int dim) {
	for(int i = 0; i < dim; i++) {
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

void COOtoCSRmap(const multimap<array<int, 2>, float>& COOmap, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	aRow.reserve(row + 1);
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

void COOtoCSRvect(vector<tuple<int, int, float>>& COOvect, const string& filename, vector<int>& aRow, vector<int>& aCol, vector<float>& aVal) {
	fstream file(filename);
	read_banner(file);
	if(storage == "general") {
		fetch_general_matrix_vector(file, COOvect);
	}
	else {
		fetch_symmetric_matrix_vector(file, COOvect);
	}
	sort(COOvect.begin(), COOvect.end(), [] (const tuple<int, int ,float>& a, const tuple<int, int ,float>& b) {
		return (get<0>(a) < get<0>(b) || get<0>(a) == get<0>(b) && get<1>(a) < get<1>(b));
	});
	aRow.reserve(row + 1);
	aCol.reserve(COOvect.size());
	aVal.reserve(COOvect.size());
	for (const auto& item: COOvect) {
		aRow.push_back(get<0>(item));
		aCol.push_back(get<1>(item));
		aVal.push_back(get<2>(item));
	}
	nnz = aRow.size();
	COOtoCSR(aRow);
}

void fetch_general_matrix_vector(fstream& file, vector<tuple<int, int, float>>& COOvect) {
	while(file.peek() == '%') {
		file.ignore(2048, '\n');
	}
	file >> row >> col >> nnz;
	COOvect.reserve(nnz);
	int r, c;
	double real, imaginary;		// double in order to avoid some problems when reading very small numbers (<1e-45)
	if(datatype == "real" || datatype == "integer") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real;
			r--; c--;
			COOvect.push_back(make_tuple(r, c, (float) real));
		}
	} else if(datatype == "pattern") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c;
			r--; c--;
			COOvect.push_back(make_tuple(r, c, 1.0f));
		}
	} else {	// it's complex
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real >> imaginary;
			r--; c--;
			COOvect.push_back(make_tuple(r, c, (float) real));
		}
	}
}

void fetch_symmetric_matrix_vector(fstream& file, vector<tuple<int, int, float>>& COOvect) {
	while(file.peek() == '%') {
		file.ignore(2048, '\n');
	}
	file >> row >> col >> nnz;
	COOvect.reserve(2 * nnz);
	int r, c;
	double real, imaginary;		// double in order to avoid some problems when reading very small numbers (<1e-45)
	if(datatype == "real" || datatype == "integer") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real;
			r--; c--;
			COOvect.push_back(make_tuple(r, c, (float) real));
			if(r != c) {
				COOvect.push_back(make_tuple(c, r, (float) real));		// insert symmetric element
			}
		}
	} else if(datatype == "pattern") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c;
			r--; c--;
			COOvect.push_back(make_tuple(r, c, 1.0f));
			if(r != c) {
				COOvect.push_back(make_tuple(c, r, 1.0f));		// insert symmetric element
			}
		}
	} else {	// it's complex
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real >> imaginary;
			r--; c--;
			COOvect.push_back(make_tuple(r, c, (float) real));
			if(r != c) {
				COOvect.push_back(make_tuple(c, r, (float) real));		// insert symmetric element
			}
		}
	}
}

// COO to CSR using vector
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

void read_banner(fstream& file) {
	string line;
	getline(file, line);
	istringstream iss(line);
	string banner;
    iss >> banner >> object >> format >> datatype >> storage;
	/*
	cout << "banner: " << banner << endl;
	cout << "object: " << object << endl;
	cout << "format: " << format << endl;
	cout << "datatype: " << datatype << endl;
	cout << "storage: " << storage << endl;
	*/
    if(banner != "%%MatrixMarket" || object != "matrix" || format != "coordinate" || storage == "hermitian" || storage == "skew-symmetric") {
    	cout << "Matrix format not supported" << endl;
    	exit(1);
    }
}

void fetch_matrix(const string& filename, multimap<array<int, 2>, float>& COOmap) {
	fstream file(filename);
	read_banner(file);
	if(storage == "general") {
		fetch_general_matrix(file, COOmap);
	}
	else {
		fetch_symmetric_matrix(file, COOmap);
	}
}

void fetch_general_matrix(fstream& file, multimap<array<int, 2>, float>& COOmap) {
	while(file.peek() == '%') {
		file.ignore(2048, '\n');
	}
	file >> row >> col >> nnz;
	int r, c;
	double real, imaginary;		// double in order to avoid some problems when reading very small numbers (<1e-45)
	if(datatype == "real" || datatype == "integer") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real;
			r--; c--;
			COOmap.insert({{r, c}, (float) real});
		}
	} else if(datatype == "pattern") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c;
			r--; c--;
			COOmap.insert({{r, c}, 1.0f});
		}
	} else {	// it's complex
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real >> imaginary;
			r--; c--;
			COOmap.insert({{r, c}, (float) real});
		}
	}
}

void fetch_symmetric_matrix(fstream& file, multimap<array<int, 2>, float>& COOmap) {
	while(file.peek() == '%') {
		file.ignore(2048, '\n');
	}
	file >> row >> col >> nnz;
	int r, c;
	double real, imaginary;		// double in order to avoid some problems when reading very small numbers (<1e-45)
	if(datatype == "real" || datatype == "integer") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real;
			r--; c--;
			COOmap.insert({{r, c}, (float) real});
			if(r != c) {
				COOmap.insert({{c, r}, (float) real});		// insert symmetric element
			}
		}
	} else if(datatype == "pattern") {
		for(int i = 0; i < nnz; i++) {
			file >> r >> c;
			r--; c--;
			COOmap.insert({{r, c}, 1.0f});
			if(r != c) {
				COOmap.insert({{c, r}, 1.0f});		// insert symmetric element
			}
		}
	} else {	// it's complex
		for(int i = 0; i < nnz; i++) {
			file >> r >> c >> real >> imaginary;
			r--; c--;
			COOmap.insert({{r, c}, (float) real});
			if(r != c) {
				COOmap.insert({{c, r}, (float) real});		// insert symmetric element
			}
		}
	}
}

float* printFullMatrix(const multimap<array<int, 2>, float>& COOmap) {
	float* mat = new float[row * col]();

	for(const auto& item: COOmap) {
		mat[item.first[0] * col + item.first[1]] = item.second;
	}

	for(int i = 0; i < row; i++) {
		cout << "[";
		for(int j = 0; j < col; j++) {
			cout << "\t" << mat[i * col + j];
		}
		cout << "\t]\n";
	}
	//delete[] mat;
	return mat;
}
