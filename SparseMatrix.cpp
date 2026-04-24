#include "SparseMatrix.h"

double SparseMatrix::getEntry(int i, int j) {
	std::vector<double> a = L[i];
	std::vector<int> ids = Ids[i];
	int n = (int)a.size();
	for (int k = 0; k < n; k++) {
		if (ids[k] == j) {
			return a[k];
		}
	}
	return 0;
}

void SparseMatrix::writeEntry(int i, int j, double value) {
	bool rewrite = false;
	std::vector<double> a = L[i];
	std::vector<int> ids = Ids[i];
	int n = a.size();
	for (int k = 0; k < n; k++) {
		if (ids[k] == j) {
			rewrite = true;
			L[i][k] = value;
		}
	}
	if (!rewrite) {
		L[i].push_back(value);
		Ids[i].push_back(j);
	}
}