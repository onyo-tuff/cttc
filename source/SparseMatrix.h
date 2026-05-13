#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>

class SparseMatrix {
public:
	// Constructors
	SparseMatrix() {} // To create empty matrix
	SparseMatrix(int N_in) : N(N_in) , L(N_in) , Ids(N_in) {} // To create NxN empty initialization

	// Public Functions
	double getEntry(int i, int j);
	void writeEntry(int i, int j, double value);
	std::vector<int> getIds(int i) { return Ids[i]; }
	int getN() { return N; }

private:
	// Members
	int N;
	std::vector<std::vector<double>> L;
	std::vector<std::vector<int>> Ids;

	// Private Member Functions

};

#endif


