#ifndef OPS_H
#define OPS_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include "SparseMatrix.h"
#include "Mesh.h"

class Laplacian {
public:
	Laplacian(Mesh& mesh);
	SparseMatrix OprA() { return A; }
	std::vector<double> OprB() { return b; }

public:
	void calc_boundary(Mesh& mesh, int n);
	int k;
	double aP, aE, aW, aN, aS, bP;
	SparseMatrix A;
	std::vector<double> b;
	std::vector<double> vBound_Value;
	std::vector<std::string> vBound_Type;
	std::ifstream input{ "Physical_BoCos.txt" };
	void getInputs();
	std::string strInput;
};


#endif
