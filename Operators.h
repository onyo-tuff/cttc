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

private:
	double kWall(Mesh& mesh, int n, int p);
	double dist(std::vector<double> a, std::vector<double> b);
	void calc_boundary(Mesh& mesh, int n);
	int k;
	double aP, aE, aW, aN, aS, bP;
	SparseMatrix A;
	std::vector<double> b;
	std::vector<double> vBoCo_Value;
	std::vector<std::string> vBoCo_Type, vBound_Name, vBound_Type;
	std::ifstream bocos_input{ "Physical_BoCos.txt" };
	std::ifstream bound_input{ "Boundary_Input.txt" };
	void getInputs();
	std::string strInput;
};


#endif
