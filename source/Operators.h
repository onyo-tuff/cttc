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

class Diffusivity {
public:
	Diffusivity(Mesh& mesh);
	SparseMatrix OprA() { return A; }
	std::vector<double> OprB() { return b; }

private:
	double kWall(Mesh& mesh, int n, int p);
	double dist(std::vector<double> a, std::vector<double> b);
	void calc_boundary(Mesh& mesh, int n);
	int k, i_Robins;
	double aP, aE, aW, aN, aS, bP;
	SparseMatrix A;
	std::vector<double> b;
	std::vector<double> vBoCo_Value, vBoCo_Value1;
	std::vector<std::string> vBoCo_Type, vBound_Name, vBound_Type;
	std::ifstream bocos_input{ "inputs/BoCos" };
	std::ifstream bound_input{ "inputs/Boundary" };
	void getInputs();
	std::string strInput;
};


#endif
