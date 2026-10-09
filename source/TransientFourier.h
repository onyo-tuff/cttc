#ifndef TRANS_F
#define TRANS_F

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include "Operators.h"
#include "print.h"

class TransientFourier {
public:
	//Constructor
	TransientFourier(Mesh& mesh);

	//Public functions
	void execute(Mesh& mesh, Diffusivity& Diff);

private:
	//Private elements
	int maxSteps, N;
	double t, t_0, t_fin, dt, dt_print, sumTime, modDiff;
	std::vector<double> T_n, T_n1;

	//Private functions
	void initialField(); 
	void readTimeInput();

};

#endif
