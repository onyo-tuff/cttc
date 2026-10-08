#include <iostream>
#include <memory>
#include <chrono>

#include "Mesh.h"
#include "SparseMatrix.h"
#include "Solver.h"
#include "Operators.h"
#include "print.h"

/*
int main() { // Stationary Fourier

	// Meshing
	Mesh mesh;
	std::cout << "Number of control volumes: " << mesh.N_vol << "\n";
	
	// Compute operators
	Diffusivity Lap(mesh);
	std::vector<double> b;
	SparseMatrix A;
	b = Lap.OprB();
	A = Lap.OprA();
	
	// Solve
	Solver sol(A, b);
	auto start = std::chrono::high_resolution_clock::now();
	sol.solve();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Total program runtime: " << elapsed.count() << " seconds" << std::endl;
	sol.printNumIterations();

	// Print to vtk
	Mesh *mesh_ptr = &mesh;
	printVTR(mesh_ptr, sol.getSolution(), "T");

	return 0;
}
*/

void initialField(std::vector<double> T_0); 
int main() { // Fourier Transitorio 
	// Meshing
	Mesh mesh;
	Mesh *mesh_ptr = &mesh;
	std::cout << "Number of control volumes: " << mesh.N_vol << "\n";

	// Generate initial field
	std::vector<double> T_n(mesh.N_vol);
	initialField(T_n);

	// Compute diffusivity operator
	Diffusivity Lap(mesh);

	// Time loop
	int k = 0;
	double t_fin = 0.1;
	double dt = 1e-6; // Cumple CFL 
	double diffSum = 0;
	std::vector<double> T_n1(mesh.N_vol);
	for (double t = 0; t <= t_fin; t+=dt) {
		k++;
		// Compute next time step
		std::vector<double> diffT = Lap.applyOp(T_n);
		for (int i = 0; i < T_n1.size(); i++) {
			T_n1[i] = T_n[i] + dt*diffT[i];
		}
		// Test
		diffSum = 0;
		for (int i = 0; i < diffT.size(); i++)
			diffSum+=diffT[i];
		std::cout << "diffSum = " << diffSum << "  |  ";
		std::cout << "t = " << t+dt << "\n";
		// Print solutions
		if (k%1000==0)
			printVTR(mesh_ptr, T_n1, std::to_string(k));
		// Update n
		T_n = T_n1;
		// Test time loop
		//std::cout << "t = " << t << "\n";
		// Break if stability
    if (std::abs(diffSum) < 1e-6)
			break;
	}
}

// Function to read input uniform T_0 field
void initialField(std::vector<double> T_0) {
	std::ifstream input{ "inputs/T0" };
	std::string strInput;
	double T;
	while (input>>strInput) {
		if (strInput == "----------------------------------------") { 
			input>>strInput;
			T = std::stod(strInput);
			// Test
			//std::cout << "T_0 = " << T << "\n";
		}
	}
	for (int i = 0; i < T_0.size(); i++) 
		T_0[i] = T;
}


