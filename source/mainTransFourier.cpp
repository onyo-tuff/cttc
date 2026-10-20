#include <iostream>
#include <memory>
#include <chrono>

#include "Mesh.h"
#include "SparseMatrix.h"
#include "Solver.h"
#include "Operators.h"
#include "print.h"
#include "TransientFourier.h"
#include "CFL.h"

//void initialField(std::vector<double>& T_0); 
int main() { // Transient Fourier  
	// Meshing
	Mesh mesh;
	std::cout << "Number of control volumes: " << mesh.N_vol << "\n";

	// CFL TESTING
	//CFL cfl; // Test creation
	//cfl.cflDiff(mesh); // Test diffussive cfl application
	// Orden de magnitud ok

	// Generate initial field
	//std::vector<double> T_n(mesh.N_vol);
	//initialField(T_n);

	// Compute diffusivity operator
	Diffusivity Lap(mesh);
	//std::cout << "ok \n";

	// Apply transient algorithm	
	TransientFourier transient(mesh);
	transient.execute(mesh,Lap);

}

// Function to read input uniform T_0 field
/*
void initialField(std::vector<double>& T_0) {
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
*/

