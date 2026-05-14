#include <iostream>
#include <memory>
#include <chrono>

#include "Mesh.h"
#include "SparseMatrix.h"
#include "Solver.h"
#include "Operators.h"
#include "print.h"

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
	print(mesh_ptr, sol.getSolution(), "T.vtk");

	return 0;
}
