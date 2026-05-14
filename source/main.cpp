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
	//bool read = false;
	//std::string meshType, strInput;
	//std::ifstream select_mesh{"inputs/SelectMesh"};
	//while (select_mesh>>strInput) {
	//	if (strInput=="----------------------------------------")
	//		read = true;
	//	if (strInput=="Mesh" && read==true)
	//		select_mesh>>meshType;
	//}
	//if (meshType=="Obstacle") {
	//	mesh = Mesh(1);
	//}
	//else if (meshType!="NoObstacle") {
	//	std::cout << "Unrecognized mesh type." << "\n";
	//	return 1;
	//}
	std::cout << "Number of control volumes: " << mesh.N_vol << "\n";
	
	// Compute operators
	Diffusivity Lap(mesh);
	std::vector<double> b;
	SparseMatrix A;
	b = Lap.OprB();
	A = Lap.OprA();
	
	//for (int i = 0; i < A.getN(); i++) {
	//	for (int j = 0; j < A.getN(); j++) {
	//		if(A.getEntry(i,j)!=0)
	//			std::cout  << "1 ";
	//		else
	//			std::cout << "0 ";
	//	}
	//	std::cout << "\n";
	//}
	
	//std::cout << b.size() << "\n";
	
	//for (int j = 0; j < A.getN(); j++) {
	//		std::cout << b[j] << "\n";
	//}
	
	
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
	//print(mesh_ptr, sol.getSolution(), "T");
	printVTR(mesh_ptr, sol.getSolution(), "T");

	return 0;
}
