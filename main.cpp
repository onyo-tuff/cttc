#include <iostream>
#include <memory>
#include <chrono>

#include "Mesh.h"
#include "SparseMatrix.h"
#include "Solver.h"
#include "Operators.h"
#include "print_vec.h"

int main() {
	//Mesh mesh = new Mesh();
	Mesh mesh;
	//bool a;
	std::cout << mesh.N_vol << "\n";
	//a = mesh.hasIdentifier(n, "Obj_Bottom");
	Laplacian Lap(mesh);
	std::cout << "laplacian ok" << "\n";


	// Extraer A, b
	std::vector<double> b;
	SparseMatrix A;
	b = Lap.OprB();
	A = Lap.OprA();
	//std::cout << Lap.vBound_Type[0] << "\n";
	//std::cout << Lap.vBound_Type[1] << "\n";
	//std::cout << Lap.vBound_Type[2] << "\n";
	//std::cout << Lap.vBound_Type[3] << "\n";
	//std::cout << Lap.vBound_Type[4] << "\n";
	//std::cout << Lap.vBound_Type[5] << "\n";
	//std::cout << Lap.vBound_Type[6] << "\n";
	//std::cout << Lap.vBound_Type[7] << "\n";

	//SparseMatrix I(mesh.N_vol);
	//std::cout << "b:" << "\n";
	//for (int i = 0; i < b.size(); i++) {
	//	I.writeEntry(i, i, 1);
	//	std::cout << b[i] << ", ";
	//	if (A.getEntry(i,i) < 0.000001)
	//		std::cout << "zero in diagonal " << i << "!" << "\n";
	//}
	//std::cout << "\n";

	// for (int i = 0; i < b.size(); i++) {
	//	std::cout << "i = " << i << " : " << A.getEntry(i, i) << "\n";
	//}
	//std::cout << "\n";
	//std::cout << "diag(A):" << "\n";
	for (int i = 0; i < b.size(); i++) {
		std::cout << A.getEntry(i,i) << ", ";
	}
	std::cout << "\n";

	// Solucionar
	Solver sol(A, b);
	// for (int i=0; i<b.size(); i++)
	//	std::cout << b[i] << "\n";

	auto start = std::chrono::high_resolution_clock::now();
	sol.solve();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Total program runtime: " << elapsed.count() << " seconds" << std::endl;
	sol.printNumIterations();
	//std::cout<<sol.solverType;

	// Extraer solucion
	std::vector<double> x;
	x = sol.getSolution();
	int i,j;
	double N, M;
	N = mesh.getN();
	M = mesh.getM();
	std::vector<std::vector<double>> x2D(N, std::vector<double>(M));
	for (int n = 0; n < b.size(); n++) {
		i = (int)(n / M);
		j = n - i * M;
		x2D[i][j] = x[n];
	}
	//std::ofstream file("DirichletNeumann.csv");
	//for (int i = x2D.size() - 1; i >= 0; i--) {
	//	for (size_t j = 0; j < x2D[i].size(); j++) {
	//		file << x2D[i][j] << (j < x2D[i].size() - 1 ? "," : "\n");
	//	}
	//}

	// Pasar a VTK
	Mesh *mesh_ptr = &mesh;
	print_vec(mesh_ptr, x2D, "test.vtk");




	return 0;
}