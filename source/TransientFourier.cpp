#include "TransientFourier.h"

TransientFourier::TransientFourier(Mesh& mesh)
	: N{mesh.N_vol}, T_n(N), T_n1(N)
{
	initialField();
	readTimeInput();
	// Temporarily hardcoded
	dt = 1e-6;
}

void TransientFourier::execute(Mesh& mesh, Diffusivity& Diff) {
	sumTime = 0;
	for ( t = 0; t < t_fin ; t += dt ) {
		modDiff = 0;
		// Find solution for next timestep
		std::vector<double> diffT = Diff.applyOp(T_n);
		for (int i = 0; i < N; i++) 
			T_n1[i] = T_n[i] + dt*diffT[i];
		// Print ?
		sumTime += dt;
		if ( sumTime > dt_print ) {
			std::string name = std::to_string(t+dt);
			printVTR(mesh, T_n1, name);
			sumTime = 0;
		}
		// Update
		T_n = T_n1;
		// Print terminal stuff
		for (int i = 0; i < N; i++) 
			modDiff += diffT[i]*diffT[i];
		std::cout << "t = " << t+dt << "  |  modDiff = " << modDiff << "\n";
		// Stop if reached sationary
		if ( modDiff < 1e-6 ) {
			std::string name = std::to_string(t+dt) + "/T";
			printVTR(mesh, T_n1, name);
			break;
		}
	}
}

void TransientFourier::initialField() {
	std::ifstream input{ "inputs/T0" };
	std::string strInput;
	double T;
	while (input>>strInput) {
		if (strInput == "----------------------------------------") { 
			input>>strInput;
			T = std::stod(strInput);
			std::cout << "T read ok \n";
		}
	}
	for (int i = 0; i < N; i++) 
		T_n[i] = T;
}

void TransientFourier::readTimeInput() {
	std::ifstream input{"inputs/TransientControl"};
	std::string strInput;
	bool read = false;
	while (input>>strInput) {
		if (strInput == "----------------------------------------")  
			read = true;
		if (strInput == "t_0" && read == true) {
			input>>strInput;
			t_0 = std::stod(strInput);
		}
		if (strInput == "t_fin" && read == true) {
			input>>strInput;
			t_fin = std::stod(strInput);
		}
		if (strInput == "dt_print" && read == true) {
			input>>strInput;
			dt_print = std::stod(strInput);
		}
	}
}


