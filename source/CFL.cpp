#include "CFL.h"

CFL::CFL() 
	: dt{1.0}
{
	//Read input file containing the selected scheme
	readScheme();
}

void CFL::readScheme() {
	//Local variables
	std::ifstream inputFile{"inputs/TransientControl"};
	bool read = false;
	while(inputFile>>integrationSch) {
		if(integrationSch=="----------------------------------------")
			read = true;
		if(read&&integrationSch=="IntegrationScheme") {
			inputFile>>integrationSch;
			break;
		}
	}
}

void CFL::cflDiff(Mesh& mesh) {
	// CHeck integrationSch is ok
	std::cout << "Integration scheme: " << integrationSch << " \n";
	//Local variables
	double f, dx, dy, Dt;
	double N = mesh.N_vol;
	double f_max = 0.0;
	// Find critical volume
	for (int k=1; k < N; k++) {
		dx = mesh.getS_n(k);
		dy = mesh.getS_e(k);
		f = mesh.getK(k)*(1.0/dx/dx + 1.0/dy/dy);
		if (f > f_max) 
			f_max = f;
	}
	// Compute Dt according to the integration scheme
	if (integrationSch=="Euler") {
		Dt = 1.0/2.0/f_max;
	} else if (integrationSch=="AB2") {
		Dt = 1.0/4.0/f_max;
	} else {
		std::cout << "No integration scheme selected. Defaulting to forward Euler. \n";
		Dt = 1.0/2.0/f_max;
	}
	// Set this Dt if its the minimum CFL condition
	if (Dt<dt)
		dt = Dt;
	// Check value (testing)
	std::cout << "dt_CFL_diff = " << dt << "\n";
}
