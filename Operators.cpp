#include "Operators.h"

Laplacian::Laplacian(Mesh& mesh)
	: A(mesh.N_vol-mesh.N_obstacle), b(mesh.N_vol-mesh.N_obstacle)
{
	std::vector<double> xn;
	int N = mesh.N;
	int M = mesh.M;
	getInputs();
	k = 0;
	for (int n = 0; n < mesh.N_vol; n++) {
		//std::cout << n << "\n";
		if (mesh.hasIdentifier(n,"Obstacle"))
			k++;
		else {
		aP = 0;
		aE = 0;
		aW = 0;
		aN = 0;
		aS = 0;
		bP = 0;
		// For all the inner volumes
		std::vector<int> ij = mesh.ij_given_n(n);
		xn = mesh.getCenter(n);
		if (ij[0]!=0)
			aW = kWall(mesh,n-M,n)*mesh.getS_w(n)/dist(xn,mesh.getCenter(n-M));
		if (ij[0]!=mesh.N-1)
			aE = kWall(mesh,n,n+M)*mesh.getS_e(n)/dist(xn,mesh.getCenter(n+M));
		if (ij[1]!=0)
			aS = kWall(mesh,n-1,n)*mesh.getS_s(n)/dist(xn,mesh.getCenter(n-1));
		if (ij[1]!=mesh.M-1)
			aN = kWall(mesh,n,n+1)*mesh.getS_n(n)/dist(xn,mesh.getCenter(n+1));
		//std::cout << "aW, aE, aS, aN ok" << "\n";
		// For boundary volumes
		calc_boundary(mesh, n);
		//std::cout << "boundary ok" << "\n";
		// 
		aP += -aE -aW -aN -aS;

		A.writeEntry(n-k,n-k,aP);
		if (aW!=0)
			A.writeEntry(n-k,n-M-k,aW);
		if (aE!=0)
			A.writeEntry(n-k,n+M-k,aE);
		if (aS!=0)
			A.writeEntry(n-k,n-1-k,aS);
		if (aN!=0)
			A.writeEntry(n-k,n+1-k,aN);
		b[n-k] = bP;
		}
		//std::cout << "A,b entries ok" << "\n";
	}
	// Check for zeros in diagonal
	for (size_t i = 0; i < b.size(); i++) {
		if (A.getEntry(i,i) == 0)
			std::cout << "zero in diagonal " << i << "!" << "\n";
	}
}

void Laplacian::calc_boundary(Mesh& mesh, int n) {
	std::vector<double> xn = mesh.getCenter(n);
	double k = mesh.getK(n);
	for (size_t i = 0; i < vBound_Name.size(); i++) {
		if (mesh.hasIdentifier(n,vBound_Name[i])) {
			if (vBound_Type[i]=="WestBoundary") {
				if (vBoCo_Type[i]=="Dirichlet") {
					aW = 0;
					aP += -k*mesh.getS_w(n)/dist(xn,mesh.getWest(n));
					bP += vBoCo_Value[i]*aP;
				}
				else if (vBoCo_Type[i]=="Neumann")
					bP += k*vBoCo_Value[i]*mesh.getS_w(n);
				else
					std::cout<<"No boundary condition type detected for" << vBound_Name[i] << " boundary." << "\n";
			}
			if (vBound_Type[i]=="EastBoundary") {
				if (vBoCo_Type[i]=="Dirichlet") {
					aE = 0;
					aP += -k*mesh.getS_e(n)/dist(xn,mesh.getEast(n));
					bP += vBoCo_Value[i]*aP;
				}
				else if (vBoCo_Type[i]=="Neumann")
					bP += -k*vBoCo_Value[i]*mesh.getS_e(n);
				else
					std::cout<<"No boundary condition type detected for" << vBound_Name[i] << " boundary." << "\n";
			}
			if (vBound_Type[i]=="SouthBoundary") {
				if (vBoCo_Type[i]=="Dirichlet") {
					aS = 0;
					aP += -k*mesh.getS_s(n)/dist(xn,mesh.getSouth(n));
					bP += vBoCo_Value[i]*aP;
				}
				else if (vBoCo_Type[i]=="Neumann")
					bP += k*vBoCo_Value[i]*mesh.getS_s(n);
				else
					std::cout<<"No boundary condition type detected for" << vBound_Name[i] << " boundary." << "\n";
			}
			if (vBound_Type[i]=="NorthBoundary") {
				if (vBoCo_Type[i]=="Dirichlet") {
					aN = 0;
					aP += -k*mesh.getS_n(n)/dist(xn,mesh.getNorth(n));
					bP += vBoCo_Value[i]*aP;
				}
				else if (vBoCo_Type[i]=="Neumann")
					bP += -k*vBoCo_Value[i]*mesh.getS_s(n);
				else
					std::cout<<"No boundary condition type detected for" << vBound_Name[i] << " boundary." << "\n";
			}
		}
	}
}

void Laplacian::getInputs() {
	bool read = false;
	while (bocos_input>>strInput) {
		//std::cout << "Reading bocos_input ok" << "\n";
		if (strInput == "----------------------------------------------------")
			read = true;
		if (strInput == ">" && read) {
			bocos_input >> strInput;
			vBound_Name.push_back(strInput);
			bocos_input >> strInput;
			vBoCo_Type.push_back(strInput);
			bocos_input >> strInput;
			vBoCo_Value.push_back(std::stod(strInput));
		}
	}
	read = false;
	vBound_Type = vBound_Name;
	while (bound_input>>strInput) {
		//std::cout << "Reading bocos_input ok" << "\n";
		if (strInput == "----------------------------------")
			read = true;
		if (strInput == ")" && read) {
			bound_input>>strInput;
			for (size_t k = 0; k < vBound_Name.size(); k++) {
				if (vBound_Name[k] == strInput)
					bound_input>>vBound_Type[k];
			}
		}
	}
	for (size_t i = 0; i < 4; i++) {
	std::cout << vBoCo_Type[i]<< " ";
	std::cout << vBoCo_Value[i] << " ";
	std::cout << vBound_Name[i] << " ";
	std::cout << vBound_Type[i] << "\n";
	}
}

double Laplacian::kWall(Mesh& mesh, int n, int p) { //Important!!!! n<p
	std::vector<double> xn = mesh.getCenter(n);
	std::vector<double> xp = mesh.getCenter(p);
	std::vector<double> xwall;
	if (xn[1]==xp[1])
		xwall = mesh.getEast(n);
	else if (xn[0]==xp[0])
		xwall = mesh.getNorth(n);
	else
		std::cout << "Unexpected boundary between x_n and x_p!!!" << "\n";
	double D = dist(xn,xp);
	double dn = dist(xn,xwall);
	double dp = dist(xp,xwall);
	//std::cout<<"Ok before getK" << "\n";
	double kn = mesh.getK(n);
	double kp = mesh.getK(p);
	//std::cout<<"Ok after getK" << "\n";
	return D/(dp/kp + dn/kn);
}

double Laplacian::dist(std::vector<double> a, std::vector<double> b) {
	std::vector<double> r(a.size());
	for (size_t i = 0; i < a.size(); i++) 
		r[i] = a[i] - b[i];
	return std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
}













