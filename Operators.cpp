#include "Operators.h"

Laplacian::Laplacian(Mesh& mesh)
	: A(mesh.N_vol-mesh.N_obstacle), b(mesh.N_vol-mesh.N_obstacle), vBound_Value(8), vBound_Type(8)
{
	getInputs();
	k = 0;
	for (int n = 0; n < mesh.N_vol; n++) {
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
		if (!mesh.hasIdentifier(n,"Front") && !mesh.hasIdentifier(n,"Obj_Back"))
			aW = mesh.getS_w(n)/(mesh.getCenter(n)[0]-mesh.getWest(n)[0]);
		if (!mesh.hasIdentifier(n,"Back") && !mesh.hasIdentifier(n,"Obj_Front"))
			aE = mesh.getS_e(n)/(mesh.getEast(n)[0]-mesh.getCenter(n)[0]);
		if (!mesh.hasIdentifier(n,"Bottom") && !mesh.hasIdentifier(n,"Obj_Top"))
			aS = mesh.getS_s(n)/(mesh.getCenter(n)[1]-mesh.getSouth(n)[1]);
		if (!mesh.hasIdentifier(n,"Top") && !mesh.hasIdentifier(n,"Obj_Bottom"))
			aN = mesh.getS_n(n)/(mesh.getNorth(n)[1]-mesh.getCenter(n)[1]);
		// For volumes inside the obstacle
		if (mesh.hasIdentifier(n,"Obstacle"))
			k++;
		// For boundary volumes
		calc_boundary(mesh, n);
		// 
		aP += -aE -aW -aN -aS;

		A.writeEntry(n-k,n-k,aP);
		if (aW!=0)
			A.writeEntry(n-k,n-mesh.M-k,aW);
		if (aE!=0)
			A.writeEntry(n-k,n+mesh.M-k,aE);
		if (aS!=0)
			A.writeEntry(n-k,n-1-k,aS);
		if (aN!=0)
			A.writeEntry(n-k,n+1-k,aN);
		b[n-k] = bP;
		}
	}
	// Check for zeros in diagonal
	for (int i = 0; i < b.size(); i++) {
		if (A.getEntry(i,i) == 0)
			std::cout << "zero in diagonal " << i << "!" << "\n";
	}
}

void Laplacian::calc_boundary(Mesh& mesh, int n) {
	if (mesh.hasIdentifier(n,"Front")) {
		if (vBound_Type[0]=="Dirichlet") {
			aP += -2*mesh.getS_w(n)/mesh.getS_n(n);
			bP += vBound_Value[0]*aP;
		}
		else if (vBound_Type[0]=="Neumann") 
			bP += vBound_Value[0]*mesh.getS_w(n);
		else 
			std::cout<<"No boundary condition detected for Front patch." << "\n";
	}
	if (mesh.hasIdentifier(n,"Back")) {
		if (vBound_Type[1]=="Dirichlet") {
			aP += -2*mesh.getS_e(n)/mesh.getS_n(n);
			bP += vBound_Value[1]*aP;
		}
		else if (vBound_Type[1]=="Neumann") 
			bP += -vBound_Value[1]*mesh.getS_e(n);
		else 
			std::cout<<"No boundary condition detected for Back patch." << "\n";
	}
	if (mesh.hasIdentifier(n,"Top")) {
		if (vBound_Type[2]=="Dirichlet") {
			aP += -2*mesh.getS_n(n)/mesh.getS_w(n);
			bP += vBound_Value[2]*aP;
		}
		else if (vBound_Type[2]=="Neumann") 
			bP += -vBound_Value[2]*mesh.getS_n(n);
		else 
			std::cout<<"No boundary condition detected for Top patch." << "\n";
	}
	if (mesh.hasIdentifier(n,"Bottom")) {
		if (vBound_Type[3]=="Dirichlet") {
			aP += -2*mesh.getS_s(n)/mesh.getS_w(n);
			bP += vBound_Value[3]*aP;
		}
		else if (vBound_Type[3]=="Neumann") 
			bP += vBound_Value[3]*mesh.getS_s(n);
		else 
			std::cout<<"No boundary condition detected for Bottom patch." << "\n";
	}
	if (mesh.hasIdentifier(n,"Obj_Front")) {
		if (vBound_Type[4]=="Dirichlet") {
			aP += -2*mesh.getS_e(n)/mesh.getS_n(n);
			bP += vBound_Value[4]*aP;
		}
		else if (vBound_Type[4]=="Neumann") 
			bP += -vBound_Value[4]*mesh.getS_e(n);
		else 
			std::cout<<"No boundary condition detected for Obj_Front patch." << "\n";
	}
	if (mesh.hasIdentifier(n,"Obj_Back")) {
		if (vBound_Type[5]=="Dirichlet") {
			aP += -2*mesh.getS_w(n)/mesh.getS_n(n);
			bP += vBound_Value[5]*aP;
		}
		else if (vBound_Type[5]=="Neumann") 
			bP += vBound_Value[5]*mesh.getS_w(n);
		else 
			std::cout<<"No boundary condition detected for Obj_Back patch." << "\n";	
	}
	if (mesh.hasIdentifier(n,"Obj_Top")) {
		if (vBound_Type[6]=="Dirichlet") {
			aP += -2*mesh.getS_s(n)/mesh.getS_w(n);
			bP += vBound_Value[6]*aP;
		}
		else if (vBound_Type[6]=="Neumann") 
			bP += vBound_Value[6]*mesh.getS_s(n);
		else 
			std::cout<<"No boundary condition detected for Obj_Top patch." << "\n";
	}
	if (mesh.hasIdentifier(n,"Obj_Bottom")) {
		if (vBound_Type[7]=="Dirichlet") {
			aP += -2*mesh.getS_n(n)/mesh.getS_w(n);
			bP += vBound_Value[7]*aP;
		}
		else if (vBound_Type[7]=="Neumann") 
			bP += -vBound_Value[7]*mesh.getS_n(n);
		else 
			std::cout<<"No boundary condition detected for Obj_Bottom patch." << "\n";
	}

}

void Laplacian::getInputs() {
	if (!input.is_open()) { std::cerr << "Error: Could not open file!" << std::endl; }
	while (std::getline(input, strInput)) {
		if (strInput == "Front:") {
			std::getline(input, strInput);
			vBound_Type[0] = strInput;
			std::getline(input, strInput);
			vBound_Value[0] = std::stod(strInput);
		}
		if (strInput == "Back:") {
			std::getline(input, strInput);
			vBound_Type[1] = strInput;
			std::getline(input, strInput);
			vBound_Value[1] = std::stod(strInput);
		}
		if (strInput == "Top:") {
			std::getline(input, strInput);
			vBound_Type[2] = strInput;
			std::getline(input, strInput);
			vBound_Value[2] = std::stod(strInput);
		}
		if (strInput == "Bottom:") {
			std::getline(input, strInput);
			vBound_Type[3] = strInput;
			std::getline(input, strInput);
			vBound_Value[3] = std::stod(strInput);
		}
		if (strInput == "Obj_Front:") {
			std::getline(input, strInput);
			vBound_Type[4] = strInput;
			std::getline(input, strInput);
			vBound_Value[4] = std::stod(strInput);
		}
		if (strInput == "Obj_Back:") {
			std::getline(input, strInput);
			vBound_Type[5] = strInput;
			std::getline(input, strInput);
			vBound_Value[5] = std::stod(strInput);
		}
		if (strInput == "Obj_Top:") {
			std::getline(input, strInput);
			vBound_Type[6] = strInput;
			std::getline(input, strInput);
			vBound_Value[6] = std::stod(strInput);
		}
		if (strInput == "Obj_Bottom:") {
			std::getline(input, strInput);
			vBound_Type[7] = strInput;
			std::getline(input, strInput);
			vBound_Value[7] = std::stod(strInput);
		}
	}
}	















