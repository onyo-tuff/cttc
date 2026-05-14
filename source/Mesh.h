#ifndef MESH_H
#define MESH_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

class Mesh {
public:
	// Constructors 
	Mesh(); // For the case of continuous domain
	Mesh(bool rectangle); // For the case with rectangular discontinuiny inside the domain

	// Functions
	int getN() { return N; }
	int getM() { return M; }
	int node_number(int i, int j);
	std::vector<int> ij_given_n(int n); // Indices de la malla colocada
	// Primary mesh
	 // - access with node number
	double getK(int n) { return MaterialK[n]; }
	std::vector<double> getCenter(int n) { return Coords[n][0]; }
	std::vector<double> getWest(int n) { return Coords[n][1]; }
	std::vector<double> getEast(int n) { return Coords[n][2]; }
	std::vector<double> getSouth(int n) { return Coords[n][3]; }
	std::vector<double> getNorth(int n) { return Coords[n][4]; }
	double getVol(int n) { return Vols_n_Surfaces[n][0]; }
	double getS_w(int n) { return Vols_n_Surfaces[n][1]; }
	double getS_e(int n) { return Vols_n_Surfaces[n][2]; }
	double getS_s(int n) { return Vols_n_Surfaces[n][3]; }
	double getS_n(int n) { return Vols_n_Surfaces[n][4]; }
	bool hasIdentifier(int n, std::string Id);
	std::vector<double> getVxWest(int n) { std::vector<int> idx; idx = ij_given_n(n); return Coords_Vx[idx[0]][idx[1]]; }
	std::vector<double> getVxEast(int n) { std::vector<int> idx; idx = ij_given_n(n); return Coords_Vx[idx[0]+1][idx[1]]; }
	std::vector<double> getVySouth(int n) { std::vector<int> idx; idx = ij_given_n(n); return Coords_Vy[idx[0]][idx[1]]; }
	std::vector<double> getVyNorth(int n) { std::vector<int> idx; idx = ij_given_n(n); return Coords_Vy[idx[0]][idx[1]+1]; }
	 // - access with (i,j)
	std::vector<double> getCenter(int i, int j) { return Coords[node_number(i,j)][0]; }
	std::vector<double> getWest(int i, int j) { return Coords[node_number(i,j)][1]; }
	std::vector<double> getEast(int i, int j) { return Coords[node_number(i,j)][2]; }
	std::vector<double> getSouth(int i, int j) { return Coords[node_number(i,j)][3]; }
	std::vector<double> getNorth(int i, int j) { return Coords[node_number(i,j)][4]; }
	double getVol(int i, int j) { return Vols_n_Surfaces[node_number(i,j)][0]; }
	double getS_w(int i, int j) { return Vols_n_Surfaces[node_number(i,j)][1]; }
	double getS_e(int i, int j) { return Vols_n_Surfaces[node_number(i,j)][2]; }
	double getS_s(int i, int j) { return Vols_n_Surfaces[node_number(i,j)][3]; }
	double getS_n(int i, int j) { return Vols_n_Surfaces[node_number(i,j)][4]; }
	bool hasIdentifier(int i, int j, std::string Id);
	// Vx mesh
	std::vector<double> getCoords_VxMesh(int i, int j) { return Coords_Vx[i][j]; }
	double getS_VxMesh(int i, int j) { return SnV_Vx[i][j][0]; }
	double getVol_VxMesh(int i, int j) { return SnV_Vx[i][j][1]; }
	std::string getId_VxMesh(int i, int j) { return VxMesh_Ids[i][j]; }
	// Vy mesh
	std::vector<double> getCoords_VyMesh(int i, int j) { return Coords_Vy[i][j]; }
	double getS_Vymesh(int i, int j) { return SnV_Vy[i][j][0]; }
	double getVol_Vymesh(int i, int j) { return SnV_Vy[i][j][1]; }
	std::string getId_VyMesh(int i, int j) { return VyMesh_Ids[i][j]; }
    // Others
	double getX(int i) { return x_coords[i]; }
    double getY(int j) { return y_coords[j]; }

	// Public members
	int N, N1, N2, N3, M, M1, M2, M3, N_vol, N_obstacle;
	
	std::vector<std::vector<std::string>> Identif;

private:
	// Members
	std::ifstream input{ "inputs/Mesh" };
	std::ifstream Boundary_Input{ "inputs/Boundary"};
	std::ifstream Materials_Input{ "inputs/Materials"};
	std::string strInput;
	double L, H, hx_min, hx_max1, hx_max2, hx_max3, hy_min, hy_max1, hy_max2, hy_max3, 
		alpha_x1, alpha_x2, alpha_x3, alpha_y1, alpha_y2, alpha_y3, x_r, y_r, L_r, H_r;
	int n_nod;
	// General 
	std::vector<double> x_coords;
	std::vector<double> y_coords;
	// For primary mesh
	std::vector<std::vector<int>> Connect; // n_nod x 2
	std::vector<std::vector<double>> Vols_n_Surfaces; // n_nod x 5
	std::vector<std::vector<std::vector<double>>> Coords; // n_nod x 5 x 2
	std::vector<double> MaterialK; // Gives material for each node
	//std::vector<std::vector<std::string>> Identif;
	// For Vx mesh
	std::vector<std::vector<std::vector<double>>> Coords_Vx;
	std::vector<std::vector<std::vector<double>>> SnV_Vx;
	std::vector<std::vector<std::string>> VxMesh_Ids;
	// For Vy mesh
	std::vector<std::vector<std::vector<double>>> Coords_Vy;
	std::vector<std::vector<std::vector<double>>> SnV_Vy;
	std::vector<std::vector<std::string>> VyMesh_Ids;
	

	// Private Functions
	bool inObstacle(int i, int j);
	bool inObstacleVx(int i, int j);
	bool inObstacleVy(int i, int j);
	void getInputs();
	void Build();
	std::vector<double> Discretize(double h_min, double h_max, double alpha, double x0, double xf);
	void Build1();
	std::vector<double> Discretize1(double h_min, double h_max, double alpha, double x0, double xf); 
	std::vector<double> Center(int i, int j);
	std::vector<double> West(int i, int j);
	std::vector<double> East(int i, int j);
	std::vector<double> South(int i, int j);
	std::vector<double> North(int i, int j);
	void buildArrays();
	void buildIdentifiers(bool yesObstacle);
	std::vector<int> findBoundaryElems(double x1, double y1, double x2, double y2, std::string type);
	void buildVxMesh();
	void buildVyMesh();
	void buildVxVyIds(double x1, double y1, double x2, double y2, std::string name);
	void AssignMaterials();
};


#endif 
