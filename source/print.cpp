#include "print.h"

void print(Mesh* mesh, std::vector<double> x, std::string name) { 

	// Write in i,j format
	int i1;
	int j1;
    int N = mesh->N;
    int M = mesh->M;
    std::vector<std::vector<double>> x2D(N, std::vector<double>(M));
    for (size_t n = 0; n < x.size(); n++) {
		i1 = (int)(n / M);
		j1 = n - i1 * M;
		x2D[i1][j1] = x[n];
	}
    
    // Write file
    std::ofstream file;
    std::string directory  = "results/";
    file.open(directory + name);
    file<<"# vtk DataFile Version 2.0"<<"\n";
    file<<name<<"\n";
    file<<"ASCII"<<"\n";
    file<<"\n";
    file<<"DATASET RECTILINEAR_GRID"<<"\n";
    file<<"DIMENSIONS"<<"   "<<N+1<<"   "<<M+1<<"   "<<1<<"\n";
    file<<"\n";
    file<<"X_COORDINATES"<<"   "<<N+1<<"   "<<"double"<<"\n";
    for(int i=0;i<=N;i++){
        file << mesh->getX(i) << "   ";
    }
    file<<"\n";
    file<<"Y_COORDINATES"<<"   "<<M+1<<"   "<<"double"<<"\n";
    for(int j=0;j<=M;j++){
        file << mesh->getY(j) << "   ";
    }
    file<<"\n";
    file<<"Z_COORDINATES"<<"   "<<1<<"   "<<"double"<<"\n";
    file<<0<<"\n";
    file<<"\n";
    file<<"CELL_DATA"<<"   "<<M*N<<"\n";
    file<<"SCALARS x double"<<"\n";
    file<<"LOOKUP_TABLE default"<<"\n";
    file<<"\n";
    for (int j = 0;j < M;j++) {
        for (int i = 0;i < N;i++) {
            file << x2D[i][j] << "\n";
        }
    }
    file.close();
}
