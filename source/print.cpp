#include "print.h"

void print(Mesh* mesh, std::vector<std::vector<double>> x, std::string name) { 
    int N = mesh->N;
    int M = mesh->M;
    std::ofstream file;
    //std::string directory  = "results/";
    file.open(name.c_str());
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
            file << x[i][j] << "\n";
        }
    }
    file.close();
}
