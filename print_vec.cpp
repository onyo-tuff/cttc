#include "print_vec.h"

void print_vec(Mesh* mesh, std::vector<std::vector<double>> x, std::string name) { 
    int N = mesh->N;
    int M = mesh->M;
    std::ofstream file;
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
    //file<<"VECTORS "<<name<<" double"<<"\n";
    //file<<"\n";
    //std::string id;
    // .stk espera los nodos en un determinado orden, usar siempre este formato j, i
    //for (int j = 0;j < M;i++) {
    //    for (int j = 0;j < N;j++) {
    //        file << mesh->getCenter(i,j)[0] << "   " << mesh->getCenter(i,j)[1] << "   " << "0.0" << "\n";
    //    }
    //}
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