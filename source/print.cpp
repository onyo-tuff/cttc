#include "print.h"

void print(Mesh* mesh, std::vector<double> x, std::string name) { 

	// Write in i,j format
	int i1;
	int j1;
    int N = mesh->N;
    int M = mesh->M;
    
    // Add obstacle filler to x
    std::vector<double> x1(mesh->N_vol);
    int n = 0;
    int n_noFiller = 0;
    if (mesh->hasObstacle) {
    	for (int i = 0; i < N; i++) {
    		for (int j = 0; j < M; j++) {
				if (mesh->inObstacle(i,j))
					x1[n] = 0;
				else {
					x1[n] = x[n_noFiller];
					n_noFiller++;
				}
				n++;
			}
		}
		x = x1; // Update vector
		//x1.erase(); // Free memory
	}


    // Write file
    std::ofstream file;
    std::string directory  = "results/";
    file.open(directory + name + ".vtk");
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
    file<<"SCALARS " << name << " double"<<"\n";
    file<<"LOOKUP_TABLE default"<<"\n";
    file<<"\n";
    // Write scalar values
    for (int j = 0;j < M;j++) {
        for (int i = 0;i < N;i++) {
        	if(!mesh->inObstacle(i,j)) {
        		n = mesh->node_number(i,j);
            	file << x[n] << "\n";
            	n++;
        	}
        	else 
        		file << "0 \n";
        	// To check boundaries
        	//if(mesh->inObstacle(i,j)) 
        	//	file << 0 << "\n";
    		//else if(mesh->hasIdentifier(i,j,"Int"))
    		//	file << 1 << "\n";
			//else if(mesh->hasIdentifier(i,j,"Ext"))
			//	file << 2 << "\n";
			//else 
			//	file << 3 << "\n";
			// Check mesh definition
			//file << mesh->getCenter(i,j)[0];
			//file << mesh->getCenter(i,j)[1];
        }
    }
    file.close();
}

void printVTR(Mesh* mesh, std::vector<double> x, std::string name) {
	// Write in i,j format
	int i1;
	int j1;
    int N = mesh->N;
    int M = mesh->M;
    
    // Add obstacle filler to x
    std::vector<double> x1(mesh->N_vol);
    int n = 0;
    int n_noFiller = 0;
    if (mesh->hasObstacle) {
    	for (int i = 0; i < N; i++) {
    		for (int j = 0; j < M; j++) {
				if (mesh->inObstacle(i,j))
					x1[n] = 0;
				else {
					x1[n] = x[n_noFiller];
					n_noFiller++;
				}
				n++;
			}
		}
		x = x1; // Update vector
		//x1.erase(); // Free memory
	}
	
	// Write XML file
	std::ofstream file;
	std::string directory  = "results/";
	file.open(directory + name + ".vtr");
	file << "<?xml version=\"1.0\"?>\n";
	file << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << "\n";
	file << "<RectilinearGrid WholeExtent=\"0 " << N << " 0 " << M << " 0 0\">\n";
	file << "<Piece Extent=\"0 " << N << " 0 " << M << " 0 0\">\n";
	file << "<Coordinates>\n";
	file << "<DataArray type=\"Float64\" " << "Name=\"XCoordinates\" " << "format=\"ascii\">\n";
    for(int i=0;i<=N;i++){
    	file << mesh->getX(i) << " ";
    }
    file<<"\n";
    file<<"</DataArray>"<<"\n";
	file << "<DataArray type=\"Float64\" " << "Name=\"YCoordinates\" " << "format=\"ascii\">\n";
	for(int j=0;j<=M;j++){
        file << mesh->getY(j) << " ";
    }
    file<<"\n";
    file<<"</DataArray>"<<"\n";
    file << "<DataArray type=\"Float64\" " << "Name=\"YCoordinates\" " << "format=\"ascii\">\n";
    file<<"0 \n";
    file<<"</DataArray>"<<"\n";
    file << "</Coordinates>\n";
    file << "<CellData Scalars=\"x\">\n";
	file << "<DataArray type=\"Float64\" " << "Name=\"" << name << "\" " << "format=\"ascii\">\n";
	    for (int j = 0;j < M;j++) {
        for (int i = 0;i < N;i++) {
        	if(!mesh->inObstacle(i,j)) {
        		n = mesh->node_number(i,j);
            	file << x[n] << "\n";
            	n++;
        	}
        	else 
        		file << "nan \n";
        }
    }
    file << "</DataArray>" << "\n";
    file << "</CellData>" << "\n";
    file << "</Piece>" << "\n";
    file << "</RectilinearGrid>" << "\n";
    file << "</VTKFile>" << "\n";
	
}

























