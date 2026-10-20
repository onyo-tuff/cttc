#ifndef CFL_CLASS
#define CFL_CLASS

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Mesh.h"

class CFL {
public:
	//Constructor
	CFL();

	//Public functions
	void cflDiff(Mesh& mesh);
	//void cflConv(); (future ampliation)
	double getDt() { return dt; }

private:
	//Private elements
	std::string integrationSch;
	double dt;
	//Private functions
	void readScheme();
};

#endif
