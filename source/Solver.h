#ifndef SOLVER_H
#define SOLVER_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include "SparseMatrix.h"

//---------- Definition of the Solver Interface ----------
class SolverBase
{
public:
	// Constructor
	SolverBase(SparseMatrix in_A, std::vector<double> in_b)
		: A {in_A}, b {in_b}
	{
		N = A.getN();
		setGuess();
	}
	// Virtual Destructor
	virtual ~SolverBase() = default;
	// Public Functions
	virtual void solve();
	void printSolution();
	void printNumIterations() { std::cout << "Converged in " << k << " iterations. " << "\n"; }
	// Public Variables
	int k = 0;
	std::vector<int> ks_POSTPROC;
	std::vector<double> rs_POSTPROC;
	std::vector<double> x;

protected:
	int N;
	SparseMatrix A;
	std::vector<double> b, x_guess;
	double error = 1E-09;
	void setGuess();
	bool compareUpdate();
	//double dotProd(std::vector<double> a, std::vector<double> b);
	//std::vector<double> MatrixVec_prod(SparseMatrix B, std::vector<double> y);
	//std::vector<double> subtractVects(std::vector<double> a, std::vector<double> b);
	//std::vector<double> ScalarVec_prod(double a, std::vector<double> b);
	//For postproc
	//double residual();

};

//---------- Particular Solvers ----------
class Jacobi : public SolverBase
{
public:
	// Constructor
	Jacobi(SparseMatrix in_A, std::vector<double> in_b)
		: SolverBase{ in_A, in_b }
	{
	}
	virtual void solve();

};

class GaussSeidel : public SolverBase
{
public:
	// Constructor
	GaussSeidel(SparseMatrix in_A, std::vector<double> in_b)
		: SolverBase{ in_A, in_b }
	{
	}
	virtual void solve();

};

class SOR : public SolverBase
{
public:
	// Constructor
	SOR(SparseMatrix in_A, std::vector<double> in_b, double in_omega)
		: SolverBase{ in_A, in_b }, omega{ in_omega }
	{
	}
	virtual void solve();

private:
	double omega;
};

class CG : public SolverBase 
{
public:
	// Constructor
	CG(SparseMatrix in_A, std::vector<double> in_b)
		: SolverBase{ in_A, in_b }
	{
	}
	virtual void solve();

private:
	std::vector<double> r, r_guess, v;
	double s_k, t_k;

};

//Preconditioned Conjugate Gradient 
class PCG : public SolverBase 
{
public:
	// Constructor
	PCG(SparseMatrix in_A, std::vector<double> in_b)
		: SolverBase{ in_A, in_b }
	{
	}
	virtual void solve();

protected:
	virtual void precondition() { z=r; }
	std::vector<double> r, r_guess, v, z, z_guess;
	double s_k, t_k;
};
 //> For diagonal preconditioner
class DiagPCG : public PCG 
{
public:
	//Constructor
	DiagPCG(SparseMatrix in_A, std::vector<double> in_b)
		: PCG{ in_A, in_b }
	{
		buildDiag();
	}

protected:
	virtual void precondition();
	void buildDiag();
	std::vector<double> diag;
};
 //> For Gauss-Seidel preconditioner
class GS_PCG : public PCG 
{
public:
	//Constructor
	GS_PCG(SparseMatrix in_A, std::vector<double> in_b)
		: PCG{ in_A, in_b }
	{
	}

protected:
	virtual void precondition();
};



//---------- Context -----------
class Solver
{
private:
	SolverBase* solv;
	std::string solverType, strInput, preconType;
	std::ifstream input{ "inputs/Solver_Input.txt" };
	// Extra inputs for specific solvers
	double omega;

public:
	Solver(SparseMatrix in_A, std::vector<double> in_b);
	void ReadInputs();
	void solve() {
		solv->solve();
	}
	void printSolution() {
		solv->printSolution();
	}
	void printNumIterations() {
		solv->printNumIterations();
	}
	std::vector<double> getSolution() {
		return solv->x;
	}
};


#endif
