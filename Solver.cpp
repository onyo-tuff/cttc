#include "Solver.h"

//----------- Functions of the Base Solver class ----------
void SolverBase::solve() {
	std::cout << "Error: no solver specified." << "\n";
}

void SolverBase::setGuess() {
	for (int i = 0; i < N; i++) {
		x_guess.push_back(3);
	}
	x = x_guess;
}

void SolverBase::printSolution() { // For testing
	std::cout << "(";
	for (int i = 0; i < N; i++) {
		std::cout << " " << x[i] << ",";
	}
	std::cout << ")" << "\n";
	std::cout << "Number of iterations: " << k << "\n";
	//k++; // Count iterations after each while
	//ks_POSTPROC.push_back(k);
	//rs_POSTPROC.push_back(residual());
}

bool SolverBase::compareUpdate() {
	// Bool true when diff > error at least once
	// Bool false when diff < error for all
	double diff;
	for (int i = 0; i < N; i++) {
		diff = std::abs(x[i] - x_guess[i]);
		if (diff > error) {
			x_guess = x;
			return true;
		}
	}
	x_guess = x;
	return false;
}

double SolverBase::dotProd(std::vector<double> a, std::vector<double> b) {
	double z = std::inner_product(
		a.begin(),
		a.end(),
		b.begin(),
		0.0
	);
	return z;
}

std::vector<double> SolverBase::MatrixVec_prod(SparseMatrix& B, std::vector<double>& y) {
	std::vector<double> z = y;
	std::vector<int> Ids;
	double sum;
	int n, j;
	for (int i = 0; i < N; i++) {
		sum = 0;
		Ids = B.getIds(i);
		n = (int)Ids.size();
		for (int p = 0; p < n; p++) {
			j = Ids[p];
			sum = sum + B.getEntry(i,j)*y[j];
		}
		z[i] = sum;
	}
	return z;
}

std::vector<double> SolverBase::ScalarVec_prod(double a, std::vector<double> b) {
	std::vector<double> z = b;
	for (int i = 0; i < N; i++) {
		z[i] = a*b[i];
	}
	return z;
}

std::vector<double> SolverBase::subtractVects(std::vector<double> a, std::vector<double> b) {
	std::vector<double> z = a;
	for (int i = 0; i < N; i++) {
		z[i] = a[i] - b[i];
	}
	return z;
}

//---------- Jacobi Solver ----------
void Jacobi::solve() {
	bool iterate = true;
	std::vector<int> Ids;
	int n, j;
	while (iterate) {
		// Recalculate x
		for (int i = 0; i < N; i++) {
			x[i] = b[i];
			Ids = A.getIds(i);
			n = (int)Ids.size();
			for (int p = 0; p < n; p++) {
				j = Ids[p];
				if (j!=i)
					x[i] = x[i] - A.getEntry(i,j)*x_guess[j];
			}
			x[i] = x[i]/A.getEntry(i,i);
		}
		// Compare x and x_guess and update
		iterate = compareUpdate();
		// Print solution after each iteration (testing)
		// printSolution();
		// Count iterations
		k++;
	}
}

//---------- Gauss-Seidel Solver ----------
void GaussSeidel::solve() {
	bool iterate = true;
	std::vector<int> Ids;
	int n, j;
	double sum = 0;
	while (iterate) {
		// Recalculate x
		for (int i = 0; i < N; i++) {
			sum = b[i];
			Ids = A.getIds(i);
			n = (int)Ids.size();
			for (int p = 0; p < n; p++) {
				j = Ids[p];
				if ( j != i ) {
					sum += - A.getEntry(i,j)*x[j];
				}
			}
			x[i] = sum/A.getEntry(i,i);
		}
		// Compare x and x_guess and update
		iterate = compareUpdate();
		// Print solution after each iteration (testing)
		// printSolution();
		// Count iterations
		k++;
	}
}

//---------- SOR Solver ----------
void SOR::solve() {
	bool iterate = true;
	std::vector<int> Ids;
	double sum, r;
	int n, j;
	while (iterate) {
		for (int i = 0; i < N; i++) {
			sum = 0;
			//sum2 = 0;
			Ids = A.getIds(i);
			n = (int)Ids.size();
			for (int p = 0; p < n; p++) {
				j = Ids[p];
				sum += A.getEntry(i,j)*x[j];
			}
			r = b[i] - sum;
			x[i] = x_guess[i] + omega*r/A.getEntry(i,i);
		}
		// Compare x vs x_guess and update
		iterate = compareUpdate();
		// Print solution after each iteration (testing)
		// printSolution();
		// Count iterations
		k++;
	}
}

//---------- Conjugate Gradient Solver ----------
void CG::solve() {
	bool iterate = true;
	// Initial setup
	r_guess = subtractVects(b,MatrixVec_prod(A,x_guess));
	v = r_guess;
	while (std::sqrt(dotProd(r_guess, r_guess)) > error) {
		t_k = dotProd(r_guess,r_guess)/dotProd(v,MatrixVec_prod(A,v));
		x = subtractVects(x_guess,ScalarVec_prod(-t_k,v));
		r = subtractVects(r_guess,ScalarVec_prod(t_k,MatrixVec_prod(A,v)));
		s_k = dotProd(r,r)/dotProd(r_guess,r_guess);
		v = subtractVects(r,ScalarVec_prod(-s_k,v));
		// update
		x_guess = x;
		r_guess = r;
		//Print solution after each iteration (testing)
		// printSolution();
		//Count iterations
		k++;
	}
}

//---------- Preconditioned Conjugate Gradient ----------
void PCG::solve() {
	bool iterate = true;
	// Initial setup
	r_guess = subtractVects(b,MatrixVec_prod(A,x_guess));
	z_guess = precondition(r_guess);
	v = z_guess;
	while (std::sqrt(dotProd(r_guess, r_guess)) > error) {
		t_k = dotProd(z_guess,r_guess)/dotProd(v,MatrixVec_prod(A,v));
		x = subtractVects(x_guess,ScalarVec_prod(-t_k,v));
		r = subtractVects(r_guess,ScalarVec_prod(t_k,MatrixVec_prod(A,v)));
		z = precondition(r);
		s_k = dotProd(z,r)/dotProd(z_guess,r_guess);
		v = subtractVects(z,ScalarVec_prod(-s_k,v));
		// Compare x vs x_guess and update
		x_guess = x;
		r_guess = r;
		z_guess = z;
		//Print solution after each iteration (testing)
		//printSolution();
		//Count iterations
		k++;
	}
}
//--> For Diagonal preconditioner
void DiagPCG::buildDiag() {
	for (int i = 0; i < N; i++) {
		diag.push_back(A.getEntry(i,i));
	}
}
std::vector<double> DiagPCG::precondition(std::vector<double> r) {
	std::vector<double> z = r;
	for (int i = 0; i < N; i++) {
		z[i] = r[i]/diag[i];
	}
	return z;
}
//--> For Gauss-Seidel preconditioner
std::vector<double> GS_PCG::precondition(std::vector<double> r) {
	std::vector<double> z = r;
	double sum;
	int n, j;
	std::vector<int> ids;
	// L + D precon
	z[0] = r[0]/A.getEntry(0,0);
	for (int i = 1; i < N; i++) {
		sum = 0;
		ids = A.getIds(i);
		n = (int)ids.size();
		for (int p = 0; p < n; p++) {
			j = ids[p];
			if (j < i) {
				sum += A.getEntry(i,j)*z[j];
			}
		}
		z[i] = (r[i] - sum) / A.getEntry(i,i);
	}
	// D + U precon
	std::vector<double> z1 = z;
	z1[N-1] = z[N-1]/A.getEntry(N-1,N-1);
	for (int i = N - 2; i >= 0; i = i - 1) {
		sum = 0;
		ids = A.getIds(i);
		n = (int)ids.size();
		for (int p = 0; p < n; p++) {
			j = ids[p];
			if ( j > i )
			sum += A.getEntry(i,j)*z1[j];
		}
		z1[i] = (z[i] - sum) / A.getEntry(i,i);
	}
	return z1;
}

//---------- Context (selector class) ----------
Solver::Solver(SparseMatrix in_A, std::vector<double> in_b) 
	: solv(nullptr)
{
	ReadInputs();
	if (solverType == "Jacobi") {
		solv = new Jacobi(in_A, in_b);
	}
	else if (solverType == "Gauss-Seidel") {
		solv = new GaussSeidel(in_A, in_b);
	}
	else if (solverType == "SOR") {
		solv = new SOR(in_A, in_b, omega);
	}
	else if (solverType == "Conjugate Gradient") {
		solv = new CG(in_A, in_b);
	}
	else if (solverType == "PCG") {
		if (preconType == "Diagonal Precon") {
			solv = new DiagPCG(in_A, in_b);
		}
		else if (preconType == "Gauss-Seidel Precon") {
			solv = new GS_PCG(in_A, in_b);
		}
		else {
			solv = new CG(in_A, in_b);
		}
	}
	else {
		solv = new SolverBase(in_A, in_b);
	}
}

void Solver::ReadInputs() {
	while (std::getline(input, strInput)) {
		if (strInput == "Solver to use:") {
			std::getline(input, solverType);
		}
		if (strInput == "Relaxation factor 'omega' (if necessary):") {
			std::getline(input, strInput);
			omega = std::stod(strInput);
		}
		if (strInput == "Preconditioning Scheme (for PCG):") {
			std::getline(input, preconType);
		}
	}
}
