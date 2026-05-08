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

//double SolverBase::dotProd(std::vector<double> a, std::vector<double> b) {
//	double z = std::inner_product(
//		a.begin(),
//		a.end(),
//		b.begin(),
//		0.0
//	);
//	return z;
//}

//std::vector<double> SolverBase::MatrixVec_prod(SparseMatrix B, std::vector<double> y) {
//	std::vector<double> z = y;
//	std::vector<int> Ids;
//	double sum;
//	int n, j;
//	for (int i = 0; i < N; i++) {
//		sum = 0;
//		Ids = B.getIds(i);
//		n = (int)Ids.size();
//		for (int p = 0; p < n; p++) {
//			j = Ids[p];
//			sum = sum + B.getEntry(i,j)*y[j];
//		}
//		z[i] = sum;
//	}
//	return z;
//}

//std::vector<double> SolverBase::ScalarVec_prod(double a, std::vector<double> b) {
//	std::vector<double> z = b;
//	for (int i = 0; i < N; i++) {
//		z[i] = a*b[i];
//	}
//	return z;
//}

//std::vector<double> SolverBase::subtractVects(std::vector<double> a, std::vector<double> b) {
//	std::vector<double> z = a;
//	for (int i = 0; i < N; i++) {
//		z[i] = a[i] - b[i];
//	}
//	return z;
//}

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
	// Local variables
	r = b;
	r_guess = r;
	int counter = 0;
	double sum, rGrG, rr;
	bool stop=false;
	std::vector<int> ids;
	std::vector<double> Av(N);
	// Initial setup
	for (int i = 0; i < N; i++) {
		sum = 0;
		ids = A.getIds(i);
		for (size_t k = 0; k<ids.size(); k++)
			sum += A.getEntry(i,ids[k])*x_guess[ids[k]];
		r_guess[i] = b[i] - sum; 
	} 
	v = r_guess;
	// Loop
	while (!stop) {
		counter++;
		if(counter>10000)
			break;
		// Calculate A*v
		for (int i = 0; i < N; i++) {
			sum = 0;
			ids = A.getIds(i);
			for (size_t k = 0; k<ids.size(); k++)
				sum += A.getEntry(i,ids[k])*v[ids[k]];
			Av[i] = sum;
		}
		// Calculate s_k
		rGrG = std::inner_product(r_guess.begin(), r_guess.end(), r_guess.begin(), 0.0);
		s_k = rGrG/std::inner_product(v.begin(), v.end(), Av.begin(), 0.0);
		// Update x and r
		for (int i = 0; i < N; i++) {
			x[i] = x_guess[i] + s_k*v[i];
			r[i] = r_guess[i] - s_k*Av[i];
		}
		// Check convergence
		rr = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
		if (std::sqrt(rr)<error)
			stop = true;
		// Update v
		t_k = rr/rGrG;
		for (int i = 0; i < N; i++) 
			v[i] = r[i] + t_k*v[i];
		// Update guesses
		r_guess = r;
		x_guess = x;
		k = counter;
	}
}

//---------- Preconditioned Conjugate Gradient ----------
void PCG::solve() {
	// Local variables
	r = b;
	z = b;
	int counter = 0;
	double sum, rGzG, rz, rr;
	bool stop=false;
	std::vector<int> ids;
	std::vector<double> Av(N);
	// Initial setup
	for (int i = 0; i < N; i++) {
		sum = 0;
		ids = A.getIds(i);
		for (size_t k = 0; k<ids.size(); k++)
			sum += A.getEntry(i,ids[k])*x_guess[ids[k]];
		r[i] = b[i] - sum; 
	} 
	precondition();
	r_guess = r;
	z_guess = z;
	v = z_guess;
	while (!stop) {
		counter++;
		if(counter>10000)
			break;
		//Calculate A*v
		for (int i = 0; i < N; i++) {
			sum = 0;
			ids = A.getIds(i);
			for (size_t k = 0; k<ids.size(); k++)
				sum += A.getEntry(i,ids[k])*v[ids[k]];
			Av[i] = sum;
		}
		//Calculate s_k
		rGzG = std::inner_product(r_guess.begin(), r_guess.end(), z_guess.begin(), 0.0);
		s_k = rGzG/std::inner_product(v.begin(), v.end(), Av.begin(), 0.0);
		//Update x, r and z
		for (int i = 0; i < N; i++) {
			x[i] = x_guess[i] + s_k*v[i];
			r[i] = r_guess[i] - s_k*Av[i];
		}
		precondition();
		//Check convergence
		rz = std::inner_product(z.begin(), z.end(), z.begin(), 0.0);
		rr = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
		if(std::sqrt(rr)<error)
			stop=true;
		//Update v
		t_k = rz/rGzG;
		for (int i = 0; i < N; i++)
			v[i] = z[i] + t_k*v[i];
		//Update guesses
		z_guess = z;
		r_guess = r;
		x_guess = x;
		k = counter;
	}
}
//--> For Diagonal preconditioner
void DiagPCG::buildDiag() {
	for (int i = 0; i < N; i++) {
		diag.push_back(A.getEntry(i,i));
	}
}
void DiagPCG::precondition() {
	for (int i = 0; i < N; i++) {
		z[i] = r[i]/diag[i];
	}
}
//--> For Gauss-Seidel preconditioner
void GS_PCG::precondition() {
	std::vector<double> z(N);
	std::vector<double> y(N);
	double sum;
	int n, j;
	std::vector<int> ids;
	// L + D precon
	y[0] = r[0]/A.getEntry(0,0);
	for (int i = 1; i < N; i++) {
		sum = 0;
		ids = A.getIds(i);
		n = (int)ids.size();
		for (int p = 0; p < n; p++) {
			j = ids[p];
			if (j < i) {
				sum += A.getEntry(i,j)*y[j];
			}
		}
		y[i] = (r[i] - sum) / A.getEntry(i,i);
	}
	// D + U precon
	z[N-1] = y[N-1];
	for (int i = N - 2; i >= 0; i = i - 1) {
		sum = 0;
		ids = A.getIds(i);
		n = (int)ids.size();
		for (int p = 0; p < n; p++) {
			j = ids[p];
			if ( j > i )
			sum += A.getEntry(i,j)*z[j];
		}
		z[i] = y[i] - sum/A.getEntry(i,i);
	}
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
