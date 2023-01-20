#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>

#include <vector>
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

using namespace std;

// solve AX=b in CG with eigen lib in CPU
 
//time info
timeval start, end, end2;
double elapsedTime, elapsedTime2;
 
int main(int argc, char** argv)
{
	//read square matrix for compte
	string matrixName(argv[1]);
	ifstream matrixFile(matrixName.c_str());
	if (!matrixFile.is_open()) {
		std::cerr << "Problem with file " << matrixName << ".\n";
		exit(1);
	}
	int M, N, NZ;
	matrixFile >> M >> N >> NZ;
	if (M != N) {
		std::cerr << "Only square matrices are accepted.\n";
    	return 0;
	}

	//compress mode for A
	int* row_ptr = new int[M+1];
	int* col_ind = new int[NZ];
	double* val = new double[NZ];

	// define the x b and A 
	Eigen::VectorXd x(M), b(M);
	Eigen::SparseMatrix<double> A(M,N);

	for(int i = 0; i < M+1; i++)
		matrixFile >> row_ptr[i];
	
	for(int i = 0; i < NZ; i++)
		matrixFile >> col_ind[i];

	for(int i = 0; i < NZ; i++)
		matrixFile >> val[i];

	for(int i = 0; i < M; i++)
		matrixFile >> b[i];  // fill in b

	for(int i = 0; i < M; i++)
		matrixFile >> x[i];  // fill in x 
    matrixFile.close();

	//CG solver

	//use the triple to store the matrix
    vector<Eigen::Triplet<double> > triple;
	for(int i = 0; i < M; i++){
        for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
            triple.push_back(Eigen::Triplet<double>(i, col_ind[j], val[j]));
        }
    }
	// fill in A
	A.setFromTriplets(triple.begin(), triple.end()); // prepare sparse matrix to go
	
	// compute Ax=b
	gettimeofday(&start, NULL);
	// CG solver
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver; // set Lower|Upper can boost performance
	// BiCGSTAB solver
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	solver.compute(A);
	x = solver.solve(b);
	gettimeofday(&end, NULL);
 
	//Output time  iterations and error
	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;

	cout<<"Time: "<<elapsedTime<<" ms."<<endl;
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error()      << std::endl;
	// update b, and solve again
	x = solver.solve(b);
 
  return 0;
}
 
 
