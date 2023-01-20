#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "Eigen/Sparse" 
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <sys/time.h>


using namespace std;
using namespace Eigen;

int main(int argc, char * argv[]){
	int npes, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double start, end;
	MPI_Barrier(MPI_COMM_WORLD); /* Timing */
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

	int* row_ptr = new int[M+1];
	int* col_ind = new int[NZ];
	double* val = new double[NZ];

	for(int i = 0; i < M+1; i++)
		matrixFile >> row_ptr[i];
	
	for(int i = 0; i < NZ; i++)
		matrixFile >> col_ind[i];

	for(int i = 0; i < NZ; i++)
		matrixFile >> val[i];
    matrixFile.close();

    //SparseMatrix<double, RowMajor> matrix(M, N);

    vector<Eigen::Triplet<double> > triple;
	for(int i = 0; i < M; i++){
        for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
            triple.push_back(Triplet<double>(i, col_ind[j], val[j]));
        }
    }
    //matrix.setFromTriplets(triple.begin(), triple.end());


	if (rank == 0) {
		start = MPI_Wtime();

	}
	int i, j, k;
    //for each column
	for (j = 0; j < N; j++) {

	// Step 0: Replace the entries above the diagonal with zeroes

		if (rank == 0) {
			for (i = 0; i < j; i++) {
				matrix.coeffRef(i,j)=0.0;
				
			}
		}

    // Step 1: Update the diagonal element

		if (j%npes == rank) {
				for (k = 0; k < j; k++) {
				matrix.coeffRef(j,j)= matrix.coeffRef(j,j)-matrix.coeffRef(j,k) * matrix.coeffRef(j,k);
			}
				matrix.coeffRef(j,j) = sqrt(matrix.coeffRef(j,j));
			}
			

		// Broadcast row with new values to other processes

		MPI_Bcast(&matrix, N, MPI_DOUBLE, j%npes, MPI_COMM_WORLD);


    // Step 2: Update the elements below the diagonal element

		// Divide the rest of the work

		for (i = j+1; i < N; i++) {
			if (i%npes == rank) {
			//if (matrix.coeffRef(i,j)!=0 && matrix.coeffRef(i,j)!=0){
				for (k = 0; k < j; k++) {
					matrix.coeffRef(i,j) =matrix.coeffRef(i,j)- matrix.coeffRef(i,k)*matrix.coeffRef(j,k);
				}
			
				matrix.coeffRef(i,j) = matrix.coeffRef(i,j) / matrix.coeffRef(j,j);
			
				
		}
		
	}
	}

	MPI_Barrier(MPI_COMM_WORLD); /* Timing */
	if (rank == 0){
		end = MPI_Wtime();
		printf("Runtime = %lf\n", end-start);

        cout<<matrix.transpose()*matrix<<endl;
		

        // Print Matrix
        // double ** LLT = matrixMultiply(L, transpose(L, n), n);
		// printf("L=\n");
		// print(L,n);
        // printf("L*L^T = \n");
        // print(LLT, n);
		// printf("A=\n");
		// print(A,n);

	}
	MPI_Finalize();
	
}