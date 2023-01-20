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



	Map<SparseMatrix<double> > sm3(M,N,NZ,row_ptr, col_ind,val); //RW 可读写

	int* row_ptr3 = new int[M+1];
	int* col_ind3 = new int[NZ];
	double* val3 = new double[NZ];

    int g = 0;
    int nnz = 0;
    row_ptr3[0] = 0;

	for(int i = 0; i < M; i++){
		for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
			if(i < col_ind[j]){
				col_ind3[g] = col_ind[j];
				val3[g] = val[j];
				g++;
				nnz++;
			}
		}
		row_ptr3[i+1] = nnz;
	}

	Map<SparseMatrix<double> > sm1(M,N,nnz,row_ptr3, col_ind3,val3); //RW 可读写

	int i, j, k;
	
	
	if (rank == 0) {
		start = MPI_Wtime();
	}
	
    //for each column
	for (j = 0; j < N; j++) {

	// Step 0: Replace the entries above the diagonal with zeroes

    // Step 1: Update the diagonal element

		if (j%npes == rank) {
			if(col_ind3[row_ptr3[j + 1] - 1] == j){
				for(int l = row_ptr3[j]; l < row_ptr3[j + 1] - 1; l++){
					val3[row_ptr3[j + 1] - 1] -= val[l] * val[l];
				}
				val3[row_ptr3[j + 1] - 1] = sqrt(val3[row_ptr3[j + 1] - 1]);
			}	
		}

		// Broadcast row with new values to other processes

		MPI_Bcast(&val3[row_ptr3[j]], row_ptr3[j+1]-row_ptr3[j], MPI_DOUBLE, j%npes, MPI_COMM_WORLD);


    // Step 2: Update the elements below the diagonal element

		// Divide the rest of the work
	
		for (i = j+1; i < N; i++) {
			if (i%npes == rank) {
				if ( sm1.coeffRef(j,j)!=0){
					for (k = 0; k < j; k++) {
						sm1.coeffRef(i,j) = sm1.coeffRef(i,j)- sm1.coeffRef(i,k)*sm1.coeffRef(j,k);
					}
				sm1.coeffRef(i,j) = sm1.coeffRef(i,j) / sm1.coeffRef(j,j);
				}
			}
		
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD); 
	
	if (rank == 0){
		end = MPI_Wtime();
		printf("Runtime = %lf\n", end-start);

        //cout<<sm1<<endl;
		

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
	return 0;
}

