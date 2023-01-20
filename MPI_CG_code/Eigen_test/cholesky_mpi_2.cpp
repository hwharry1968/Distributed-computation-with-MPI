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
	int i, j, k, r;
	if (rank == 0) {
		start = MPI_Wtime();

	}
	
    //for each column
	for (j = 0; j < N; j++) {
	// Step 0: Replace the entries above the diagonal with zeroes
		if (rank == 0) {
			for (i = col_ind[row_ptr[j]]; i<col_ind[row_ptr[j+1]]; i++) {
				if (i < j)
                {
                    val[i]=0.0;
                }
                

				
			}
		}	
    // Step 1: Update the diagonal element
		if (j%npes == rank) {
            for (i = col_ind[row_ptr[j]]; i<col_ind[row_ptr[j+1]]; i++) {
				if (i == j)
                {
                    for(k=0; k<row_ptr[j+1]-row_ptr[j];k++){
                    if (col_ind[row_ptr[i]+k]<i)
                    {
                        val[i]=val[i]-val[row_ptr[j-1]+k]*val[row_ptr[j-1]+k]
                    }
                    }
                    
                }
                val[i]=sqrt(val[i]);
				
			}
        }

		// Broadcast row with new values to other processes

		MPI_Bcast(&val[row_ptr[j]:row_ptr[j+1]], N, MPI_DOUBLE, j%npes, MPI_COMM_WORLD);


    // Step 2: Update the elements below the diagonal element

		// Divide the rest of the work

		for (i=0; i < row_ptr[j+1]-row_ptr[j]; i++) {
			if (i%npes == rank) {
				for (k = 0; k <i; k++) {
                    val[row_ptr[j]+i]=val[row_ptr[j]+i]-val[row_ptr[i]+k]*val[row_ptr[j]+k]

                }
                
			val[row_ptr[j]+i]=val[row_ptr[j]+i]/val[row_ptr[j]+r]
		}
		
        }
	}
	
	

	MPI_Barrier(MPI_COMM_WORLD); /* Timing */
	if (rank == 0){
		end = MPI_Wtime();
		printf("Runtime = %lf\n", end-start);

        cout<<matrix<<endl;
		

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