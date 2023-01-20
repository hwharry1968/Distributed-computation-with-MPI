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


	//A矩阵
	Map<SparseMatrix<double,RowMajor> > sm1(M,N,NZ,row_ptr, col_ind,val); //RW 可读写


	int* row_ptr4 = new int[M+1];
	int* col_ind4 = new int[NZ];
	double* val4 = new double[NZ];

    int g = 0;
    int nnz = 0;
    row_ptr4[0] = 0;
	if(rank == 0)cout<<sm1<<endl<<"*************"<<endl;
	for(int i = 0; i < M; i++){
		for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
			if(i <= col_ind[j]){
				col_ind4[g] = col_ind[j];
				val4[g] = val[j];
				g++;
				nnz++;
			}
		}
		row_ptr4[i+1] = nnz;
	}
	
	Map<SparseMatrix<double,RowMajor> > sm4(M,N,NZ,row_ptr4, col_ind4,val4); //RW 可读写


	//if(rank ==0)cout<<sm4<<endl<<"*************"<<endl;
	int i, j, k;
	
	
	if (rank == 0) {
		start = MPI_Wtime();
	}
	
    //for each row
	for (j = 0; j < M; j++) {

	// Step 0: Replace the entries above the diagonal with zeroes

    // Step 1: Update the diagonal element

		if (j%npes == rank &&col_ind4[row_ptr4[j]]==j) {
			if (j==0)
			{
				val4[row_ptr4[j]] = sqrt(val4[row_ptr4[j]]);
			}
			else{
					val4[row_ptr4[j]] = val4[row_ptr4[j]]-sm4.col(j).cwiseProduct(sm4.col(j)).middleRows(0,j).sum();
					val4[row_ptr4[j]] = sqrt(val4[row_ptr4[j]]);
			}

		}

		// Broadcast row with new values to other processes

		//MPI_Bcast(&val3[row_ptr3[j]], row_ptr3[j+1]-row_ptr3[j], MPI_DOUBLE, j%npes, MPI_COMM_WORLD);
		MPI_Bcast(&val4[row_ptr4[j]], row_ptr4[j+1]-row_ptr4[j], MPI_DOUBLE, j%npes, MPI_COMM_WORLD);
		

    // Step 2: Update the elements below the diagonal element

		// Divide the rest of the work
	
		for (i = row_ptr4[j]+1; i < row_ptr4[j+1]; i++) {
			if (i%npes == rank &&col_ind4[row_ptr4[j]]==j) {
			if (j==0)
			{
				val4[i] = val4[i] / val4[row_ptr4[j]];

			}
			else
			{
				val4[i] = val4[i]- sm4.col(col_ind[i]).cwiseProduct(sm4.col(j)).middleRows(0,j).sum();
				val4[i] = val4[i] / val4[row_ptr4[j]];
			}
					
			}
		
		}
	
	
	}

	
	MPI_Barrier(MPI_COMM_WORLD); 
	
	if (rank == 0){
		end = MPI_Wtime();
		printf("Runtime = %lf\n", end-start);
		// Matrix<double, 5000, 1> x,b,solution;
		// x.setRandom(5000);
		// b = sm4 * x;
		// double solvet = MPI_Wtime();
		// solution = sm4.triangularView<Upper>().solve(b);
		// double solveend = MPI_Wtime();
		//printf("Runtime = %lf\n", solveend-solvet);
		// cout<<sm4<<endl;
    	 cout<<sm4.transpose()*sm4<<endl;
		// for(i=0;i<NZ;i++)
		// {
		// 	cout<<val4[i]<<endl;
		// }
		
		
		


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

