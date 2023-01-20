#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include <mpi.h>
 
#include "Eigen/Sparse"

using namespace std;

const double eps = 1e-6;
int myid, numprocs;
MPI_Status status;

double VecMulVec(double *a, double *b, int N) {
	double ans = 0;
	for (int i = 0; i < N; i++) {
		ans += a[i] * b[i];
	}
	return ans;
}
//可以优化为内敛函数
bool check(double *r, int N) {
	double res = 0;
	for (int i = 0; i < N; i++) {
		res += fabs(r[i]);
	}
	if (res < eps) return true;
	else return false;
}

int main (int argc, char* argv[]){
	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//read square matrix
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
	double *b = new double [M];

	for(int i = 0; i < M+1; i++)
		matrixFile >> row_ptr[i];

  	for(int i = 0; i < NZ; i++)
		matrixFile >> col_ind[i];

	for(int i = 0; i < NZ; i++)
		matrixFile >> val[i];

	for(int i = 0; i < M; i++)
		matrixFile >> b[i];
	matrixFile.close();

//time info
	timeval start, end, end2;
	double elapsedTime, elapsedTime2;
//
	double *x = new double [M];
	double *r = new double [M];
	double *p = new double [M];

	double *tmp, *nr;
	if (myid == 0) {    
		gettimeofday(&start, NULL);           //记录开始时间戳
		tmp = new double [M];         //保留中间迭代结果
		nr = new double [M];          //残差
		for (int i = 0; i < M; i++) {
			x[i] = 0;     //向量x随机赋初值
			r[i] = b[i];
			p[i] = b[i];
		}
	} 

	double *sub_vec_y = new double[M / numprocs];
	int done = 0;
	//矩阵分块不均衡，同时进程0需要计算剩余部分，需要优化?
	while(true){
		MPI_Bcast(p, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for(int i = 0; i < (M / numprocs); ++i){
			sub_vec_y[i] = 0;
			int idx = myid * (M / numprocs) + i;
			for(int j = row_ptr[idx]; j < row_ptr[idx + 1]; ++j){
				sub_vec_y[i] += val[j] * p[col_ind[j]];
	 		}
		}
		MPI_Gather(sub_vec_y, M / numprocs, MPI_DOUBLE, tmp, M / numprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(myid == 0){
			for(int i = M / numprocs * numprocs; i < M ; i++){
				for(int j = row_ptr[i]; j < row_ptr[i + 1]; j++){
					tmp[i] += val[j] * p[col_ind[j]]; 
				}
			}
			double alpha = VecMulVec(r, r, M) / VecMulVec(p, tmp, M);
			for (int i = 0; i < N; i++)
				x[i] += alpha * p[i];
			for (int i = 0; i < M; i++)
				nr[i] = r[i] - alpha * tmp[i];
			if (check(nr, M))
				done = 1;
			double beta = VecMulVec(nr, nr, M) / VecMulVec(r, r, M);
			for (int i = 0; i < M; i++)
				p[i] = nr[i] + beta * p[i];
			for (int i = 0; i < M; i++)
				r[i] = nr[i];
		}
		MPI_Bcast(&done, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(done != 0 )break;
	}
	if(myid == 0){
	
		gettimeofday(&end, NULL);				//记录结束时间戳

		//Output time and iterations
		elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
		elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
		cout<<"Time: "<<elapsedTime<<" ms."<<endl;
		
		// Output the matrixes for debug
		/*
		for (int i = 0; i < N; i++) 
			cout<<x[i]<<endl;
		*/
		//Generate files for debug purpouse
		ofstream Af;
		Af.open("mpi");
		for (int i = 0; i < N; i++) 
			Af<<x[i]<<endl;                 
		Af.close();
	}
	delete []row_ptr;
	delete []col_ind;
	delete []val;
	delete []b;
	delete []x;
	delete []r;
	delete []p;
	delete []sub_vec_y;
	MPI_Finalize();
	return 0;
}
