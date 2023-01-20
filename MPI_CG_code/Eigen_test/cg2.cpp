#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cg2.h"

// Solve Ax = b for x, using the Conjugate Gradient method.
// Terminates once the maximum number of steps or tolerance has been reached
double *conjugate_gradient_serial(double *A, double *b, int N, int max_steps, double tol) {
  double *x, *r, *p, *a, *z;
  double gamma, gamma_new, alpha, beta;
  int k = 0;
  
  x = malloc(N*sizeof(double));
  r = malloc(N*sizeof(double));
  p = malloc(N*sizeof(double));
  z = malloc(N*sizeof(double));

  malloc_test(x);
  malloc_test(r);
  malloc_test(p);
  malloc_test(z);

  // x = [0 ... 0]
  // r = b - A * x
  // p = r
  // gamma = r' * r
  gamma = 0.0;
  for (int i = 0; i < N; ++i) {
    x[i] = 0.0;
    r[i] = b[i];
    p[i] = r[i];
    gamma += r[i] * r[i];
  }


  for (int n = 0; n < max_steps; ++n) {
    // z = A * p
    for (int i = 0; i < N; ++i) {
      a = A + (i*N);
      z[i] = 0.0;
      for (int j = 0; j < N; ++j)
	      z[i] += a[j] * p[j]; 
    }

    // alpha = gamma / (p' * z)
    alpha = 0.0;
    for (int i = 0; i < N; ++i)
      alpha += p[i] * z[i];
    alpha = gamma / alpha;
    
    // x = x + alpha * p
    // r = r - alpha * z
    // gamma_new = r' * r
    gamma_new = 0.0;
    for (int i = 0; i < N; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * z[i];
      gamma_new += r[i] * r[i];
    }

    if (sqrt(gamma_new) < tol)
      break;
    
    beta = gamma_new / gamma;

    // p = r + (gamma_new / gamma) * p;
    for (int i = 0; i < N; ++i)
      p[i] = r[i] + beta * p[i];
    
    // gamma = gamma_new
    gamma = gamma_new;
  }

  free(r);
  free(p);
  free(z);

  return x;
}


void conjugate_gradient_parallel(process_data row, equation_data equation, int N, int max_steps, double tol) {
  double *a, *p, *p_recv, *p_send, *p_tmp, *r, *z, *x;
  double alpha, alpha_tmp, beta, gamma, gamma_new, gamma_tmp;
  int rank_up, rank_down, rank_block, displ_send, count_send;
  MPI_Request send_req; //用来查询发送是否完成
  
  x = &(equation.x_star[0]);

  p = malloc(row.count*sizeof(double));
  r = malloc(row.count*sizeof(double));
  z = malloc(row.count*sizeof(double));

  malloc_test(p);
  malloc_test(r);
  malloc_test(z);

  p_recv = malloc(row.count_max*sizeof(double));
  p_send = malloc(row.count_max*sizeof(double));
  malloc_test(p_recv);
  malloc_test(p_send);

  rank_up = row.ranks[(row.Np + row.coord - 1)%row.Np];
  rank_down = row.ranks[(row.coord + 1)%row.Np];

  // x = [0 ... 0] - initial guess
  // r = b
  // p = r
  // gamma = r' * r
  gamma_tmp = 0.0;
  for (int i = 0; i < row.count; ++i) {
    x[i] = 0.0;
    r[i] = equation.b[i];
    p[i] = r[i]; 
    p_recv[i] = p[i];
    z[i] = 0.0;
    gamma_tmp += r[i] * r[i];
  }
  MPI_Allreduce(&gamma_tmp, &gamma, 1, MPI_DOUBLE, MPI_SUM, row.comm);

  for (int n = 0; n < max_steps; ++n) {
    // z = A * p
    for (int m = 0; m < row.Np; ++m) {
      p_tmp = p_recv;
      p_recv = p_send;
      p_send = p_tmp;   //流水线发送

      if (m < row.Np-1)                             //目的进程号
	      MPI_Isend(p_send, row.count_max, MPI_DOUBLE, rank_up, 222, row.comm, &send_req);

      rank_block = row.ranks[(row.coord + m)%row.Np];
      displ_send = row.displs[rank_block];       //取出每个block的偏移，大小
      count_send = row.counts[rank_block];

      for (int i = 0; i < row.count; ++i) {
	      a = equation.A + (i*row.N + displ_send);  //分块矩阵A的首地址
	      for (int j = 0; j < count_send; ++j) 
	        z[i] += p_send[j]*a[j];
      }

      if (m < row.Np-1) {
	      MPI_Recv(p_recv, row.count_max, MPI_DOUBLE, rank_down, 222, row.comm, MPI_STATUS_IGNORE);
	      MPI_Wait(&send_req, MPI_STATUS_IGNORE);
      }
    }

    // alpha = gamma / (p' * z)
    alpha_tmp = 0.0;
    for (int i = 0; i < row.count; ++i)
      alpha_tmp += p[i] * z[i];

    MPI_Allreduce(&alpha_tmp, &alpha, 1, MPI_DOUBLE, MPI_SUM, row.comm);

    alpha = gamma / alpha;

    // x = x + alpha * p
    // r = r - alpha * z
    // gamma_new = r' * r
    gamma_tmp = 0.0;
    for (int i = 0; i < row.count; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * z[i];
      gamma_tmp += r[i] * r[i];
    }

    MPI_Allreduce(&gamma_tmp, &gamma_new, 1, MPI_DOUBLE, MPI_SUM, row.comm);

    if (sqrt(gamma_new) < tol)
      break;
    
    beta = gamma_new / gamma;

    // p = r + beta * p;
    for (int i = 0; i < row.count; ++i) {
      p[i] = r[i] + beta * p[i];
      p_recv[i] = p[i];
      z[i] = 0.0;
    }
    
    // gamma = gamma_new
    gamma = gamma_new;
  }

  free(p);
  free(r);
  free(z);
  free(p_send);
  free(p_recv);

  return;
}
