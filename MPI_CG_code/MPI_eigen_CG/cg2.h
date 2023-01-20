#include <mpi.h>

typedef struct {
  // Struct representing local data in a system of linear equations Ax = b, 
  // where A is symmetric and positive definite
  int N;

  double *A;             //分块的矩阵
  double *x;
  double *x_star; // vector x* for approximated solution
  double *b;             //分块的b
} equation_data;


typedef struct {
  int N;
  int Np;

  int rank;
  int coord; // Coordinate of this row in the 1D grid
  int displ; // The displacement of this row from the top
 
  int count; // Height of this row
  int count_max; // The maximum height of a row 
  int count_min; // The minimum height of a row

  MPI_Datatype block_trans_t; // Type for a block of columns in a row
  MPI_Datatype row_t; // Type for a full row

  int *ranks; // Rank of all processes, indexed by coordinate
  int *counts; // Dimension of all rows, indexed by rank
  int *displs; // Displacement of all rows, indexed by rank

  MPI_Comm comm;
} process_data;

process_data set_up_world(int Np, int N);
double *random_matrix(int N, int M);
equation_data random_linear_system(process_data row);
double max_error(double *real_x, double *approx_x, int N);
void malloc_test(void *ptr);

double *conjugate_gradient_serial(double *A, double *b, int N, int max_steps, double tol);

void conjugate_gradient_parallel(process_data row, equation_data equation, int N, int max_steps, double tol);
