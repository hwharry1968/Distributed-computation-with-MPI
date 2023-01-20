#include <mpi.h>
#include <stdio.h>

// compile command
// mpicxx hello_MPI.cpp -o hello_world
// mpiexec -n myid_number ./hello_world
// //more than one node will be spycial, myid_number canbe you want to test

int main(int argc,char *argv[])
{
    int myid,numprocs;

    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    printf("Hello world! Process %d of %d\n", myid, numprocs);

    MPI_Finalize();
}