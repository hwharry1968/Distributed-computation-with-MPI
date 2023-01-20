parallel conjugate gradient solver
===================================
This is a conjugate gradient solver for large sparse matrix. ie, to get x in Ax = b.
Read data from "csrmatrix", save result in "ans".

### generate N*N matrix
		python3 Random_matrix_copy.py N  

### serial algorithm
		g++ -I ./Eigen CPU_Eigen_CG.cpp -o CPU_Eigen_CG
		./cpu_conjugate csrmatrix -t -d   

### parallel algorithm
		mpicxx MPI_iterate_CG.cpp -o mpi_conjugate
		mpiexec -n 4 ./mpi_conjugate csrmatrix

### validation on multi-nodes
		scp mpi_conjugate  eda220618@ecs-dc50-0002:/home/eda220618               //copy the binary file to remote node
		scp csrmatrix eda220618@ecs-dc50-0002:/home/eda220618                    //copy the matrix file to remote node
		mpiexec -n 16 --hostfile hostfile  ./mpi_conjugate csrmatrix             //run on different node

### run the CPU_Eigen_CG: with the algorithm in eigen lib
		python3 Random_matrix_copy.py N
		g++ -I ./Eigen CPU_Eigen_CG.cpp -o CPU_Eigen_CG
		./CPU_Eigen_CG csrmatrix2 -t -d

### run the EIgen_CG：MPI and Eigen
		mpicxx -I ./Eigen EigenCG_copy.cpp -o cg
		mpiexec -n 4 ./cg 200 400
`mpirun -np <p> ./cg <N> <max_steps>`, where `<p>` is the desired number of processes,
    `<N>` is the problem size - i.e. the N in the NxN matrix generated,
    `<max_steps>` is the max number of steps of conjugate gradient method.

