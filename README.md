# Distributed-computation-with-MPI
Parallel Distributed-computation-with-MPI

CONFIGURE THE ENVIRONMENT
=========================

STEP 1: Load Environment Variables

    The oepnmpi and intel_mkl_lib have been installed in /opt,
    you should edit your ~/.bashrc to add following contents: 

        export PATH=$PATH:/opt/openmpi/bin
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/openmpi/lib
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/lib/intel64:/opt/intel/lib/intel64 (if you need it)

    Then run "source ~/.bashrc" and "mpiexec --vertion" get the infomation about mpi. 


STEP 2: Configure an MPI Cluster within a LAN

    You need to communicate between the four computers,
    and the /etc/hosts is used to map hostnames to IP addresses,
    you can view the contents with command "cat /etc/hosts".

    You need to do the following commands on each machine:

        ssh-keygen -t rsa -b 4096

        ssh-copy-id username@esc-dc50-0001
        ssh-copy-id username@esc-dc50-0002
        ssh-copy-id username@esc-dc50-0003
        ssh-copy-id username@esc-dc50-0004

    Then you must be able to login to other machines without any password prompt.

NOTE:

    If other third-party software is required, install it by yourself, please.

    We have provided matrix data and their description in ./data.

