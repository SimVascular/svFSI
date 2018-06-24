MPI_NAME    = mpich
MPI_INCDIR  = $(wordlist 2,99,$(shell /opt/local/bin/mpif90-mpich-devel-clang -link_info))
MPI_LIBS    = $(wordlist 2,99,$(shell /opt/local/bin/mpicxx-mpich-devel-clang -link_info)) $(wordlist 2,99,$(shell /opt/local/bin/mpif90-mpich-devel-clang -link_info))
MPI_SO_PATH = 
MPIEXEC_PATH  = 
MPIEXEC       = /opt/local/bin/mpiexec-mpich-devel-clang

