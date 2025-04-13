#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    int N = atoi(argv[1]);
    printf("eee: %d\n", N);

    MPI_Finalize();
}