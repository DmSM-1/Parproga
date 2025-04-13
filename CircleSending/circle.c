#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char **argv){
    int val = 0;
    int commsize, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    if (!rank){
        if (argc > 1)
            val = atoi(argv[1]);

        printf("Proc: %2d Val(int): %6d Send -> %2d\n", rank, val, (rank+1)%commsize);

        MPI_Send((void*)&val, 1, MPI_DOUBLE, (rank+1)%commsize, 0, MPI_COMM_WORLD);
        MPI_Recv((void*)&val, 1, MPI_DOUBLE, (rank-1)%commsize, 0, MPI_COMM_WORLD, &status);

        printf("Proc: %2d Val(int): %6d (result)\n", rank, val);
    }
    else{
        MPI_Recv((void*)&val, 1, MPI_DOUBLE, (rank-1)%commsize, 0, MPI_COMM_WORLD, &status);

        val *= 2;
        printf("Proc: %2d Val(int): %6d Send -> %2d\n", rank, val, (rank+1)%commsize);

        MPI_Send((void*)&val, 1, MPI_DOUBLE, (rank+1)%commsize, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}