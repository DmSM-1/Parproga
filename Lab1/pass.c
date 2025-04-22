#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define TEST_CONST 10000

int main(int argc, char** argv){
    double total_time, mean_time, buf = 0;
    int commsize, rank;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0){
        for (int i = 0; i < TEST_CONST; i++){
            total_time = -MPI_Wtime();

            MPI_Send((void*)&buf, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv((void*)&buf, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);

            total_time += MPI_Wtime();
            mean_time += total_time;
        }
        mean_time /= TEST_CONST;
        printf("ping: %lf\n", mean_time);
    }

    if (rank == 1){
        for(int i = 0; i < TEST_CONST; i++){
            MPI_Recv((void*)&buf, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Send((void*)&buf, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
}