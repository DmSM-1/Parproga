#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>


#ifdef DEBUG
#endif


int main(int argc, char **argv){
    int commsize, rank;
    int N, load, overload = 0;
    double sum = 0;
    double proc_time = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc == 1){
        if (!rank) printf("There is no N val\n");
    }
    else if (!(N = atoi(argv[1]))){
        if (!rank) printf("N val is ZERO\n");
    }
    else{
        proc_time -= MPI_Wtime();
        overload = N%commsize;
        load = N/commsize;
        
        int start = rank*load + 1;
        int end = start + load;
        for(int i = start; i < end; i++)
            sum += 1.0/i;

        if (rank+1 == commsize)
            for(int i = end; i <= N; i++)
                sum += 1.0/i;

        #ifdef DEBUG
        printf("Sum: %lf Proc: %d\n", sum, rank);
        #endif

        if (rank){
            MPI_Send((void*)&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        else{
            double buf = 0;
            MPI_Status status;
            for (int i = 1; i < commsize; i++){
                MPI_Recv((void*)&buf, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                sum += buf;
            }
            proc_time += MPI_Wtime();
            printf("Result: %lf Time: %lf\n", sum, proc_time);
        }
        
    }
    
    MPI_Finalize();
}