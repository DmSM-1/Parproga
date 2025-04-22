#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//answer = sin(x) + cos(t)
//f = cos(x) - sin(t)


int K, M;
double X,T, t, h;

typedef struct State 
{
    int j;
    int k;
    int s;
    double u[3];
}State;



double psi(double x){
    double psi_val = sin(x) + 1;

    return psi_val;
}


double phi(double t){
    double phi_val = cos(t);

    return phi_val;
}


double fanc(double t, double x){
    double f_val = cos(x) - sin(t);

    return f_val;
}


double main_method(double** net, int j, int m){
    return 0;
}


double print_res(double* f){
    printf("\n");
    for(int i = 0, pos = 0; i < M; i++){
        for(int j = 0; j < K; j++, pos++){
            printf("%3.3lf ", f[pos]);
        }
        printf("\n");
    }
}


int main(int argc, char** argv){

    if (argc < 5)
    return EXIT_FAILURE;
    
    if (!(T = atof(argv[1])))
        return EXIT_FAILURE;

    if (!(X = atof(argv[2])))
        return EXIT_FAILURE;

    if (!(K = atoi(argv[3])))
        return EXIT_FAILURE;

    if (!(M = atoi(argv[4])))
        return EXIT_FAILURE;

    if(!M || !K)
        return EXIT_FAILURE;

    t = T/(double)M;
    h = X/(double)K;

    int commsize, rank;
    double proc_time = 0;

    K++;
    M++;
    int total = K*M;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Status status;

    int executors = commsize-1;
    if (executors > M)
        executors = M;

    if (!rank){
        //stage 1
        proc_time -= MPI_Wtime();
        double *f = (double*)calloc(total, sizeof(double));

        int cycles = M / executors;
        int overload = M % executors;
        int executors_num = 0;

        for(int i = 0; i < total; i++)
            f[i] = 0;

        for(int i = 0, T = 0, X = 0; i < cycles + (overload>0); i++){
            executors_num = (i==cycles)?overload:executors;
            for(int j = 0; j < executors_num; j++, T++){
                MPI_Send((void*)&T, 1, MPI_INTEGER, (T)%executors+1, 0, MPI_COMM_WORLD);
                MPI_Send((void*)&X, 1, MPI_INTEGER, (T)%executors+1, 0, MPI_COMM_WORLD);
            }
            T -= executors_num;
            for(int j = 0; j < executors_num; j++, T++){
                MPI_Recv((void*)&(f[T]), 1, MPI_DOUBLE, (T)%executors+1, 0, MPI_COMM_WORLD, &status);
                //print_res(f);
            }
                
        }
        
        free(f);
        proc_time += MPI_Wtime();
        printf("execution time:%lf\n", proc_time);
    }

    else if(rank < M+2){
        int cycles = M / executors + (rank <= M%executors);
        int X = 0, T = 0;
        double f_val;
        for(int i = 0; i < cycles; i++){
            MPI_Recv((void*)&(T), 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv((void*)&(X), 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
            f_val = phi(T/t);
            MPI_Send((void*)&f_val, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

    }
    

    MPI_Finalize();

    return 0;
}