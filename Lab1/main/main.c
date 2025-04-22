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


int main(int argc, char** argv){

    if (argc < 5)
        return EXIT_FAILURE;

    if (!(K = atoi(argv[1]) && M == atoi(argv[2])))
        return EXIT_FAILURE;

    if(!M || !K)
        return EXIT_FAILURE;

    if (!(K = atof(argv[3]) && M == atof(argv[4])))
        return EXIT_FAILURE;

    t = T/(double)M;
    h = X/(double)K;

    int commsize, rank;
    double proc_time = 0;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank){
        int  *planner    = (State*)calloc(M*K, sizeof(int*));
        int  *result   = (int*)calloc(M*K, sizeof(int*));

        for()

    }
    



    return 0;
}