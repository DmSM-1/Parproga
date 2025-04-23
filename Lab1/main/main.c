#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//answer = sin(x) + cos(t)
//f = cos(x) - sin(t)
//xt
//x+t

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


double error(double *f){
    double err = 0;
    for(int i = 0, pos = 0; i < K; i++){
        for(int j = 0; j < M; j++, pos++){
            err += pow(f[pos]-(cos(i/(K-1))+sin(j/(M-1))),2);
        }
        err /= M*K;
        return sqrt(err);
    }
}


double print_res(double* f){
    for(int i = 0, pos = 0; i < K; i++){
        for(int j = 0; j < M; j++, pos++){
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

    t = T/(double)(M-1);
    h = X/(double)(K-1);


    int commsize, rank;
    double proc_time = 0;

    int total = K*M;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Status status;

    int executors = commsize-1;
    if (executors > M)
        executors = M;

    double *data = (double*)calloc((M/executors+1)*4, sizeof(double));
    double *res = (double*)calloc((M/executors+1), sizeof(double));

    if (!rank){
        //stage 1
        proc_time -= MPI_Wtime();
        double *f = (double*)calloc(total, sizeof(double));

        int cycles = M / executors;
        int overload = M % executors;
        int executors_num = 0;
        int vals = 0;
        for(int i = 0; i < total; i++)
            f[i] = 0;

        
        for(int l = 0; l < K; l++){
            for(int i = 0; i < executors; i++){
                int vals = cycles + (i < overload);
                for(int j = 0; j < vals; j++){
                    data[4*j]   = (double)(executors*j+i);
                    data[4*j+1] = (double)l;
                    int pos = executors*j+i;
                    if (l > 0 && pos!=0 && pos!=M){
                        data[4*j+2] = f[(l-1)*M+j*executors+i-1];
                        data[4*j+3] = f[(l-1)*M+j*executors+i+1];
                    }
                }
                MPI_Send((void*)data, vals*4, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
            }
            for(int i = 0; i < executors; i++){
                int vals = cycles + (i < overload);
                MPI_Recv((void*)res, vals, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &status);
                for (int j = 0; j < vals; j++){
                    f[l*M+j*executors+i] = res[j];
                }
            }
        }
        double err = error(f);
        free(f);
        proc_time += MPI_Wtime();
        printf("execution time:%lf error %lf\n", proc_time, err);
    }

    else if(rank < M+2){
        int cycles = M / executors + (rank <= M%executors);
        for(int l = 0; l < K; l++){
            MPI_Recv((void*)data, cycles*4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            for(int i = 0; i < cycles; i++){
                if (data[4*i+1] < 1.0){
                    res[i] = phi(data[4*i]/(M-1));
                }else if(data[4*i] < 1.0){
                    res[i] = psi(data[4*i+1]/(K-1));
                }else if(data[4*i] < M-1){
                    res[i] = (fanc(data[4*i]/(M-1), data[4*i+1]/(K-1)));
                    res[i] -= (data[4*i+3]-data[4*i+2])/2/h;
                    res[i] *= t;
                    res[i] += (data[4*i+3]+data[4*i+2])/2;
                }else{
                    res[i] = 0.0;
                }
            }
            MPI_Send((void*)res, cycles, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

    }
    

    MPI_Finalize();
    free(res);
    free(data);
    return 0;
}