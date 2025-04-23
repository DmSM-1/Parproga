#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

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

    double time, x_cor = 0;
    double answer = 0;

    if(argc == 7){
        time = atof(argv[5]);
        x_cor = atof(argv[6]);
    }

    MPI_Status status;
    MPI_Request recv_req;

    int proc_len = (int)(M/(commsize-1));
    int proc_tail = M%(commsize-1);
    int exec_num = (M>commsize-1)?(commsize-1):M;

    if (rank==commsize-1){
        double proc_time = -MPI_Wtime();
        double *res = (double*)calloc(total, sizeof(double));
        int num = 0;

        for(int k = 0; k < K; k++){
            for(int m = 0, proc = 0; proc < exec_num; proc++){
                num = (proc<proc_tail)?(proc_len+1):(proc_len);
                MPI_Irecv((void*)(res+k*M+m), num, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &recv_req);
                MPI_Wait(&recv_req, &status);
                m += num;
            }
        }

        proc_time += MPI_Wtime();
        double err = error(res);
        printf("time: %lf error: %lf\n", proc_time, err);
        
        free(res);

    }else if(rank < exec_num){
        int fix = (rank < proc_tail)?rank:proc_tail;
        int begin = rank*proc_len+fix;
        int end = begin+proc_len+(rank<proc_tail);
        bool first, last, odd;
        
        first = (rank == 0);
        last = (rank == exec_num - 1);
        odd = rank % 2;
        
        //if (rank == exec_num - 1)
        //end+=proc_tail;
        
        int iter = end-begin;
        //printf("rank %d %d %d\n", rank, begin, end);
        
        double *data = (double*)calloc(iter, sizeof(double));
        double *previos = (double*)calloc(iter+2, sizeof(double));
        int m = 0;

        for(int k = 0; k < K; k++){
            for(m = 0; m < iter; m++){
                if(k==0){
                    data[m] = phi((double)(m+begin)/(M-1));
                }else if(m+begin==0){
                    data[m] = psi((double)(k)/(K-1));
                }else if(m+begin==K-1){
                    data[m] = (fanc((double)m/(M-1),(double)(k-1)/(K-1))-(previos[m+1]-previos[m])/h)*t+previos[m+1];
                }else{
                    data[m] = (fanc((double)m/(M-1),(double)(k-1)/(K-1))-(previos[m+2]-previos[m])/2/h)*t+0.5*(previos[m+2]+previos[m]);
                }
            }
            memcpy((void*)(previos+1), (void*)(data), iter*sizeof(double));
            MPI_Isend((void*)(data), iter, MPI_DOUBLE, commsize-1, 0, MPI_COMM_WORLD, &recv_req);

                
            if(!odd){
                if(!first)
                    MPI_Send((void*)(data), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                if(!last){
                    MPI_Send((void*)(data+iter-1), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
                    MPI_Recv((void*)(previos+iter+1), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
                }
                if(!first)
                    MPI_Recv((void*)(previos), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);

            }else{
                MPI_Recv((void*)(previos), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
                if(!last)
                    MPI_Recv((void*)(previos+iter+1), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
                MPI_Send((void*)(data), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                if(!last)
                    MPI_Send((void*)(data+iter-1), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            }

            MPI_Wait(&recv_req, &status);
        }

        //printf("rank:%d ", rank);
        //for(int i = 0; i < iter+2;i++){
        //    printf("%lf ", previos[i]);
        //}
        //printf("\n");
    
        free(previos);
        free(data);
        
    }

    MPI_Finalize();
    return 0;
}