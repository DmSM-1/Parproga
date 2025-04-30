#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>


/* 
FUNCTION DESCRIPTION 
U(t,x) = (x^2+t^2)/2
U(t,0) = phi(t) = t^2/2
U(0,x) = psi(x) = x^2/2

MAIN TASK
du/dt + du/dx = f

0 <= t <= T,    t = k*tau,  0 <= k <= K
0 <= x <= X,    x = m*h,    0 <= m <= M

MAIN APPROACH
U(k+1,m) =  f(k,m)*tau 
            - (U(k,m+1) - U(k, m-1))*tau/(2*h) 
            + (U(k,m+1) + U(k, m-1))/2

HELPER
U(k+1,m) =  f(k,m)*tau
            - (U(k,m) - U(k, m-1))*tau/h
            + U(k,m) 
*/


#define MIN(a,b) ((a<b)?a:b)


double U(double t, double x){
    return (x*x+t*t)/2;
    // return cos(x) + sin(t);
}

double phi(double t){
    return t*t/2;
    // return 1 + sin(t);
}

double psi(double x){
    return x*x/2;
    // return cos(x);
}

double f(double t, double x){
    return x + t;
    // return -sin(x) + cos(t);
}

void dump(double *data, int K, int M){
    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            printf("%3.3lf ", data[k*M+m]);
        }
        printf("\n");
    }
}


double deviation(double *data_1, double *data_2, int K, int M){
    double dev = 0.0;
    for (int k = 0; k < K; k++){
        for (int m = 0; m < M; m++){
            dev +=  pow(data_1[k*M+m] - data_2[k*M+m], 2);
        }
    }
    dev /= M*K;
    dev = sqrt(dev);

    return dev;
}


int main(int argc, char* argv[]){

    double T = atof(argv[1]);
    double X = atof(argv[2]);

    int K = atoi(argv[3]);
    int M = atoi(argv[4]);

    double tau  = T/(double)(K-1);
    double h    = X/(double)(M-1);

    int commsize, rank;

    double time = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    commsize = MIN(commsize, M);

    
    if (commsize == 1){
        double *U_km_answer = (double*)calloc(K*M, sizeof(double));
        double *U_km = (double*)calloc(K*M, sizeof(double));

        for (int k = 0; k < K; k++){
            for (int m = 0; m < M; m++){
                U_km_answer[k*M+m] = U(k*tau, m*h);
            }
        }
        time = -MPI_Wtime();

        for (int k = 0; k < K; k++){
            for(int m = 0; m < M; m++){
                if (k == 0){
                    U_km[k*M+m] = psi(m*h);
                }else if (m == 0){
                    U_km[k*M] = phi(k*tau);
                }else if (m != M-1){
                    U_km[k*M+m] = f((k-1)*tau, m*h)*tau
                                - (U_km[(k-1)*M+m+1] - U_km[(k-1)*M+m-1])*tau/(2*h)
                                + (U_km[(k-1)*M+m+1] + U_km[(k-1)*M+m-1])/2;             
                }else{
                    U_km[k*M+m] = f((k-1)*tau, m*h)*tau
                                - (U_km[(k-1)*M+m] - U_km[(k-1)*M+m-1])*tau/h
                                + U_km[(k-1)*M+m];
                }
            }
        }
        time += MPI_Wtime();

        // dump(U_km, K, M);
        // printf("\n");
        // dump(U_km_answer, K, M);
        // printf("\n");
        printf("Deviatoin: %lf\n", deviation(U_km_answer, U_km, K, M));
        printf("time:%lf\n", time);

        free(U_km);
        free(U_km_answer);

    } else if(rank < commsize){

        if (!rank)
            time = -MPI_Wtime();

        int rem = M%commsize;
        int div = M/commsize;
        int M_local = div + (rank < rem);

        int m_begin = div*rank + MIN(rank, rem);
        int m_end = m_begin + M_local;

        double *U_km_local = (double*)calloc((M_local)*K, sizeof(double));
        double *U_buf = (double*)calloc((M_local + 2), sizeof(double));

        MPI_Request recv_req_send_left;
        MPI_Request recv_req_send_right;
        MPI_Request recv_req_recv_left;
        MPI_Request recv_req_recv_right;
        MPI_Status status;

        for (int k = 0; k < K; k++){
            for (int m = 0; m < M_local; m++){
                if (k == 0){
                    U_km_local[k*M_local+m] = psi((m+m_begin)*h);
                }
                else if (m == 0 && rank==0){
                    U_km_local[k*M_local] = phi(k*tau);

                }
                else if (m == M_local-1 && rank == commsize - 1){
                    U_km_local[k*M_local+m] = f((k-1)*tau, (m+m_begin)*h)*tau
                                            - (U_buf[m+1] - U_buf[m])*tau/h
                                            + U_buf[m+1];

                // U_km[k*M+m] = f((k-1)*tau, m*h)*tau
                // - (U_km[(k-1)*M+m] - U_km[(k-1)*M+m-1])*tau/h
                // + U_km[(k-1)*M+m];
                }
                else{
                    U_km_local[k*M_local+m] = f((k-1)*tau, (m+m_begin)*h)*tau
                                - (U_buf[m+2] - U_buf[m])*tau/(2*h)
                                + (U_buf[m+2] + U_buf[m])/2;             
                }
                
                if (rank && !m){
                    MPI_Isend((void*)(U_km_local+k*M_local+m), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &recv_req_send_left);
                }

                if (rank < commsize - 1 && m == M_local-1){
                    MPI_Isend((void*)(U_km_local+k*M_local+m), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_req_send_right);
                }
            }
            if (rank)
                MPI_Irecv((void*)(U_buf), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &recv_req_recv_left);
            if (rank < commsize - 1)
                MPI_Irecv((void*)(U_buf+M_local+1), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &recv_req_recv_right);

            memcpy((void*)(U_buf+1) ,(void*)(U_km_local+M_local*k), M_local*sizeof(double));

            if (rank){
                MPI_Wait(&recv_req_recv_left, &status);
                // MPI_Wait(&recv_req_send_left, &status);
            }
            if (rank < commsize - 1){
                MPI_Wait(&recv_req_recv_right, &status);
                // MPI_Wait(&recv_req_send_right, &status);
            }

            // if (rank == 2){
            //     printf("inter %d\n", k);
            //     dump(U_buf, 1, M_local+2);
            //     printf("\n");
            //     dump(U_km_local, K, M_local);
            //     printf("\n");
            // }
        }

        if (rank){
            MPI_Isend((void*)(U_km_local), K*M_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recv_req_send_right);
        }
        else{
            double *U_km_answer = (double*)calloc(K*M, sizeof(double));
    
            for (int k = 0; k < K; k++){
                for (int m = 0; m < M; m++){
                    U_km_answer[k*M+m] = U(k*tau, m*h);
                }
            }

            double *U_km = (double*)calloc(K*M, sizeof(double));
            for (int i = 0; i < commsize; i++){
                M_local = div + (i < rem);
                m_begin = div*i + MIN(i, rem);
                if(i){
                    MPI_Irecv((void*)(U_km_local), K*M_local, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &recv_req_recv_right);
                    MPI_Wait(&recv_req_recv_right, &status);
                }
                for (int k = 0; k < K; k++){
                    memcpy((void*)(U_km+k*M+m_begin), (void*)(U_km_local+k*M_local), M_local*sizeof(double));
                }
            }
            time += MPI_Wtime();

            printf("Deviatoin: %lf\n", deviation(U_km_answer, U_km, K, M));
            printf("time: %lf\n", time);
            // dump(U_km, K, M); 
            free(U_km_answer);
            free(U_km);
        }
        // printf("rank:%d ", rank);
        // dump(U_km_local, K, M_local);
        // printf("\n");
        // if (!rank)
            // printf("time: %lf\n", time);

        // dump(U_km_local, K, M_local);

        free(U_km_local);
        free(U_buf);

    }


    MPI_Finalize();

    return 0;
}