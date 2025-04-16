#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define NUM_THREADS 4

double result = 0.0;
pthread_mutex_t lock_r;


typedef struct _thread_data_t{
    int tid;
    int begin;
    int load;
}thread_data_t;


void* thr_func(void* arg){
    thread_data_t *thr_data = (thread_data_t*)arg;
    double sum = 0.0;

    for(int i = 0; i < thr_data->load; i++)
        sum += 1.0/(double)(i+thr_data->begin+1);

    printf("Thr%d: sum:%lf load:%d begin:%d\n", thr_data->tid, sum, thr_data->load, thr_data->begin);
    
    // pthread_mutex_lock(&lock_r);
    //     result += sum;
    //     printf("Thr%d: sum:%lf load:%d begin:%d\n", thr_data->tid, result, thr_data->load, thr_data->begin);
    // pthread_mutex_unlock(&lock_r);

    double *r = (double*)malloc(sizeof(double));
    *r = sum;
    pthread_exit((void*)r);
}


int main(int argc, char **argv){
    int N, num_thr = 0;
    
    if (argc < 2)
        return EXIT_FAILURE;
    
    if (!(N = atoi(argv[1])))
        return EXIT_FAILURE;

    if (argc < 3){
        num_thr = NUM_THREADS;
    }
    else if (!(num_thr = atoi(argv[2]))){
        num_thr = NUM_THREADS;
    }

    pthread_t *thr          = calloc(num_thr, sizeof(pthread_t));
    thread_data_t *thr_data = calloc(num_thr, sizeof(thread_data_t));

    int load = N/num_thr;
    int overload = N%num_thr;
    
    pthread_mutex_init(&lock_r, NULL);

    for(int i = 0, localN = 0; i < NUM_THREADS; i++){
        thr_data[i].tid     = i;
        thr_data[i].load    = load+(overload>i);
        thr_data[i].begin   = localN;

        pthread_create(&thr[i], NULL, thr_func, &thr_data[i]);
        localN += load+(overload>i);
    }

    for(int i = 0; i < NUM_THREADS; i++){
        void* ret;
        pthread_join(thr[i], &ret);
        result += *(double*)ret;
        free((double*)ret);
    }

    printf("Thr MAIN: sum:%lf\n", result);

    pthread_mutex_destroy(&lock_r);

    free(thr);
    free(thr_data);
}