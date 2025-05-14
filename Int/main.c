#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>


/*
FOR THIS TASK 
f(x) = cos(1/(x+5))
x = [-4,99; 2] 
*/


#define abslf(x) (((x)>0.0)?(x):(-(x)))
// #define DEBUG


typedef struct Task{
    double A;
    double B;
    double res;
    int stat;
}task;


typedef struct TaskQueue tQ;



typedef struct TaskQueue{
    size_t cap;
    size_t size;
    int exec;
    task** tasks;
    double result;
    int (*add)(tQ* ,double, double);
    task* (*get)(tQ*);
    void (*remove)(tQ*, task);
}tQ;

pthread_mutex_t mutex_queue;
pthread_mutex_t mutex_exec;
tQ *queue;
double acc;
int exec_task = 0;

int add_task(tQ *queue, double A, double B){
    if (queue->size >= queue->cap)
        return 0;

    queue->tasks[queue->size] = (task*)malloc(sizeof(task));
    queue->tasks[queue->size]->A = A;
    queue->tasks[queue->size]->B = B;
    queue->tasks[queue->size]->res = 0;
    queue->tasks[queue->size]->stat = 0;

    queue->size++;

    return 1;
}


task *get_task(tQ *queue){
    if (queue->size == 0)
        return NULL;
    
    queue->size--;
    queue->exec++;

    return queue->tasks[queue->size];
}


tQ *init_Queue(int cap){
    tQ* queue = (tQ*)malloc(sizeof(tQ));

    queue->cap = cap;
    queue->size = 0;
    queue->exec = 0;
    queue->result = 0.0;
    queue->tasks = (task**)calloc(queue->cap, sizeof(task*));
    queue->add = add_task;
    queue->get = get_task;
    return queue;
}


void *destroy_Queue(tQ *queue){
    for (int i = 0; i < queue->size; i++)
        free(queue->tasks[i]);
    free(queue);
}


double fanc(double x){
    return cos(1/(x+5));
}


double integ(double A, double B, double acc, void* external_f, int deep){

    double (*f)(double) = (double (*)(double))external_f;
    
    double J =  0;
    double C = (A+B)/2;
    
    double fa = f(A);
    double fb = f(B);
    double fc = f(C);
    
    double Sab = (fa+fb)*(B-A)/2;
    double Sac = (fa+fc)*(C-A)/2;
    double Scb = (fb+fc)*(B-C)/2;

    #ifdef DEBUG
        printf("\rA=%3.3lf deep:%d", A, deep);
    #endif

    if ((abslf(Sab-Sac-Scb) >= acc*abslf(Sac+Scb)) && deep < 25000)
        return integ(A, C, acc, external_f, deep+1) + integ(C, B, acc, external_f, deep+1);

    return Sac+Scb;

}


void* pthread_fanc(void* args){
    double (*f)(double) = (double (*)(double))fanc;

    int exec, size, new_task;
    task* task;

    int local_buf_cap = 16;
    int local_buf_size = 0;

    double A_buf[16];
    double B_buf[16];

    double A, B, C;
    double fa, fb, fc, Sab, Sac, Scb;
    int external = 0;
    double buf = 0;
    double sm_int = pow(2, -50);

    while (1){
        if (local_buf_size){
            local_buf_size--;
            A = A_buf[local_buf_size];
            B = B_buf[local_buf_size];
            external = 0;
        }else{


            pthread_mutex_lock(&mutex_queue);
                exec = queue->exec;
                size = queue->size;
                task = queue->get(queue);
            pthread_mutex_unlock(&mutex_queue);
            external = 1;
            
            if (exec == 0 && size == 0)
            break;
            
            if (task == NULL){
                usleep(1);
                continue;
            }

            A = task->A;
            B = task->B;

            free(task);
        }

        while (1){
            C = (A+B)/2;
            fa = f(A);
            fb = f(B);
            fc = f(C);
            
            Sab = (fa+fb)*(B-A)/2;
            Sac = (fa+fc)*(C-A)/2;
            Scb = (fb+fc)*(B-C)/2;

            if (abslf(Sab-Sac-Scb) >= acc*abslf(Sac+Scb) && B-A > sm_int){

                if (local_buf_size < local_buf_cap){
                    A_buf[local_buf_size] = C;
                    B_buf[local_buf_size] = B;
                    local_buf_size++;
                }else{
                    pthread_mutex_lock(&mutex_queue);
                        new_task = (queue->add(queue, C, B));
                    pthread_mutex_unlock(&mutex_queue);

                    if (!new_task){
                        buf += integ(C, B, acc, fanc, 0);
                    }
                }
                B = C;
            }
            else{
                buf += Sac+Scb;
                if (external){
                    pthread_mutex_lock(&mutex_queue);
                        queue->exec -= external;
                    pthread_mutex_unlock(&mutex_queue);
                }
                break;
            }
            
        }
        
    }
    pthread_mutex_lock(&mutex_queue);
        queue->result += buf;
    pthread_mutex_unlock(&mutex_queue);
    return NULL;
}


int main(int argc, char* argv[]){ 

    if (argc < 5) return EXIT_FAILURE;

    double A, B;
    int threads_num;
    
    if(!(A = atof(argv[1]))) return EXIT_FAILURE;
    if(!(B = atof(argv[2]))) return EXIT_FAILURE;
    if(!(acc = atof(argv[3]))) return EXIT_FAILURE;
    if(!(threads_num = atoi(argv[4]))) return EXIT_FAILURE;

    struct timespec start, end;

    clock_gettime(CLOCK_MONOTONIC, &start);

        
        queue = init_Queue(102400);

        for (int i = 0; i < threads_num; i++){
            queue->add(queue, A+(B-A)/threads_num*i, A+(B-A)/threads_num*(i+1));
        }

        pthread_mutex_init(&mutex_queue, NULL);
        pthread_mutex_init(&mutex_exec, NULL);
        pthread_t *threads = (pthread_t*)malloc(threads_num * sizeof(pthread_t));

        for (int i = 0; i < threads_num; i++) {
            pthread_create(&threads[i], NULL, pthread_fanc, NULL);
        }

        for (int i = 0; i < threads_num; i++) {
            pthread_join(threads[i], NULL);
        }

    clock_gettime(CLOCK_MONOTONIC, &end);


    printf("%lf\n", queue->result);
    printf("%lf\n", (end.tv_sec-start.tv_sec)+1e-9*(end.tv_nsec-start.tv_nsec));

    pthread_mutex_destroy(&mutex_queue);
    pthread_mutex_destroy(&mutex_exec);

    destroy_Queue(queue);

    return 0;
}