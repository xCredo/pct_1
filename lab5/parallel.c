#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 100000000
#define THRESHOLD 1000
#define THREADS 2

void swap(int *x, int *y)
{
    int tmp = *x;
    *x = *y;
    *y = tmp;
}

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void partition(int *v, int *i, int *j, int low, int high)
{
    *i = low;
    *j = high;
    int pivot = v[(low + high) / 2];
    do
    {
        while (v[*i] < pivot)
            (*i)++;
        while (v[*j] > pivot)
            (*j)--;
        if (*i <= *j)
        {
            swap(&(v[*i]), &(v[*j]));
            (*i)++;
            (*j)--;
        }
    } while (*i <= *j);
}

void quicksort_tasks(int *v, int low, int high)
{
    int i, j;
    partition(v, &i, &j, low, high);
    if (high - low < THRESHOLD || (j - low < THRESHOLD || high - i < THRESHOLD))
    {
        if (low < j)
            quicksort_tasks(v, low, j);
        if (i < high)
            quicksort_tasks(v, i, high);
    }
    else
    {
#pragma omp task untied
        {
            quicksort_tasks(v, low, j);
        }
        quicksort_tasks(v, i, high);
    }
}

void quicksort(int *v, int low, int high)
{
    int i, j;
    // print_arr(v);
    partition(v, &i, &j, low, high);
    if (low < j)
        quicksort(v, low, j);
    if (i < high)
        quicksort(v, i, high);
}

void init(int **arr)
{
    for (int i = 0; i < N; i++)
        (*arr)[i] = rand() % 100;
}

void print_arr(int *arr)
{
    for (int i = 0; i < N; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

int main()
{
    int *arr = malloc(sizeof(int) * N);
    init(&arr);
    // print_arr(arr);
    double t = wtime();
    quicksort(arr, 0, N - 1);
    t = wtime() - t;
    printf("%lf - время последовательной программы\n", t);
// quicksort(arr, 0, N - 1);
    for (int i = 2 ;i <= 8; i+=2)
    {
        init(&arr);
        double time = wtime();
    #pragma omp parallel num_threads(i)
        {
        #pragma omp single
            quicksort_tasks(arr, 0, N - 1);
        }
    time = wtime() - time;
    printf("время работы паралел прог - %lf, потоков - %d  speedup:%lf\n", time,i, t/time);
    }
    return 0;
}
