#include <stdio.h>
#include <math.h>
#include <omp.h>
#define _POSIX_C_SOURCE 1
#include <stdlib.h>

double getrand()
{
    return (double)rand() / RAND_MAX;
}

double getrand_p(unsigned int *seed)
{
    return (double)rand_r(seed) / RAND_MAX;
}

double func(double x, double y)
{
    return pow(exp(x + y), 2);
}

const double PI = 3.14159265358979323846;
const int n = 100000000;
int main(int argc, char **argv)
{
    int in = 0;
    double s = 0;
    double t = omp_get_wtime();
    for (int i = 0; i < n; i++)
    {
        double x = getrand() * PI; /* x in [0, pi] */
        double y = getrand();      /* y in [0, sin(x)] */
        if (y <= sin(x))
        {
            in++;
            s += func(x, y);
        }
    }
    double v = PI * in / n;
    double res = v * s / in;
    t = omp_get_wtime() - t;
    printf("Elapsed time (sec.): %.6f\n", t);
    printf("Result: %.12f, n %d\n\n", res, n);



    printf("Numerical (parallel) integration by Monte Carlo method: n = %d\n", n);
    t = omp_get_wtime();
    in = 0;
    s = 0;
    #pragma omp parallel num_threads(2)
    {
        double s_loc = 0;
        int in_loc = 0;
        unsigned int seed = omp_get_thread_num();
        #pragma omp for nowait
        for (int i = 0; i < n; i++)
        {
            double x = getrand_p(&seed) * PI; /* x in [0, pi] */
            double y = getrand_p(&seed);      /* y in [0, sin(x)] */
            if (y <= sin(x))
            {
                in_loc++;
                s_loc += func(x, y);
            }
        }
        #pragma omp atomic
        s += s_loc;
        #pragma omp atomic
        in += in_loc;
    }
    v = PI * in / n;
    res = v * s / in;
    t = omp_get_wtime() - t;
    printf("Elapsed time (sec.): %.6f\n", t);
    printf("Result: %.12f, n %d\n", res, n);
    return 0;
}
