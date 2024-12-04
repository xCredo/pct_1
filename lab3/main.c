#include <stdio.h>
#include <math.h>
#include <omp.h>

double func(double x)
{
  double result;
    
    if (x == 0) {
        result = 1;
    } else {
        result = log(1 + x) / x;
    }
  
    return result;
}


int main()
{
  double t = omp_get_wtime();
  const double eps = 1E-5;
  const double a = -0.1;
  const double b = 1;
  const int n0 = 100000000;
  printf("Numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);
  int n = n0, k;
  double sq[2], delta = 1;
  for (k = 0; delta > eps; n *= 2, k ^= 1)
  {
    double h = (b - a) / n;
    double s = 0.0;
    for (int i = 0; i < n; i++)
      s += func(a + h * (i + 0.5));
    sq[k] = s * h;
    if (n > n0)
      delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
  }
  printf("Result Pi: %.12f; Runge rule: EPS %e, n %d\n", sq[k] * sq[k], eps, n / 2);
  t = omp_get_wtime() - t;
  printf("Elapsed time (sec.): %.6f\n", t);



  t = omp_get_wtime();
  printf("Numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);

  #pragma omp parallel num_threads(8)
  {
    int n = n0, k;
    double delta = 1;
    for (k = 0; delta > eps; n *= 2, k ^= 1)
    {
      double h = (b - a) / n;
      double s = 0.0;
      sq[k] = 0;
    // Ждем пока все потоки закончат обнуление sq[k], s
    #pragma omp barrier

    #pragma omp for nowait
      for (int i = 0; i < n; i++)
        s += func(a + h * (i + 0.5));
    #pragma omp atomic
      sq[k] += s * h;
    // Ждем пока все потоки обновят sq[k]
    #pragma omp barrier
      if (n > n0)
        delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
    #if 0
    printf("n=%d i=%d sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
    #endif
    }

    #pragma omp master
    printf("Result Pi: %.12f; Runge rule: EPS %e, n %d\n", sq[k] * sq[k], eps, n / 2);
  }
  t = omp_get_wtime() - t;
  printf("Elapsed time (sec.): %.6f\n", t);

  return 0;
}
