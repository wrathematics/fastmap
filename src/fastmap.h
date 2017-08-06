#ifndef FASTMAP_FASTMAP_
#define FASTMAP_FASTMAP_


#include "safeomp.h"


static inline int sample(const int min, const int max)
{
  return min + (int) ((double)(max-min + 1) * unif_rand());
}



struct CustomMax
{
  double value;
  int index;
};
#pragma omp declare reduction(maximum : struct CustomMax : omp_out = (omp_in.value>omp_out.value ? omp_in : omp_out))

static inline void find_most_distant(const int m, const int n, const double *const restrict x, double *restrict a, double *restrict b, double *restrict work)
{
  struct CustomMax max;
  max.value = 0.0;
  max.index = 0;
  
  #pragma omp parallel for default(shared) reduction(maximum:max) if(m*n>OMP_MIN_SIZE)
  for (int i=0; i<m; i++)
  {
    const int tid = omp_get_thread_num();
    SAFE_SIMD
    for (int j=0; j<n; j++)
      work[j + tid*n] = -x[i + m*j];
    
    daxpy_(&n, &(double){1.0}, b, &(int){1}, work + tid*n, &(int){1});
    double tmp = dnrm3(n, work + tid*n, 2);
    
    if (tmp > max.value)
    {
      max.value = tmp;
      max.index = i;
    }
  }
  
#ifdef OMP_VER_4
  #pragma omp parallel for simd if(n>OMP_MIN_SIZE)
#else
  #pragma omp parallel for      if(n>OMP_MIN_SIZE)
#endif
  for (int j=0; j<n; j++)
    a[j] = x[max.index + m*j];
}



static inline void fastmap(const int m, const int n, const double *const restrict x, double *const restrict a, double *const restrict b, double *restrict work)
{
  // Take random row b in x;
  const int index = sample(0, m-1);
  
  for (int j=0; j<n; j++)
    b[j] = x[index + m*j];
  
  // a = most distant point in x from b
  find_most_distant(m, n, x, a, b, work);
  
  // b = the most distant point in x from a
  find_most_distant(m, n, x, b, a, work);
}


#endif
