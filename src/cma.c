#include <R.h>
#include <Rinternals.h>

#include "blas.h"
#include "dnrm3.h"
#include "fastmap.h"

#define SIGN(x) ((x)>0?1:((x)==0?0:-1))

#define FREE(x) {if(x)free(x);}

#define BAD_DIMS -12
#define BAD_K -2
#define OOM 1


static inline int sweep(const int m, const int n, double *restrict x, double *restrict vec)
{
  if (m == 0 || n == 0)
    return 0;
  
  #pragma omp parallel for if(m*n>OMP_MIN_SIZE)
  for (int j=0; j<n; j++)
  {
    SAFE_SIMD
    for (int i=0; i<m; i++)
      x[i + m*j] -= vec[j];
  }
  
  return 0;
}



static inline int cma(const int m, const int n, const int k, double *const restrict x)
{
  char trans = 'N';
  int ncol = n;
  
  if (m < 1 || n < 1)
    return BAD_DIMS;
  if (k > n)
    return BAD_K;
  
  double *a = malloc(n * sizeof(*a));
  double *b = malloc(n * sizeof(*b));
  double *work = malloc(nthreads()*n * sizeof(*work));
  double *y = malloc(m * sizeof(*y));
  if (a == NULL || b == NULL || work == NULL || y == NULL)
  {
    FREE(a);FREE(b);FREE(work);FREE(y);
    return OOM;
  }
  
  
  for (int i=0; i<k; i++)
  {
    double *const restrict x_i = x + (m*i);
    
    // Select pivot row pair
    fastmap(m, ncol, x_i, a, b, work);
    
    // Translate rows to pivot origin
    sweep(m, ncol, x_i, a);
    
    // Apply householder reflection on right
    daxpy_(&ncol, &(double){-1.0}, a, &(int){1}, b, &(int){1});
    
    const double wnorm = dnrm3(ncol, b, 2);
    b[0] += SIGN(b[0]) * wnorm;
    double vnorm = dnrm3(ncol, b, 2);
    
    SAFE_FOR_SIMD
    for (int j=0; j<ncol; j++)
      b[j] /= vnorm;
    
    dgemv_(&trans, &m, &ncol, &(double){1.0}, x_i, &m, b, &(int){1}, &(double){0.0}, y, &(int){1});
    
    dger_(&m, &ncol, &(double){-2.0}, y, &(int){1}, b, &(int){1}, x_i, &m);
    
    ncol--;
  }
  
  
  free(y);
  free(work);
  free(b);
  free(a);
  
  return 0;
}



SEXP R_cma(SEXP x, SEXP k_)
{
  SEXP ret;
  const int m = nrows(x);
  const int n = ncols(x);
  const int k = INTEGER(k_)[0];
  
  PROTECT(ret = allocMatrix(REALSXP, m, k));
  
  double *cpy = malloc((size_t)m*n * sizeof(*cpy));
  if (cpy == NULL)
    error("OOM");
  
  memcpy(cpy, REAL(x), (size_t)m*n*sizeof(*cpy));
  
  GetRNGstate();
  int info = cma(m, n, k, cpy);
  PutRNGstate();
  
  if (info != 0)
    error("TODO\n");
  
  memcpy(REAL(ret), cpy, (size_t)m*k * sizeof(*cpy));
  free(cpy);
  
  UNPROTECT(1);
  return ret;
}
