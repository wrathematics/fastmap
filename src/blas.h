#ifndef FASTMAP_BLAS_
#define FASTMAP_BLAS_


void daxpy_(const int *const restrict n, const double *const restrict a,
  const double *const restrict x, const int *const restrict incx,
  const double *const restrict y, const int *const restrict incy);

void dgemv_(const char *trans, const int *m, const int *n, const double *restrict alpha, 
  const double *restrict a, const int *lda, const double *restrict x, const int *incx, 
  const double *restrict beta, double *restrict y, const int *incy);

void dger_(const int *m, const int *n, const double *restrict alpha, const double *restrict x, 
  const int *incx, const double *restrict y, const int *incy, double *restrict a, 
  const int *lda);


#endif
