#include <R.h>
#include <Rinternals.h>

#include "blas.h"
#include "dnrm3.h"
#include "fastmap.h"


SEXP R_fastmap(SEXP x)
{
  SEXP a, b;
  SEXP ret, ret_names;
  const int m = nrows(x);
  const int n = ncols(x);
  
  double *work = malloc(nthreads()*n * sizeof(*work));
  if (work == NULL)
    error("OOM");
  
  PROTECT(a = allocVector(REALSXP, n));
  PROTECT(b = allocVector(REALSXP, n));
  
  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(ret_names = allocVector(STRSXP, 2));
  
  GetRNGstate();
  fastmap(m, n, REAL(x), REAL(a), REAL(b), work);
  PutRNGstate();
  
  SET_VECTOR_ELT(ret, 0, a);
  SET_VECTOR_ELT(ret, 1, b);
  
  SET_STRING_ELT(ret_names, 0, mkChar("a"));
  SET_STRING_ELT(ret_names, 1, mkChar("b"));
  
  setAttrib(ret, R_NamesSymbol, ret_names);
  
  UNPROTECT(4);
  return ret;
}
