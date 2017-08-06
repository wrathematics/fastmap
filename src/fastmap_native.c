/* Automatically generated. Do not edit by hand. */
  
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP R_cma(SEXP x, SEXP k_);
extern SEXP R_fastmap(SEXP x);

static const R_CallMethodDef CallEntries[] = {
  {"R_cma", (DL_FUNC) &R_cma, 2},
  {"R_fastmap", (DL_FUNC) &R_fastmap, 1},
  {NULL, NULL, 0}
};

void R_init_fastmap(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
