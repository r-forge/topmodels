#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP c_CRPS_numeric(SEXP, SEXP, SEXP, SEXP, SEXP);

const R_CallMethodDef CallEntries[] = {
  {"c_CRPS_numeric", (DL_FUNC) &c_CRPS_numeric, 5},
  {NULL, NULL, 0}
};

void R_init_topmodels(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


