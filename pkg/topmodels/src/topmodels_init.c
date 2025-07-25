#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP c_CRPS_numeric(SEXP, SEXP, SEXP, SEXP, SEXP);

/* Functions to retrieve CDF and quantiles from the PDF */
SEXP c_d2pq_numeric(SEXP at, SEXP p, SEXP dim, SEXP y, SEXP continuous, SEXP what);
void c_d2p_calculate(SEXP dim, SEXP x, double* at, double* pdf, double* res, bool do_cdf);
void c_d2q_calculate(SEXP dim, SEXP x, double* at, double* pdf, double* res, bool do_cdf);

/* Helper functions used for 2d interpolation and numeric integration */
double interpolate_linear(double xlo, double xhi, double ylo, double yhi, double x);
double integrate_2d(double xlo, double xhi, double ylo, double yhi, double x);

/* Custom type: stuctured object with ...
 * a double vector and length.
 * NOTE: .values must be freed by the user! */
typedef struct {
    double* values;
    int length;
} doubleVec;

const R_CallMethodDef CallEntries[] = {
  {"c_CRPS_numeric", (DL_FUNC) &c_CRPS_numeric, 5},
  {"c_d2pq_numeric", (DL_FUNC) &c_d2pq_numeric, 6},
  {NULL, NULL, 0}
};

void R_init_topmodels(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


