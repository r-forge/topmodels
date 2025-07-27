
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>
#include "topmodels.h"


/* Calculating weights
 * 
 * @param i   index of distribution (row index for 'at').
 * @param n   number of distributions (rows).
 * @param pdf numeric matrix of dimension (n x w->length) with the PDFs.
 *
 * No return, overwrites w->values, where
 * the sum over w->values equals 1.0 */
void c_d2moments_calculate_weights(int i, int n, double* pdf, doubleVec* w) {
    double sump = 0.0;
    int k;
    // Sum up PDF
    for (k = 0; k < w->length; k++) { sump += pdf[i + n * k]; }
    // Calculating weights
    for (k = 0; k < w->length; k++) { w->values[k] = pdf[i + n * k] / sump; }
}

/* Calculating weighted mean based on PDF.
 *
 * @param i    index of distribution (row index for 'at').
 * @param n    number of distributions (rows).
 * @param at   numeric matrix of grid position (n x w->length) where the PDF was
 *             provided. The weights are stored on 'w'.
 * @param w    doubleVec with weights, also contains the number of columns in 'at' on
 *             its attribute w->length.
 * @param what 1 = mean, 2 = variance, 3 = skewness, 4 = kurtosis
 * Returns single numeric value.
 */
double c_d2moments_calculate_moment(int i, int n, double* at, doubleVec* w, int what) {

    // Calculating mean
    double mean = 0.0;
    for (int k = 0; k < w->length; k++) {
       mean += w->values[k] * at[i + n * k];
    }
    if (what == 1) { return mean; }

    // Calculating variance
    double variance = 0.0;
    for (int k = 0; k < w->length; k++) {
       variance += w->values[k] * pow(at[i + n * k] - mean, 2.0);
    }
    if (what == 2) { return variance; }

    // If the variance is ~0 we return NA for skewness/kurtosis
    Rprintf(" variance = %.10f\n", variance);
    if (variance < 1e-9) {
        return R_NaN;
    }

    // Calculating skewness or kurtosis
    double p = (double)what;   // Define power parameter (3 for skewness, 4 = kurtosis)
    variance = sqrt(variance); // Convert to standard deviation
    double res = 0.0;
    Rprintf("  p = %.3f\n", p);
    for (int k = 0; k < w->length; k++) {
       res += w->values[k] * pow((at[i + n * k] - mean) / variance, p);
    }
    // Calculate excess kurtosis if needed
    if (what == 4) { res = res - 3.0; }

    return res;
}

/*
 * Numerically integrating the density to evaluate the cummulative distribution function (CDF)
 * given we know nothing else than the PDF of a distribution.
 * 
 * at      numeric vector, grid points along the x-axis of the distribution.
 * p       the density at grid points 'x'
 * int     integer of length 2, dimension of matrices at/p
 * x       numeric, where to evaluate the CDF.
 * whatint integer which defines what to calculate (hardcoded).
 *         1 = mean, 2 = variance, ... TODO(R) extend.
 * 
 * Dimensions
 * ----------
 * Input 'dim' contains the dimension of the matrices, whereof we will
 * use dim[0] = n and dim[1] = k, i.e., assuming matrix dimension (n x k).
 */
SEXP c_d2moments_numeric(SEXP at, SEXP p, SEXP dim, SEXP continuous, SEXP whatint) {

    // Setting up pointers for at, pat, dim, x
    double *atptr = REAL(at);
    double *pptr  = REAL(p);

    // Getting matrix dimension
    int *dimptr = INTEGER(dim);
    int n = dimptr[0], k = dimptr[1];

    // What to calculate?
    int what = INTEGER(whatint)[0];

    // Loop index
    int i;
    int nprotected = 0;

    // Vector to store weights
    doubleVec weights;
    weights.values = (double*)malloc(k * sizeof(double));
    weights.length = k;

    // Setting up results vector (length n)
    SEXP rval;
    PROTECT(rval = allocVector(REALSXP, n));
    double *rvalptr = REAL(rval);
    nprotected += 1;

    // Looping trough all observations i = 1, ..., n;
    for (i = 0; i < n; i++) {
        // Calculating weights (overwrites weights.values)
        c_d2moments_calculate_weights(i, n, pptr, &weights);

        // Calculate requested moment
        rvalptr[i] = c_d2moments_calculate_moment(i, n, atptr, &weights, what);
    }

    // Releasing allocated memory
    free(weights.values);

    // Releasing protected object(s) and return
    UNPROTECT(nprotected);
    return rval;
}

