
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>
#include "topmodels.h"


/* Calculating moments based on a PDF vector.
 *
 * @param i    index of distribution (row index for 'at').
 * @param dim  dimension of matrices 'at' and 'p'.
 * @param at   numeric matrix of grid position (n x w->length) where the PDF was
 *             provided. The weights are stored on 'w'.
 * @param p    numeric matrix with PDF evaluated at 'at', same dimension as 'at'.
 * @param what 1 = mean, 2 = variance, 3 = skewness, 4 = kurtosis
 * Returns single numeric value.
 */
double c_d2moments_calculate_moment(int i, int* dim, double* at, double* p, int what) {

    // Extracting matrix dimensions for convenience
    int n = dim[0], nk = dim[1];

    // Temporary vector to store width of intervals
    doubleVec width;
    width.values = (double*)malloc((nk - 1) * sizeof(double));
    width.length = nk - 1;

    // Normalization using Trapezoidal Rule
    double cdfint = 0.0, avgp;
    for (int k = 0; k < (nk - 1); k++) {
        width.values[k] = at[i + n * (k + 1)] - at[i + n * k]; // Width of the trapezoid
        avgp            = (p[i + n * k]  + p[i + n * (k + 1)])  * 0.5; // Average pdf
        cdfint += width.values[k] * avgp;
    }
    if (cdfint <= 1e-9) { error("Error: The area under the 'pdf' curve is zero or negative."); }

    // Calculating mean
    double y1, y2;
    double mean = 0.0;
    for (int k = 0; k < (nk - 1); k++) {
        y1 = at[i + n * k]       * p[i + n * k];
        y2 = at[i + n * (k + 1)] * p[i + n * (k + 1)];
        mean += width.values[k] * (y1 + y2) * 0.5;
    }
    mean = mean / cdfint; // Normalize
    if (what == 1) { return mean / cdfint; }


    // Calculating variance
    double variance = 0.0;
    for (int k = 0; k < (nk - 1); k++) {
        y1 = pow(at[i + n * k]       - mean, 2.0) * p[i + n * k];
        y2 = pow(at[i + n * (k + 1)] - mean, 2.0) * p[i + n * (k + 1)];
        variance += width.values[k] * (y1 + y2) * 0.5;
    }
    variance = variance / cdfint; // Normalize
    if (what == 2) { return variance / cdfint; }

    // Calculating skewness or kurtosis
    double powexp = (double)what;   // Define power parameter (3 for skewness, 4 = kurtosis)
    variance = sqrt(variance); // Convert to standard deviation
    double res = 0.0;

    for (int k = 0; k < (nk - 1); k++) {
       y1 = pow((at[i + n * k]       - mean) / variance, powexp) * p[i + n * k];
       y2 = pow((at[i + n * (k + 1)] - mean) / variance, powexp) * p[i + n * (k + 1)];
       res += width.values[k] * (y1 + y2) * 0.5;
    }
    res = res / cdfint; // Normalize
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
        // Calculate requested moment
        rvalptr[i] = c_d2moments_calculate_moment(i, dimptr, atptr, pptr, what);
    }

    // Releasing allocated memory
    free(weights.values);

    // Releasing protected object(s) and return
    UNPROTECT(nprotected);
    return rval;
}

