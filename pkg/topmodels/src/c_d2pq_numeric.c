
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>

/* Two-dimensional linear interpolation
 *
 * xlo/xhi: Lower and upper value on the first dimension ('x axis')
 * ylo/yhi: Corresponding values on the second dimension ('y axis')
 * x: Value at which the corresponding y should be evaluated trough linear interpolation.
 *
 * Returns a sinle real value.
 *
 * Reto, July 2025.
 */
double interpolate_2d(double xlo, double xhi, double ylo, double yhi, double x) {
    // Interval width 0 - return average of ylo, yhi.
    if (xlo == xhi) { return (ylo + yhi) / 2.; }
    // Performing interpolation
    return ylo + (x - xlo) * (yhi - ylo) / (xhi - xlo);
}

/* Simple numeric integration using the trapezodial rule.
 *
 * If 'x > xhi' the return value is the area under 'xlo/xhi given ylo/yhi', if
 * 'x' is within the interval [xlo, xhi] we first calculate 'y' at 'x' using
 * two-dimensional linear interpolation, and then calculate the area under
 * 'xlo/x given ylo/y'.
 *
 * @param xlo/xhi: Lower and upper value on the first dimension ('x axis')
 * @param ylo/yhi: Corresponding values on the second dimension ('y axis')
 * @param x: Value at which the corresponding y should be evaluated trough linear interpolation.
 *
 * @return Single real value with the area.
 *
 * Reto, July 2025.
 */
double integrate_2d(double xlo, double xhi, double ylo, double yhi, double x) {
    if (ISNAN(xlo) || ISNAN(yhi) || ISNAN(ylo) || ISNAN(yhi) || ISNAN(x)) { return NA_REAL; }
    // Sanity check
    if (x < xlo || xlo > xhi) { error("Problematic arguments - can't perform integration."); }
    // If 'x < xlo' we must integrate "xlo to x", wherefore we re-calculate 'yhi' by performing
    // linear interpolation first and overwrite xhi by x, so that we can use the same formula.
    if (x < xhi) {
        yhi = interpolate_2d(xlo, xhi, ylo, yhi, x);
        xhi = x;
    }
    // Numerical integration w/ trapezoidal rule.
    return (xhi - xlo) * (ylo + yhi) / 2.0;
}


/* Calculate CDF from PDF
 *
 * @param dim  Integer vector with the dimension of the matrices 'at' and 'pdf'.
 * @param x    Real vector defining where to evaluate the CDF.
 * @param at   Real pointer; matrix of dimension (dim[0] x dim[1]) where each row
 *             corresponds to a distribution; contains the points of a grid for
 *             which the PDF (pdf) is given.
 * @param pdf  Real pointer; matrix of dimension (dim[0] x dim[1]) containing the PDF
 *             for each distribution (rows) at position 'at'.
 * @param res  Real pointer; matrix of dimension (dim[0] x length(x)) to store the result (CDF).
 *
 * @return Void function, updates 'res'.
 *
 * Reto, July 2025.
 */
void c_d2p_calculate(SEXP dim, SEXP x, double* at, double* pdf, double* res) {

    // Element index
    int idx, residx; // Element indices for matrices
    // Dimensions
    int *dimptr = INTEGER(dim), n = dimptr[0], k = dimptr[1];
    // Vector x
    int nx = length(x);
    double *xptr = REAL(x);

    // Integer counter for current position in 'x'
    for (int i = 0; i < n; i++) {
        // Looping over elements in 'x' (at which the CDF is evaluated)
        for (int ix = 0; ix < nx; ix++) {
            // Index pointer for results matrix
            residx = i * n + ix;

            // Catching the simple case where x[ix] <= lowest 'at' (below defined grid)
            // Or where x[idx] is >= highest 'at' value (above defined grid)
            // - Assign 0.0 if below
            // - Assign 1.0 if above (assuming we cover the entire range)
            if (xptr[ix] <= at[i]) {
                res[residx] = 0.0; continue;
            } else if (xptr[nx] >= at[i + (k - 1) * n]) {
                res[residx] = 1.0;  continue;
            }

            // TODO(R): Inefficiently solved as I go trough the entire matrix j=0, ... every time;
            //          should take a shortcut and use the result from the previous 'ix'.
            // Else perform numeric integration
            res[residx] = pdf[i];  // Initialization
            for (int j = 1; j < k; j++) {
                idx = i * n + j; // Current element index for matrices
                // If x is below the current interval we already reached our result; breaking loop
                if (xptr[ix] < at[idx - n]) { break; }
                // Numerical integration
                res[residx] += integrate_2d(at[idx - n], at[idx], pdf[idx - n], pdf[idx], xptr[ix]);
            }
        }
    }
}

/* Calculate Quantile from PDF
 *
 * @param dim: Integer vector with the dimension of the matrices 'at' and 'pdf'.
 * @param x:   Real vector defining where to evaluate the Quantiles (in (0, 1)).
 * @param at:  Real pointer; matrix of dimension (dim[0] x dim[1]) where each row
 *             corresponds to a distribution; contains the points of a grid for
 *             which the PDF (pdf) is given.
 * @param pdf: Real pointer; matrix of dimension (dim[0] x dim[1]) containing the PDF
 *             for each distribution (rows) at position 'at'.
 * @param res: Real pointer; matrix of dimension (dim[0] x length(x)) to store the result (quantile).
 *
 * @return Void function, updates 'res'.
 *
 * Reto, July 2025.
 */
void c_d2q_calculate(SEXP dim, SEXP x, double* at, double* pdf, double* res) {

    // Element index
    int idx, residx; // Element indices for matrices
    // Dimensions
    int *dimptr = INTEGER(dim), n = dimptr[0], k = dimptr[1];
    // Vector x
    int nx = length(x);
    double *xptr = REAL(x);
    double tmp;
    double infinity = 1./0.;

    // Flag used to check if we found the quantile we were looking for; required
    // to be able to store the highest 'grid point' if the quantile falls above
    // the grid we use for evaluation.
    bool found;

    // Integer counter for current position in 'x'
    for (int i = 0; i < n; i++) {
        // Looping over elements in 'x' (at which the CDF is evaluated)
        for (int ix = 0; ix < nx; ix++) {
            // Flag used to check if we have fou
            // Index pointer for results matrix
            residx = i * n + ix;
            found = false; // Default

            // Carching the simplest case: if the current quantile (xptr[ix])
            // falls below the lowest pdf[i] we assign x[i] as our result.
            if (xptr[ix] <= pdf[i]) { res[residx] = at[i]; continue; }

            // TODO(R): Inefficiently solved as I go trough the entire matrix j=0, ... every time;
            //          should take a shortcut and use the result from the previous 'ix'.
            // Else perform numeric integration
            res[residx] = pdf[i];  // Initialization
            for (int j = 1; j < k; j++) {
                idx = i * n + j; // Current element index for matrices

                // Numerical integration of the current interval
                tmp = integrate_2d(at[idx - n], at[idx], pdf[idx - n], pdf[idx], infinity);

                // If xptr[ix] < (res[residx] + tmp) the quantile we are looking for falls into
                // the current interval. Interpolate, store the result and continue the loop.
                if (xptr[ix] < (res[residx] + tmp)) {
                    res[residx] = interpolate_2d(res[residx], res[residx] + tmp, at[idx - n], at[idx], xptr[ix]);
                    found = true; // Result found
                    break;
                }

                // Else adding 'tmp' to the current numeric integral.
                res[residx] += tmp;
            }

            // Not found? Seems we are above the grid; store highest grid value
            if (!found) { res[residx] = at[idx]; }
        }
    }
}

/*
 * Numerically integrating the density to evaluate the cummulative distribution function (CDF)
 * given we know nothing else than the PDF of a distribution.
 * 
 * at:  Numeric vector, grid points along the x-axis of the distribution.
 * p:   The density at grid points 'x'
 * int: Integer of length 2, dimension of matrices at/p
 * x:   Numeric, where to evaluate the CDF.
 * what: Character, either "cdf" or "quantile".
 * 
 * Dimensions
 * ----------
 * Input 'dim' contains the dimension of the matrices, whereof we will
 * use dim[0] = n and dim[1] = k, i.e., assuming matrix dimension (n x k).
 */
SEXP c_d2pq_numeric(SEXP at, SEXP p, SEXP dim, SEXP x, SEXP continuous, SEXP what) {
    // Setting up pointers for at, pat, dim, x
    double *atptr = REAL(at);
    double *pptr  = REAL(p);
    int *dimptr = INTEGER(dim);

    // What do we need to calculate?
    const char* dowhat = CHAR(STRING_ELT(what, 0));
    bool do_cdf      = strcmp(dowhat, "cdf")      == 0;
    //bool do_quantile = strcmp(dowhat, "quantile") == 0;

    int n = dimptr[0], nx = length(x);

    // Setting up results matrix (return value)
    SEXP rval;
    PROTECT(rval = allocVector(REALSXP, n * nx));
    double *rvalptr = REAL(rval);

    // Looping trough all observations i = 1, ..., n;
    if (do_cdf) {
        c_d2p_calculate(dim, x, atptr, pptr, rvalptr);
    } else {
        c_d2q_calculate(dim, x, atptr, pptr, rvalptr);
    }


    // Releasing protected object(s) and return
    UNPROTECT(1);
    return rval;
}

