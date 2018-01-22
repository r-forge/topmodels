#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>



void wobsmu(double *x, double *z, double *y, double *par, double *left, double *right, int p, int q, int n, double *wobs, double *hess) 
{
  int i, j;



  double ddist, a, pdist, mills, mu, sigma;

  for(i = 0; i < n; i++) {
    mu = 0;
    sigma = 0;
    for(j = 0; j < p; j++) {
      mu = mu + x[n*j+i]*par[j];
    }
    for(j = 0; j < q; j++) {
      sigma = sigma + z[n*j+i]*par[p+j];
    }

    sigma = exp(sigma);
    if(y[i] <= left[i]) {

      ddist = dnorm((left[i] - mu) / sigma, 0.0, 1.0, 0) / sigma;
      pdist = pnorm5((left[i] - mu) / sigma, 0.0, 1.0, 1, 0);
      mills = ddist / pdist;
      a = (left[i] - mu)/ pow(sigma, 2.0) + mills;
      hess[i] =  a * mills;
      wobs[i] = mu - 1/a;


    } else {
      if(y[i] >= right[i]){
        ddist = dnorm((right[i] - mu) / sigma, 0.0, 1.0, 0);
        pdist = pnorm5((right[i] - mu) / sigma, 0.0, 1.0, 0, 0);
        mills = ddist / pdist;
        a = (-right[i] + mu)/ pow(sigma, 2.0) + mills;
        hess[i] = a*mills;
        wobs[i] = mu + 1/a;
      } else {
        hess[i] = 1 / pow(sigma, 2.0);
        wobs[i] = y[i];
      }
    }

    
  }
}




void wobssigma(double *x, double *z, double *y, double *par, double *left, double *right, int p, int q, int n, double *wobs, double *hess) 
{
  int i, j;



  double ddist, pdist, mills, mu, sigma, score, dcm, dcm2, sd2;

  for(i = 0; i < n; i++) {
    mu = 0;
    sigma = 0;
    for(j = 0; j < p; j++) {
      mu = mu + x[n*j+i]*par[j];
    }
    for(j = 0; j < q; j++) {
      sigma = sigma + z[n*j+i]*par[p+j];
    }
    sigma = exp(sigma);
    if(y[i] <= left[i]) {
      ddist = dnorm((left[i] - mu) / sigma, 0.0, 1.0, 0) / sigma;
      pdist = pnorm5((left[i] - mu) / sigma, 0.0, 1.0, 1, 0);
      mills = ddist / pdist;
      sd2 = pow(sigma, 2.0);
      dcm = left[i] - mu;
      dcm2 = pow(dcm, 2.0);
      score = -1 * mills * (left[i] - mu);
      hess[i] =  (-2 * dcm + dcm2*dcm/sd2)*mills + pow(mills, 2.0)*dcm2 - score*sigma;
      wobs[i] = log(sigma) + 1/hess[i]*score;
    } else {
      if(y[i] >= right[i]){
        ddist = dnorm((right[i] - mu) / sigma, 0.0, 1.0, 0) / sigma;
        pdist = pnorm5((right[i] - mu) / sigma, 0.0, 1.0, 0, 0);
        mills = ddist / pdist;
        sd2 = pow(sigma, 2.0);
        dcm = left[i] - mu;
        dcm2 = pow(dcm, 2.0);
        score = mills * (right[i] - mu);
        hess[i] = ( 2 * dcm - dcm2*dcm/sd2)*mills + pow(mills, 2.0)*dcm2- score*sigma;
        wobs[i] = log(sigma) + 1/hess[i]*score;
      } else {
        sd2 = pow(sigma, 2.0);
        dcm2 = pow((y[i] - mu), 2.0);
        score = (dcm2 - sd2) / pow(sigma, 2.0);
        hess[i] = -(sd2 - 3 * dcm2)/sd2- score*sigma;
        wobs[i] = log(sigma) + 1/hess[i]*score;
      }
    }
    
    
  }

}


double dcnorm2(double *x, double *z, double *y, double *par, double *left, double *right, int p, int q, int n)
{
  int i, j;

  double rval, mu, sigma, dens;

  rval = 0;
  for(i = 0; i < n; i++) {
    mu = 0;
    sigma = 0;
    for(j = 0; j < p; j++) {
      mu = mu + x[n*j+i]*par[j];
    }
    for(j = 0; j < q; j++) {
      sigma = sigma + z[n*j+i]*par[p+j];
    }
    sigma = exp(sigma);
    if(y[i] <= left[i]) {
      dens = pnorm5((left[i] - mu) / sigma, 0.0, 1.0, 1, 1);
    } else {
      if(y[i] >= right[i]){
        dens = pnorm5((right[i] - mu) / sigma, 0.0, 1.0, 0, 1);
      } else {
        dens = dnorm((y[i] - mu) / sigma, 0.0, 1.0, 1) - log(sigma);
      }
    }
    rval = rval + dens;
  }

  return rval;
}


SEXP crchglmnet(SEXP x, SEXP z, SEXP y, SEXP left, SEXP right, SEXP lambdaseq, SEXP maxit, SEXP reltol)
{
  int n = INTEGER(GET_DIM(x))[0];
  int p = INTEGER(GET_DIM(x))[1];
  int q = INTEGER(GET_DIM(z))[1];
  int nlambda = length(lambdaseq);

  SEXP coefpath = PROTECT(allocMatrix(REALSXP, nlambda, p+q));
  SEXP llpath = PROTECT(allocVector(REALSXP, nlambda));
  SEXP rval = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(rval, 0, coefpath);
  SET_VECTOR_ELT(rval, 1, llpath);

  double *coefpathptr = REAL(coefpath);
  double *llpathptr = REAL(llpath);

  double *yptr = REAL(y);
  double *xptr = REAL(x);
  double *zptr = REAL(z);
  double *leftptr = REAL(left);
  double *rightptr = REAL(right);
  double *lambdaseqptr = REAL(lambdaseq);
  int *maxitptr = INTEGER(maxit);
  double *reltolptr = REAL(reltol);

  int i, lambdaind, it1, it2, colind, rowind, colind2;
  int activex[p], activez[q];
  double a, xbeta, coefold, zgamma, dev, devold, ll, llold, lambda, coef[p+q];
  double wxsum[p], wzsum[q], resid[n], w[n], wobs[n];


  for(i = 0; i < p+q; i++) {
    coef[i] = 0;
  }

  for(lambdaind = 0; lambdaind < nlambda; lambdaind++) {
    lambda = lambdaseqptr[lambdaind];
    for(it1 = 0; it1 < *maxitptr; it1++) {
      /* location */
      wobsmu(xptr, zptr, yptr, coef, leftptr, rightptr, p, q, n, wobs, w);
      for(i = 0; i < p; i++) {
        activex[i] = 1;
      }
      for(colind = 0; colind < p; colind++) {
        wxsum[colind] = 0;
        for(rowind = 0; rowind < n; rowind++) {
          xbeta = 0;
          for(colind2 = 0; colind2 < p; colind2++) { 
            xbeta = xbeta + xptr[rowind+colind2*n]*coef[colind2];
          }
          resid[rowind] = w[rowind]*(wobs[rowind] - xbeta);
          wxsum[colind] = wxsum[colind] + w[rowind]*pow(xptr[rowind+colind*n], 2);
        }
      }
      for(it2 = 0; it2 < *maxitptr; it2++) {
        for(colind = 0; colind < p; colind++) {
          if(activex[colind] == 1) {
            a = 0;
            for(rowind = 0; rowind < n; rowind++) {
              a = a + xptr[rowind + colind*n] * resid[rowind];
            }
            a=(a + wxsum[colind]*coef[colind])/wxsum[colind];
            coefold = coef[colind];
            if(lambda < fabs(a)){
              if(a>0) { 
                coef[colind] = (a - lambda);
              } else {
                coef[colind] = (a + lambda);
              }

            } else {

              coef[colind] = 0; 
              activex[colind] = 0;
            }
            if(coef[colind] != coefold) {
              for(rowind = 0; rowind < n; rowind++) {
                resid[rowind] = resid[rowind] + w[rowind] * xptr[rowind + colind*n]*(coefold - coef[colind]);
              }
            }
          }
        } 
        dev = 0;
        for(rowind = 0; rowind < n; rowind++) {
          dev = dev + pow(resid[rowind], 2.0);
        }

        if(fabs((dev-devold)/devold) < *reltolptr) {
          for(colind2 = 0; colind2 < p; colind2++) { 
            coefpathptr[colind2*nlambda + lambdaind] = coef[colind2];
          }

          break;
        }
        devold = dev;
      }
 
      /* log-lik */
      ll = dcnorm2(xptr, zptr, yptr, coef, leftptr, rightptr, p, q, n);
      if(fabs((ll-llold)/llold) < *reltolptr) {
        break;
      }
      llold = ll;
    }

    /* scale */
    for(it1 = 0; it1 < *maxitptr; it1++) {

      wobssigma(xptr, zptr, yptr, coef, leftptr, rightptr, p, q, n, wobs, w);
      activez[0] = 2;   /* exclude intercept from regularization TODO: what to do if no intercept */
      for(i = 1; i < q; i++) {
        activez[i] = 1;
      }
      
      for(colind = 0; colind < q; colind++) {
        wzsum[colind] = 0;
        for(rowind = 0; rowind < n; rowind++) {
          xbeta = 0;
          for(colind2 = 0; colind2 < q; colind2++) { 
            xbeta = xbeta + zptr[rowind+colind2*n]*coef[colind2+p];
          }
          resid[rowind] = w[rowind]*(wobs[rowind] - xbeta);
          wzsum[colind] = wzsum[colind] + w[rowind]*pow(zptr[rowind+colind*n], 2); 
        }
      }
      for(it2 = 0; it2 < *maxitptr; it2++) {
        for(colind = 0; colind < q; colind++) {
          if(activez[colind] != 0) {
            a = 0;
            for(rowind = 0; rowind < n; rowind++) {
              a = a + zptr[rowind + colind*n] * resid[rowind];
            }
            a=(a + wzsum[colind]*coef[colind+p])/wzsum[colind];
            coefold = coef[colind+p];
            if(activez[colind] == 2) {
              coef[colind + p] = a;
            } else {
              if(lambda < fabs(a)){
                if(a>0) { 
                  coef[colind + p] = (a - lambda);
                } else {
                  coef[colind + p] = (a + lambda);
                }
              } else {
                coef[colind + p] = 0;
                activez[colind] = 0; 
              }
            }
            if(coef[colind] != coefold) {
              for(rowind = 0; rowind < n; rowind++) {
                resid[rowind] = resid[rowind] + w[rowind] * zptr[rowind + colind*n]*(coefold - coef[colind+p]);
              }
            }
          }
        } 
        dev = 0;
        for(rowind = 0; rowind < n; rowind++) {
          dev = dev + pow(resid[rowind], 2.0);
        }

        if(fabs((dev-devold)/devold) < *reltolptr) {
          for(colind2 = 0; colind2 < q; colind2++) { 
            coefpathptr[(colind2 + p)*nlambda + lambdaind] = coef[colind2 + p];
          }

          break;
        }
        devold = dev;
      }

      ll = dcnorm2(xptr, zptr, yptr, coef, leftptr, rightptr, p, q, n);
      if(fabs((ll-llold)/llold) < *reltolptr) {
        break;
      }
      llold = ll;
    }
    llpathptr[lambdaind] = ll;
  }

 UNPROTECT(3);
 return rval;
}

