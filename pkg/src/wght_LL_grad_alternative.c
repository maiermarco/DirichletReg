/*============================================================================*\
 * WARNING:
 * These routines are intended for internal usage and optimized for speed.
 * Therefore, checks, coercions, etc. are only minimally done and R will most
 * likely crash if you do not prepare input properly.
\*============================================================================*/

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>



SEXP wght_LL_grad_alternative(
  SEXP logY,            // [matrix] log of the response Y
  SEXP alpha,           // [matrix] alpha
  SEXP mu,              // [matrix] mu
  SEXP phi,             // [double] phi = row-sums of alpha
  SEXP digamma_alpha,   // [matrix] digamma of alpha
  SEXP digamma_phi,     // [double] digamma of phi (row sums of alpha)
  SEXP X,               // [matrix] X
  SEXP X_cols,          //    [int] number of columns of X
  SEXP Z,               // [matrix] Z
  SEXP Z_cols,          //    [int] number of columns of Z
  SEXP r,               //    [int] number of rows
  SEXP c,               //    [int] number of columns
  SEXP base,            //    [int] number of the base category
  SEXP npar,            //    [int] number of parameters
  SEXP wghts            // [double] probability/frequency weights
  ){



  // input and pointers
  SEXP real_logY   = PROTECT(coerceVector(logY, REALSXP));
  double *p_logY   = REAL(real_logY);

  SEXP real_alpha  = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha  = REAL(real_alpha);

  SEXP real_mu     = PROTECT(coerceVector(mu, REALSXP));
  double *p_mu     = REAL(real_mu);

  SEXP real_phi    = PROTECT(coerceVector(phi, REALSXP));
  double *p_phi    = REAL(real_phi);

  SEXP real_di_al  = PROTECT(coerceVector(digamma_alpha, REALSXP));
  double *p_di_al  = REAL(real_di_al);

  SEXP real_di_phi = PROTECT(coerceVector(digamma_phi, REALSXP));
  double *p_di_phi = REAL(real_di_phi);

  SEXP real_X      = PROTECT(coerceVector(X, REALSXP));
  double *p_X      = REAL(real_X);

  int int_X_cols   = asInteger(X_cols);

  SEXP real_Z      = PROTECT(coerceVector(Z, REALSXP));
  double *p_Z      = REAL(real_Z);

  int int_Z_cols   = asInteger(Z_cols);
  int v_r          = asInteger(r);
  int v_c          = asInteger(c);
  int base_m1      = asInteger(base) - 1;
  int int_npar     = asInteger(npar);

  SEXP real_wghts  = PROTECT(coerceVector(wghts, REALSXP));
  double *p_wghts  = REAL(real_wghts);



// output vector and pointer - LL ----------------------------------------------
  SEXP result = PROTECT(allocVector(REALSXP, v_r));
  double *p_result = REAL(result);



// output vector and pointer - gradient ----------------------------------------
  SEXP gradient      = PROTECT(allocMatrix(REALSXP, v_r, int_npar));
  double *p_gradient = REAL(gradient);

  double aux_term[v_r];

  int deriv_par = 0;



// weighted log-likelihood -----------------------------------------------------
  for(int row = 0; row < v_r; ++row){
    p_result[row] = 0.0;
    for(int col = 0; col < v_c; ++col){
      p_result[row] +=
        ( p_alpha[row + col * v_r] - 1.0 ) * p_logY[row + col * v_r]
        - lgammafn( p_alpha[row + col * v_r] );
    }
    p_result[row] = p_wghts[row] * (p_result[row] + lgammafn(p_phi[row]));
  }



// weighted gradient -----------------------------------------------------------
  // beta parameters -----------------------------------------------------------
  for(int comp = 0; comp < v_c; ++comp){ // big loop for mu-variables
    if(comp == base_m1) continue; // skip if base category

    memset(aux_term, 0.0, v_r * sizeof(double));

    for(int row = 0; row < v_r; ++row){
      for(int co = 0; co < v_c; ++co){
        if(co == comp){ continue; } else { // skip the current component
          aux_term[row] +=
            p_mu[row + co*v_r] * (
              p_di_al[row + co*v_r] -
              p_logY[row + co*v_r]
            );
        }
      }
      aux_term[row] =
        p_wghts[row] *
        p_alpha[row + comp*v_r] * (
          (1.0 - p_mu[row + comp*v_r]) *
          (p_logY[row + comp*v_r] - p_di_al[row + comp*v_r]) +
          aux_term[row]
        );
    }

    for(int iv = 0; iv < int_X_cols; ++iv){
      for(int row = 0; row < v_r; ++row){
        p_gradient[row + deriv_par*v_r] =
          p_X[row + iv*v_r] *
          aux_term[row];
      }
    ++deriv_par;
    }
  }

  // gamma parameters ----------------------------------------------------------
  memset(aux_term, 0.0, v_r * sizeof(double));

  for(int row = 0; row < v_r; ++row){
    for(int comp = 0; comp < v_c; ++comp){
      aux_term[row] +=
        p_alpha[row + comp*v_r] * (
          p_logY[row + comp*v_r] -
          p_di_al[row + comp*v_r]
        );
    }
    aux_term[row] = p_wghts[row] * (
      p_phi[row] *
      p_di_phi[row] +
      aux_term[row]
    );
  }

  for(int iv = 0; iv < int_Z_cols; ++iv){
    for(int row = 0; row < v_r; ++row){
      p_gradient[row + deriv_par*v_r] =
        p_Z[row + iv*v_r] *
        aux_term[row];
    }
    ++deriv_par;
  }



// finish ----------------------------------------------------------------------

  setAttrib(result, install("gradient"), gradient);

  UNPROTECT(11);

  return result;

}
