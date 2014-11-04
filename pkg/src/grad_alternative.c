/*============================================================================*\
 * WARNING:
 * These routines are intended for internal usage and optimized for speed.
 * Therefore, checks, coercions, etc. are only minimally done and R will most
 * likely crash if you do not prepare input properly.
\*============================================================================*/

/*============================================================================*\
 * grad_alternative ... computes individual weighted log-likelihoods for the
 *                      alternative param.
 *
 * y ........ dirichlet-distributed matrix
 * X ........ list of design matrices X
 * alpha .... vector or matrix of alpha-values
 * Aplus .... vector of row sums of alpha
 * npar ..... number of parameters
 * X_dims ... list of vectors with the dimensions of elements in X
 * A_dims ... vector with the dimension of y and alpha
 * wghts .... vector of non-negative weights
\*============================================================================*/

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

SEXP grad_alternative(SEXP Y, SEXP X, SEXP Z, SEXP alpha, SEXP eps, SEXP esum, SEXP f, SEXP npar, SEXP X_cols, SEXP Z_cols, SEXP A_dims, SEXP base, SEXP wghts){

  // input and pointers
  SEXP real_Y     = PROTECT(coerceVector(Y, REALSXP));
  double *p_Y     = REAL(real_Y);

  SEXP real_X     = PROTECT(coerceVector(X, REALSXP));
  double *p_X     = REAL(real_X);

  SEXP real_Z     = PROTECT(coerceVector(Z, REALSXP));
  double *p_Z     = REAL(real_Z);

  SEXP real_alpha = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha = REAL(real_alpha);

  SEXP real_eps   = PROTECT(coerceVector(eps, REALSXP));
  double *p_eps   = REAL(real_eps);

  SEXP real_esum  = PROTECT(coerceVector(esum, REALSXP));
  double *p_esum  = REAL(real_esum);

  SEXP real_f     = PROTECT(coerceVector(f, REALSXP));
  double *p_f     = REAL(real_f);

  SEXP int_npar   = PROTECT(coerceVector(npar, INTSXP));

  SEXP int_X_cols = PROTECT(coerceVector(X_cols, INTSXP));

  SEXP int_Z_cols = PROTECT(coerceVector(Z_cols, INTSXP));

  SEXP int_A_dims = PROTECT(coerceVector(A_dims, INTSXP));
  int v_r = INTEGER(int_A_dims)[0];
  int v_c = INTEGER(int_A_dims)[1];

  SEXP int_base   = PROTECT(coerceVector(base, INTSXP));
  int base_m1     = INTEGER(int_base)[0] - 1;

  SEXP real_wghts = PROTECT(coerceVector(wghts, REALSXP));
  double *p_wghts = REAL(real_wghts);

  // output matrix and pointer
  SEXP gradient = PROTECT(allocMatrix(REALSXP, v_r, INTEGER(int_npar)[0]));
  double *p_gradient = REAL(gradient);

  SEXP t_1 = PROTECT(allocVector(REALSXP, v_r));
  double *term1 = REAL(t_1);
  SEXP t_2 = PROTECT(allocVector(REALSXP, v_r));
  double *term2 = REAL(t_2);
  
  SEXP esum2quared = PROTECT(allocVector(REALSXP, v_r));
  double *esum2    = REAL(esum2quared);

  for(int row = 0; row < v_r; row++){
    esum2[row] = R_pow_di(p_esum[row], 2);
  }

  int deriv_par = 0;

//------------------------------------------------------------------------------
  for(int comp = 0; comp < v_c; comp++){ // big loop for mu-variables
    if(comp == base_m1) continue; // skip if base category

    for(int row = 0; row < v_r; row++){ // initialize terms for the helper terms
      term1[row] = term2[row] = 0.0;
    }

    for(int co = 0; co < v_c; co++){
      if(co == comp){
        continue; // skip the current component
      } else {
        for(int row = 0; row < v_r; row++){
          term1[row] += p_eps[row + co*v_r];
          term2[row] += p_eps[row + co*v_r] * (digamma(p_alpha[row + co*v_r]) - log(p_Y[row + co*v_r]));
        }
      }
    }

    for(int row = 0; row < v_r; row++){
      term1[row] = p_wghts[row] * p_eps[row + comp*v_r] * p_f[row] * (
            term1[row] * (log(p_Y[row + comp*v_r]) - digamma(p_alpha[row + comp*v_r])) + term2[row]
          ) / esum2[row];
    }
    
    //--------------------------------------------------------------------------
    for(int iv = 0; iv < INTEGER(int_X_cols)[0]; iv++){
      for(int row = 0; row < v_r; row++){
        p_gradient[row + deriv_par*v_r] = term1[row] * p_X[row + iv*v_r];
      }
    deriv_par++;
    }
  }

//------------------------------------------------------------------------------
  for(int row = 0; row < v_r; row++){
    term1[row] = term2[row] = 0.0;
  }
  
  for(int comp = 0; comp < v_c; comp++){
    for(int row = 0; row < v_r; row++){
      term1[row] += p_eps[row + comp*v_r] * (log(p_Y[row + comp*v_r]) - digamma(p_alpha[row + comp*v_r]));
    }
  }
  for(int row = 0; row < v_r; row++){
    term2[row] += p_wghts[row] * p_f[row] * (p_esum[row] * digamma(p_f[row]) + term1[row]) / p_esum[row];
  }

  for(int iv = 0; iv < INTEGER(int_Z_cols)[0]; iv++){ // big loop for phi-variables
    for(int row = 0; row < v_r; row++){
      p_gradient[row + deriv_par*v_r] = p_Z[row + iv*v_r] * term2[row];        
    }
    deriv_par++;
  }

//------------------------------------------------------------------------------
  UNPROTECT(17);
  return gradient;

}
