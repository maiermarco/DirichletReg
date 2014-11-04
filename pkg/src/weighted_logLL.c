/*============================================================================*\
 * WARNING:
 * These routines are intended for internal usage and optimized for speed.
 * Therefore, checks, coercions, etc. are only minimally done and R will most
 * likely crash if you do not prepare input properly.
\*============================================================================*/

/*============================================================================*\
 * weighted_logLL ... computes individual weighted log-likelihoods
 * 
 * y ....... dirichlet-distributed matrix
 * alpha ... vector or matrix of alpha-values
 * rc ...... vector of number of rows and columns
 * wghts ... vector of non-negative weights
\*============================================================================*/



#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>





SEXP weighted_logLL(SEXP y, SEXP alpha, SEXP Aplus, SEXP rc, SEXP wghts){
  // input and pointers
  SEXP real_y       = PROTECT(coerceVector(y, REALSXP));
  double *p_y       = REAL(real_y);

  SEXP real_alpha   = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha   = REAL(real_alpha);

  SEXP real_Aplus   = PROTECT(coerceVector(Aplus, REALSXP));
  double *p_Aplus   = REAL(real_Aplus);

  SEXP int_rc = PROTECT(coerceVector(rc, INTSXP));
  int v_r = INTEGER(int_rc)[0];
  int v_c = INTEGER(int_rc)[1];

  SEXP real_wghts = PROTECT(coerceVector(wghts, REALSXP));
  double *p_wghts = REAL(real_wghts);
  
  // output vector and pointer
  SEXP result = PROTECT(allocVector(REALSXP, v_r));
  double *p_result = REAL(result);

  for(int row = 0; row < v_r; row++){
    p_result[row] = 0.0;
    for(int col = 0; col < v_c; col++){
      p_result[row] += ( p_alpha[row + col * v_r] - 1.0 ) * log( p_y[row + col * v_r] ) - lgammafn( p_alpha[row + col * v_r] );
    }
    p_result[row] = p_wghts[row] * (p_result[row] + lgammafn(p_Aplus[row]));
  }
  
  UNPROTECT(6);
  return result;

}
