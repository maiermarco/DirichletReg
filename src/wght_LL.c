/*============================================================================*\
 * WARNING:
 * These routines are intended for internal usage and optimized for speed.
 * Therefore, checks, coercions, etc. are only minimally done and R will most
 * likely crash if you do not prepare input properly.
\*============================================================================*/

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>



SEXP wght_LL(
  SEXP logY,            // [matrix] log of the response Y
  SEXP alpha,           // [matrix] alpha
  SEXP Aplus,           // [double] alpha+ ... row-sums of alpha
  SEXP rc,              //    [int] number of rows and columns
  SEXP wghts            // [double] probability/frequency weights
  ){

  
  
  // input and pointers
  SEXP real_logY   = PROTECT(coerceVector(logY, REALSXP));
  double *p_logY   = REAL(real_logY);

  SEXP real_alpha  = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha  = REAL(real_alpha);

  SEXP real_Aplus  = PROTECT(coerceVector(Aplus, REALSXP));
  double *p_Aplus  = REAL(real_Aplus);
  
  SEXP int_rc      = PROTECT(coerceVector(rc, INTSXP));
  int v_r          = INTEGER(int_rc)[0];
  int v_c          = INTEGER(int_rc)[1];

  SEXP real_wghts  = PROTECT(coerceVector(wghts, REALSXP));
  double *p_wghts  = REAL(real_wghts);
  

  
// output vector and pointer - LL ----------------------------------------------
  double result = 0.0;
  
//  SEXP result = PROTECT(allocVector(REALSXP, v_r));
//  double *p_result = REAL(result);

  
  
// weighted log-likelihood -----------------------------------------------------
  for(int row = 0; row < v_r; ++row){
    for(int col = 0; col < v_c; ++col){
      result +=
        ( p_alpha[row + col * v_r] - 1.0 ) * p_logY[row + col * v_r]
        - lgammafn( p_alpha[row + col * v_r] );
    }
    result = p_wghts[row] * (result + lgammafn(p_Aplus[row]));
  }



// finish ----------------------------------------------------------------------
  
  UNPROTECT(5);

  return ScalarReal(result);

}
