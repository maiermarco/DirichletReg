/*============================================================================*\
 * WARNING:
 * These routines are intended for internal usage and optimized for speed.
 * Therefore, checks, coercions, etc. are only minimally done and R will most
 * likely crash if you do not prepare input properly.
\*============================================================================*/

/*============================================================================*\
 * grad_common ... computes individual weighted log-likelihoods for the common
 *                 param.
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

SEXP grad_common(SEXP Y, SEXP X, SEXP alpha, SEXP Aplus, SEXP npar, SEXP X_dims, SEXP A_dims, SEXP wghts){

  // input and pointers
  SEXP real_Y     = PROTECT(coerceVector(Y, REALSXP));
  double *p_Y     = REAL(real_Y);

  SEXP X_list     = PROTECT(coerceVector(X, VECSXP));

  SEXP real_alpha = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha = REAL(real_alpha);

  SEXP real_Aplus = PROTECT(coerceVector(Aplus, REALSXP));
  double *p_Aplus = REAL(real_Aplus);

  SEXP X_dimsList = PROTECT(coerceVector(X_dims, VECSXP));

  SEXP int_npar   = PROTECT(coerceVector(npar, INTSXP));

  SEXP int_A_dims = PROTECT(coerceVector(A_dims, INTSXP));
  int v_r = INTEGER(int_A_dims)[0];
  int v_c = INTEGER(int_A_dims)[1];

  SEXP real_wghts = PROTECT(coerceVector(wghts, REALSXP));
  double *p_wghts = REAL(real_wghts);

  // output matrix and pointer
  SEXP gradient = PROTECT(allocMatrix(REALSXP, v_r, INTEGER(int_npar)[0]));
  double *p_gradient = REAL(gradient);

  int deriv_par = 0;

  for(int comp = 0; comp < v_c; comp++){
    for(int iv = 0; iv < INTEGER(coerceVector(VECTOR_ELT(X_dimsList, comp), INTSXP))[0]; iv++){
      for(int row = 0; row < v_r; row++){
        p_gradient[row + deriv_par*v_r] =
          p_wghts[row] * REAL(coerceVector(VECTOR_ELT(X_list, comp), REALSXP))[row + iv*v_r] * p_alpha[row + comp*v_r]  * (
            log(p_Y[row + comp*v_r]) - digamma(p_alpha[row + comp*v_r]) + digamma(p_Aplus[row])
          );
      }
    deriv_par++;
    }
  }

  UNPROTECT(9);
  return gradient;

}
