/*============================================================================*\
 * WARNING:
 * These routines are intended for internal usage and optimized for speed.
 * Therefore, checks, coercions, etc. are only minimally done and R will most
 * likely crash if you do not prepare input properly.
\*============================================================================*/

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>



SEXP wght_LL_grad_common(
  SEXP logY,            // [matrix] log of the response Y
  SEXP alpha,           // [matrix] alpha
  SEXP Aplus,           // [double] row-sums of alpha
  SEXP digamma_alpha,   // [matrix] digamma of alpha
  SEXP digamma_Aplus,   // [double] digamma of the row sums of alpha
  SEXP X,               //   [list] list of X matrices
  SEXP X_dims,          //    [int] ncols of X
  SEXP rc,              //    [int] number of rows and columns
  SEXP npar,            //    [int] number of parameters
  SEXP wghts            // [double] probability/frequency weights
  ){



  // input and pointers
  SEXP real_logY  = PROTECT(coerceVector(logY, REALSXP));
  double *p_logY  = REAL(real_logY);

  SEXP real_alpha = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha = REAL(real_alpha);

  SEXP real_Aplus = PROTECT(coerceVector(Aplus, REALSXP));
  double *p_Aplus = REAL(real_Aplus);

  SEXP real_di_al = PROTECT(coerceVector(digamma_alpha, REALSXP));
  double *p_di_al = REAL(real_di_al);

  SEXP real_di_Ap = PROTECT(coerceVector(digamma_Aplus, REALSXP));
  double *p_di_Ap = REAL(real_di_Ap);

  SEXP X_list     = PROTECT(coerceVector(X, VECSXP));

  SEXP int_X_dims = PROTECT(coerceVector(X_dims, INTSXP));
  int *p_X_dims   = INTEGER(int_X_dims);

  SEXP int_rc     = PROTECT(coerceVector(rc, INTSXP));
  int v_r         = INTEGER(int_rc)[0];
  int v_c         = INTEGER(int_rc)[1];

  SEXP int_npar   = PROTECT(coerceVector(npar, INTSXP));

  SEXP real_wghts = PROTECT(coerceVector(wghts, REALSXP));
  double *p_wghts = REAL(real_wghts);



// output vector and pointer - LL ----------------------------------------------
  SEXP result = PROTECT(allocVector(REALSXP, v_r));
  double *p_result = REAL(result);



// output vector and pointer - gradient ----------------------------------------
  SEXP gradient = PROTECT(allocMatrix(REALSXP, v_r, INTEGER(int_npar)[0]));
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
    p_result[row] = p_wghts[row] * (p_result[row] + lgammafn(p_Aplus[row]));
  }



// weighted gradient -----------------------------------------------------------
  memset(aux_term, 0.0, v_r * sizeof(double));
  for(int comp = 0; comp < v_c; ++comp){
    for(int row = 0; row < v_r; ++row){
      aux_term[row] =
        p_wghts[row] *
        p_alpha[row + comp*v_r] * (
          p_logY[row + comp*v_r] +
          p_di_Ap[row] -
          p_di_al[row + comp*v_r]
        );
    }
    for(int iv = 0; iv < p_X_dims[comp]; ++iv){
      for(int row = 0; row < v_r; ++row){
        p_gradient[row + deriv_par*v_r] =
          REAL(coerceVector(VECTOR_ELT(X_list, comp), REALSXP))[row + iv*v_r] * // this is X_m for component "d"
          aux_term[row];
      }
    ++deriv_par;
    }
  }



// finish ----------------------------------------------------------------------

  setAttrib(result, install("gradient"), gradient);

  UNPROTECT(12);

  return result;

}
