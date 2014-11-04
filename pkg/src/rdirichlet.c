/*============================================================================*\
 *
 * ddirichlet_log_vector ... RNG with a vector of alphas
 * ddirichlet_log_matrix ... RNG with a matrix of alphas
 * 
 * y ........ dirichlet-distributed matrix
 * alpha .... vector or matrix of alpha-values
 * rc ....... vector of number of rows and columns
 *
\*============================================================================*/



#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>





SEXP rdirichlet_vector(SEXP n, SEXP alpha){
  // helpers
  SEXP v_n = PROTECT(coerceVector(n, INTSXP));
  int p_n = INTEGER(v_n)[0];
  
  int dims = length(alpha);
  double norming_constant;
  
  // pointers
  SEXP real_alpha = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha = REAL(real_alpha);
  
  // check alphas ...
  for(int i = 0; i < dims; i++){
    if(p_alpha[i] <= 0.0) error("alphas must be > 0");
  }
  
  // output vector and pointer
  SEXP result = PROTECT(allocMatrix(REALSXP, p_n, dims));
  double *p_result = REAL(result);

  GetRNGstate();

  for(int row = 0; row < p_n; row++){
    R_CheckUserInterrupt();
    norming_constant = 0.0;
    p_result[row] = 0.0;
    for(int col = 0; col < dims; col++){
      norming_constant += p_result[row + col * p_n] = rgamma(p_alpha[col], 1.0);
    }
    for(int col = 0; col < dims; col++){
      p_result[row + col * p_n] = p_result[row + col * p_n] / norming_constant;
    }
  }

  PutRNGstate();

  // unprotect output vector and return
  UNPROTECT(3);
  
  return result;
  
}





SEXP rdirichlet_matrix(SEXP n, SEXP alpha, SEXP alpha_rc){
  // helpers
  SEXP v_n = PROTECT(coerceVector(n, INTSXP));
  int p_n = INTEGER(v_n)[0];
  SEXP int_alpha_rc = PROTECT(coerceVector(alpha_rc, INTSXP));
  int dims = INTEGER(int_alpha_rc)[1];
  
  // pointers
  SEXP real_alpha = PROTECT(coerceVector(alpha, REALSXP));
  double *p_alpha = REAL(real_alpha);
  
  // check alphas ...
  if(p_n != INTEGER(alpha_rc)[0]) error("n and alpha do not match");
  for(int i = 0; i < length(real_alpha); i++){
    if(p_alpha[i] <= 0.0) error("alphas must be > 0");
  }
  
  // output vector and pointer
  SEXP result = PROTECT(allocMatrix(REALSXP, p_n, dims));
  double *p_result = REAL(result);

  double norming_constant[p_n];
  for(int i = 0; i < p_n; i++) norming_constant[i] = 0;
  
  GetRNGstate();

  R_CheckUserInterrupt();
  for(int col = 0; col < dims; col++){
    for(int row = 0; row < p_n; row++){
      norming_constant[row] += p_result[row + col * p_n] = rgamma(p_alpha[row + col * p_n], 1.0);
    }
  }

  PutRNGstate();

  for(int col = 0; col < dims; col++){
    for(int row = 0; row < p_n; row++){
      p_result[row + col * p_n] = p_result[row + col * p_n] / norming_constant[row];
    }
  }
  
  // unprotect output vector and return
  UNPROTECT(4);
  return result;
  
}
