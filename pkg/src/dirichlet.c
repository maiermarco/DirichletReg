/*============================================================================*\

  ddirichlet_log_vector ... computes the log-densities for a vector of alphas
  ddirichlet_log_matrix ... computes the log-densities for a matrix of alphas
  
  y ........ dirichlet-distributed matrix
  alpha .... vector or matrix of alpha-values
  r ........ number of rows
  c ........ number of columns
  result ... vector with log-likelihoods

\*============================================================================*/


#include<R.h>
#include<Rmath.h>


void ddirichlet_log_vector(double *y, double *alpha, int *r, int *c, double *result){
  int row, col;
  double a_sum = 0.0, aux = 0.0, norm_const;
  
  for(col = 0; col < *c; col++){
    a_sum += alpha[col];
    aux   += lgammafn(alpha[col]);
  }
  norm_const = lgammafn(a_sum) - aux;

  for(row = 0; row < *r; row++){
    for(col = 0; col < *c; col++){
      result[row] += ( alpha[col] - 1.0 ) * log( y[row + col * *r] );
    }
    result[row] += norm_const;
  }
}


void ddirichlet_log_matrix(double *y, double *alpha, int *r, int *c, double *result){
  int row, col;
  double a_sum;

  for(row = 0; row < *r; row++){
    a_sum = 0.0;
    for(col = 0; col < *c; col++){
      a_sum += alpha[row + col * *r];
      result[row] += ( alpha[row + col * *r] - 1.0 ) * log( y[row + col * *r] ) - lgammafn( alpha[row + col * *r] );
    }
    result[row] += lgammafn(a_sum);
  }
}
