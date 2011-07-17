AIC.DirichletRegMixture <- function(object, ..., k = 2){
  -object$logLik/2 + k*object$npar
}

BIC.DirichletRegMixture <- function(object, ...){
  -object$logLik/2 + log(object$nobs)*object$npar
}
