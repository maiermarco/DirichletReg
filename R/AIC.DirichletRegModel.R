AIC.DirichletRegModel <- function(object, ..., k = 2){
  -object$logLik/2 + k*object$npar
}

BIC.DirichletRegModel <- function(object, ...){
  -object$logLik/2 + log(nrow(object$X[[1]]))*object$npar
}
