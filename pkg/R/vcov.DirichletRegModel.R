vcov.DirichletRegModel <- function(object, ...){
  solve(-object$hessian)
}