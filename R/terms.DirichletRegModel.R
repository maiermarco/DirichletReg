terms.DirichletRegModel <- function(x, ...){
  y <- terms(x$formula)
  if(is.null(y)){
    y <- attr(x, "terms")
    if(is.null(y)) stop("no terms component nor attribute")
  }
  return(y)
}

extractAIC.DirichletRegModel <- function(fit, scale, k = 2, ...) 
{
#  if(!missing(scale)) warning("scale argument not used.")
  return(AIC(object = fit, k = k))
}
