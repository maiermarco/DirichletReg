fitted.DirichletRegModel <- function(object, mu = TRUE, alpha = FALSE, phi = FALSE, ...){
  
  if((mu + alpha + phi) == 0L) stop("Either \"mu\", \"alpha\", or \"phi\" has to be requested.")
  
  if((mu + alpha + phi) == 1L){
    if(mu)    return(object$fitted.values$mu)
    if(alpha) return(object$fitted.values$alpha)
    if(phi)   return(object$fitted.values$phi)
  } else {
    res <- list()
    if(mu)    res[["mu"]]    <- object$fitted.values$mu
    if(alpha) res[["alpha"]] <- object$fitted.values$alpha
    if(phi)   res[["phi"]]   <- object$fitted.values$phi
    return(res)
  }
}
