predict.DirichletRegModel <- function(object, newdata, alpha=FALSE, ...){

  if(missing(newdata)) return(fitted(object, alpha))

  
  
  res <- matrix(NA, nrow=400, ncol=object$Y$dims)
  
  coef.ind <- cumsum(object$n.vars)


  for(i in 1:ncol(res)){
    res[,i] <- newdata[[i]] %*% object$coefficients[ifelse(i==1,1,coef.ind[i-1]+1):coef.ind[i]]
  }
  
  res <- exp(res)
  
  if(!alpha) res <- res / rowSums(res)
  
  return(res)

}
