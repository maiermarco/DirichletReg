get_starting_values <- function(Y, X.mats, Z.mat, repar, base){

  if(!repar){

    beta.LL <- function(x, y, X){
      b <- matrix(x, ncol=2)
      if(ncol(X) > 1){
        dbeta(y, exp(X%*%b[,1]), exp(X%*%b[,2]), log=TRUE)
      } else {
        dbeta(y, unlist(exp(X*x[1])), unlist(exp(X*x[2])), log=TRUE)
      }
    }
    
    beta.LL.deriv <- function(x, y, X){
      b <- matrix(x, ncol=2)
      grad <- matrix(0, nrow=nrow(X), ncol=prod(dim(b)))
      element <- 1
      if(ncol(X) > 1){
        for(cc in 1:ncol(b))for(rr in 1:nrow(b)){
          grad[,element] <- X[,rr]*(psigamma(exp(X%*%b[,1]+X%*%b[,2]))-psigamma(exp(X%*%b[,cc]))+log(y))
          element <- element + 1
        }
      } else {
        for(cc in 1:ncol(b)){
          grad[,element] <- X*(psigamma(exp(X%*%x[1]+X%*%x[2]))-psigamma(exp(X%*%x[cc]))+log(y))
          element <- element + 1
        }
      }
    }
 
    unidim.fit <- sapply(1:ncol(Y), function(i){
      suppressWarnings(
        maxBFGS(beta.LL, beta.LL.deriv,
          start=rep(0, 2*ncol(X.mats[[i]])),
          tol=1e-05,
          finalHessian=FALSE,
          X=X.mats[[i]],
          y=Y[,i])$estimate[1:ncol(X.mats[[i]])]
      )
    })

  } else {
  
    Y.logit <- log(Y/(1-Y))
    unidim.fit <- sapply((1:ncol(Y))[-base], function(i){ lm(Y[,i]~X.mats[[i]]-1)$coefficients })
    unidim.fit <- c(unidim.fit, lm(I(rep(.5,nrow(Y)))~Z.mat-1)$coefficients)
  
  }

  return(as.numeric(unlist(unidim.fit)))

}
                                                                               