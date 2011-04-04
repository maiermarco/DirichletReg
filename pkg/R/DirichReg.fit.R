DirichReg.fit <- function(Y, X, Z, sv, d, k, w, ctls, repar, base){

  if(repar){ 
  res <- suppressWarnings(
           maxBHHH(fn=DReg.repar, grad=DRegGrad.repar,
                  start=sv,
                  finalHessian=TRUE, iterlim=ctls$iterlim,
                  tol=1e-15, reltol=1e-15, print.level=ctls$trace,
                  Y=Y, X=X[[1]], Z=Z, d=d, k=k[1], w=w, base=base)
         )






  } else { 
  res <- suppressWarnings(
           maxBHHH(fn=DReg, grad=DRegGrad, start=sv,
                  finalHessian=TRUE, iterlim=ctls$iterlim,
                  Y=Y, X=X, d=d, k=k, w=w, print.level=ctls$trace)
         )
  }

  return(res)

}





DReg <- function(x, Y, X, d, k, w){
  B <- list(x[1:k[1]])
  for(i in 2:d)B[[i]] <- x[(cumsum(k)[i-1]+1):cumsum(k)[i]]

  A <- sapply(1:d, function(i){ exp(X[[i]]%*%B[[i]]) })



  res <- lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A-1)*log(Y))
  res <- w * res

  return(res)
}





DRegGrad <- function(x, Y, X, d, k, w){
  B <- list(x[1:k[1]])
  for(i in 2:d)B[[i]] <- x[(cumsum(k)[i-1]+1):cumsum(k)[i]]

  A <- sapply(1:d, function(i){ exp(X[[i]]%*%B[[i]]) })




  res <- matrix(NA, ncol=length(x), nrow=nrow(Y))
  ind <- 1

  for(dimm in 1:d){
    for(varr in 1:k[dimm]){
      res[,ind] <- 
                     X[[dimm]][,varr]*A[,dimm]*(  psigamma(rowSums(A))
                                                 -psigamma(A[,dimm])
                                                 +log(Y[,dimm]) )
                     

      ind <- ind + 1
    }
  }
  
  return(res)
}



























DReg.repar <- function(x, Y, X, Z, d, k, w, base){
  B <- matrix(0, nrow=k, ncol=d)
  B[cbind(rep(1:k, (d-1)), rep((1:d)[-base], each=k))] <- x[1:((d-1)*k)]
  
  g <- matrix(x[((d-1)*k+1):length(x)], ncol=1)

  XB <- exp(apply(B, 2, function(b){ X%*%b }))
  MU <- apply(XB, 2, function(x){ x /rowSums(XB) })

  f <- exp(Z%*%g)
  
  A <- apply(MU, 2, function(x){ x * f })
  
  lgamma(rowSums(A))-rowSums(lgamma(A))+rowSums((A-1)*log(Y))
}



DRegGrad.repar <- function(x, Y, X, Z, d, k, w, base){
  B <- matrix(0, nrow=k, ncol=d)
  B[cbind(rep(1:k, (d-1)), rep((1:d)[-base], each=k))] <- x[1:((d-1)*k)]
  
  g <- matrix(x[((d-1)*k+1):length(x)], ncol=1)

  XB <- exp(apply(B, 2, function(b){ X%*%b }))
  MU <- apply(XB, 2, function(x){ x /rowSums(XB) })

  f <- exp(Z%*%g)
  
  denomB <- rowSums(XB)^2
  denomF <- rowSums(XB)

  res <- matrix(NA,nrow=nrow(Y),ncol=length(x))

  i <- 0

  for(co in (1:d)[-base]){
    for(iv in 1:k){
      i <- i + 1
      res[,i] <- -X[,iv]*w*XB[,co]*f*rowSums(
                   sapply((1:d)[-co], function(p){
                     XB[,p]*(log(Y[,p]) - log(Y[,co]) + psigamma(MU[,co]*f) - psigamma(MU[,p]*f))
                   })
                 )/denomB
    }
  }
  
  for(iv in 1:ncol(Z)){
    i <- i + 1
    res[,i] <- Z[,iv]*w*f*rowSums(
                 sapply(1:d, function(p){
                   XB[,p]*(log(Y[,p]) + psigamma(f) - psigamma(MU[,p]*f))
                 })
               )/denomF
  }
  
  return(res)
  
}
