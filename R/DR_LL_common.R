DReg <- function(x, Y, X, d, k, w, NR){
  npar <- length(x)

  B <- sapply(1:d, function(i){ x[(cumsum(c(0,k))[i]+1) : cumsum(k)[i]] }, simplify=FALSE)

  A <- sapply(1:d, function(i){ exp(X[[i]] %*% B[[i]]) })



################################################################################
### LOG-LIKELIHOOD #############################################################
################################################################################

  LL <- w * (lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A-1)*log(Y)))



################################################################################
### GRADIENT ###################################################################
################################################################################

  gradient <- matrix(NA, nrow=nrow(Y), ncol=npar)

  i <- 0
  for(comp in 1:d){
    for(iv in 1:k[comp]){
      i <- i + 1
      gradient[,i] <- X[[comp]][,iv] * A[,comp] * ( log(Y[,comp]) - psigamma(A[,comp]) + psigamma(rowSums(A)) )
    }
  }

  attr(LL, "gradient") <- gradient



################################################################################
### HESSIAN ####################################################################
################################################################################

  if(NR){

  hessian <- matrix(NA, nrow=npar, ncol=npar)

  hessian.ind <- cbind(rep(1:d, k), unlist(sapply(k, function(i) 1:i, simplify=FALSE)))

  for(hess.j in 1:npar){
    for(hess.i in 1:npar){
      if(hess.i < hess.j) next   # skip symmetric elements

      derv <- hessian.ind[c(hess.i, hess.j),1]
      
      vars <- hessian.ind[c(hess.i, hess.j),2]

      ##########################################################################
      ####################################################### SAME RESPONSES ###
      if(derv[1] == derv[2]) {
        derv <- derv[1]

        hessian[hess.i, hess.j] <-
        sum(X[[derv]][,vars[1]] * X[[derv]][,vars[2]] * A[,derv] * (
          log(Y[,derv]) + psigamma(rowSums(A)) - psigamma(A[,derv]) + A[,derv] * (
            psigamma(rowSums(A), 1) - psigamma(A[,derv], 1)
          )
        ))
      ##########################################################################
      ################################################## DIFFERENT RESPONSES ###
      } else {
        hessian[hess.i, hess.j] <-
        sum(
          X[[derv[1]]][,vars[1]]*X[[derv[2]]][,vars[2]]*A[,derv[1]]*A[,derv[2]]*psigamma(rowSums(A), 1)
        )
      }
    }
  }

  hessian <- make.symmetric(hessian)

  attr(LL, "hessian") <- hessian

  }

  return(LL)

}
