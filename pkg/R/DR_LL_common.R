DReg <- function(x, Y, X, d, k, w, NR){
  npar <- length(x)
  seq_along_d <- seq_len(d)

  B <- lapply(seq_along_d, function(i){ x[ (cumsum(c(0L, k))[i] + 1L) : cumsum(k)[i] ] })

  A <- matrix(unlist(lapply(seq_along_d, function(i){ exp(X[[i]] %*% B[[i]]) })), nrow=nrow(Y))
  rsumA <- rowSums(A)


################################################################################
### LOG-LIKELIHOOD
################################################################################

  LL <- w * (lgamma(rsumA) - rowSums(lgamma(A)) + rowSums((A-1.0)*log(Y)))



################################################################################
### GRADIENT
################################################################################

  gradient <- matrix(NA, nrow=nrow(Y), ncol=npar)

  i <- 0L
  for(comp in seq_along_d){
    for(iv in seq_len(k[comp])){
      i <- i + 1L
      gradient[,i] <- w * X[[comp]][,iv] * A[,comp] * ( log(Y[,comp]) - psigamma(A[,comp]) + psigamma(rsumA) )
    }
  }

  attr(LL, "gradient") <- gradient



################################################################################
### HESSIAN
################################################################################

  if(NR){

    hessian <- matrix(NA, nrow=npar, ncol=npar)

    hessian.ind <- cbind(rep(seq_along_d, k), unlist(sapply(k, function(i) seq_len(i), simplify=FALSE)))

    for(hess.j in seq_len(npar)){
      for(hess.i in seq_len(npar)){
        if(hess.i < hess.j) next

        derv <- hessian.ind[c(hess.i, hess.j), 1L]
        
        vars <- hessian.ind[c(hess.i, hess.j), 2L]



        if(derv[1L] == derv[2L]) {
          derv <- derv[1L]

          hessian[hess.i, hess.j] <-
          sum(w*(X[[derv]][,vars[1L]] * X[[derv]][,vars[2L]] * A[,derv] * (
            log(Y[,derv]) + psigamma(rsumA) - psigamma(A[,derv]) + A[,derv] * (
              psigamma(rsumA, 1L) - psigamma(A[,derv], 1L)
            )
          )))


        } else {
          hessian[hess.i, hess.j] <-
          sum(w*(
            X[[derv[1L]]][,vars[1L]]*X[[derv[2L]]][,vars[2L]]*A[,derv[1L]]*A[,derv[2L]]*psigamma(rsumA, 1L)
          ))
        }
      }
    }

    hessian <- make.symmetric(hessian)

    attr(LL, "hessian") <- hessian

  }

  return(LL)

}
