DReg <- function(x, Y, X, d, k, w, NR){
  npar <- length(x)
  seq_along_d <- seq_len(d)

  B <- lapply(seq_along_d, function(i){ x[ (cumsum(c(0L, k))[i] + 1L) : cumsum(k)[i] ] })

  A <- matrix(unlist(lapply(seq_along_d, function(i){ exp(X[[i]] %*% B[[i]]) })), nrow=nrow(Y))
  Aplus <- rowSums(A)






  
  
  LL <- .Call("weighted_logLL", Y, A, Aplus, dim(Y), w)







  
  
  
  
  
  
  
  
  
  
  

  attr(LL, "gradient") <- .Call("grad_common", Y, X, A, Aplus, npar, lapply(X, ncol), dim(A), w)







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
            log(Y[,derv]) + psigamma(Aplus) - psigamma(A[,derv]) + A[,derv] * (
              psigamma(Aplus, 1L) - psigamma(A[,derv], 1L)
            )
          )))
        
        
        } else {
          hessian[hess.i, hess.j] <-
          sum(w*(
            X[[derv[1L]]][,vars[1L]]*X[[derv[2L]]][,vars[2L]]*A[,derv[1L]]*A[,derv[2L]]*psigamma(Aplus, 1L)
          ))
        }
      }
    }

    hessian <- make.symmetric(hessian)

    attr(LL, "hessian") <- hessian

  }

  return(LL)

}
