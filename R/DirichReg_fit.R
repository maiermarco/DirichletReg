DirichReg_fit <- function(Y, X, Z, sv, d, k, w, ctls, repar, base, vrb){

  n <- nrow(Y)
  npar <- length(sv)

  ops <- options(warn = -1L)
  on.exit(options(ops))

 
 
  if(repar){
    k <- k[1L]
    beta_ind <- as.integer(matrix(seq_len(d*k), nrow = k, ncol = d)[,-base])
    beta_x_ind <- seq_len((d-1L)*k)
    gamma_ind <- seq.int((d-1L)*k + 1L, npar)
    ncolX <- ncol(X[[1L]])
    ncolZ <- ncol(Z)

   





      bfgs <- maxBFGS(fn=DReg.repar,
        start=sv,
        finalHessian=FALSE, iterlim=ctls$iterlim, tol=ctls$tol1, reltol=ctls$tol1, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=FALSE)

      res <- maxNR(fn=DReg.repar,
        start=bfgs$estimate,
        iterlim=ctls$iterlim, tol=ctls$tol2, reltol=ctls$tol2, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=TRUE)

 
 
  } else {
    seq_along_d <- seq_len(d)
    ncolX <- unlist(lapply(X, ncol))
    beta_x_ind <- lapply(seq_along_d, function(i){ seq.int(cumsum(c(0L, k))[i] + 1L, cumsum(k)[i]) })

      bfgs <- maxBFGS(fn=DReg,
        start=sv,
        finalHessian=FALSE, iterlim=ctls$iterlim, tol=ctls$tol1, reltol=ctls$tol1, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X, ncolX=ncolX, n=n, d=d, k=k, w=w, npar=npar, seq_along_d=seq_along_d, bx=beta_x_ind, NR=FALSE)

      res <- maxNR(fn=DReg,
        start=bfgs$estimate,
        iterlim=ctls$iterlim, tol=ctls$tol2, reltol=ctls$tol2, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X, ncolX=ncolX, n=n, d=d, k=k, w=w, npar=npar, seq_along_d=seq_along_d, bx=beta_x_ind, NR=TRUE)
  }

  res$bfgs.it <- bfgs$iterations

  return(res)

}
