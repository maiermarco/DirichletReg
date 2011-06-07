DirichReg_fit <- function(Y, X, Z, sv, d, k, w, ctls, repar, base){



  if(repar){
    bfgs <- suppressWarnings(maxBFGS(fn=DReg.repar,
              start=sv,
              finalHessian=FALSE, iterlim=ctls$iterlim, tol=1e-10, reltol=1e-10, print.level=ctls$trace,
              Y=Y, X=X[[1]], Z=Z, d=d, k=k[1], w=w, base=base, NR=FALSE)
            )
    res  <- suppressWarnings(maxNR(fn=DReg.repar,
              start=bfgs$estimate,
              iterlim=ctls$iterlim, tol=1e-13, reltol=1e-13, print.level=ctls$trace,
              Y=Y, X=X[[1]], Z=Z, d=d, k=k[1], w=w, base=base, NR=TRUE)
            )


  } else {
    bfgs <- suppressWarnings(maxBFGS(fn=DReg, 
              start=sv,
              finalHessian=FALSE, iterlim=ctls$iterlim, tol=1e-10, reltol=1e-10, print.level=ctls$trace,
              Y=Y, X=X, d=d, k=k, w=w, NR=FALSE)
            )
    res  <- suppressWarnings(maxNR(fn=DReg, 
              start=bfgs$estimate,
              iterlim=ctls$iterlim, tol=1e-13, reltol=1e-13, print.level=ctls$trace,
              Y=Y, X=X, d=d, k=k, w=w, NR=TRUE)
            )
  }

  res$bfgs.it <- bfgs$iterations

  return(res)

}
