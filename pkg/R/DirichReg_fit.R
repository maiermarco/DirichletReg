DirichReg_fit <- function(Y, X, Z, sv, d, k, w, ctls, repar, base){# BEGIN DirichReg.fit

################################################################################
################################################ alternative parametrization ###
  if(repar){
    ## first a couple of iterations only for precisions to accelerate
    bfgs1 <- suppressWarnings(maxBFGS(fn=DReg.repar,
               start=sv, fixed=1:((d-1)*ncol(X[[1]])),
               finalHessian=FALSE, iterlim=10, tol=1e-2, reltol=1e-2, print.level=ctls$trace,
               Y=Y, X=X[[1]], Z=Z, d=d, k=k[1], w=w, base=base, NR=FALSE)
             )
    bfgs <- suppressWarnings(maxBFGS(fn=DReg.repar,
               start=bfgs1$estimate,
               finalHessian=FALSE, iterlim=ctls$iterlim, tol=1e-7, reltol=1e-7, print.level=ctls$trace,
               Y=Y, X=X[[1]], Z=Z, d=d, k=k[1], w=w, base=base, NR=FALSE)
             )
    res  <- suppressWarnings(maxNR(fn=DReg.repar,
              start=bfgs$estimate,
              iterlim=ctls$iterlim, tol=1e-13, reltol=1e-13, print.level=ctls$trace,
              Y=Y, X=X[[1]], Z=Z, d=d, k=k[1], w=w, base=base, NR=TRUE)
            )
################################################################################
##################################################### common parametrization ###
  } else {
    bfgs <- suppressWarnings(maxBFGS(fn=DReg, 
              start=sv,
              finalHessian=FALSE, iterlim=ctls$iterlim, tol=1e-7, reltol=1e-7, print.level=ctls$trace,
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

}### END DirichReg.fit
