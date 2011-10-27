DirichMixture <- function(Y, classes=2, EMs=5, verbosity=0, ctl=list(em.tol=1e-8, op.tol=1e-8)){

  if(class(Y) != "DirichletRegData") stop("Y must be prepared by DR_data first")
  
  X <- Y$Y
  dims <- Y$dims

  LLs <- numeric(EMs)
  
  params <- list()
  
  
  em.cycle <- 0
  for(i in 1:EMs){
    cat("."); flush.console()
    if(((em.cycle <- em.cycle + 1) %% 10) == 0) cat("\n")

    theta <- sapply(1:classes, function(i) rnorm(dims), simplify=FALSE)
    ll_diff <- -1
    ll <- ll_old <- rep(0, classes)
  
    em.iter <- 0
  
    while(abs(ll_diff) > ctl$em.tol){
      em.iter <- em.iter + 1
  
      
      dens <- as.data.frame(lapply(theta, function(th) ddirichlet(X, exp(th))))
      post.prob <- dens/rowSums(dens)
      
      
      opti <- sapply(1:classes, function(i){
                maxBFGS(fn = function(x){
                  log(post.prob[,i]*ddirichlet(X, exp(x)))
                },
                grad = function(x){
                  sapply(1:length(x), function(a){
                    post.prob[,i] * exp(x[a]) * (
                      psigamma(sum(exp(x))) - psigamma(exp(x[a])) + log(X[,a])
                    )
                  }, simplify=TRUE)
                },
                start=theta[[i]], tol=ctl$op.tol, finalHessian=FALSE)
              }, simplify=FALSE)
      ll <- unlist(lapply(opti, `[[`, "maximum"))
      ll_diff <- sum(abs(ll_old - ll))
      ll_old <- ll
      theta <- lapply(opti, `[[`, "estimate")
    }
    
    params[[em.cycle]] <- theta
    LLs[em.cycle] <- sum(ll)
  }

  optima <- which(LLs == max(LLs))
  
  theta <- params[[optima[1]]]
  ll <- LLs[[optima[1]]]

  dens <- as.data.frame(lapply(theta, function(th) ddirichlet(X, exp(th))))
  post.prob <- dens/rowSums(dens)
  colnames(post.prob) <- paste("Class", 1:classes)

  .e <- -post.prob*log(post.prob)
  if(any(is.na(.e))) .e[is.na(.e)] <- 1
  entropy <- 1 - rowSums(.e)/log(classes)
  
  res <- list(data=Y,
              classes=classes,
              groups=colSums(post.prob),
              coefficients=theta,
              logLik=ll,
              npar=sum(unlist(lapply(theta, length))),
              nobs=nrow(X),
              post.prob=post.prob,
              entropy=ifelse(classes == 1, NA, mean(entropy)),
              entropy.obs=ifelse(classes == 1, NA, unname(entropy)),
              EMs=EMs,
              LLs=LLs)

  class(res) <- "DirichletRegMixture"

  return(res)
  
}
