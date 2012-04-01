predict.DirichletRegModel <- function(object, newdata, mu=TRUE, alpha=FALSE, phi=FALSE, ...){

  if(missing(newdata)) return(fitted(object, mu, alpha, phi))

  repar <- object$parametrization == "alternative"
  dims  <- ncol(object$Y)





  model_formula <- object$mf_formula
  model_formula$data <- as.name("newdata")
  model_formula$lhs <- 0
  mf <- eval(model_formula)
  if(!repar){ 
    if(length(model_formula$formula)[2]){
      X <- lapply(rep(1,ncol(object$Y)), function(i) model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=i)), mf) )
    } else {             
      X <- lapply(1:ncol(object$Y), function(i) model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=i)), mf) )
    }
    Z <- NULL
  } else { 
    X <- lapply(1:ncol(object$Y), function(i) model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=1)), mf) )
    Z <- model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=2)), mf)
  }
  
  cc <- coef(object)
  if(repar) cc[[1]] <- unlist(cc[[1]])
  if(!repar) cc <- unlist(cc)
  
  from <- 1
  base <- object$base

  if(repar){
    ETA <- matrix(0, nrow=nrow(newdata), ncol=ncol(object$Y))
    for(i in (1:dims)[-base]){
      ETA[,i] <- X[[i]] %*% cc[[1]][from:(from+ncol(X[[i]])-1)]
      from <- from + ncol(X[[i]])
    }
    MU  <- exp(ETA)/rowSums(exp(ETA))
    PHI <- exp(Z %*% unlist(cc[[2]]))
    ALPHA <- MU*as.numeric(PHI)
  } else {
    ALPHA <- matrix(0, nrow=nrow(newdata), ncol=ncol(object$Y))
    for(i in (1:dims)){
      ALPHA[,i] <- exp(X[[i]] %*% cc[from:(from+ncol(X[[i]])-1)])
      from <- from + ncol(X[[i]])
    }
    PHI <- rowSums(ALPHA)
    MU <- ALPHA/PHI
  }
  
  if(!any(mu | alpha | phi)) stop("Either mu, alpha or phi has to be requested.")

  if(sum(mu + alpha + phi) == 1){
    if(mu)    return(MU)
    if(alpha) return(ALPHA)
    if(phi)   return(PHI)
  } else {
    res <- list()
    if(mu)    res[["mu"]]    <- MU
    if(alpha) res[["alpha"]] <- ALPHA
    if(phi)   res[["phi"]]   <- PHI
    
    return(res)
  }
  
}
