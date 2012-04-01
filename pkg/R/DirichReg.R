DirichReg <- function(formula,
                      data,
                      model = c("common", "alternative"),
                      subset,
                      sub.comp,                                                 
                      base,
                      weights,
                      control,
                      verbosity=0
                      ){

  this.call <- match.call()
  
  if(!(verbosity %in% 0:4)){
    verbosity <- 0
    warning("invalid value for verbosity.")
  }

if(verbosity > 0){
  cat("- PREPARING DATA\n")
  if(interactive()) flush.console()
}

  
  if(missing(data)) data <- environment(formula)

  if(missing(formula)){
    stop("specification of \"formula\" is necessary.")
  } else {
    oformula <- formula
    formula <- as.Formula(formula)
  }
  model <- match.arg(model)
  if(missing(control)){
    control <- list(sv=NULL, iterlim=1000, tol1=1e-5, tol2=1e-10)
  } else {
    if(is.null(control$sv))      control$sv       <-  NULL
    if(is.null(control$iterlim)) control$iterlim  <-  1000
    if(is.null(control$tol1))    control$tol1     <-  1e-5
    if(is.null(control$tol2))    control$tol2     <- 1e-10
  }

  if(ifelse(!missing(data), deparse(oformula[[2]]) %in% names(data), FALSE)){
    Y_full <- data[[deparse(oformula[[2]])]]                                    
  } else {                                                                      
    Y_full <- get(deparse(oformula[[2]]), environment(oformula))
    if(missing(data)){
      data[[deparse(oformula[[2]])]] <- Y_full
    } else {
      assign(deparse(oformula[[2]]), Y_full)
    }
  } 
  if(class(Y_full) != "DirichletRegData") stop("the response must be prepared by DR_data")

  repar <- ifelse(model == "common", FALSE, TRUE)
  
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1, match(c("formula", "data", "subset", "weights"), names(mf), 0))]
  mf$formula <- as.Formula(mf$formula)
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf_formula <- mf
  d <- mf <- eval(mf, parent.frame())
  
  if("(weights)" %in% names(mf)) weights <- mf[["(weights)"]] else weights <- rep(1, nrow(mf))

  Y <- model.response(mf, "numeric")

  if(missing(sub.comp)){
    sub.comp <- 1:ncol(Y)
  } else {
    if(length(sub.comp) == ncol(Y)) warning("no subcomposition made, because all variables were selected")
    if(length(sub.comp) == (ncol(Y) - 1)) stop("no subcomposition made, because all variables except one were selected")
    if(any((sub.comp < 1) | (sub.comp > ncol(Y)))) stop("subcompositions must contain indices of variables of the Dirichlet data object")
    y_in  <- (1:ncol(Y))[sub.comp]
    y_out <- (1:ncol(Y))[-sub.comp]
    y_in_labels <- colnames(Y)[y_in]
    y_out_labels <- paste(colnames(Y)[y_out], sep="", collapse=" + ")
    Y <- cbind(rowSums(Y[,y_out]), Y[,y_in])
    colnames(Y) <- c(y_out_labels, y_in_labels)
  }
  
  base <- ifelse(missing(base), attr(Y_full, "base"), base)
  if(!(base %in% (1:ncol(Y)))) stop("the base variable lies outside the number of variables")

  n.dim <- ncol(Y)


  if(length(formula)[1] != 1) stop("the left hand side of the model must contain one object prepared by DR_data()")

  if(!repar){   
    if(length(formula)[2] == 1) for(i in 2:ncol(Y)) attr(formula, "rhs") <- lapply(1:ncol(Y), function(i) attr(formula, "rhs")[[1]])
    if(length(formula)[2] > ncol(Y)) stop("the right hand side must contain specifications for either one or all variables")
  } else {   
    if(length(formula)[2] == 1) formula <- as.Formula(formula(formula), ~ 1)
    if(length(formula)[2] > 2) stop("the right hand side can only contain one or two specifications in the alternative parametrization")
  }



  if(!repar){
    X.mats <- lapply(1:ncol(Y), function(i) model.matrix(terms(formula, data=data, rhs=i), mf) )
    Z.mat <- NULL
    n.vars <- unlist(lapply(X.mats, ncol))
  } else {
    X.mats <- lapply(1:ncol(Y), function(i) model.matrix(terms(formula, data=data, rhs=1), mf) )
    Z.mat  <- model.matrix(terms(formula, data=data, rhs=2), mf)
    n.vars <- c(unlist(lapply(X.mats, ncol))[-1], ncol(Z.mat))
  }
  
  

if(verbosity > 0){
  cat("- COMPUTING STARTING VALUES\n")
  if(interactive()) flush.console()
}

  
  if(is.null(control$sv)){
    starting.vals <- get_starting_values(Y=Y, X.mats=lapply(X.mats, as.matrix),
                       Z.mat={if(repar) as.matrix(Z.mat) else Z.mat},
                       repar=repar, base=base, weights=weights) * if(repar){ 1 } else { 1/n.dim }
  } else {
    if(length(control$sv) != n.vars) stop("wrong number of starting values supplied.")
    starting.vals <- control$sv
  }

  parametrization <- ifelse(repar, "alternative", "common")

if(verbosity > 0){
  cat("- ESTIMATING PARAMETERS\n")
  if(interactive()) flush.console()
}

  
  fit.res <- DirichReg_fit(Y     = Y,
                           X     = lapply(X.mats, as.matrix),
                           Z     = as.matrix(Z.mat),
                           sv    = starting.vals,
                           d     = n.dim,
                           k     = n.vars,
                           w     = as.vector(weights),
                           ctls  = control,
                           repar = repar,
                           base  = base,
                           vrb   = verbosity)
  

  varnames <- colnames(Y)
  
  coefs <- fit.res$estimate

  if(repar){
    names(coefs) <- unlist(as.vector(c(rep(colnames(X.mats[[1]]),n.dim-1),colnames(Z.mat))))
  } else {
    names(coefs) <- unlist(as.vector(sapply(X.mats, colnames)))
  }


  
  if(repar){

    B <- matrix(0, nrow=n.vars[1], ncol=n.dim)
    B[cbind(rep(1:n.vars[1], (n.dim-1)), rep((1:n.dim)[-base], each=n.vars[1]))] <- coefs[1:((n.dim-1)*n.vars[1])]
    
    g <- matrix(coefs[((n.dim-1)*n.vars[1]+1):length(coefs)], ncol=1)
  
    XB <- exp(apply(B, 2, function(b){ as.matrix(X.mats[[1]]) %*% b }))
    MU <- apply(XB, 2, function(x){ x /rowSums(XB) })
  
    PHI <- exp(as.matrix(Z.mat) %*% g)
    
    ALPHA <- apply(MU, 2, "*", PHI)
  
  } else {

    B <- sapply(1:n.dim, function(i){ coefs[(cumsum(c(0,n.vars))[i]+1) : cumsum(n.vars)[i]] }, simplify=FALSE)
    
    ALPHA <- sapply(1:n.dim, function(i){ exp(as.matrix(X.mats[[i]]) %*% matrix(B[[i]], ncol=1)) })

    PHI <- rowSums(ALPHA)
    MU  <- apply(ALPHA, 2, "/", PHI)

  }
  
  colnames(ALPHA) <- varnames
  colnames(MU) <- varnames

  hessian <- fit.res$hessian

  vcov <- tryCatch(solve(-fit.res$hessian),
                   error=function(x){ return(matrix(NA, nrow=nrow(hessian), ncol=ncol(hessian))) },
                   silent=TRUE)
                   
  if(!repar){   
    coefnames <- apply(cbind(rep(varnames, n.vars), unlist(lapply(X.mats, colnames))), 1, paste, collapse=":")
  } else {   
    coefnames <- apply(cbind(rep(c(varnames[-base], "(phi)"), n.vars), c(unlist(lapply(X.mats, colnames)[-base]), colnames(Z.mat))), 1, paste, collapse=":")
  }

  dimnames(hessian) <- list(coefnames, coefnames)
  dimnames(vcov) <- list(coefnames, coefnames)
  shortnames <- names(coefs)
  names(coefs) <- coefnames

  se <- if(!any(is.na(vcov))) sqrt(diag(vcov)) else rep(NA,length(coefs))
  
  res <- list(call=this.call,
              parametrization=parametrization,
              varnames=varnames,
              n.vars=n.vars,
              dims=length(varnames),
              Y=Y,
              X=X.mats,
              Z=Z.mat,
              subset=if(missing(subset)) NULL else subset,
              sub.comp=sub.comp,
              base=base,
              weights=weights,
              orig.resp=Y_full,
              data=data,
              d=d,
              formula=formula,
              mf_formula=mf_formula,
              npar=length(coefs),
              coefficients=coefs,
              coefnames=shortnames,
              fitted.values=list(mu=MU,phi=PHI,alpha=ALPHA),
              logLik=fit.res$maximum,
              vcov=vcov,
              hessian=hessian,
              se=se,
              optimization=list(convergence=fit.res$code,
                                iterations=fit.res$iterations,
                                bfgs.it=fit.res$bfgs.it,
                                message=fit.res$message))
              
  class(res) <- "DirichletRegModel"
              
  return(res)


}
