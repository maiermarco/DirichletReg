DirichReg <- function(
  formula,
  data,
  model = c("common", "alternative"),
  subset,
  sub.comp,                           # subcompositions
  base,
  weights,
  control,
  verbosity = getOption("verbose")
){
  
  if(!(verbosity %in% (0L:4L))){ # will become (verbosity %notin% (0L:4L)) in future releases
    verbosity <- 0L
    warning("invalid value for \"verbosity\".")
  } else verbosity <- as.integer(verbosity)
  
  if(verbosity > 0L){
    cat("- PREPARING DATA\n")
    if(interactive()) flush.console()
  }
  
  this_call <- match.call()
#  this_call$formula <- formula # TODO this should no longer be necessary -- test
  
  if(missing(formula)) stop("specification of \"formula\" is necessary.")
  original_formula <- as.formula(formula)
  
  model <- match.arg(model, c("common", "alternative"))
  repar <- model == "alternative"
  
  # checks and preliminary work
  if(missing(data)) data <- environment(formula)
  
  if(missing(control)){
    control <- list(sv = NULL, iterlim = 10000L, tol1 = .Machine$double.eps^(1.0 / 2.0), tol2 = .Machine$double.eps^(3.0 / 4.0))
  } else {
    if(!all(names(control) %in% c("sv", "iterlim", "tol1", "tol2")) || any(duplicated(names(control)))) stop("duplicated or inadmissible elements in \"control\".")
    if(is.null(control$sv))      control$sv       <-   NULL
    if(is.null(control$iterlim)) control$iterlim  <- 10000L
    if(is.null(control$tol1))    control$tol1     <- .Machine$double.eps^(1.0 / 2.0)
    if(is.null(control$tol2))    control$tol2     <- .Machine$double.eps^(3.0 / 4.0)
  }
  
  
  
  formula <- as.Formula(formula)
  formula_length <- length(formula) # as Formula now!
  if(length(formula_length) != 2L) stop("unexpected error.")
  if(formula_length[1L] != 1L) stop("problem on the left hand side of formula (check for | expressions)")
  
  
  
  call_model_frame <- match.call(expand.dots = FALSE)
  call_model_frame <- call_model_frame[c(1L, which(names(call_model_frame) %in% c("formula", "data", "subset", "weights")))] # maybe add "na.action" and "offset"?
  call_model_frame[[1L]] <- as.name("model.frame") # .default
  call_model_frame[["formula"]] <- formula
  call_model_frame[["drop.unused.levels"]] <- TRUE
  
  mf_formula <- call_model_frame # for compatibility w/old code
  
  
  model_frame <- eval(call_model_frame, parent.frame())
  
  
  Y_full <- model.response(model_frame, type = "numeric")
  if(!inherits(Y_full, "DirichletRegData")) stop("the response must be prepared by DR_data")
  n_obs <- nrow(Y_full)
  if(is.null(dim(Y_full))) stop('\"Y_full\" must be matrix-like.')
  if(n_obs < 1L) stop("empty model") # TODO better description!
  if(ncol(Y_full) < 2L) stop("\"Y_full\" has less than 2 components.")
  if(anyNA(Y_full)) stop("\"Y_full\" must noch contain missing values.")
  
  
  
  
  Y <- Y_full # TODO check
  attributes(Y)[c("Y.original", "dims", "dim.names", "obs", "valid_obs", "normalized", "transformed", "base", "class")] <- NULL
  Y <- as.matrix(Y)
  
  
  
  
## SUBCOMPOSITIONS # TODO CHECK
  if(missing(sub.comp)){
    sub.comp <- seq_len(ncol(Y))
  } else {
    if(length(sub.comp) == ncol(Y)) warning("no subcomposition made, because all variables were selected")
    if(length(sub.comp) == (ncol(Y) - 1L)) stop("no subcomposition made, because all variables except one were selected")
    if(any((sub.comp < 1L) | (sub.comp > ncol(Y)))) stop("subcompositions must contain indices of variables of the Dirichlet data object")
    y_in         <- seq_len(ncol(Y))[sub.comp]
    y_out        <- seq_len(ncol(Y))[-sub.comp]
    y_in_labels  <- colnames(Y)[y_in]
    y_out_labels <- paste(colnames(Y)[y_out], sep = "", collapse = " + ")
    Y            <- cbind(rowSums(Y[, y_out]), Y[, y_in])
    colnames(Y)  <- c(y_out_labels, y_in_labels)
  }
  
  if(missing(base)) base <- attr(Y_full, "base")
  if(!(base %in% seq_len(ncol(Y)))) stop("the base variable lies outside the number of variables")
  
  # effective components (possibly fewer than in the original Y_full due to subcompositions)
  n.dim <- ncol(Y)
  
  
  
  # SANITY CHECKS FOR THE FORMULA
  def_formula <- formula
  
  len_def_frm <- length(def_formula)
  if(length(len_def_frm) != 2L) stop("is \"def_formula\" not a \"Formula\" object?") # should not happen
  
  if(len_def_frm[1L] != 1L) stop("the left hand side of the model must contain one object prepared by DR_data()")
  
  if(!repar && (                                                     # if common model AND
      (len_def_frm[2L] > ncol(Y)) ||                                 #   more predictor sets than components in Y OR
      (ncol(Y) > 2L && (len_def_frm[2L] %in% (2L:(ncol(Y) - 1L))))   #   for Y w/3 or more components: predictor sets between 2:(n_components - 1)
    )) stop("the right hand side must contain specifications for either one or all variables")
  
  if( repar && (len_def_frm[2L] > 2L)) stop("the right hand side can only contain one or two specifications in the alternative parametrization")
  
  # expand "shorthands"
  if(len_def_frm[2L] == 1L){                # if there is only one RHS part in def_formula
    frml_parts <- as.character(def_formula) # split def_formula into c("~", <LHS>, <RHS>)
    if(length(frml_parts) != 3L) stop("unexpected error (more than 3 formlua parts as character).") # should not happen
    if(repar){                              # if an alternative model only has predictors in the mean model, add an intercept for the precision model
      message("Only mean model specified, assuming precision model to be ~ 1.")
      def_formula <- sprintf("%s ~ %s | 1", frml_parts[2L], frml_parts[3L])
    } else {                                # if a common model has only one set of predictors, use the same for every component
      message("Only one set of predictors specified, assuming the same for all components.")
      def_formula <- sprintf("%s ~ %s", frml_parts[2L], paste(rep(frml_parts[3L], ncol(Y)), collapse = " | "))
    }
  }
  def_formula <- as.Formula(def_formula)
  
  
  
  
  
  
  
  model_terms   <- terms(def_formula, data = data)
  model_terms_X <- if(repar){
      replicate(ncol(Y), simplify = FALSE, delete.response(terms(def_formula, data = data, rhs = 1L)))
    } else {
      lapply(seq_len(length(def_formula)[2L]), function(rhs) delete.response(terms(def_formula, data = data, rhs = rhs)))
    }
  model_terms_Z <- if(repar) delete.response(terms(def_formula, data = data, rhs = 2L)) else NULL
  
  # set up predictor variables for both models
  X.mats <- lapply(model_terms_X, model.matrix, data = model_frame)
  Z.mat  <- if(repar) model.matrix(model_terms_Z, data = model_frame) else NULL
  n.vars <- if(repar) c(unlist(lapply(X.mats, ncol))[-1L], ncol(Z.mat)) else unlist(lapply(X.mats, ncol))
  
  model_terms   <- .add_predvars_and_dataClasses(model_terms, model_frame)
  model_terms_X <- lapply(model_terms_X, function(this_X){ .add_predvars_and_dataClasses(this_X, model_frame) })
  model_terms_Z <- if(repar) .add_predvars_and_dataClasses(model_terms_Z, model_frame) else NULL
  
  weights <- model.weights(model_frame)
  if(is.null(weights)) weights <- 1.0
  if(length(weights) == 1L) weights <- rep(weights, n_obs)
  if(any(weights < 0.0)) stop("weights < 0 found.")
  weights <- as.vector(weights)
  names(weights) <- rownames(model_frame)
  
  
  
  
  ##############################################################################
  ### simplify and typecast variables ##########################################
  ##############################################################################
  
  # simplify dependent variable
  Y_fit <- unclass(Y)
  attributes(Y_fit) <- NULL
  dim(Y_fit) <- dim(Y)
  storage.mode(Y_fit) <- "double"
  
### TODO start check
  ### typecasting
  X_fit <- lapply(X.mats, function(this_mat){
    attr(this_mat, "dimnames") <- NULL
    attr(this_mat, "assign") <- NULL
    return(this_mat)
  })
  for(i in seq_along(X_fit)) storage.mode(X_fit[[i]]) <- "double"
  
  if(!is.null(Z.mat)) storage.mode(Z.mat) <- "double"
  
  storage.mode(n.dim)  <- "integer"
  storage.mode(n.vars) <- "integer"
  storage.mode(base)   <- "integer"
### TODO end check
  
  
  
  
  
  
  
  
  if(verbosity > 0L){
    cat("- COMPUTING STARTING VALUES\n")
    if(interactive()) flush.console()
  }
  
  # compute starting values
  if(is.null(control$sv)){
    starting.vals <- get_starting_values(Y = Y_fit, X.mats = X_fit,
                       Z.mat = if(repar) as.matrix(Z.mat) else Z.mat,
                       repar = repar, base = base, weights = weights) * if(repar){ 1.0 } else { 1.0 / n.dim }
  } else {
    if(length(control$sv) != sum(n.vars)) stop("wrong number of starting values supplied.")
    starting.vals <- control$sv
  }
  
  parametrization <- if(repar) "alternative" else  "common"
  
  
  
  
  
  
  
  
  
  
  if(verbosity > 0){
    cat("- ESTIMATING PARAMETERS\n")
    if(interactive()) flush.console()
  }
  
  # fit and store the results
  fit.res <- DirichReg_fit(Y     = Y_fit,
                           X     = X_fit,
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
    names(coefs) <- unlist(as.vector(c(rep(colnames(X.mats[[1L]]), n.dim - 1L), colnames(Z.mat))))
  } else {
    names(coefs) <- unlist(lapply(X.mats, colnames))
  }
  
  
  # FITTED VALUES
  if(repar){
    
    B <- matrix(0.0, nrow = n.vars[1L], ncol = n.dim)
    B[cbind(rep(seq_len(n.vars[1L]), (n.dim - 1L)), rep(seq_len(n.dim)[-base], each = n.vars[1]))] <- coefs[1L:((n.dim - 1L) * n.vars[1L])]
    
    g <- matrix(coefs[((n.dim - 1L) * n.vars[1L] + 1L):length(coefs)], ncol = 1L)
    
    XB <- exp(apply(B, 2L, function(b){ as.matrix(X.mats[[1L]]) %*% b }))
    MU <- apply(XB, 2L, function(x){ x / rowSums(XB) })
    
    PHI <- exp(as.matrix(Z.mat) %*% g)
    
    ALPHA <- apply(MU, 2L, "*", PHI)
    
  } else {
    
    B <- sapply(seq_len(n.dim), function(i){ coefs[(cumsum(c(0L, n.vars))[i] + 1L) : cumsum(n.vars)[i]] }, simplify = FALSE)
    
    ALPHA <- sapply(seq_len(n.dim), function(i){ exp(as.matrix(X.mats[[i]]) %*% matrix(B[[i]], ncol = 1L)) })
    
    PHI <- rowSums(ALPHA)
    MU  <- apply(ALPHA, 2L, "/", PHI)
    
  }
  
  colnames(ALPHA) <- varnames
  colnames(MU) <- varnames
  
  hessian <- fit.res$hessian
  
  vcov <- tryCatch(solve(-fit.res$hessian),
                   error = function(x){ return(matrix(NA_real_, nrow = nrow(hessian), ncol = ncol(hessian))) },
                   silent = TRUE)
  
  if(!repar){   ## COMMON                              #nolint
    coefnames <- apply(cbind(rep(varnames, n.vars), unlist(lapply(X.mats, colnames))), 1L, paste, collapse = ":")
  } else {   ## ALTERNATIVE
    coefnames <- apply(cbind(rep(c(varnames[-base], "(phi)"), n.vars), c(unlist(lapply(X.mats, colnames)[-base]), colnames(Z.mat))), 1L, paste, collapse = ":")
  }
  
  dimnames(hessian) <- list(coefnames, coefnames)
  dimnames(vcov)    <- list(coefnames, coefnames)
  shortnames        <- names(coefs)
  names(coefs)      <- coefnames
  
  se <- if(anyNA(vcov)) rep(NA_real_, length(coefs)) else sqrt(diag(vcov))
  
  res <- structure(list(
      call            = this_call,
      parametrization = parametrization,
      varnames        = varnames,
      n.vars          = n.vars,
      dims            = length(varnames),
      Y               = Y,
      X               = X.mats,
      Z               = Z.mat,
      sub.comp        = sub.comp,
      base            = base,
      weights         = weights,
      orig.resp       = Y_full,
      data            = data,
      d               = model_frame,                                  # the model.frame
      formula         = def_formula,                   # new def_formula instead of formula
      mf_formula      = mf_formula,
      npar            = length(coefs),
      coefficients    = coefs,
      coefnames       = shortnames,
      fitted.values   = list(mu = MU, phi = PHI, alpha = ALPHA),
      logLik          = fit.res$maximum,
      vcov            = vcov,
      hessian         = hessian,
      se              = se,
      optimization    = list(
                          convergence = fit.res$code,
                          iterations  = fit.res$iterations,
                          bfgs.it     = fit.res$bfgs.it,
                          message     = fit.res$message
                        ),
      # new components
      terms           = list(X = model_terms_X, Z = model_terms_Z, full = model_terms),
      levels          = list(
                          X    = lapply(model_terms_X, function(this_X){ .getXlevels(this_X, model_frame) }),
                          Z    = if(repar) .getXlevels(model_terms_Z, model_frame) else NULL,
                          full = .getXlevels(model_terms, model_frame)
                        ),
      contrasts       = list(
                          X = lapply(model_terms_X, function(this_X){ attr(this_X, "contrasts") }),
                          Z = if(repar) attr(Z.mat, "contrasts") else NULL
                        )
      ),
      class = "DirichletRegModel"
    )
  
  # remove stuff from maxLik in the parent frame...
  for(maxLik_ob in c("lastFuncGrad", "lastFuncParam")){
    if(exists(maxLik_ob, envir = parent.frame(), inherits = FALSE)) rm(list = maxLik_ob, envir = parent.frame(), inherits = FALSE)
  }
  # ...and make some space
  on.exit(gc(verbose = FALSE, reset = TRUE, full = TRUE), add = TRUE)
  
  res
}





# Adapted from package "betareg" (version 3.2-6, 2026-08-26, Achim Zeileis) https://doi.org/10.32614/CRAN.package.betareg
# https://codeberg.org/zeileis/betareg/src/commit/dc6c031917006426f7e82afbe83b2de26aa4c045/R/betareg.R#L44

# obtain correct subset of predvars/dataClasses to terms
.add_predvars_and_dataClasses <- function(terms, model.frame){
  ## original terms
  rval <- terms
  ## terms from model.frame
  nval <- if(inherits(model.frame, "terms")) model.frame else terms(model.frame)
  
  ## associated variable labels
  ovar <- sapply(as.list(attr(rval, "variables")), deparse)[-1L]
  nvar <- sapply(as.list(attr(nval, "variables")), deparse)[-1L]
  if(!all(ovar %in% nvar)) stop(gettextf(
      "The following terms variables are not part of the model.frame: %s",
      pretty_list(ovar[!(ovar %in% nvar)])
    ), domain = NA)
  ix <- match(ovar, nvar)
  
  ## subset predvars
  if(!is.null(attr(rval, "predvars"))) warning("terms already had 'predvars' attribute, now replaced")
  attr(rval, "predvars") <- attr(nval, "predvars")[1L + c(0L, ix)]
  
  ## subset dataClasses
  if(!is.null(attr(rval, "dataClasses"))) warning("terms already had 'dataClasses' attribute, now replaced")
  attr(rval, "dataClasses") <- attr(nval, "dataClasses")[ix]
  
  rval
}
