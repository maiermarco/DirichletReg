print.DirichletRegData <- function(x, type = c("processed", "original"), ...){
  allowed_dddots <- c("digits", "na.print", "print.gap", "max", "width") # "quote", "useSource", "right"
  chkDots(..., allowed = allowed_dddots)
  
  type <- match.arg(type)
  
  if(type == "processed"){
    out <- unclass(x)
    attr(out, "Y.original" ) <- NULL
    attr(out, "dims"       ) <- NULL
    attr(out, "dim.names"  ) <- NULL
    attr(out, "obs"        ) <- NULL
    attr(out, "valid_obs"  ) <- NULL
    attr(out, "normalized" ) <- NULL
    attr(out, "transformed") <- NULL
    attr(out, "base"       ) <- NULL
    attr(out, "class"      ) <- NULL
  } else {
    out <- attr(x, "Y.original")
  }
  
  do.call("print", c(list("x" = out), list(...)[...names() %in% allowed_dddots]))
  
  invisible(x)
}



`[.DirichletRegData` <- function (x, i, j, drop){
  has_i <- !missing(i)
  has_j <- !missing(j)
  
  if(has_i && !is.matrix(i)){
    if(is.logical(i) && length(i) != nrow(x)) stop("if logical, i must exactly have nrow(x) elements.")
  }
  if(has_j && is.logical(j) && length(j) != ncol(x)) stop("if logical, j must have exactly ncol(x) elements.")
  if(has_j && is.character(j)){
    if(any(duplicated(j))){ j <- j[!duplicated(j)]; warning("duplicates in j were removed.") }
    if(!all(j %in% colnames(x))) stop("some names in j were not found in x.")
  }
  
  if(has_i && has_j && is.matrix(i)){ warning("???"); browser() }
  
  o <- attr(x, "Y.original") # TODO this should never not be a matrix...
  
  if( has_i && is.matrix(i)){
    res <- .subset(x, i)
    ooo <- .subset(o, i)
  } else if(!has_i && !has_j){
    res <- .subset(x)
    ooo <- .subset(o)
  } else if( has_i && !has_j){
    res <- .subset(x, i = i, j = rep_len(TRUE, ncol(x)))
    ooo <- .subset(o, i = i, j = rep_len(TRUE, ncol(o)))
  } else if(!has_i &&  has_j){
    res <- .subset(x, i = rep_len(TRUE, nrow(x)), j = j, drop = if(missing(drop)) TRUE else drop)
    ooo <- .subset(o, i = rep_len(TRUE, nrow(o)), j = j, drop = if(missing(drop)) TRUE else drop)
  } else if( has_i &&  has_j){
    res <- .subset(x, i = i, j = j, drop = if(missing(drop)) TRUE else drop)
    ooo <- .subset(o, i = i, j = j, drop = if(missing(drop)) TRUE else drop)
  } else {
    warning("unexpected error!")
    browser() # this should not happen ...
  }
  
  # if x results in a vector, return the object without classes or attributes
  if(is.null(dim(res))){ return(unclass(res)) }
  
  # restore attributes from x
  attr(res, "Y.original" ) <- ooo
  attr(res, "dims"       ) <- ncol(res)
  attr(res, "dim.names"  ) <- colnames(res)
  attr(res, "obs"        ) <- nrow(res)
  attr(res, "valid_obs"  ) <- sum(complete.cases(res))
  attr(res, "normalized" ) <- attr(x, "normalized" )
  attr(res, "transformed") <- attr(x, "transformed")
  attr(res, "base"       ) <- attr(x, "base"       )
  
  class(res) <- "DirichletRegData"
  
  return(res)
}



str.DirichletRegData <- function(object, ...){
  NextMethod("str")
}
