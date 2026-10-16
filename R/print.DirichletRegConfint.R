print.DirichletRegConfint <- function(x, digits = 3L, ...){
  
  e     <- x$e
  repar <- x$repar
  
  if(interactive()) cat("\n")
  
  cat(
    sprintf("%s%% ", paste(round(100.0 * x$level, 3L), collapse = "%, ")), gettext("Confidence Intervals"), " (",
    if(e) gettext("exponentiated") else gettext("original form"), ")\n\n", sep = ""
    )
  
  ctab <- x$coefficients
  
  
  if(repar){
  ##############################################################################
  ################################################################## ALTERNATIVE
    if(e){
      for(v in seq_along(ctab)){
        ctab[[v]] <- lapply(ctab[[v]], function(l) if(is.null(l)) NULL else exp(l))
      }
    }
    
    iter <- 1L
    for(cc in rev(seq_along(x$ci[[1L]]))[-1L]){
        iter <- iter + 1L
      for(ll in seq_along(x$ci)){
        if(is.null(x$ci[[1L]][[cc]])) next
        ctab[[1L]][[cc]] <- cbind(ctab[[1L]][[cc]], x$ci[[ll]][[cc]])
      }
    }
    for(ll in seq_along(x$ci)){
      ctab[[2L]][[1L]] <- cbind(ctab[[2L]][[1L]], x$ci[[ll]][[iter]])
    }
    
    for(cc in seq_along(ctab[[1L]])){
      if(is.null(ctab[[1L]][[cc]])) next
      ll <- ncol(ctab[[1L]][[cc]])
      ind <- c(rev(seq(2L, ll, by = 2L)), seq(1L, ll, by = 2L))
      
      ctab[[1L]][[cc]] <- ctab[[1L]][[cc]][, ind, drop = FALSE]
    }
    ll <- ncol(ctab[[2L]][[1L]])
    ind <- c(rev(seq(2L, ll, by = 2L)), seq(1L, ll, by = 2L))
    ctab[[2L]][[1L]] <- ctab[[2L]][[1L]][, ind, drop = FALSE]
    
  } else {
  ##############################################################################
  ####################################################################### COMMON
    if(e){
      ctab <- lapply(ctab, exp)
    }
    
    for(cc in rev(seq_along(x$ci[[1L]]))){
      for(ll in seq_along(x$ci)){
        ctab[[cc]] <- cbind(ctab[[cc]], x$ci[[ll]][[cc]])
      }
    }
    for(cc in seq_along(x$ci[[1L]])){
      ll <- ncol(ctab[[cc]])
      ind <- c(rev(seq(2L, ll, by = 2L)), seq(1L, ll, by = 2L))
      
      ctab[[cc]] <- ctab[[cc]][, ind, drop = FALSE]
    }
    
  }
  
  
  lo_lab <- paste0(format(100.0 * (1.0 - rev(x$level)) / 2.0       , format = "f"), "%")
  hi_lab <- paste0(format(100.0 * (x$level + (1.0 - x$level) / 2.0), format = "f"), "%")
  
  
  if(repar){
  ##############################################################################
  ################################################################## ALTERNATIVE
    for(tt in seq_len(2L)){
    cat("- ", if(tt == 1L) gettext("Beta-Parameters") else gettext("Gamma-Parameters"), ":\n", sep = "")
      for(i in seq_along(ctab[[tt]])){
        if(tt == 1L) cat(gettext("Variable"), ": ", names(ctab[[1L]])[i], "\n", sep = "")
        if(is.null(ctab[[tt]][[i]])){ cat("  ", gettext("variable omitted"), "\n\n", sep = ""); next }
        ttab <- round(ctab[[tt]][[i]], digits)
        colnames(ttab) <- c(lo_lab, if(e) "exp(Est.)" else "Est.", hi_lab)
        print(ttab, digits = digits, print.gap = 2L)
      }
    }
  } else {
  ##############################################################################
  ####################################################################### COMMON
    for(i in seq_along(ctab)){
      cat(gettext("Variable"), ": ", names(x$coefficients)[i], "\n", sep = "")
      ttab <- round(ctab[[i]], digits)
      colnames(ttab) <- c(lo_lab, if(e) "exp(Est.)" else "Est.", hi_lab)
      print(ttab, digits = digits, print.gap = 2L)
    }
  }
  
  if(interactive()) cat("\n")
  
  invisible(x)
}
