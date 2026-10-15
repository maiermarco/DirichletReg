print.DirichletRegModel <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  
  .wd <- getOption("width")
  
  names(x$coefficients) <- x$coefnames
  
  maxLik_code_msg(x$optimization$convergence)
  
  if(interactive()) cat("\n")
  
  cat(gettext("Call"), ":\n", sep = "")
  cat(strwrap(deparse_nocutoff(x$call), .wd), "\n")
  cat(gettextf("using the %s parametrization", dQuote(x$parametrization)), "\n", sep = "")
  
  cat("\n")
  
  cat(gettextf("Log-likelihood: %s on %d df (%d BFGS + %d NR Iterations)", format(x$logLik, digits = digits), x$npar, x$optimization$bfgs.it, x$optimization$iterations), "\n", sep = "")
  
  cat("\n")
  
  coef.ind <- cumsum(x$n.vars)
  
  if(x$parametrization == "common"){
  
    for(i in seq_along(x$varnames)){
      cat(strrep("-", min(41L, .wd)), "\n", sep = "")
      cat(gettextf("Coefficients for variable no. %d: %s", i, x$varnames[i]), "\n", sep = "")
      print.default(
        x         = format(x$coefficients[(if(i == 1L) 1L else (coef.ind[i - 1L] + 1L)):coef.ind[i]], digits = digits),
        print.gap = 2L,
        quote     = FALSE
      )
    }
    cat(strrep("-", min(41L, .wd)), "\n", sep = "")
    
  } else {
    
    printed.var <- 1L
    set.size    <- ncol(x$X[[1L]])
    
    cat(toupper(gettext("Mean models")), ":\n", sep = "")
    
    for(i in seq_along(x$varnames)){
      if(i == x$base){
        cat(strrep("-", min(41L, .wd)), "\n", sep = "")
        cat(gettextf("Coefficients for variable no. %d: %s", i, x$varnames[i]), "\n", sep = "")
        cat("- ", gettext("variable omitted (reference category)"), " -\n", sep = "")
      } else {
        cat(strrep("-", min(41L, .wd)), "\n", sep = "")
        cat(gettextf("Coefficients for variable no. %d: %s", i, x$varnames[i]), "\n", sep = "")
        print.default(
          x         = format(x$coefficients[printed.var:(printed.var + set.size - 1L)], digits = digits),
          print.gap = 2L,
          quote     = FALSE
        )
        
        printed.var <- printed.var + set.size
      }
    }
    
    cat(strrep("-", min(41L, .wd)), "\n", sep = "")
    
    cat("\n")
    
    cat(toupper(gettext("Precision model")), ":\n", sep = "")
    
    cat(strrep("-", min(41L, .wd)), "\n", sep = "")
    print.default(
      x         = format(x$coefficients[printed.var:length(x$coefficients)], digits = digits),
      print.gap = 2L,
      quote     = FALSE
    )
    cat(strrep("-", min(41L, .wd)), "\n", sep = "")
  }
  
  if(interactive()) cat("\n")
  
  invisible(x)
}
