summary.DirichletRegModel <- function(object, ...) {
  
  maxLik_code_msg(object$optimization$convergence)
  
  res <- structure(
    .Data = list(
      call            = object[["call"]],
      terms           = terms(object$formula),
      logLik          = object$logLik,
      df              = object$npar,
      deviance        = -2.0 * object$logLik,
      aic             = AIC(object),
      bic             = BIC(object),
      residuals       = residuals(object, type = "standardized"),
      coefficients    = object$coefficients,
      varnames        = object$varnames,
      base            = object$base,
      n.vars          = object$n.vars,
      npar            = object$npar,
      coef.ind        = cumsum(object$n.vars),
      nobs            = nobs(object),
      parametrization = object$parametrization,
      optimization    = object$optimization
    ),
    class = "summary_DirichletRegModel"
  )
  
  names(object$coefficients) <- object$coefnames
  
  resid.mat <- round(t(apply(residuals(object, type = "standardized"), 2L, quantile)), 4L)
  colnames(resid.mat) <- c("Min", "1Q", "Median", "3Q", "Max")
  
  res$resid.mat <- resid.mat
  
  z.values <- object$coefficients / object$se
  p.values <- 2.0 * pnorm(-abs(z.values))
  coef.mat <- cbind(object$coefficients, object$se, z.values, p.values)
  colnames(coef.mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  
  res$coef.mat <- coef.mat
  
  res
}



print.summary_DirichletRegModel <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  .wd <- getOption("width")
  
  signif_codes <- gettextf("Significance codes: %s", attr(symnum(0.0, corr = FALSE, na = FALSE,
    cutpoints = c(0.0, 0.001, 0.01, 0.05, 0.1, 1.0),
    symbols   = c("***", "**", "*", ".", " ")), "legend"))
  sig_code_len <- nchar(signif_codes)
  
  if(interactive()) cat("\n")
  
  cat(gettext("Call:\n"))
  cat(strwrap(deparse_nocutoff(getCall(x)), .wd, exdent = 4L), sep = "\n")
  
  cat(gettext("\nStandardized Residuals:\n"))
  print(x$resid.mat, digits = digits, print.gap = 2L)
  cat("\n")
  
  ##############################################################################
  ################################################################# COMMON PARAM
  
  if(x$parametrization == "common"){
    
    for(i in seq_along(x$varnames)){
      cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
      cat(gettextf("Beta-Coefficients for variable no. %d: %s\n", i, x$varnames[i]))
      
      printCoefmat(
        x             = x$coef.mat[(if(i == 1L) 1L else (x$coef.ind[i - 1L] + 1L)):x$coef.ind[i], , drop = FALSE],
        digits        = digits,
        cs.ind        = 1L:2L,
        tst.ind       = 3L,
        P.values      = TRUE,
        signif.legend = FALSE
      )
    }
    cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
    cat(signif_codes, "\n", sep = "")
    
  } else {
  
  ##############################################################################
  ############################################################ ALTERNATIVE PARAM
  
    printed.var <- 1L
    set.size    <- x$n.vars[1L]
    
    cat(toupper(gettext("Mean models:\n")))
    
    for(i in seq_along(x$varnames)){
      if(i == x$base){
        cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
        cat(gettextf("Coefficients for variable no. %d: %s", i, x$varnames[i]), "\n", sep = "")
        cat("- ", gettext("variable omitted (reference category)"), " -\n", sep = "")
      } else {
        cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
        cat(gettextf("Coefficients for variable no. %d: %s", i, x$varnames[i]), "\n", sep = "")
        
        printCoefmat(
          x             = x$coef.mat[printed.var:(printed.var + set.size - 1L), , drop = FALSE],
          digits        = digits,
          cs.ind        = 1L:2L,
          tst.ind       = 3L,
          P.values      = TRUE,
          signif.legend = FALSE
        )
        
        printed.var <- printed.var + set.size
      }
    }
    
    cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
    
    cat("\n")
    
    cat(toupper(gettext("Precision model")), ":\n", sep = "")
    cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
    printCoefmat(
      x             = x$coef.mat[printed.var:length(x$coefficients), , drop = FALSE],
      digits        = digits,
      cs.ind        = 1L:2L,
      tst.ind       = 3L,
      P.values      = TRUE,
      signif.legend = FALSE
    )
    cat(strrep("-", min(sig_code_len, .wd)), "\n", sep = "")
    cat(signif_codes, "\n", sep = "")
    
  }
  
  ##############################################################################
  ########################################################################## FIN
  
  cat("\n")
  cat(gettextf("Log-likelihood: %s on %d df (%d BFGS + %d NR Iterations)", format(x$logLik, digits = digits), x$npar, x$optimization$bfgs.it, x$optimization$iterations), "\n", sep = "")
  cat(sprintf("AIC: %s, BIC: %s", format(x$aic, digits = digits), format(x$bic, digits = digits)), "\n", sep = "")
  cat(gettextf("Number of Observations: %d", x$nobs), "\n", sep = "")
  if(x$parametrization == "common"){
    cat(gettextf("Link: Log\nParametrization: %s", dQuote(first_cap(x$parametrization, all_lower = TRUE))), "\n", sep = "")
  } else {
    cat(gettextf("Links: Logit (Means) and Log (Precision)\nParametrization: %s", dQuote(first_cap(x$parametrization, all_lower = TRUE))), "\n", sep = "")
  }
  
  if(interactive()) cat("\n")
  
  invisible(x)

}



coef.summary_DirichletRegModel <- function(object, ...){
  
  coefmat <- as.matrix(object$coef.mat)
  rownames(coefmat) <- names(object$coefficients)
  
  coefmat

}
