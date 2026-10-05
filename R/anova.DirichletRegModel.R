anova.DirichletRegModel <- function(object, ..., sorted = FALSE) {

  comp.objs <- list(object, ...)

  if(!all("DirichletRegModel" %in% unlist(lapply(comp.objs, class)))){
    stop('only models fitted using "DirichReg()" can be compared.')
  }
  for(i in seq_along(comp.objs)[-1L]){
    if(!identical(comp.objs[[i-1L]][["Y"]], comp.objs[[i]][["Y"]])){
      stop("models appear not to be nested.")
    }
  }

  n.mods  <- length(comp.objs)
  n.pars  <- unlist(lapply(comp.objs, `[[`, "npar"))
  calls   <- unlist(lapply(comp.objs, `[[`, "call"))
  sorting <- if(sorted){
               order(n.pars, decreasing = TRUE)
             } else {
               seq_len(n.mods)
             }
  comp.objs <- comp.objs[sorting]

  deviances <- -2.0 * unlist(lapply(comp.objs, `[[`, "logLik"))
  dev_diffs <- abs(c(NA_real_, deviances[1L] - deviances[-1L]))
  d_f       <- abs(c(NA_real_, n.pars[1L] - n.pars[-1L]))

  res <- structure(list(
    "Deviance"   = deviances,
    "N. par"     = n.pars,
    "Difference" = dev_diffs,
    "df"         = d_f,
    "Pr(>Chi)"   = pchisq(dev_diffs, ifelse(d_f <= 0.0, NA_real_, d_f), lower.tail = FALSE),
    "sorting"    = sorting,
    "n.mods"     = n.mods,
    "calls"      = calls
  ), class = "anova_DirichletRegModel")

  return(res)

}



print.anova_DirichletRegModel <- function(x, ...){
  
  if(interactive()) writeLines("")
  
  writeLines("Analysis of Deviance Table\n")
  
  model_numbers <- sprintf("Model %s", format(x$sorting))
  model_names   <- unlist(lapply(x$calls[x$sorting], deparse_nocutoff))
  
  for(i in seq_len(x$n.mods)){
    writeLines(strwrap(paste0(model_numbers[[i]], ": ", model_names[[i]]),
      width = getOption("width") - 2L, exdent = 2L))
  }
  
  writeLines("")
  
  res <- as.data.frame(x[seq_len(5L)])
  colnames(res) <- names(x[seq_len(5L)])
  rownames(res) <- model_numbers
  
  printCoefmat(res, cs.ind = 1L, tst.ind = 3L, P.values = TRUE, na.print = "")
  
  if(interactive()) writeLines("")
  
  invisible(x)
}
