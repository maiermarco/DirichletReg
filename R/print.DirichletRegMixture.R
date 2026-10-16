print.DirichletRegMixture <- function(x, ...){
  
  if(interactive()) cat("\n")
  
  cat(gettextf("Dirichlet Mixture of %d Latent Classes", x$classes), "\n\n", gettext("Group Sizes and Probabilities"), ":\n", sep = "")
  
  sp.tab <- as.table(rbind(format(round(x$groups, 0L)),
                           format(round(x$groups / sum(x$groups), 3L))))
  dimnames(sp.tab) <- list(c(gettext("Size"), gettext("Probability")), gettextf("Class %d", seq_len(x$classes)))
  print(sp.tab, print.gap = 2L, justify = "right")
  
  cat("\n\n", gettext("Log-Likelihood"), ": ", round(x$logLik, 3L), " with ", x$npar, " Parameters", sep = "")
  cat("\nBIC: ", round(BIC(x), 1L), sep = "")
  cat("\n", gettext("Entropy"), ": ", round(x$entropy, 3L), "\n", sep = "")
  
  if(interactive()) cat("\n")
  
  invisible(x)
}
