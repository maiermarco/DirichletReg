print.DirichletRegMixture <- function(x, ...){

  cat("\n\nDirichlet Mixture of", x$classes, "Latent Classes\n\n\nGroup Sizes and Probabilities:\n")

  sp.tab <- as.table(rbind(format(round(x$groups, 0L)),
                           format(round(x$groups / sum(x$groups), 3L))))
  dimnames(sp.tab) <- list(c("Size", "Probability"), paste("Class", 1:x$classes))
  print(sp.tab, print.gap = 2L, justify = "right")

  cat("\n\nLog-Likelihood:", round(x$logLik, 3L), "with", x$npar, "Parameters")
  cat("\nBIC:", round(BIC(x), 1L))
  cat("\nEntropy:", round(x$entropy, 3L), "\n\n")
  
  invisible(x)
}
