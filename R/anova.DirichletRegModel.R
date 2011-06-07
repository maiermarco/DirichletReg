anova.DirichletRegModel <- function(object, ..., sorted = FALSE){

  comp.objs <- list(object, ...)
  
  if(any(!(sapply(comp.objs, class, simplify=TRUE) == "DirichletRegModel"))) stop("only models fitted using DirichReg can be compared.")
  
  n.mods <- length(comp.objs)
  
  for(i in 2:n.mods){
    if(!all.equal(comp.objs[[i-1]]$Y, comp.objs[[i]]$Y)) stop("models appear not to be nested.")
  }
  
  n.pars <- rep(NA, n.mods)
  
  for(i in 1:n.mods){
    n.pars[i] <- comp.objs[[i]]$npar
  }

  if(sorted){
    sorting <- order(n.pars, decreasing=TRUE)
  } else {
    sorting <- 1:n.mods
  }

  cat("\nAnalysis of Deviance Table\n\n")
  
  for(i in 1:n.mods){
    cat("Model ",i,":\n",
      paste(strsplit(deparse(comp.objs[[i]]$call), split=getOption("width")),sep="\n",collapse="\n"),
      "\n",sep="",collapse="")
  }
  
  cat("\n")
  
  res <- matrix(NA, ncol=5, nrow=n.mods)
  
  colnames(res) <- c("Deviance", "N. par", "Difference", "df", "p-value")
  rownames(res) <- paste("Model",sorting)
  
  r <- 1
  for(i in sorting){
    res[r,1:2] <- c(-2*comp.objs[[i]]$logLik, comp.objs[[i]]$npar)
    r <- r + 1
  }

  for(i in 2:n.mods){
    res[i,3] <- abs(res[i-1, 1] - res[i, 1])
    res[i,4] <- abs(res[i-1, 2] - res[i, 2])
  }

  res[-1,5] <- 1 - pchisq(res[-1,3], res[-1,4])
      
  print.default(res, quote=F, print.gap=3,na.print="-")
  
  cat("\n")
  
}
