summary.DirichletRegModel <- function(object, digits=max(3, getOption("digits") - 3), ...) 
{
    .wd <- getOption("width")

    if(object$optimization$convergence > 2) stop("\n",paste(strwrap(paste("\nOptimization did not converge in",object$optimization$bfgs.it,"+",object$optimization$iterations,"iterations and exited with code",object$optimization$convergence),.wd),sep="\n",collapse="\n"))
    
    resid.mat <- round(t(apply(residuals(object, type="standardized"), 2, quantile)), 4)
    colnames(resid.mat) <- c("Min", "1Q", "Median", "3Q", "Max")

    cat("\nCall:\n",
        paste(strwrap(deparse(object$call), .wd), sep="\n", collapse="\n"),
        "\n\n", sep = "")

    cat("\nStandardized Residuals:\n")
    print(resid.mat, print.gap=2)
    cat("\n\n")
    
    coef.ind <- cumsum(object$n.vars)
    
    z.values <- object$coefficients / object$se
    p.values <- 2 * pnorm(-abs(z.values))
    coef.mat <- cbind(object$coefficients, object$se, z.values, p.values)
    colnames(coef.mat) <- c("Estimate","Std. Error","z-Value","p-Value")

################################################################################
################################################################### COMMON PARAM

    if(object$parametrization == "common"){
    
      for(i in 1:length(object$varnames)){
        cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
        cat("\nBeta-Coefficients for variable no. ",i,": ",object$varnames[i],"\n",sep="",collapse="")


        printCoefmat(coef.mat[ ifelse(i==1,1,coef.ind[i-1]+1):coef.ind[i] , , drop=F],
                     digits = digits, cs.ind=1:2, tst.ind=3, has.Pvalue=T, signif.legend = F)


      }
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n",sep="",collapse="")
      cat(paste(strwrap("Signif. codes:  \u60***' < .001, \u60**' < 0.01, \u60*' < 0.05, \u60.' < 0.1",.wd),sep="\n",collapse="\n"),sep="")

    } else {

################################################################################
############################################################## ALTERNATIVE PARAM

      printed.var <- 1
      set.size    <- ncol(object$X[[1]])

      cat("\nMEAN MODELS:\n",sep="",collapse="")

      for(i in 1:length(object$varnames)){
        if(i == object$orig.resp$base){
          cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
          cat("\nCoefficients for variable no. ",i,": ",object$varnames[i],"\n",sep="",collapse="")
          cat("- variable omitted (reference category) -\n")
        } else {
          cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
          cat("\nCoefficients for variable no. ",i,": ",object$varnames[i],"\n",sep="",collapse="")

          coef.mat <- cbind(object$coefficients[printed.var:(printed.var+set.size-1)],
                            object$se[printed.var:(printed.var+set.size-1)],
                            z.values[printed.var:(printed.var+set.size-1)],
                            p.values[printed.var:(printed.var+set.size-1)])
          colnames(coef.mat) <- c("Estimate","Std. Error","z-Value","p-Value")
          
          printCoefmat(coef.mat, digits = digits, cs.ind=1:2, tst.ind=3, has.Pvalue=T, signif.legend = F)
          
          printed.var <- printed.var + set.size
        }
      }
      
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n\n",sep="",collapse="")

      cat("PRECISION MODEL:\n",sep="",collapse="")

      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n",sep="",collapse="")

      coef.mat <- cbind(object$coefficients[printed.var:length(object$coefficients)],
                        object$se[printed.var:length(object$coefficients)],
                        z.values[printed.var:length(object$coefficients)],
                        p.values[printed.var:length(object$coefficients)])
      colnames(coef.mat) <- c("Estimate","Std. Error","z-Value","p-Value")
      
      printCoefmat(coef.mat, digits = digits, cs.ind=1:2, tst.ind=3, has.Pvalue=T, signif.legend = F)

      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n",sep="",collapse="")
      cat(paste(strwrap("Signif. codes:  \u60***' < .001, \u60**' < 0.01, \u60*' < 0.05, \u60.' < 0.1",.wd),sep="\n",collapse="\n"),sep="")

    }
    
################################################################################
############################################################################ FIN

    cat("\n\n\nLog-likelihood: ",format(object$logLik,digits=digits)," on ",object$npar," df (",
        object$optimization$bfgs.it,"+",object$optimization$iterations," iterations)\n",sep="",collapse="")
    if(object$parametrization == "common"){
      cat("Link: Log\nParametrization: ", object$parametrization, "\n\n",sep="")
    } else {
      cat("Links: Logit (Means) and Log (Precision)\nParametrization: ", object$parametrization, "\n\n",sep="")
    }
    
}

#summary(rrr)
