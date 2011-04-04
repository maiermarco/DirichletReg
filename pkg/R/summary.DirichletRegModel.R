summary.DirichletRegModel <- function(object, digits=max(3, getOption("digits") - 3), ...) 
{
    .wd <- getOption("width")

    if(object$optimization$convergence > 2) stop("\n",paste(strwrap(paste("\nOptimization did not converge in",object$optimization$counts,"iterations and exited with code",object$optimization$convergence),.wd),sep="\n",collapse="\n"))

    cat("\nCall:\n",
        paste(strwrap(deparse(object$call), .wd), sep="\n", collapse="\n"),
        "\n\n", sep = "")

    cat("RESIDUALS!\n\n")
    
    coef.ind <- cumsum(object$n.vars)
    
    z.values <- object$coefficients / object$se
    p.values <- 2 * pnorm(-abs(z.values))
    coef.mat <- cbind(object$coefficients, object$se, z.values, p.values)
    colnames(coef.mat) <- c("Estimate","Std. Error","z-Value","p-Value")
    
    if(object$parameterization == "common"){
    
      for(i in 1:length(object$varnames)){
        cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
        cat("\nCoefficients for variable no. ",i,": ",object$varnames[i],"\n",sep="",collapse="")


        printCoefmat(coef.mat[ ifelse(i==1,1,coef.ind[i-1]+1):coef.ind[i] , , drop=F],
                     digits = digits, cs.ind=1:2, tst.ind=3, has.Pvalue=T, signif.legend = F)


      }
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n",sep="",collapse="")
      cat(paste(strwrap("Signif. codes:  `***' < .001, `**' < 0.01, `*' < 0.05, `.' < 0.1",.wd),sep="\n",collapse="\n"),sep="")

    } else {

      printed.var <- 1
      set.size    <- ncol(object$X[[1]])

      cat("MEAN MODELS:\n",sep="",collapse="")

      for(i in 1:length(object$varnames)){
        if(i == object$orig.resp$base){
          cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
          cat("\nCoefficients for variable no. ",i,": ",object$varnames[i],"\n",sep="",collapse="")
          cat("- variable omitted (reference category) -\n")
        } else {
          cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
          cat("\nCoefficients for variable no. ",i,": ",object$varnames[i],"\n",sep="",collapse="")
          print.default(format(object$coefficients[printed.var:(printed.var+set.size-1)], digits=digits), print.gap=2, quote=F)
          
          printed.var <- printed.var + set.size
        }
      }
      
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n\n",sep="",collapse="")

      cat("PRECISION MODEL:\n",sep="",collapse="")

      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n",sep="",collapse="")
      print.default(format(object$coefficients[printed.var:length(object$coefficients)], digits=digits), print.gap=2, quote=F)
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n\n",sep="",collapse="")

    }

    cat("\n\nLog-likelihood: ",format(object$logLik,digits=digits)," on ",object$npar," df (",
        object$optimization$counts," iterations)\n",sep="",collapse="")
    cat("Link: Log\nParameterization: ", object$parameterization, "\n\n",sep="")
}


