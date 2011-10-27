print.DirichletRegModel <- function(x, digits=max(3, getOption("digits") - 3), ...) 
{
    .wd <- getOption("width")

    if(x$optimization$convergence == 3) cat("\n",strwrap("CAUTION! Possible convergence problems!",.wd),"\n",sep="")
    if(x$optimization$convergence > 3) stop("\n",paste(strwrap(paste("\nOptimization did not converge in",x$optimization$bfgs.it,"+",x$optimization$iterations,"iterations and exited with code",x$optimization$convergence),.wd),sep="\n",collapse="\n"))

    cat("\nCall:\n",
        paste(strwrap(deparse(x$call), .wd), sep="\n", collapse="\n"),
        "\nusing the ", x$parametrization, " parametrization\n\n",
        sep = "")
    
    cat("Log-likelihood: ",format(x$logLik,digits=digits)," on ",x$npar," df (",
        x$optimization$bfgs.it,"+",x$optimization$iterations," iterations)\n\n",sep="",collapse="")

    coef.ind <- cumsum(x$n.vars)
    
    if(x$parametrization == "common"){

      for(i in 1:length(x$varnames)){
        cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
        cat("\nCoefficients for variable no. ",i,": ",x$varnames[i],"\n",sep="",collapse="")
        print.default(format(x$coefficients[ifelse(i==1,1,coef.ind[i-1]+1):coef.ind[i]], digits=digits), print.gap=2, quote=F)
      }
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n\n",sep="",collapse="")

    } else {

      printed.var <- 1
      set.size    <- ncol(x$X[[1]])

      cat("MEAN MODELS:\n",sep="",collapse="")

      for(i in 1:length(x$varnames)){
        if(i == x$orig.resp$base){
          cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
          cat("\nCoefficients for variable no. ",i,": ",x$varnames[i],"\n",sep="",collapse="")
          cat("- variable omitted (reference category) -\n")
        } else {
          cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),sep="",collapse="")
          cat("\nCoefficients for variable no. ",i,": ",x$varnames[i],"\n",sep="",collapse="")
          print.default(format(x$coefficients[printed.var:(printed.var+set.size-1)], digits=digits), print.gap=2, quote=F)
          
          printed.var <- printed.var + set.size
        }
      }
      
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n\n",sep="",collapse="")

      cat("PRECISION MODEL:\n",sep="",collapse="")

      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n",sep="",collapse="")
      print.default(format(x$coefficients[printed.var:length(x$coefficients)], digits=digits), print.gap=2, quote=F)
      cat(paste(rep("-", min(80, .wd)),sep="",collapse=""),"\n\n",sep="",collapse="")

    }
}


#print.DirichletRegModel(res2)

#res2 <- DirichReg(DR_data(Rocks[,1:5])~type, Rocks)
