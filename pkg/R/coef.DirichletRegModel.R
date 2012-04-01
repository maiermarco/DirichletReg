coef.DirichletRegModel <- function(object, type=c("both", "beta", "gamma"), ...){

  type <- match.arg(type, c("both", "beta", "gamma"))

  cc <- object$coefficients
  names(cc) <- object$coefnames

  if(object$parametrization == "common"){

    ind <- cumsum(c(1, object$n.vars))
    cc <- sapply(seq_along(ind)[-1], function(i) cc[ind[i-1]:(ind[i]-1)], simplify=FALSE)
    names(cc) <- object$varnames
    return(cc)

  } else {

    ind <- cumsum(c(1, object$n.vars))
    bb <- list()
    
    bb_par <- cc[-((length(cc)-ncol(object$Z)+1):length(cc))]
    iter <- 0
    for(i in 1:object$dims){
      bb[[i]] <- if(i == object$base) NULL else {
        iter <- iter + 1
        bb_par[ind[iter]:(ind[iter+1]-1)]
      }
    }
    names(bb) <- c(object$varnames)
    gg <- list(cc[ ((length(cc)-ncol(object$Z)+1):length(cc))])

    if(type == "both"){
      return(list("beta"=bb, "gamma"=gg))
    } else if(type == "beta"){
      return(bb)
    } else if(type == "gamma"){
      return(gg)
    }

  }

}
