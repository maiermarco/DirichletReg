coef.DirichletRegModel <- function(object, type=c("both", "beta", "gamma"), ...){

  type <- match.arg(type, c("both", "beta", "gamma"))

  cc <- object$coefficients

  if(object$parametrization == "common"){

    ind <- cumsum(c(1, object$n.vars))
    cc <- sapply(seq_along(ind)[-1], function(i) cc[ind[i-1]:(ind[i]-1)], simplify=FALSE)
    return(cc)

  } else {

    ind <- cumsum(c(1, object$n.vars))
    bb <- list()
    
    bb_par <- cc[-((length(cc)-ncol(object$Z)+1):length(cc))]
    iter <- 0
    for(i in 1:object$orig.resp$dims){
      bb[[i]] <- if(i == object$base) NULL else {
        iter <- iter + 1
        bb_par[ind[iter]:(ind[iter+1]-1)]
      }
    }
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
