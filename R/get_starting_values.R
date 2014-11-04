get_starting_values <- function(Y, X.mats, Z.mat, repar, base, weights){

  if(!repar){

    
    exclude_par <- lapply(X.mats, function(list_el){
      lin_coef <- lm.fit(x = list_el, y = runif(nrow(list_el)))[["coefficients"]]
      if(any(na_pos <- is.na(lin_coef))){
        temp_names <- names(lin_coef)
        lin_coef <- seq_along(lin_coef)
        names(lin_coef) <- temp_names
        return(lin_coef[na_pos])  
      } else {
        return(NULL)
      }
    })
    

    beta.LL <- function(x, y, X, w){
      b <- matrix(x, ncol = 2L)
      if(ncol(X) > 1L){
        LL <- w * dbeta(y, exp(X%*%b[,1]), exp(X%*%b[,2]), log=TRUE)
      } else {
        LL <- w * dbeta(y, unlist(exp(X*x[1])), unlist(exp(X*x[2])), log=TRUE)
      }
      return(LL)
    }
    
    beta.LL.deriv <- function(x, y, X, w){
      b <- matrix(x, ncol=2L)
      grad <- matrix(0.0, nrow=nrow(X), ncol=prod(dim(b)))
      element <- 1L
      if(ncol(X) > 1L){
        for(cc in seq_len(ncol(b)))for(rr in seq_len(nrow(b))){
          grad[,element] <- w * X[,rr]*(psigamma(exp(X%*%b[,1L]+X%*%b[,2L]))-psigamma(exp(X%*%b[,cc]))+log(y))
          element <- element + 1L
        }
      } else {
        for(cc in seq_len(ncol(b))){
          grad[,element] <- w * X*(psigamma(exp(X%*%x[1L]+X%*%x[2L]))-psigamma(exp(X%*%x[cc]))+log(y))
          element <- element + 1L
        }
      }
      return(grad)
    }

    unidim.fit <- lapply(seq_len(ncol(Y)), function(i){
      if(is.null(exclude_par[[i]])){
        correctX <- X.mats[[i]]
      } else {
        correctX <- X.mats[[i]][ , -exclude_par[[i]], drop = FALSE]
      }
      suppressWarnings(
        maxBFGS(beta.LL, beta.LL.deriv,
          start        = rep(0, 2*ncol(correctX)),
          tol          = 1e-05,
          finalHessian = FALSE,
          X            = correctX,
          y            = Y[,i],
          w            = weights)$estimate[seq_len(ncol(correctX))]
      )
    })

    for(cmp in seq_len(ncol(Y))){
      if(is.null(exclude_par[[cmp]])){
        break
      } else {
        for(NA_vars in seq_along(exclude_par[[cmp]])) unidim.fit[[cmp]] <- append(unidim.fit[[cmp]], NA, exclude_par[[cmp]][NA_vars] - 1L)
      }
    }
    
    unidim.fit <- unlist(unidim.fit)

  } else {
  
    Y.logit <- log(Y/(1-Y))
    unidim.fit <- sapply((1:ncol(Y))[-base], function(i){ lm(Y[,i]~X.mats[[i]]-1, weights=weights)$coefficients })
    unidim.fit <- c(unidim.fit, lm(I(rep(.5,nrow(Y)))~Z.mat-1, weights=weights)$coefficients)
  
  }

  return(as.numeric(unlist(unidim.fit)))

}
