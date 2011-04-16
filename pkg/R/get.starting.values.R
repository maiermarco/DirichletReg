get_starting_values <- function(Y, X.mats, Z.mat, repar, base){

  if(!repar){

    beta.LL <- function(x, y, X){
      b <- matrix(x, ncol=2)
      if(nrow(b) > 2){
        dbeta(y, exp(X%*%b[,1]), exp(X%*%b[,2]), log=TRUE)
      } else {
        dbeta(y, unlist(exp(X*x[1])), unlist(exp(X*x[2])), log=TRUE)
      }
    }
 
    unidim.fit <- sapply(1:ncol(Y), function(i){
      maxBHHH(beta.LL,
              start=rep(0, 2*ncol(X.mats[[i]])),
              X=X.mats[[i]],
              y=Y[,i])$estimate[1:ncol(X.mats[[i]])]
    })

  } else {
  
    Y.logit <- log(Y/(1-Y))
    unidim.fit <- sapply((1:ncol(Y))[-base], function(i){ lm(Y[,i]~X.mats[[i]]-1)$coefficients })
    unidim.fit <- c(unidim.fit, lm(I(rep(.5,nrow(Y)))~Z.mat-1)$coefficients)

  
  }

  return(as.numeric(unlist(unidim.fit)))




  


}
