blank.trim <- function(x){ paste(unlist(strsplit(x, "^\ +|\ +$")),sep="",collapse="") }



make.symmetric <- function(x){
  if(nrow(x) != ncol(x)) stop("x must be a square matrix")

  cell.ind <- which(is.na(x), arr.ind=TRUE)
  
  x[cell.ind] <- x[cell.ind[,2:1]]
  
  return(x)
}



inv.logit <- function(x){
  log(x) - log(1-x)
}



rdirichlet <- function(n,      # a single integer specifying the sample size
                       alpha   # alpha which can either be a single vector or a matrix (rows = n)
                      ){###rdirichlet

  # check if the sample size is an interger > 0
  if( ((n %% 1) != 0) | (n <= 0)) stop("n must be an integer > 0")
  # check if any value in alpha is <= 0
  if( any(alpha <= 0) ) stop("all values in alpha must be > 0")

  .vec <- is.vector(alpha)
  .mat <- is.matrix(alpha)
  
  if(!.vec & !.mat){ # alpha is neither a vector nor a matrix

    stop("alpha must be a vector or a matrix")

  } else if(.vec & !.mat){ # alpha is a vector

    dims <- length(alpha)
    G <- matrix(rgamma(n*dims, alpha, 1), byrow=TRUE, ncol=dims)
    X <- G / rowSums(G)

  } else { # alpha is a matrix

    dims <- ncol(alpha)
    if(n != nrow(alpha)) stop("when alpha is a matrix, the number of its rows must be equal to n")
  
    G <- matrix(rgamma(n*dims, as.vector(alpha), 1), nrow=n)
    X <- G / rowSums(G)
  
  }
  
  return(X)
  

}###END OF rdirichlet



ddirichlet <- function(x, alpha, log=FALSE, sum.up=FALSE){
  # some checking!
  if(is.null(dim(x))) stop("x must be a matrix")
  if(is.vector(alpha)){
    if(ncol(x) != length(alpha)) stop("alpha must be a vector/matrix fitting to the data in x")
    alpha <- matrix(rep(alpha,nrow(x)),nrow(x),byrow=T)
  }
  if(any(dim(alpha) != dim(x))) stop("check if x and alpha are correctly specified")
                                                                             
  res <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) + rowSums((alpha-1)*log(x))

  if(sum.up) res <- sum(res)

  if(log){ return(res) } else { return(exp(res)) }
}
