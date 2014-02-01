DReg.repar <- function(x, Y, X, Z, d, k, w, base, NR){
################################################################################
### PREPARATION ################################################################
################################################################################

  npar <- length(x)
  seq_along_d <- seq_len(d)
  seq_along_k <- seq_len(k)

  B <- matrix(0, nrow=k, ncol=d)
  B[cbind(rep(seq_along_k, (d-1)), rep(seq_along_d[-base], each=k))] <- x[1:((d-1)*k)]
  
  g <- matrix(x[((d-1)*k+1):npar], ncol=1L)

  eps  <- exp(apply(B, 2L, function(b){ X%*%b }))
  esum <- rowSums(eps)
  
  mu <- apply(eps, 2L, function(x){ x / esum })

  f <- exp(Z%*%g)
  
  A <- apply(mu, 2, function(MU){ MU*f })
  


################################################################################
### LOG-LIKELIHOOD #############################################################
################################################################################

  LL <- w*(lgamma(f)-rowSums(lgamma(A))+rowSums((A-1)*log(Y)))
  


################################################################################
### GRADIENT ###################################################################
################################################################################

  gradient <- matrix(NA, nrow=nrow(Y), ncol=npar)

  i <- 0
  for(co in seq_along_d[-base]){
    for(iv in seq_along_k){
      i <- i + 1
      gradient[,i] <- w * X[,iv] * eps[,co] * f * (
                      + rowSums(eps[,-co,drop=F]) * (log(Y[,co]) - psigamma(A[,co]))
                      + rowSums(eps[,-co,drop=F] * (psigamma(A[,-co,drop=F]) - log(Y[,-co,drop=F])))
                      ) / esum^2
    }
  }
  
  for(iv in 1:ncol(Z)){
    i <- i + 1
    gradient[,i] <- w * Z[,iv] * f * (
                    + esum * psigamma(f)
                    + rowSums( eps * (log(Y) - psigamma(A)) )
                    ) / esum
  }
  
  attr(LL, "gradient") <- gradient



################################################################################
### HESSIAN ####################################################################
################################################################################

  if(NR){

  if(d == 2){
    eEps <- list(0, 0)
  } else {
    eEps <- sapply(seq_along_d, function(i){
              rowSums(combn(seq_along_d[-i], 2, function(a){
                eps[,a[1]]*eps[,a[2]]
              }))
            }, simplify=F)
  }

  hessian <- matrix(NA, nrow=npar, ncol=npar)

  hessian.ind <- rbind(as.matrix(expand.grid(seq_along_k, seq_along_d[-base])[,2:1]), cbind(-1, 1:ncol(Z)))
  
  for(hess.j in 1:npar){
    for(hess.i in 1:npar){
      if(hess.i < hess.j) next   # skip symmetric elements
    
      v1 <- hessian.ind[hess.i,2]
      v2 <- hessian.ind[hess.j,2]
      
      derv <- hessian.ind[c(hess.i, hess.j),1]
      
      ##########################################################################
      ############################################### BETAs - SAME RESPONSES ###
      if((derv[1] == derv[2]) & all(derv != -1)) {
        derv <- derv[1]
        nder <- seq_along_d[-derv]

        hessian[hess.i, hess.j] <-
        sum(w*( -X[,v1]*X[,v2]*f*eps[,derv]*(( 
          rowSums(eps[,nder,drop=F]*log(Y[,nder,drop=F])) - log(Y[,derv]) * rowSums(eps[,nder,drop=F])  
        )*(
          2 * eEps[[derv]] + rowSums(eps[,nder,drop=F]^2) - eps[,derv]^2
        )
      + rowSums(sapply(nder, function(i){ eps[,i,drop=F] * (
          eps[,i,drop=F] * eps[,derv,drop=F] * f * psigamma(A[,i,drop=F], 1)
        - psigamma(A[,i,drop=F]) * ( 2 * eEps[[derv]] + rowSums(eps[,nder,drop=F]^2) - eps[,derv,drop=F]^2 )
        ) }))
      + (rowSums(eps[,nder,drop=F]) * (
          rowSums(eps[,nder,drop=F]) * eps[,derv] * f * psigamma(A[,derv], 1) + psigamma(A[,derv]) * (
            2 * eEps[[derv]]
          + rowSums(eps[,nder,drop=F]^2)
          - eps[,derv]^2) 
        ))
        ) / esum^4))
      ##########################################################################
      ########################################## BETAs - DIFFERENT RESPONSES ###
      } else if((derv[1] != derv[2]) & all(derv != -1)) {
        nder <- seq_along_d[-derv]

        hessian[hess.i, hess.j] <-
        sum(w*(X[,v1]*X[,v2]*f*eps[,derv[1]]*eps[,derv[2]]*(
          rowSums(sapply(derv, function(i){
            f*eps[,i]*psigamma(A[,i], 1)*rowSums(eps[,-i])
          }))
        + rowSums(sapply(derv, function(i){
            psigamma(A[,i])*(2*eEps[[i]]+rowSums(eps[,-i]^2)-eps[,i]^2)
          }))
        + esum * (
            rowSums(sapply(nder, function(i){
              -f*mu[,i]*eps[,i]*psigamma(A[,i],1) - 2*eps[,i]*psigamma(A[,i]) + 2*log(Y[,i])*eps[,i]
            }))
          + rowSums(sapply(derv, function(i){ log(Y[,i])*(eps[,i]-rowSums(eps[,-i])) }))
          )    
        ) / esum^4))
      ##########################################################################
      ######################################################### BETA / GAMMA ###
      } else if(any(derv != -1) & any(derv == -1)) {
        derv <- derv[which(derv != -1)]
        nder <- seq_along_d[-derv]
        
        hessian[hess.i, hess.j] <-
        sum(w*(Z[,v1] * X[,v2] * f * eps[,derv] * (
          esum * (rowSums(sapply(nder, function(i){
            eps[,i] * (psigamma(A[,i]) - log(Y[,i]))
            }) ))
        + (log(Y[,derv]) - psigamma(A[,derv])) * (
            rowSums(eps[,-derv,drop=F]^2)
          + eps[,derv] * rowSums(eps[,-derv,drop=F])
          + 2 * eEps[[derv]]
          )
        + f * ( rowSums(sapply(nder, function(i){
            eps[,i]^2*psigamma(A[,i], 1)
            })) - eps[,derv] * psigamma(A[,derv], 1) * rowSums(eps[,nder,drop=F])  
          )
        ) / esum^3))
      ##########################################################################
      ############################################################### GAMMAs ###
      } else if(all(derv == -1)){
        hessian[hess.i, hess.j] <-
        sum(w*(Z[,v1]*Z[,v2]*f*(
          rowSums(sapply(seq_along_d, function(i){
            eps[,i] * ( log(Y[,i]) - psigamma(A[,i]) - A[,i]*psigamma(A[,i],1) )
          }))
        + esum * (psigamma(f) + f*psigamma(f,1)) 
        )/esum))
      }
    }
  }

  hessian <- make.symmetric(hessian)

  attr(LL, "hessian") <- hessian

  }
  
  return(LL)
  
}
