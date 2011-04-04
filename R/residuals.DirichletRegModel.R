residuals.DirichletRegModel <- function(object, type=c("score",
                                                       "raw",
                                                       "standardized",
                                                       "composite"), ...){

  Y  <- object$Y$Y
  MU <- fitted(object)
  raw.res <- Y - MU

  type <- match.arg(type)
  
  if(type == "raw") return(raw.res)
  
  wghts <- object$weights
  
  A  <- fitted(object, alphas=TRUE)
  A0 <- rowSums(A)
  
  varY <- (MU*(1-MU)) / (1+A0)
  
  switch(type,
  
    "score" = {
      res <- digamma(colSums(A)) - digamma(A) + log(Y)
    },

    "standardized" = {          
      res <- (Y-MU)/sqrt(varY)
    },

    "composite" = {
      res <- rowSums(((Y-MU)/sqrt(varY))^2)
    }

  )
  
  return(res)
}
