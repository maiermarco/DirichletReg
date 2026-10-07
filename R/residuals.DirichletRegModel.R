residuals.DirichletRegModel <- function(
  object,
  type = c(
    "standardized",
#   "score",
    "composite",
    "raw"
  ),
  ...
){

  Y           <- object$Y
  fitted.vals <- fitted(object, mu = TRUE, alpha = TRUE, phi = TRUE)
  M           <- fitted.vals[["mu"]]
  A           <- fitted.vals[["alpha"]]
  f           <- fitted.vals[["phi"]]

  raw.res <- Y - M

  type <- match.arg(type)

  wghts <- object$weights

  V <- apply(M, 2L, function(m){ (m * (1.0 - m)) / (1.0 + f) })

  switch(type,

    "standardized" = {
      res <- raw.res / sqrt(V)
    },

#    "score" = {
#      res <- digamma(colSums(A)) - digamma(A) + log(Y)
#    },

    "composite" = {
      res <- rowSums((raw.res / sqrt(V))^2.0)
    },

    "raw" = {
      res <- raw.res
    }


  )

  return(res)

}
