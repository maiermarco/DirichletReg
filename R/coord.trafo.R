coord.trafo <- function(Y){
  
  .Deprecated(new = "toSimplex(x)", old = "coord.trafo(Y)")
  
  toSimplex(Y[, c(2L, 3L, 1L)])
  
}
