coord.trafo <- function(Y){
  
  .Deprecated(msg = "coord.trafo() is now deprecated.\nTo plot in the new (qua)ternary coordinate systems, use toSimplex()")
  
  return(toSimplex(Y[, c(2L, 3L, 1L)]))
  
}
