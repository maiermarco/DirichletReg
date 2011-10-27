CIspline3d <- function(Cints, E, ...){

  if(length(Cints) != 6) stop() 
  if(length(E) != 3) stop() 
  
  aa <- list(...)
  .n <- ifelse(!is.null(aa[["n"]]), aa[["n"]], 200)
  .method <- ifelse(!is.null(aa[["method"]]), aa[["method"]], "periodic")

  y <- c(Cints, Cints[1])
  x <- seq(0,2*pi,length.out=length(y))
    
  radia <- as.data.frame(spline(x, y, n=.n, method=.method))
  
  xy <- cbind(cos(-radia[,1]+pi/2)*radia[,2],
              sin( radia[,1]+pi/2)*radia[,2]) + coord.trafo(sapply(E, rep, times=nrow(radia)))
  return(xy)
  
}




