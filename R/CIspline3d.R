CIspline3d <- function(Cints, E, ...){

  if(length(Cints) != 6) stop() # these can only be plotted for alpha > 1
  if(length(E) != 3) stop() # can only be plotted for 3d data
  
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

# plot(0,0,xlim=c(-4,4), ylim=c(-2,1), asp=1)
# for(i in 1:6) points(Cints[i]*cos(-pi/2+(i-1)*pi/3),Cints[i]*sin(pi/2+(i-1)*pi/3), pch=as.character(i))
# lines()
