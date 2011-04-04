cmyk2rgb <- function(cmyk){
  R <- (cmyk[,1] - 1) * (cmyk[,4] - 1)
  G <- (cmyk[,2] - 1) * (cmyk[,4] - 1)
  B <- (cmyk[,3] - 1) * (cmyk[,4] - 1)
  
  return(rgb(R,G,B))
}
