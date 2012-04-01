checkerboard <- function(object, maximize=TRUE){

  if(class(object) != "DirichletRegData") stop("data must be prepared by DR_data() first.")
  
  if(any(is.na(object))) Y <- object[which(rowSums(is.na(object)) == 0),] else Y <- object
  if(nrow(Y) == 0) stop("no observations left.")
  colnames(Y) <- colnames(object)
  
  n <- nrow(Y)
  k <- ncol(Y)
  
  plot(NULL, xlim=c(0, n), ylim=c(-k, 0)-0.5, main="Checkerboard Plot", xlab="Index", ylab="", axes=FALSE)
  
  if(maximize){
    y_norm <- Y - min(Y)
    y_norm <- y_norm / max(y_norm)
  } else {
    y_norm <- Y
  }
  
  for(i in -(1:ncol(Y))){
    rect(0:(n-1), i-0.5, 1:n, i+0.5, col=gray(1-y_norm[,abs(i)]), lty=0)
  }
  
  rect(0, -0.5, n, -k-.5)
  segments(rep(0, k-1), -(2:k)+0.5, rep(n, k-1), -(2:k)+0.5)
  
  axis(2, -(1:k), labels=colnames(Y), tick = FALSE, pos=2, las=1)
  
}
