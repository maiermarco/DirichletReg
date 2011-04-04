getdata <- function(object, orig = FALSE){

  if(orig){
    res <- object$Y.orig
  } else {
    res <- object$Y
  }
  
  return(res)

}