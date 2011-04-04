print.DirichletRegData <- function(x, ...){

  cat("\nThis object contains compositional data with",x$dims,"dimensions.\n")

  cat("Number of observations:", x$obs,"\n")
  
  if(!is.null(x$exclude)){
    cat("Valid obs:", sum(!x$exclude),"- excluded obs:", sum(x$exclude),"\n")
  }

  if(x$normalized | x$transformed){
    cat("* The data were ")
    if(x$normalized) cat("normalized")
    if(x$normalized & x$transformed) cat(" and ")
    if(x$transformed) cat("transformed")
    cat(".\n")
  }

  cat("\nTo access the data, use the function getdata().\n\n")

}
