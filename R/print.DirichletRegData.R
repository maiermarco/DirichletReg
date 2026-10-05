print.DirichletRegData <- function(x, type=c("processed", "original"), ...){
  
  print_data <- if(match.arg(type) == "processed") x else attr(x, "Y.original")
  
  print(print_data[,]) # use the data.frame method to avoid recursion
  
  invisible(x)
  
}
