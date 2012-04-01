summary.DirichletRegData <- function(object, ...){

  cat("\nThis object contains compositional data with", attr(object, "dims"), "dimensions.\n")

  cat("Number of observations:", attr(object, "obs"))
  
  cat(" of which", attr(object, "valid_obs"), "(", round(100*attr(object, "valid_obs")/attr(object, "obs"),2), "% ) are valid.\n")
  
  if(attr(object, "normalized") | attr(object, "transformed")){
    cat("\nNote: The data were ")
    if(attr(object, "normalized")) cat("normalized")
    if(attr(object, "normalized") & attr(object, "transformed")) cat(" and ")
    if(attr(object, "transformed")) cat("transformed")
    cat(".\n")
  }
  
  cat("\n")

}
