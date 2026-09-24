summary.DirichletRegData <- function(object, ...){
  
  if(interactive()) message()
  
  message(sprintf("This object contains compositional data with %d dimensions.", attr(object, "dims")))
  
  message(sprintf(
    "Number of observations: %d of which %d (%0.2f%%) are valid.",
    attr(object, "obs"), attr(object, "valid_obs"), 100 * attr(object, "valid_obs") / attr(object, "obs")
  ))

  message()
  
  if((is_normalized <- attr(object, "normalized")) | (is_transformed <- attr(object, "transformed"))){
    message(sprintf("Note: The data were %s.",
      paste(c("normalized", "transformed")[c(is_normalized, is_transformed)], collapse = " and ")
    ))
  }
  
  if(interactive()) message()
  
}
