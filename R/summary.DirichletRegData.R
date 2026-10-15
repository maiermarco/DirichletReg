summary.DirichletRegData <- function(object, ...){
  
  if(interactive()) cat("\n")
  
  cat(gettextf("This object contains compositional data with %d dimensions.\n", attr(object, "dims")))
  
  cat(gettextf(
    "Number of observations: %d of which %d (%0.2f%%) are valid.\n",
    attr(object, "obs"),
    attr(object, "valid_obs"),
    100.0 * attr(object, "valid_obs") / attr(object, "obs")
  ))
  
  is_normalized  <- attr(object, "normalized")
  is_transformed <- attr(object, "transformed")
  
  if(is_normalized || is_transformed){
    cat(gettextf("\nNote: The data were %s.\n",
      paste(
        c(gettext("normalized"), gettext("transformed"))[c(is_normalized, is_transformed)],
        collapse = sprintf(" %s ", gettext("and"))
      )
    ))
  }
  
  if(interactive()) cat("\n")
  
  invisible(object)  
}
