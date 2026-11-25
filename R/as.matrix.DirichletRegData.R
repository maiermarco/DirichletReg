as.matrix.DirichletRegData <- function(x, ...){
  
  attr(x, "Y.original" ) <- NULL
  attr(x, "dims"       ) <- NULL
  attr(x, "dim.names"  ) <- NULL
  attr(x, "obs"        ) <- NULL
  attr(x, "valid_obs"  ) <- NULL
  attr(x, "normalized" ) <- NULL
  attr(x, "transformed") <- NULL
  attr(x, "base"       ) <- NULL
  attr(x, "class"      ) <- NULL
  
  x

}
