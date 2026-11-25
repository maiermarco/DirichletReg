terms.DirichletRegModel_terms <- function(x, type = c("full", "separate"), ...){
  
  type <- match.arg(type)
  
  if(type == "full"){
    x[["full"]]
  } else if(type == "separate"){
    if(attr(x, "parametrization") == "common"){
      structure(x[["X"]], names = x$varnames)
    } else {
      list(
        mean      = x[["X"]][[1L]],
        precision = x[["Z"]]
      )
    }
  }

}

terms.DirichletRegModel <- function(x, type = c("full", "separate"), ...){
  
  terms.DirichletRegModel_terms(x[["terms"]], type = type)

}

terms.summary_DirichletRegModel <- function(x, type = c("full", "separate"), ...){
  
  terms.DirichletRegModel_terms(x[["terms"]], type = type)

}





model.matrix.DirichletRegModel <- function(object, ...){
  
  object[c("X", "Z")]

}
