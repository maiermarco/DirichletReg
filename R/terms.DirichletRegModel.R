terms.DirichletRegModel <- function(x, ...){
  terms(x[["formula"]])
}



model.matrix.DirichletRegModel <- function(object, ...){
  object[c("X", "Z")]
}
