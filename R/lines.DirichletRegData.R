lines.DirichletRegData <- function(x, indep, groups, ..., cumul=FALSE, orig.scale=FALSE){

  if(class(x) != "DirichletRegData") stop("x must be of class 'DirichletRegData'")

  if(missing(groups)) groups <- 1

  if(is.list(indep)){
    indep   <- indep[[1]]
    indep.f <- indep[[2]]
  }

  if( (class(groups) == "formula") | (class(indep) == "formula") ){
    if( (class(groups) == "formula") & (class(indep) == "formula") ){

    } else if(class(groups) == "formula"){
      groups <- model.matrix(groups)
    } else if(class(indep) == "formula"){
      indep <- model.matrix(indep)
    }
  }




  plot(NULL, xlim=range(indep), ylim=0:1)
  for(i in 1:x$dims) points(indep, x$Y[,i], col=i, pch=16)


}



