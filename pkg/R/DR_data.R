DR_data <- function(Y,               
                    trafo = FALSE,   
                    base  = 1        
                   ){



  state.norm <- FALSE   
  state.tran <- FALSE   


  Y.original <- Y
  

  

  if((!is.matrix(Y) & !is.data.frame(Y)) | ifelse(is.null(ncol(Y)), FALSE, ncol(Y) == 1)){
    if(any((Y < 0) | (Y > 1))) stop("only one variable supplied with values outside [0, 1]. beta distribution cannot safely be assumed. prepare your data first.")
    Y <- cbind(1-Y, Y)
    
    .name <- deparse(match.call()$Y)
    .name <- gsub(".*\\$", "", .name)   
    .name <- gsub("\\[.*", "", .name)   
    if(length(.name) == 0) .name <- "Y" 
    colnames(Y) <- c(paste("1 -",.name), .name)

    message("only one variable in [0, 1] supplied - beta-distribution assumed. check this assumption.")
  }


  

  if(!is.matrix(Y) & !is.data.frame(Y)) stop('"Y" must be either a matrix or a data.frame')
  if(ncol(Y) <= 1) stop('"Y" must at least have two columns')
  if((base < 1) | (base > ncol(Y))) stop('"base" must be in the range of variables')
  
  if(is.null(colnames(Y))) colnames(Y) <- paste("v", 1:ncol(Y), sep="")
  







  
  


  row.sums <- rowSums(Y)
  
  if( !all(na.delete(row.sums) == rep(1, length(na.delete(row.sums))) ) ){
    Y <- Y/row.sums
    state.norm <- TRUE
  }
  
  


  if(trafo | any(Y == 0, na.rm=TRUE) | any(Y == 1, na.rm=TRUE)){
    n.obs <- length(na.delete(row.sums))
    Y <- (Y * (n.obs - 1) + 1/ncol(Y)) / n.obs
    if(!trafo) state.tran <- TRUE
  }
  


  res <- as.matrix(Y)
  attr(res, "Y.original") <- as.data.frame(Y.original)
  attr(res, "dims") <- ncol(Y)
  attr(res, "dim.names") <- colnames(Y)
  attr(res, "obs") <- nrow(Y)
  attr(res, "valid_obs") <- length(na.delete(row.sums))
  attr(res, "normalized") <- state.norm
  attr(res, "transformed") <- state.tran
  attr(res, "base") <- base
  class(res) <- "DirichletRegData"

  
  
  if(state.norm & state.tran){
    warning("not all rows sum up to 1 => normalization forced\nsome entries are 0 or 1 => transformation forced", immediate.=TRUE)
  } else if(state.norm){
    warning("not all rows sum up to 1 => normalization forced", immediate.=TRUE)
  } else if(state.tran){
    warning("some entries are 0 or 1 => transformation forced", immediate.=TRUE)
  }


  
  return(res)
  
}
