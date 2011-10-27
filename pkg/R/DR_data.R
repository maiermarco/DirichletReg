DR_data <- function(data,            # the data to be 
                    trafo = FALSE,   # transform (compress) the data?
                    base  = 1        # base variable for the reparametrized model
                   ){

  state.norm <- FALSE   # was normalization forced?
  state.tran <- FALSE   # was transformation forced?


### CONVENIENTLY HANDLE BETA-DISTRIBUTED VARIABLES
  if((!is.matrix(data) & !is.data.frame(data)) | ifelse(is.null(ncol(data)), FALSE, ncol(data) == 1)){
    if(any((data < 0) | (data > 1))) stop("only one variable supplied with values outside [0, 1]. beta distribution cannot safely be assumed. prepare your data first.")
    data <- cbind(1-data, data)
                       
    .names <- colnames(data)
    colnames(data) <- c(paste("1 -",.names[2]), .names[2])

    message("only one variable supplied, beta-distribution assumed")
  }

  

  if(!is.matrix(data) & !is.data.frame(data)) stop('"data" must be either a matrix or a data.frame')
  if(ncol(data) <= 1) stop('"data" must at least have two columns')
  
  if((base < 1) | (base > ncol(data))) stop('"base" must be in the range of variables')
  
  if(any(is.na(data))){
    exclude <- as.logical(apply(is.na(data), 1, any))
    warning("missing values are not allowed. the following observation(s) will be excluded:\n",
            paste(which(exclude), collapse=", "))
  } else {
    exclude <- NULL
  }
  
  
  
  row.sums <- na.exclude(rowSums(data))
  
  if( !isTRUE( all.equal(row.sums, rep(1, length(row.sums)) ) ) ){
    data <- data/rowSums(data)
    state.norm <- TRUE
    warning("not all rows sum up to 1 => normalization forced")
  }
  
  # save the original data for reference
  data.original <- data
  
  if(trafo | any(data == 0, na.rm=TRUE) | any(data == 1, na.rm=TRUE)){
    n.obs <- ifelse(is.null(exclude), nrow(data), sum(!exclude))
    data <- (data * (n.obs - 1) + 1/ncol(data)) / n.obs
    state.tran <- TRUE
    warning("some entries are 0 or 1 => transformation forced")
  }
  
  if(is.null(colnames(data))){
    colnames(data) <- paste("v", 1:ncol(data), sep="")
  }
  
  res <- list(Y = data,
              Y.orig = data.original,
              dims = ncol(data),
              dim.names = colnames(data),
              obs = nrow(data),
              exclude = exclude,
              normalized = state.norm,
              transformed = state.tran,
              base = base)
  
  class(res) <- "DirichletRegData"
  
  return(res)
  
}
