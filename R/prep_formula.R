prep_formula <- function(f, d){

  fs <- strsplit(paste(deparse(f), collapse=""), "|", fixed=T)[[1]]
  fs <- sapply(fs, blank.trim, USE.NAMES=F)
  fl <- length(fs)

  response <- get(blank.trim(strsplit(fs[1],"~",fixed=TRUE)[[1]][1]))
  if(class(response) != "DirichletRegData") stop("the response variable must be prepared by \"DR_data()\"")
  ndim <- response$dims
  
  Y <- response$Y

  if(is.null(d)){
    d <- data.frame(rep(1,nrow(Y)),row.names=1:nrow(Y))
    names(d) <- "(Intercept)"
  }
  
  # backup rownames of Y and the data
  Y.labels <- rownames(Y)
  rownames(Y) <- 1:nrow(Y)
  
  d.labels <- rownames(d)
  rownames(d) <- 1:nrow(d)

  # exclude observations from the data if they have been excluded from the response matrix
#  if(!is.null(response$exclude)){
#    d <- d[which(!response$exclude),]
#    Y <- Y[which(!response$exclude),]
#  }
  
  
  
  ### extract formulae

  predictors <- list(RHS=NULL, PHI=NULL)
  
  if(fl == 1){ # if the formula has one element
  
    fcomm <- strsplit(fs, "~", fixed=TRUE)[[1]]
    if(length(fcomm) != 2) stop("error in formula")
    predictors$RHS <- sapply(1:ndim, function(x) blank.trim(fcomm[2]), simplify=F)
  
  } else if (fl == 2) { # if the formula has two elements
  
    fcomm <- strsplit(fs, "~", fixed=TRUE)

    if((length(fcomm[[1]]) == 2) & (length(fcomm[[2]]) == 1)){

      predictors$RHS <- list(blank.trim(fcomm[[1]][2]),
                             blank.trim(fcomm[[2]]))

    } else if((length(fcomm[[1]]) == 2) & (length(fcomm[[2]]) == 2)){

      MU  <- fcomm[[1]]
      PHI <- fcomm[[2]]
      if(toupper(blank.trim(PHI[1])) != "PHI") stop("error in formula")
      
      predictors$RHS <- sapply(1:ndim, function(x) blank.trim(MU[2]), simplify=F)
      predictors$PHI <- blank.trim(PHI[2])

    } else {

      stop("error in formula")

    }
  
  } else if (fl == ndim) { # if the formula has ndim elements
  
    fcomm <- strsplit(fs, "~", fixed=TRUE)
    if(!all(sapply(2:ndim, function(i){ length(fcomm[[i]]) == 1 }))) stop("error in formula")
    
    predictors$RHS <- sapply(1:ndim, function(i){
                        if(i==1){
                          blank.trim(fcomm[[i]][2])
                        } else {
                          blank.trim(fcomm[[i]])
                        }
                      }, simplify=F)
  
  } else {
  
    stop("error in formula")
  
  }
  
  ### extract terms and matrices

  if(!is.null(predictors$RHS) & is.null(predictors$PHI)){

    repar <- FALSE

    X.matrices <- lapply(predictors$RHS, function(fstr){
                    model.matrix(as.formula(paste("~",fstr)),
                      data=model.frame(as.formula(paste("~",fstr)),
                                       data=d,
                                       drop.unused.levels=TRUE)
                    )
                  })
                  
    Z.matrix <- NULL
                  
  } else if(!is.null(predictors$RHS) & !is.null(predictors$PHI)){

    repar <- TRUE

    X.matrices <- lapply(predictors$RHS, function(fstr){
                    model.matrix(as.formula(paste("~",fstr)),
                      data=model.frame(as.formula(paste("~",fstr)),
                                       data=d,
                                       drop.unused.levels=TRUE)
                    )
                  })

    Z.matrix <- model.matrix(as.formula(paste("~",predictors$PHI)),
                      data=model.frame(as.formula(paste("~",predictors$PHI)),
                                       data=d,
                                       drop.unused.levels=TRUE))

  } else {

    stop("error in formula")

  }

  # check for "intercept-only" terms
  X_1 <- unlist(lapply(predictors$RHS, `==`, "1"))
  if(any(X_1)) for(i in which(X_1)){
    X.matrices[[i]] <- data.frame(rep(1,nrow(Y)), row.names=1:nrow(Y))
    names(X.matrices[[i]]) <- "(Intercept)"
  }
  if(!is.null(predictors$PHI)) Z_1 <- predictors$PHI == "1" else Z_1 <- FALSE
  if(Z_1){
    Z.matrix <- data.frame(rep(1,nrow(Y)), row.names=1:nrow(Y))
    names(Z.matrix) <- "(Intercept)"
  }



  ### drop invalid observations from Y
  
  valid.obs <- sapply(X.matrices, rownames, simplify=FALSE)
  if(exists("Z.matrix")) valid.obs[[length(valid.obs)+1]] <- rownames(Z.matrix)
  valid.obs <- lapply(valid.obs, as.numeric)

  vobs <- if(is.null(response$exclude)) valid.obs[[1]] else which(!response$exclude)
  for(i in 1:length(valid.obs)){
    vobs <- intersect(valid.obs[[i]], vobs)
  }
  
  Y <- Y[vobs,]
  rownames(Y) <- Y.labels[vobs]
  
  X.matrices <- lapply(X.matrices, function(L){ L[which(as.numeric(rownames(L)) %in% vobs),,drop=FALSE] })
  for(i in 1:length(X.matrices)) rownames(X.matrices[[i]]) <- d.labels[vobs]
  if(!is.null(Z.matrix)){
    Z.matrix <- Z.matrix[which(as.numeric(rownames(Z.matrix)) %in% vobs),,drop=FALSE]
    rownames(Z.matrix) <- d.labels[vobs]
  }
  


  ### wrap up and return
  
  nvars <- sapply(X.matrices, ncol, simplify=T)

  return(list(response = response,
              Y        = Y,
              n.dim    = ndim,
              n.vars   = nvars,
              X.mats   = X.matrices,
              Z.mat    = Z.matrix,
              preds    = predictors,
              repar    = repar,
              d        = d,
              base     = response$base))

}
