get_or_else <- function(what, or_else, where) if(any(what %in% names(where))) where[[what]] else or_else



plot.DirichletRegData <- function(x,
                                  dims,   
                                  explore,   
                                  model,   
                                  c.grid=TRUE,   
                                  ticks=TRUE,   
                                  colored=TRUE,   
                                  col.scheme=c("dims", "entropy"),   
                                  entropy.contours=FALSE,   
                                  entropy.colors=FALSE,   
                                  dim.labels,
                                  args.3d=list(rgl=TRUE, theta=NULL, phi=NULL, ref.lines=NULL, ...),  
                                  rug=T,
                                  ...){

  if(missing(dims))    dims    <- NULL
  if(missing(explore)) explore <- NULL
  if(missing(model))   model   <- NULL
  
  if(class(x) != "DirichletRegData") stop("data must be prepared by 'DR.data()'")

  col.scheme <- match.arg(col.scheme)

  dotlist <- list(...)
  
  .main <- get_or_else("main", NULL, dotlist)
  .xlim <- get_or_else("xlim", c(NULL, NULL), dotlist)
  .ylim <- get_or_else("ylim", c(NULL, NULL), dotlist)
   .col <- get_or_else("col", NULL, dotlist)
   .pch <- get_or_else("pch", 16, dotlist)
   .cex <- get_or_else("cex", 1, dotlist)
   .lwd <- get_or_else("lwd", 1, dotlist)
   .lty <- get_or_else("lty", 1, dotlist)

  .marginal <- FALSE
  
  if(is.null(dims)){
    if(x$dims > 4){
      x$Y         <- x$Y[, 1:4]
      x$dims      <- 4
      x$dim.names <- x$dim.names[1:4]
      warning("data contains > 4 variables. the first four are being used. to change this, set the 'dims' argument appropriately.")
      .marginal <- TRUE
    } else {
      dims <- 1:x$dims
    }
  } else {
    if((length(dims) < 2) | (length(dims) > 4)) stop("the argument 'dims' must have 2, 3 or 4 elements")
    x$Y         <- x$Y[, dims]
    x$dims      <- length(dims)
    x$dim.names <- x$dim.names[dims]
    .marginal <- TRUE
  }
  

 
  if(missing(dim.labels)){
    dim.labels <- x$dim.names
  }


  
  if(!is.null(explore) & !is.null(model)) stop("model and explore cannot be specified at the same time.")
  
  if(!is.null(explore)){
    frmla <- explore$formula
    dta   <- explore$data

    ex.vars <- model.frame(frmla, dta)
    if( (length(ex.vars) < 1) | (length(ex.vars) > 3) ) stop("only 1 - 3 variables allowed in 'explore'.")
    if( (sum(sapply(ex.vars, is.factor)) > 2) | (sum(sapply(ex.vars, is.numeric)) > 1) ) stop("check your explore-variables")
    
    logit.Y <- as.matrix(inv.logit(x$Y))
    
    ex.model <- lm(update(logit.Y~1, frmla), dta)
    
    .to.plot <- predict(ex.model, data.frame(depth=seq(0,100,length.out=100)))
    .to.plot <- exp(.to.plot) / (1 + exp(.to.plot))
    .to.plot <- .to.plot / rowSums(.to.plot)
  } else {
    .to.plot <- NULL
  }
  
  if(!is.null(model)){
    cat("\nnot implemented yet!\n")
  }


  if(x$dims == 2){
    if(is.null(.main)) .main <- "Density Plot of a Beta-Distributed Variable"
    plot_DRdata_2d(y = x$Y[,2], rug=rug, main=.main, ylim=.ylim, colr=.col, lwd=.lwd, lty=.lty)
  } else if(x$dims == 3) {
    if(!all(is.null(c(.xlim,.ylim)))) warning("xlim and ylim not useable in a ternary plot. arguments ignored.")
    if(is.null(.main)) .main <- "Ternary Plot"

    plot_DRdata_3d(x=x, entropy.contours=entropy.contours, colored=colored, c.grid=c.grid, ticks=ticks, dim.labels=dim.labels, col.scheme=col.scheme,
    main=.main, col=.col, pch=.pch, cex=.cex, lwd=.lwd, lty=.lty,
    .to.plot=.to.plot)
  } else if(x$dims == 4){
    plot_DRdata_4d(x=x, dim.labels=dim.labels,
    main=.main,cex=.cex,
    args.3d=args.3d
    )
  }

}

















