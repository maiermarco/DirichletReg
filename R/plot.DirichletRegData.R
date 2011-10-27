get_or_else <- function(what, or_else, where) if(any(what %in% names(where))) where[[what]] else or_else



plot.DirichletRegData <- function(x,
                                  dims,   # which dimensions to plot
                                  c.grid=TRUE,   # plot a grid?
                                  ticks=TRUE,   # plot ternary ticks?
                                  colored=TRUE,   # colors?
                                  ref.lines=NULL, # reference lines for 2d and 3d plots?
                                  col.scheme=c("dims", "entropy"),   # if colors: which scheme?
                                  entropy.contours=FALSE,   # plot entropy-contour lines?
                                  entropy.colors=FALSE,   # if entropy-contours: plot colored regions?
                                  dim.labels,
                                  args.3d=list(rgl=TRUE, ...),  # theta and phi for the viewport
                                  rug=T,
                                  reset_par=TRUE,
                                  ...){

  if(reset_par){ # reset the current pars after plotting
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
  }

  if(missing(dims))    dims    <- NULL
  
  if(class(x) != "DirichletRegData") stop("data must be prepared by 'DR_data()'")

  col.scheme <- match.arg(col.scheme)

  dotlist <- list(...)
  
  .main <- get_or_else("main", NULL, dotlist)
  .xlim <- get_or_else("xlim", NULL, dotlist)
  .ylim <- get_or_else("ylim", NULL, dotlist)
   .col <- get_or_else("col", NULL, dotlist)
   .pch <- get_or_else("pch", 16, dotlist)
   .cex <- get_or_else("cex", 1, dotlist)
   .lwd <- get_or_else("lwd", 1, dotlist)
   .lty <- get_or_else("lty", 1, dotlist)

      theta <- get_or_else("theta", NULL, dotlist)
        phi <- get_or_else("phi", NULL, dotlist)
  ref.lines <- get_or_else("ref.lines", NULL, dotlist)
   
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


  


  if(x$dims == 2){
    if(is.null(.main)) .main <- "Density Plot of a Beta-Distributed Variable"
    plot_DRdata_2d(y = x$Y[,2], rug=rug, main=.main, ylim=.ylim, colr=.col, lwd=.lwd, lty=.lty)

  } else if(x$dims == 3) {
    if(!all(is.null(c(.xlim,.ylim)))) warning("xlim and ylim not useable in a ternary plot. arguments ignored.")
    if(is.null(.main)) .main <- "Ternary Plot"

    plot_DRdata_3d(x=x, entropy.contours=entropy.contours, colored=colored, c.grid=c.grid, ticks=ticks, dim.labels=dim.labels, col.scheme=col.scheme,
    .main=.main, .col=.col, .pch=.pch, .cex=.cex, .lwd=.lwd, .lty=.lty)

  } else if(x$dims == 4){
    plot_DRdata_4d(x=x, dim.labels=dim.labels, ref.lines=ref.lines,
    main=.main,cex=.cex,
    args.3d=args.3d, theta=theta, phi=phi
    )
  }

}
